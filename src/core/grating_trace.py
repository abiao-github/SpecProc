
"""
Grating order tracing for spectral reduction pipeline.

Handles order detection and tracing for grating spectrographs,
where orders run horizontally (left to right) across the detector.
"""

import numpy as np
import logging
from typing import List, Tuple, Dict, Optional
from scipy import signal, ndimage
from scipy.interpolate import InterpolatedUnivariateSpline
from numpy.polynomial import Chebyshev
from src.core.data_structures import ApertureSet, ApertureLocation
from src.config.config_manager import ConfigManager

logger = logging.getLogger(__name__)


def find_local_minima(array: np.ndarray, window: int = 10) -> Tuple[np.ndarray, np.ndarray]:
    """
    Find local minima in a 1D array.

    Args:
        array: Input array
        window: Minimum distance between peaks

    Returns:
        Tuple of (indices, values) of local minima
    """
    # Use scipy.signal.argrelextrema to find local minima
    from scipy.signal import argrelextrema
    minima_indices = argrelextrema(array, np.less, order=window)[0]
    minima_values = array[minima_indices]
    return minima_indices, minima_values


def find_grating_orders(data: np.ndarray, mask: Optional[np.ndarray] = None,
                      config: Optional[ConfigManager] = None,
                      scan_step: int = 50, minimum: float = 1e-3,
                      separation: float = 20, align_deg: int = 2,
                      filling: float = 0.3, degree: int = 3,
                      conv_core: str = 'auto', fill: bool = False,
                      fill_tol: int = 10, recenter: str = 'parabola') -> ApertureSet:
    """
    Find the positions of grating orders on a CCD image.

    This function implements a robust order tracing algorithm similar to gamse,
    adapted for grating spectrographs where orders run horizontally.

    Args:
        data: Image data (2D array)
        mask: Optional bad pixel mask (same shape as data)
        config: Configuration manager
        scan_step: Steps of pixels used to scan along the dispersion direction
        minimum: Minimum value to filter the input image
        separation: Estimated order separation (in pixels) along cross-dispersion direction
        align_deg: Degree of polynomial used to find the alignment of cross-sections
        filling: Fraction of detected pixels to total step of scanning
        degree: Degree of polynomials to fit aperture locations
        conv_core: Width of convolution core ('auto', int, or None)
        fill: If True, use interpolation to find missing orders
        fill_tol: Tolerance in pixels for finding missing orders
        recenter: Method for recentering ('no', 'parabola', 'threshold', or 'nearest')

    Returns:
        ApertureSet with detected orders
    """
    if mask is None:
        mask = np.zeros_like(data, dtype=np.int32)

    h, w = data.shape
    allx = np.arange(w)
    ally = np.arange(h)

    # Filter pixels smaller than the input "minimum" value
    logdata = np.log10(np.maximum(data, minimum))

    # Initialize order separation as a function of Y position
    try:
        separation = float(separation)
        fsep = lambda y: separation
    except:
        # If separation is a string like "1000:20,2000:25,3000:30"
        # parse it to create a polynomial
        x_lst, y_lst = [], []
        for substring in separation.split(','):
            g = substring.split(':')
            x, y = float(g[0]), float(g[1])
            x_lst.append(x)
            y_lst.append(y)
        x_lst = np.array(x_lst)
        y_lst = np.array(y_lst)
        # Sort x and y
        index = x_lst.argsort()
        x_lst = x_lst[index]
        y_lst = y_lst[index]
        from numpy.polynomial import Polynomial
        separation_poly = Polynomial.fit(x_lst, y_lst, deg=x_lst.size-1)
        fsep = separation_poly

    # Convolution core for cross-sections
    if conv_core is None:
        core = None
    elif isinstance(conv_core, str) and conv_core == 'auto':
        core = np.hanning(int(fsep(h/2)))
    elif isinstance(conv_core, int):
        if conv_core > 0:
            core = np.hanning(conv_core)
        elif conv_core == 0:
            core = None
        else:
            logger.warning(f'Warning: convolution core ({conv_core}) must be positive.')
            core = None
    elif isinstance(conv_core, float):
        core = np.hanning(abs(conv_core))
    else:
        logger.warning(f'Warning: Unknown datatype of conv_core: {conv_core}')
        core = None

    # Normalize the convolution core
    if core is not None:
        core /= core.sum()

    # Scan the image along Y axis starting from the middle row
    y0 = h // 2
    y_lst = {-1: [], 1: []}
    y1 = y0
    direction = -1
    density = 10
    irow = 0
    peak_lst = []

    # Generate a window list according to separations
    dense_x = np.linspace(0, w-1, (w-1)*density+1)
    separation_lst = fsep(dense_x)
    separation_lst = np.int32(np.round(separation_lst))
    window = 2*separation_lst*density+1

    # Two-directional parameter list
    param_lst = {-1: [], 1: []}
    nodes_lst = {}

    # Forward and backward functions for polynomial fitting
    def forward(x, p):
        deg = len(p) - 1
        res = p[0]
        for i in range(deg):
            res = res*x + p[i+1]
        return res

    def forward_der(x, p):
        deg = len(p) - 1
        p_der = [(deg-i)*p[i] for i in range(deg)]
        return forward(x, p_der)

    def backward(y, p):
        x = y
        for ite in range(20):
            dy = forward(x, p) - y
            y_der = forward_der(x, p)
            dx = dy/y_der
            x = x - dx
            if (abs(dx) < 1e-7).all():
                break
        return x

    def fitfunc(p, interfunc, n):
        return interfunc(forward(np.arange(n), p[0:-1])) + p[-1]

    def resfunc(p, interfunc, flux0, mask=None):
        res_lst = flux0 - fitfunc(p, interfunc, flux0.size)
        if mask is None:
            mask = np.ones_like(flux0, dtype=bool)
        return res_lst[mask]

    def find_shift(flux0, flux1, deg):
        p0 = [0.0 for i in range(deg+1)]
        p0[-3] = 1.0

        interfunc = InterpolatedUnivariateSpline(
                    np.arange(flux1.size), flux1, k=3, ext=3)
        mask = np.ones_like(flux0, dtype=bool)
        clipping = 5.
        for i in range(10):
            from scipy.optimize import leastsq
            p, _ = leastsq(resfunc, p0, args=(interfunc, flux0, mask))
            res_lst = resfunc(p, interfunc, flux0)
            std = res_lst.std()
            mask1 = res_lst < clipping*std
            mask2 = res_lst > -clipping*std
            new_mask = mask1*mask2
            if new_mask.sum() == mask.sum():
                break
            mask = new_mask
        return p, mask

    def find_local_peak(xdata, ydata, mask, smooth=None):
        if smooth is not None:
            core = np.hanning(min(smooth, ydata.size))
            core = core/core.sum()
            ydata = np.convolve(ydata, core, mode='same')
        argmax = ydata.argmax()
        xmax = xdata[argmax]
        if argmax < 2 or argmax > ydata.size-2:
            return xdata[xdata.size//2]
        coeff = np.polyfit(xdata[argmax-1:argmax+2], ydata[argmax-1:argmax+2], deg=2)
        a, b, c = coeff
        ypeak = -b/2/a
        return ypeak

    # Main scanning loop
    while True:
        # Get the slice of data in logarithm scale
        flux1 = logdata[y1, :]
        mask1 = (mask[y1, :] == 0)

        # Interpolate to fill bad pixels
        fixfunc = InterpolatedUnivariateSpline(
                np.arange(w)[mask1], flux1[mask1], k=3, ext=3)
        flux1 = fixfunc(np.arange(w))

        # Get the slice of data in linear scale by local median
        linflux1 = np.median(data[y1-2:y1+3, :], axis=0)

        # Convolve the cross-section if core is not None
        if core is None:
            convflux1 = flux1.copy()
        else:
            convflux1 = np.convolve(flux1, core, mode='same')

        if irow == 0:
            convflux1_center = convflux1

        # Find peaks with X precision of 1./density pixels
        n = convflux1.size
        f = InterpolatedUnivariateSpline(np.arange(n), convflux1, k=3)
        convflux2 = f(dense_x)
        imax, fmax = find_local_minima(-convflux2, window=window)
        xmax = dense_x[imax]
        fmax = -fmax

        if irow == 0:
            for x, f in zip(xmax, fmax):
                peak_lst.append((x, f))
                nodes_lst[y1] = []
                for x_val in xmax:
                    nodes_lst[y1].append(x_val)

        else:
            # Aperture alignment of each selected row
            param, _ = find_shift(convflux0, convflux1, deg=align_deg)
            param_lst[direction].append(param[0:-1])

            for x, f in zip(xmax, fmax):
                xstep = x
                for param in param_lst[direction][::-1]:
                    xstep = backward(xstep, param)
                peak_lst.append((xstep, f))
                nodes_lst[y1] = []
                for x_val in xmax:
                    xstep_val = x_val
                    for param in param_lst[direction][::-1]:
                        xstep_val = backward(xstep_val, param)
                    nodes_lst[y1].append(xstep_val)

        y1 += direction*scan_step
        if y1 <= 10:
            direction = +1
            y1 = y0 + direction*scan_step
            y_lst[direction].append(y1)
            convflux0 = convflux1_center
            irow += 1
            continue
        elif y1 >= h - 10:
            break
        else:
            y_lst[direction].append(y1)
            convflux0 = convflux1
            irow += 1
            continue

    # Create ApertureSet
    aperture_set = ApertureSet()

    # Generate a 2-D mesh grid
    yy, xx = np.mgrid[:h, :w]

    # Find all unique x positions from nodes
    all_x_positions = sorted(nodes_lst.keys())
    if not all_x_positions:
        logger.warning("No orders detected")
        return aperture_set

    # For each detected order, trace it across the image
    for order_idx, (x_start, f_start) in enumerate(peak_lst):
        if order_idx >= len(all_x_positions):
            break

        xfit, yfit = [], []

        # Trace in both directions from the starting position
        for direction in [-1, 1]:
            xstep = x_start
            for iy, param in enumerate(param_lst[direction]):
                xstep = forward(xstep, param)
                if xstep < 0 or xstep > w - 1:
                    continue
                y1 = y_lst[direction][iy]

                # Recenter if needed
                if recenter == 'parabola':
                    local_sep = fsep(y1)
                    x1 = int(max(0, xstep - local_sep/2))
                    x2 = int(min(w, xstep + local_sep/2))
                    if x2 - x1 <= 5:
                        continue
                    xdata = np.arange(x1, x2)
                    ydata = logdata[y1, x1:x2].mean(axis=0)
                    m = mask[y1, x1:x2]
                    xpeak = find_local_peak(xdata, ydata, m, smooth=15)
                    xfit.append(xpeak)
                    yfit.append(y1)
                elif recenter == 'no':
                    xfit.append(xstep)
                    yfit.append(y1)
                else:
                    xfit.append(xstep)
                    yfit.append(y1)

        # Sort xfit and yfit
        xfit, yfit = np.array(xfit), np.array(yfit)
        if len(xfit) > 0:
            argsort = xfit.argsort()
            xfit, yfit = xfit[argsort], yfit[argsort]

            # Fit Chebyshev polynomial
            if xfit[0] <= scan_step + 10:
                left_domain = 0
            else:
                left_domain = xfit[0] - scan_step

            if xfit[-1] >= w - scan_step - 10:
                right_domain = w - 1
            else:
                right_domain = xfit[-1] + scan_step

            domain = (left_domain, right_domain)

            poly = Chebyshev.fit(xfit, yfit, domain=domain, deg=degree)

            # Initialize aperture position instance
            aperture_loc = ApertureLocation(
                aperture=order_idx + 1,
                order=order_idx + 1,
                center_coef=poly.coef,
                lower_coef=None,  # Will be set later
                upper_coef=None,  # Will be set later
                width=separation if isinstance(separation, (int, float)) else 30.0
            )
            aperture_set.add_aperture(aperture_loc)

    # Set lower and upper boundaries for each aperture
    for aperture_id, aperture in aperture_set.apertures.items():
        # Get center positions
        center_positions = aperture.get_position(np.arange(w))

        # Set lower and upper boundaries
        profile_width = aperture.width
        half_width = profile_width / 2.0

        # For grating spectrographs, boundaries are along Y direction
        aperture.lower_coef = aperture.center_coef.copy()
        aperture.lower_coef[-1] -= half_width
        aperture.upper_coef = aperture.center_coef.copy()
        aperture.upper_coef[-1] += half_width

    logger.info(f"Detected {aperture_set.norders} orders (grating tracing)")

    return aperture_set
