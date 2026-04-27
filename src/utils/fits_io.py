"""
FITS file I/O utilities for reading and writing spectral data.

Provides functions for loading raw FITS images, writing intermediate
results, and managing FITS headers.
"""

import numpy as np
from astropy.io import fits
from pathlib import Path
from typing import Tuple, Dict, Optional
import logging

logger = logging.getLogger(__name__)


def read_fits_image(filepath: str) -> Tuple[np.ndarray, dict]:
    """
    Read a FITS image file.

    Args:
        filepath: Path to FITS file

    Returns:
        Tuple of (image_data, header_dict)
    """
    try:
        with fits.open(filepath) as hdul:
            if len(hdul) == 0:
                raise ValueError(f"FITS file {filepath} has no HDUs")

            # Get primary data (image)
            image_data = hdul[0].data
            if image_data is None:
                # Try first non-empty HDU
                for hdu in hdul:
                    if hdu.data is not None:
                        image_data = hdu.data
                        break

            # Get header
            header = dict(hdul[0].header)

            return image_data, header
    except Exception as e:
        logger.error(f"Error reading FITS file {filepath}: {e}")
        raise


def read_fits_table(filepath: str, hdu: int = 1) -> Tuple[np.ndarray, dict]:
    """
    Read structured table from FITS file.

    Args:
        filepath: Path to FITS file
        hdu: HDU index (default 1 for binary table)

    Returns:
        Tuple of (table_data, header_dict)
    """
    try:
        with fits.open(filepath) as hdul:
            if len(hdul) <= hdu:
                raise ValueError(f"FITS file {filepath} has no HDU at index {hdu}")

            table_hdu = hdul[hdu]
            data = table_hdu.data
            header = dict(table_hdu.header)

            return data, header
    except Exception as e:
        logger.error(f"Error reading FITS table from {filepath}: {e}")
        raise


def write_fits_image(filepath: str, data: np.ndarray, header: Optional[Dict] = None,
                     overwrite: bool = True, dtype: Optional[str] = None):
    """
    Write image to FITS file.

    Args:
        filepath: Output file path
        data: Image array
        header: Optional FITS header dictionary
        overwrite: Overwrite existing file
        dtype: Optional data type for output (e.g., 'float32', 'float64')
    """
    try:
        # Convert data to specified dtype if provided
        if dtype is not None:
            original_dtype = str(data.dtype)
            data = data.astype(dtype)
            logger.debug(f"Converted data from {original_dtype} to {dtype}")
        
        # Create primary HDU
        hdu = fits.PrimaryHDU(data=data)

        # If writing to a float format, remove BZERO and BSCALE to prevent incorrect scaling
        # by FITS readers. These keywords are for scaling integers.
        if data.dtype.kind == 'f' and header:
            if 'BZERO' in header: del header['BZERO']
            if 'BSCALE' in header: del header['BSCALE']

        # Add header keywords
        if header:
            for key, value in header.items():
                try:
                    # Skip COMMENT and HISTORY - they often contain special formatting
                    # that causes issues with astropy, and are not critical for processing
                    if key.upper() in ('COMMENT', 'HISTORY'):
                        continue

                    # astropy.io.fits considers a single dot as an illegal value.
                    # Replace it with an empty string, which is valid.
                    if isinstance(value, str) and value.strip() == '.':
                        value = ''
                    
                    # Clean up non-ASCII characters for FITS header values
                    if isinstance(value, str):
                        # Remove or replace non-ASCII characters
                        cleaned_value = value.encode('ascii', 'ignore').decode('ascii')
                        # Truncate if too long (FITS headers have max 80 chars per record)
                        if len(cleaned_value) > 70:
                            cleaned_value = cleaned_value[:67] + '...'
                        hdu.header[key] = cleaned_value
                    else:
                        hdu.header[key] = value
                except Exception as e:
                    logger.warning(f"Could not add header keyword {key}: {e}")

        # Write to file
        hdul = fits.HDUList([hdu])
        hdul.writeto(filepath, overwrite=overwrite)

        logger.info(f"Wrote image to {filepath}")
    except Exception as e:
        logger.error(f"Error writing FITS image {filepath}: {e}")
        raise


def write_fits_table(filepath: str, data: np.ndarray, header: Optional[Dict] = None,
                     overwrite: bool = True):
    """
    Write structured array as binary table to FITS file.

    Args:
        filepath: Output file path
        data: Structured numpy array with named fields
        header: Optional FITS header
        overwrite: Overwrite existing file
    """
    try:
        # Create binary table HDU from structured array
        cols = []
        for name in data.dtype.names:
            col_data = data[name]
            if col_data.dtype == np.object_:
                # Handle variable-length columns
                format_str = f'P{len(col_data[0])}'
                col = fits.Column(name=name, format=format_str, array=col_data)
            else:
                col = fits.Column(name=name, format=f'{len(col_data[0]) if col_data.ndim > 1 else 1}D',
                                array=col_data)
            cols.append(col)

        table_hdu = fits.BinTableHDU.from_columns(cols)

        # Add header keywords
        if header:
            for key, value in header.items():
                try:
                    table_hdu.header[key] = value
                except Exception as e:
                    logger.warning(f"Could not add header keyword {key}: {e}")

        # Create HDU list with primary + table
        primary_hdu = fits.PrimaryHDU()
        hdul = fits.HDUList([primary_hdu, table_hdu])
        hdul.writeto(filepath, overwrite=overwrite)

        logger.info(f"Wrote table to {filepath}")
    except Exception as e:
        logger.error(f"Error writing FITS table {filepath}: {e}")
        raise


def get_header_value(header: dict, keyword: str, default: any = None) -> any:
    """
    Get header value with case-insensitive lookup.

    Args:
        header: FITS header dictionary
        keyword: Keyword name
        default: Default value if not found

    Returns:
        Header value or default
    """
    # Try exact match
    if keyword in header:
        return header[keyword]

    # Try case-insensitive
    for key, value in header.items():
        if key.upper() == keyword.upper():
            return value

    return default


def get_image_info(filepath: str) -> Dict:
    """
    Get basic information about a FITS image.

    Args:
        filepath: Path to FITS file

    Returns:
        Dictionary with image info (shape, datatype, header keywords)
    """
    try:
        img, header = read_fits_image(filepath)
        return {
            'filepath': filepath,
            'shape': img.shape,
            'dtype': str(img.dtype),
            'exptime': get_header_value(header, 'EXPTIME', 0.0),
            'date_obs': get_header_value(header, 'DATE-OBS', ''),
            'object': get_header_value(header, 'OBJECT', ''),
        }
    except Exception as e:
        logger.error(f"Error getting info from {filepath}: {e}")
        return {}


def combine_fits_images(*filepaths: str, method: str = 'mean') -> np.ndarray:
    """
    Combine multiple FITS images.

    Args:
        filepaths: Paths to FITS files
        method: 'mean', 'median', or 'sum'

    Returns:
        Combined image array
    """
    images = []
    for filepath in filepaths:
        img, _ = read_fits_image(filepath)
        images.append(img)

    images = np.array(images)

    if method == 'mean':
        return np.mean(images, axis=0)
    elif method == 'median':
        return np.median(images, axis=0)
    elif method == 'sum':
        return np.sum(images, axis=0)
    else:
        raise ValueError(f"Unknown combine method: {method}")
