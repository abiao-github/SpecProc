import re

with open('src/core/order_tracing.py', 'r') as f:
    lines = f.readlines()

def get_method_lines(lines, method_name):
    start = -1
    for i, line in enumerate(lines):
        if line.strip().startswith(f"def {method_name}(self,"):
            start = i
            break
    if start == -1:
        return []
    
    end = -1
    for i in range(start + 1, len(lines)):
        if lines[i].startswith("    def ") and not lines[i].startswith("        "):
            end = i
            break
    if end == -1:
        end = len(lines)
    
    return lines[start:end]

m1 = get_method_lines(lines, '_smooth_order_response')
m2 = get_method_lines(lines, '_smooth_blaze_sparse_bspline')
m3 = get_method_lines(lines, 'build_order_response_map')
m4 = get_method_lines(lines, 'save_step4_diagnostics')

# Remove them from order_tracing.py
new_order_tracing_lines = []
skip = False
for line in lines:
    if line.strip().startswith("def _smooth_order_response(self,"): skip = True
    if line.strip().startswith("def _smooth_blaze_sparse_bspline(self,"): skip = True
    if line.strip().startswith("def build_order_response_map(self,"): skip = True
    if line.strip().startswith("def save_step4_diagnostics(self,"): skip = True
    
    if skip and line.startswith("    def ") and not line.startswith("        "):
        # Not the start of the ones we skip if they don't match
        if not (line.strip().startswith("def _smooth_order_response(self,") or
                line.strip().startswith("def _smooth_blaze_sparse_bspline(self,") or
                line.strip().startswith("def build_order_response_map(self,") or
                line.strip().startswith("def save_step4_diagnostics(self,")):
            skip = False
            
    if not skip:
        new_order_tracing_lines.append(line)

with open('src/core/order_tracing.py', 'w') as f:
    f.writelines(new_order_tracing_lines)

# Now add them to flat_correction.py
with open('src/core/flat_correction.py', 'r') as f:
    fc_lines = f.readlines()

class_code = """
class FlatCorrectionModelBuilder:
    def __init__(self):
        self.flat_data = None
        self.flat_mask = None
        self.response_map = None
        self.smoothed_model = None
        self.pixel_flat = None
        self.illumination_flat = None
        self.flat_corr_2d = None
        self.flat_sens = None
        self.order_diagnostics = {}

"""

import_lines = """
from scipy import ndimage
from scipy.signal import savgol_filter
from scipy.interpolate import UnivariateSpline
from scipy.ndimage import gaussian_filter1d
from typing import Dict, Tuple, List
"""

# Insert imports after existing imports
import_idx = 0
for i, line in enumerate(fc_lines):
    if line.startswith("from src.core.data_structures"):
        import_idx = i
        break

fc_lines.insert(import_idx, import_lines)

# Append class at the end or before process_flat_correction_stage
process_idx = len(fc_lines)
for i, line in enumerate(fc_lines):
    if line.startswith("def process_flat_correction_stage"):
        process_idx = i
        break

fc_lines.insert(process_idx, class_code + "".join(m1) + "".join(m2) + "".join(m3) + "".join(m4) + "\n\n")

with open('src/core/flat_correction.py', 'w') as f:
    f.writelines(fc_lines)

print("Done")
