import os
import re

# File paths
order_tracing_path = "src/core/order_tracing.py"
simple_trace_path = "src/core/grating_trace_simple.py"
legacy_trace_path = "src/core/grating_trace.py"
gui_path = "src/gui/main_window.py"

# Read grating_trace_simple
with open(simple_trace_path, 'r') as f:
    simple_content = f.read()

# Strip imports from the beginning of simple_content to avoid duplication
# or just leave them, module level imports mid-file are harmless in Python,
# but to be clean, let's keep all.

# Read order_tracing
with open(order_tracing_path, 'r') as f:
    ot_content = f.read()

# Remove the import statements from order_tracing
ot_content = re.sub(r'from src\.core\.grating_trace_simple import.*?\n', '', ot_content)

# Append simple_content
ot_content = ot_content + "\n\n" + "# " + "="*70 + "\n"
ot_content = ot_content + "# MERGED ALGORITHM KERNEL FROM grating_trace_simple.py\n"
ot_content = ot_content + "# " + "="*70 + "\n\n"
ot_content = ot_content + simple_content

with open(order_tracing_path, 'w') as f:
    f.write(ot_content)

# Fix gui file
with open(gui_path, 'r') as f:
    gui_content = f.read()

gui_content = gui_content.replace(
    "from src.core.grating_trace_simple import extend_apertures_for_background",
    "from src.core.order_tracing import extend_apertures_for_background"
)

with open(gui_path, 'w') as f:
    f.write(gui_content)

# Remove the files
os.remove(simple_trace_path)
os.remove(legacy_trace_path)

print("Merge complete.")
