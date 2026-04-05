#!/usr/bin/env python3
"""
Generate a static pipeline diagram inspired by gamse documentation.
Uses a clear vertical flow layout with numbered steps.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patheffects as PathEffects

# Create figure
fig, ax = plt.subplots(figsize=(12, 20))
ax.set_xlim(0, 12)
ax.set_ylim(0, 20)
ax.axis('off')

# Define colors for each stage
stage_colors = [
    '#e1f5ff',  # Stage 0 - light blue
    '#fff4e1',  # Stage 1 - light yellow
    '#e8f5e9',  # Stage 2 - light green
    '#fce4ec',  # Stage 3 - light pink
    '#f3e5f5',  # Stage 4 - light purple
    '#fff9c4',  # Stage 5 - light yellow
    '#ffccbc',  # Stage 6 - light orange
    '#d1c4e9',  # Stage 7 - light purple
]

# Stage configurations with detailed steps
stages = [
    {
        'name': 'Stage 0',
        'title': 'Overscan Correction',
        'y': 19,
        'color': stage_colors[0],
        'steps': ['Read raw FITS', 'Extract overscan', 'Subtract bias']
    },
    {
        'name': 'Stage 1',
        'title': 'Bias Subtraction',
        'y': 17,
        'color': stage_colors[1],
        'steps': ['Combine bias frames', 'Generate master bias', 'Subtract from all images']
    },
    {
        'name': 'Stage 2',
        'title': 'Flat Fielding & Order Tracing',
        'y': 15,
        'color': stage_colors[2],
        'steps': ['Combine flat frames', 'Detect orders', 'Fit traces', 'Extract blaze']
    },
    {
        'name': 'Stage 3',
        'title': 'Background Subtraction',
        'y': 13,
        'color': stage_colors[3],
        'steps': ['Estimate background', 'Fit 2D polynomial', 'Subtract background']
    },
    {
        'name': 'Stage 4',
        'title': 'Cosmic Ray Correction',
        'y': 11,
        'color': stage_colors[4],
        'steps': ['Detect cosmic rays', 'Identify pixels', 'Interpolate (median)']
    },
    {
        'name': 'Stage 5',
        'title': '1D Spectrum Extraction',
        'y': 9,
        'color': stage_colors[5],
        'steps': ['Select method', 'Extract spectra', 'Calculate errors']
    },
    {
        'name': 'Stage 6',
        'title': 'Wavelength Calibration',
        'y': 7,
        'color': stage_colors[6],
        'steps': ['Fit 2D polynomial', 'Pixel to wavelength', 'Calibrate all spectra']
    },
    {
        'name': 'Stage 7',
        'title': 'De-blazing',
        'y': 5,
        'color': stage_colors[7],
        'steps': ['Match orders', 'Divide by blaze', 'Normalize continuum']
    }
]

# Draw stages and steps
for i, stage in enumerate(stages):
    # Stage header box
    header_box = mpatches.FancyBboxPatch((0.5, stage['y'] - 0.5), 11, 0.6,
                                             boxstyle="round,pad=0.1",
                                             facecolor=stage['color'],
                                             alpha=0.3,
                                             edgecolor=stage['color'],
                                             linewidth=2)
    ax.add_patch(header_box)
    
    # Stage number and title
    ax.text(1.25, stage['y'], f"{stage['name']}\n{stage['title']}", 
            ha='left', va='center', fontsize=10, fontweight='bold')
    
    # Draw steps as numbered list
    steps = stage['steps']
    for j, step in enumerate(steps, 1):
        y_pos = stage['y'] - 0.4 - j * 0.5
        
        # Step number circle
        circle = plt.Circle((0.8, y_pos), 0.12, 
                           facecolor='white', edgecolor=stage['color'], linewidth=2)
        ax.add_patch(circle)
        
        # Step number
        ax.text(0.8, y_pos, str(j), ha='center', va='center', 
                fontsize=7, fontweight='bold', color=stage['color'])
        
        # Step text
        ax.text(1.2, y_pos, f"• {step}", ha='left', va='center', 
                fontsize=8)

# Draw arrows between stages
for i in range(len(stages) - 1):
    y_top = stages[i]['y'] - 0.5
    y_bottom = stages[i + 1]['y'] + 0.1
    ax.arrow(6, y_top, 0, -1.7,
               head_width=0.25, head_length=0.2, 
               fc='darkgray', ec='darkgray', linewidth=2.5)

# Add title at top
ax.text(6, 19.8, 'SpecProc Spectral Reduction Pipeline', 
        ha='center', va='center', fontsize=16, fontweight='bold')

# Add output arrow and final result
ax.arrow(6, 5.3, 0, -2.3,
           head_width=0.25, head_length=0.2, 
           fc='darkgray', ec='darkgray', linewidth=2.5)

ax.text(6, 2.5, 'Final Output: Calibrated 1D Spectra', 
        ha='center', va='center', fontsize=11, fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.5', facecolor='#4CAF50', alpha=0.8))

# Add color legend at bottom
legend_y = 1.5
legend_items = [
    (stage_colors[0], 'Stage 0: Overscan'),
    (stage_colors[1], 'Stage 1: Bias'),
    (stage_colors[2], 'Stage 2: Flat Fielding'),
    (stage_colors[3], 'Stage 3: Background'),
    (stage_colors[4], 'Stage 4: Cosmic Ray'),
    (stage_colors[5], 'Stage 5: Extraction'),
    (stage_colors[6], 'Stage 6: Wavelength'),
    (stage_colors[7], 'Stage 7: De-blazing'),
]

for i, (color, label) in enumerate(legend_items):
    y_pos = legend_y - i * 0.4
    # Color box
    ax.add_patch(plt.Rectangle((0.3, y_pos - 0.25), 0.25, 0.35,
                            facecolor=color, alpha=0.7, edgecolor=color, linewidth=1.5))
    ax.text(0.6, y_pos, label, ha='left', va='center', fontsize=8)

plt.tight_layout()
plt.savefig('docs/processing_pipeline_diagram.png', dpi=150, bbox_inches='tight', 
            facecolor='white', edgecolor='none')
print("Pipeline diagram saved to docs/processing_pipeline_diagram.png")
