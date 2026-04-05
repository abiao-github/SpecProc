#!/usr/bin/env python3
"""
Generate a static image version of the processing pipeline diagram.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Create figure
fig, ax = plt.subplots(figsize=(12, 16))
ax.set_xlim(0, 12)
ax.set_ylim(0, 16)
ax.axis('off')

# Define colors for each stage
stage_colors = [
    '#e1f5ff',  # Stage 0 - blue
    '#fff4e1',  # Stage 1 - yellow
    '#e8f5e9',  # Stage 2 - green
    '#fce4ec',  # Stage 3 - pink
    '#f3e5f5',  # Stage 4 - purple
    '#fff9c4',  # Stage 5 - yellow
    '#ffccbc',  # Stage 6 - orange
    '#d1c4e9',  # Stage 7 - purple
]

# Stage configurations
stages = [
    {
        'name': 'STAGE 0\nOverscan Correction',
        'y': 14.5,
        'color': stage_colors[0],
        'steps': [
            'Read raw FITS files',
            'Extract overscan region',
            'Calculate median/polynomial',
            'Subtract overscan bias',
            'Output: Overscan corrected'
        ]
    },
    {
        'name': 'STAGE 1\nBias Subtraction',
        'y': 12.5,
        'color': stage_colors[1],
        'steps': [
            'Read overscan corrected',
            'Combine bias frames',
            'Generate master bias',
            'Subtract master bias',
            'Output: Bias corrected'
        ]
    },
    {
        'name': 'STAGE 2\nFlat Fielding',
        'y': 10.5,
        'color': stage_colors[2],
        'steps': [
            'Read bias corrected flat',
            'Combine flat frames',
            'Generate master flat',
            'Detect echelle orders',
            'Fit polynomial traces',
            'Extract blaze profiles',
            'Output: Master flat'
        ]
    },
    {
        'name': 'STAGE 3\nBackground Subtraction',
        'y': 8.5,
        'color': stage_colors[3],
        'steps': [
            'Read bias corrected image',
            'Estimate background',
            'Fit 2D polynomial/Median',
            'Subtract background',
            'Output: Background subtracted'
        ]
    },
    {
        'name': 'STAGE 4\nCosmic Ray Correction',
        'y': 6.5,
        'color': stage_colors[4],
        'steps': [
            'Read background subtracted',
            'Detect cosmic rays',
            'Identify cosmic pixels',
            'Interpolate (median filter)',
            'Output: Cosmic corrected'
        ]
    },
    {
        'name': 'STAGE 5\n1D Extraction',
        'y': 4.5,
        'color': stage_colors[5],
        'steps': [
            'Read cosmic corrected image',
            'Select extraction method',
            'Sum/Optimal extraction',
            'Extract 1D spectra',
            'Calculate errors',
            'Output: SpectraSet'
        ]
    },
    {
        'name': 'STAGE 6\nWavelength Calibration',
        'y': 2.5,
        'color': stage_colors[6],
        'steps': [
            'Step 1: ThAr calibration',
            'Extract ThAr 1D spectrum',
            'Identify emission lines',
            'Fit 2D wavelength polynomial',
            'Step 2: Apply to science',
            'Convert pixels to wavelength',
            'Output: Wavelength calibrated'
        ]
    },
    {
        'name': 'STAGE 7\nDe-blazing',
        'y': 0.5,
        'color': stage_colors[7],
        'steps': [
            'Read wavelength calibrated',
            'Read flat blaze function',
            'Match orders',
            'Divide by blaze function',
            'Normalize to continuum',
            'Output: Final calibrated'
        ]
    }
]

# Draw stages
for i, stage in enumerate(stages):
    # Stage box
    box = mpatches.FancyBboxPatch((1, stage['y'] - 0.7), 10, 1.4,
                                     boxstyle="round,pad=0.1",
                                     facecolor=stage['color'],
                                     alpha=0.3,
                                     edgecolor=stage['color'],
                                     linewidth=2)
    ax.add_patch(box)
    
    # Stage title
    ax.text(6, stage['y'], stage['name'], ha='center', va='center',
            fontsize=11, fontweight='bold')
    
    # Draw steps
    steps = stage['steps']
    step_width = 10 / (len(steps) + 1)
    for j, step in enumerate(steps):
        x_pos = 1 + (j + 1) * step_width
        if j == 0:
            # Draw input arrow
            ax.arrow(6, stage['y'] + 0.6, 0, -1.2,
                   head_width=0.1, head_length=0.1, fc='gray', ec='gray')
        
        # Draw step indicator
        if j < len(steps) - 1:
            # Regular step
            ax.text(x_pos, stage['y'], step, ha='center', va='center',
                    fontsize=7, rotation=45, alpha=0.7)
        else:
            # Output step - make it stand out
            step_box = mpatches.FancyBboxPatch((x_pos - step_width/2 - 0.2, stage['y'] - 0.12),
                                                   step_width + 0.4, 0.24,
                                                   boxstyle="round,pad=0.02",
                                                   facecolor='white',
                                                   edgecolor=stage['color'],
                                                   linewidth=2,
                                                   alpha=0.9)
            ax.add_patch(step_box)
            ax.text(x_pos, stage['y'], step, ha='center', va='center',
                    fontsize=8, fontweight='bold')

# Draw arrows between stages
for i in range(len(stages) - 1):
    ax.arrow(6, stages[i]['y'] - 0.7, 0, -1.3,
               head_width=0.15, head_length=0.1, fc='darkgray', ec='darkgray')

# Add title
ax.text(6, 15.5, 'SpecProc 8-Stage Spectral Reduction Pipeline', 
        ha='center', va='center', fontsize=16, fontweight='bold')

# Add legend
legend_text = (
    "Pipeline Overview:\n\n"
    "Stage 0: Overscan correction (all images)\n"
    "Stage 1: Bias subtraction (median/mean combine)\n"
    "Stage 2: Flat fielding & order tracing\n"
    "Stage 3: Background subtraction\n"
    "Stage 4: Cosmic ray removal (science only)\n"
    "Stage 5: 1D spectrum extraction (sum/optimal)\n"
    "Stage 6: Wavelength calibration (pixel → wavelength)\n"
    "Stage 7: De-blazing (divide by blaze function)"
)

ax.text(12, 8, legend_text, ha='right', va='center', fontsize=8,
        bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.8))

plt.tight_layout()
plt.savefig('docs/processing_pipeline_diagram.png', dpi=150, bbox_inches='tight')
print("Pipeline diagram saved to docs/processing_pipeline_diagram.png")
