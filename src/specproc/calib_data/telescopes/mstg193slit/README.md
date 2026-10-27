# Muztagh-Ata 1.93 m slit spectrograph — FeAr anchors

`mstg193slit_R_ids.csv` is the automatic FeAr anchor file for the **RR2000V**
setting (0.8″ slit, 4096-pixel dispersion axis, VPH grating with a visibly
curved dispersion).

- 30 strong, isolated Fe/Ar lines (pixel + wavelength + ion + S/N + residual).
- Cubic pixel → wavelength solution fitted to 37 matched lamp peaks; overall
  rms ≈ 0.28 Å, and the selected anchors have |residual| ≤ 0.5 Å.
- The anchors span ≈ 5912–9784 Å.  The cubic extrapolates to ≈ 3976–11677 Å
  over the 4096 pixels (≈ 1.88 Å/pix).  The lamp-line FWHM is ≈ 2.1 pix, i.e.
  R ≈ 2000 near the 8000 Å blaze.
- The blue extrapolation is checked against the night-sky spectrum of the
  science frame: [OI] 5577.3 Å at pixel 976, Na I 5890/5896 at 1153,
  [OI] 6300.3 at 1380, [OI] 6363.8 at 1414 and Hα 6562.8 at 1522 are all
  reproduced within 0.4 Å.  (The very strong sky peak at pixel 2450 is
  OH 8345 Å, not [OI] 6300.)  The red end beyond ≈ 9800 Å remains an
  extrapolation.
- The first two columns (0-based pixel, air wavelength in Å) are exactly the
  `[slit] identification_file` format used by `SpecProc`.

The pipeline loads these anchors, extends them by matching the full
`calib_data/linelists/fear.txt` catalogue and refits the final polynomial.  On
the 2026-10-10 RR2000V FeAr frame this gives 46 matched lines with an rms of
≈ 0.25 Å over 3970–11680 Å; the final cubic still reproduces the sky lines
above to better than 0.5 Å.

## Regenerate

From the repository root (the reference atlas is 0.2 Å/pixel KPNO FeAr):

```bash
specproc-match-anchors \
    --lamp      output/step5_extraction/<FeAr frame>_spec.fits \
    --reference src/specproc/calib_data/linelists/fear.fits \
    --linelist  src/specproc/calib_data/linelists/fear.txt \
    --config    src/specproc/config/default_config_MSTG193_slit.cfg \
    --output    src/specproc/calib_data/telescopes/mstg193slit/mstg193slit_R_ids.csv \
    --plot      docs/mstg193slit_R_anchor_match.png
```

`--config` supplies the approximate `wave_start`/`wave_end` and restricts the
reference catalogue to that range; the blind curvature scan works too, but the
hint makes the match faster and less prone to a red alias.  Always inspect the
diagnostic figure.  Re-run the tool after changing the grating, slit width or
focus.

## 中文说明

`mstg193slit_R_ids.csv` 是 RR2000V 模式（0.8″ 狭缝、4096 像素、光栅色散
明显弯曲）的 FeAr 自动锚点文件，含 30 条强而孤立的 Fe/Ar 线。前两列
（0 基像素、空气波长 Å）就是 `[slit] identification_file` 的格式。

三次多项式由 37 条匹配灯谱线拟合，整体残差约 0.28 Å，入选锚点的残差
均 ≤0.5 Å。锚点覆盖约 5912–9784 Å；多项式外推到整帧约 3976–11677 Å
（约 1.88 Å/像素）。灯谱线半高全宽约 2.1 像素，在 8000 Å 附近对应
R≈2000。蓝端外推已用科学帧的夜天光验证：[OI] 5577.3（像素 976）、
Na I 5890/5896（1153）、[OI] 6300.3（1380）、[OI] 6363.8（1414）、
Hα 6562.8（1522）均在 0.4 Å 内；2450 处的强天空峰是 OH 8345 Å，不是
[OI] 6300。9800 Å 以外的红端仍属外推。

流水线会用 `fear.txt` 全量线表扩展这些锚点并重新拟合；2026-10-10 的测试帧
得到 46 条匹配线，3970–11680 Å 范围内残差约 0.25 Å，最终解对上述天空线
的误差小于 0.5 Å。更换光栅、狭缝宽度或重新调焦后，请用上面的命令重新
生成，并先检查诊断图。
