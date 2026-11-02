# 校准数据目录说明

本目录包含 SpecProc 所需的所有校准数据文件。

## 目录结构

```
calib_data/
├── linelists/              # 灯谱发射线目录
│   ├── thar-noao.dat      # ThAr 灯谱线（兴隆216 HRS 推荐）
│   ├── thar.dat           # 标准 ThAr 灯谱线
│   ├── fear.dat           # FeAr 灯谱线表（Willmarth）
│   ├── fear.fits          # FeAr 高分辨率图集光谱（NOIRLAB/KPNO）
│   ├── henear*.dat        # He-Ne-Ar 灯谱线
│   └── ...                # 其他灯谱线和夜空谱线
└── telescopes/             # 望远镜特定校准文件
    ├── generic/           # 通用配置模板
    └── xinglong216hrs/    # 兴隆216望远镜
        ├── wlcalib_20141103049.fits
        ├── wlcalib_20171202012.fits
        ├── wlcalib_20190905028_A.fits
        └── wlcalib_20211123011_A.fits
```

## 灯谱线文件

灯谱线文件存储在 `linelists/` 目录中。

**可用的灯谱线文件**：
- `thar-noao.dat` - Th-Ar 灯谱线；推荐用于兴隆 2.16 米 HRS。
- `thar.dat` - Th-Ar 灯谱线（3180-9665 Å）。
- `vacthar.dat` - `thar.dat` 转换得到的真空波长谱线。
- `fear.dat` - Fe-Ar 空心阴极灯谱线表（3000-10,000 Å，Willmarth）。
- `fear.fits` - Fe-Ar 图集光谱（NOIRLAB/KPNO，0.2 Å 采样，带线性波长 WCS）。
- `cuar.dat` - Cu-Ar 灯谱线。
- `henear.dat`、`henearhres.dat`、`idhenear.dat` - He-Ne-Ar 灯谱线，包含高分辨率和扩展谱线表。
- `vacidhenear.dat` - `idhenear.dat` 转换得到的真空波长谱线。
- `ctiohenear.dat`、`ctiohear.dat` - CTIO 使用的 He-Ne-Ar 和 He-Ar 谱线表。
- `argon.dat`、`krypton.dat`、`xenon.dat` - Ar、Kr、Xe 的红外发射谱线。
- `skylines.dat`、`lowskylines.dat`、`ohlines.dat` - 夜空谱线及 OH 发射谱线。
- `thorium.dat`、`vacthorium.dat` - 已弃用的旧文件，分别请使用 `thar.dat` 和 `vacthar.dat`。

这些目录文件为纯文本，列内容因谱线表而异。`airvac.cl` 是空气/真空波长转换工具，不是谱线表。请根据灯类型和波长范围选择合适的谱线表；目录中没有 CSV 格式的谱线文件。

### 空气波长与真空波长转换

SpecProc 提供了基于 IRAF `airvac.cl` 的波长转换命令，公式参考 Allen 的 *Astrophysical Quantities*（1973）。波长单位为 Å。在仓库根目录运行：

```bash
specproc air2vac input.dat output.dat
specproc vac2air input.dat output.dat
```

同样的功能也保留为独立命令 `specproc.air2vac` 和 `specproc.vac2air`。安装或更新 SpecProc 后，这些命令会添加到当前 Python 环境中（在仓库本地以可编辑模式安装可运行 `python -m pip install -e .`）。默认环境参数为温度 `-t 15` ℃、气压 `-p 760` mmHg、水蒸气分压 `-f 4` mmHg，也可以通过这些参数覆盖。工具只转换每行第一个数值列，并保留注释和后续文本。请使用不同的输出文件名，以免覆盖原始线表。

如需不安装、直接从仓库运行，可使用 `PYTHONPATH=src python -m specproc.utils.airvac input.dat output.dat --type air2vac` 或 `--type vac2air`。

### `thar-noao.dat` 的波长类型

`thar-noao.dat` 文件没有注明波长为空气波长还是真空波长。与 `thar.dat` 比较后，证据支持它使用空气波长：在 `thar.dat` 的波长范围（3184.9486-9664.6983 Å）内，有 2,948 个波长数值与 `thar.dat` 完全一致。IRAF 线表说明指出 `vacthar.dat` 是由 `thar.dat` 转换得到的真空波长表，因此这些完全相同的值与两份原始线表均使用空气波长的判断相符。这是有较强依据的推断，并非 `thar-noao.dat` 文件中的明确声明；如能取得原始来源元数据，应以其为准。

### 谱线证认与锚点生成（狭缝光谱，推荐）

慕士塔格 1.93 米狭缝光谱仪没有预先计算的定标文件。`specproc identify` 把实测灯谱与参考
图集（`fear.fits`）或线表对照：窗口上方把两者叠加在同一个像素轴上（参考谱线按当前色散
模型映射到像素坐标），下方给出锚点线的 Legendre/Chebyshev 拟合残差（基函数、阶数和残差裁剪可调）。点
**Auto identify** 会对强灯谱峰做色散曲率扫描与形态匹配，结果还可以人工增删；保存时同时
写出锚点 CSV 和诊断对比图。锚点文件的前两列（0 基像素、波长 Å）就是 `[slit]
identification_file` 所需格式，后续流水线会用 `full_linelist` 把锚点扩展到全量谱线并重新
拟合。波长随像素递减时可勾选 **Flip x** 翻转横向显示（自动匹配会自行判断方向）。

```bash
# 交互式：实测灯谱 + 参考图集
specproc identify \
    --lamp      output/step5_extraction/<灯谱文件>_pixel.fits \
    --ref src/specproc/calib_data/linelists/fear.fits \
    --linelist  src/specproc/calib_data/linelists/fear.dat \
    --output    src/specproc/calib_data/telescopes/mstg193slit/mstg193slit_anchors_R.csv \
    --plot      docs/mstg193slit_R_anchor_match.pdf
# 慕峰红/蓝端由 FITS 头的 TELESCOP=MOST + INSTRUME=Camera Spectrum-R/B
# 自动判断；--arm B 仅用于显式覆盖。默认文件名自动变为
# mstg193slit_anchors_B.csv。

# 全自动，不开窗口
specproc identify --auto --lamp <灯谱.fits> --ref <图集.fits> \
    --output <锚点.csv> --plot <诊断图.pdf>
```

工具会扫描色散曲率（光栅色散不是严格线性），用最强的一批参考谱线做形态匹配，再用三次
多项式迭代精修；不需要人工点选，也不需要 `wave_start`/`wave_end` 提示（加 `--config`
可给出近似范围提示；工具不会按该范围裁剪参考线表，因此范围写偏也不会被逼到假解）。两个
像素方向都会尝试，波长随像素递减的灯谱同样能自动匹配，窗口中会自动翻转横轴，使波长从左
到右增加。RR2000V 的一帧，三次多项式（默认 Legendre 基，可配置为 Chebyshev/5 阶）
外推到约 3976–11677 Å，锚点本身覆盖约 5912–9784 Å；锚点残差约 0.3 Å。蓝端外推已用
夜天光验证（[OI] 5577/6300/6363、Na I、Hα 均在 0.4 Å 内），9800 Å 以外的红端仍是
外推值。**使用前务必看一眼诊断图**（上：灯谱与锚点；中：观测/参考局部叠合；下：残差），
确认匹配无误。若某帧锚点太少，可调低 `--catalog-frac`（如 0.005）或
`--catalog-min-snr` 以保留更多参考线，再在图中核对新增的谱线。若更换光栅、狭缝或
重新调焦，应重新生成锚点文件。

**人工调整**：若自动匹配因灯谱质量差而失败，可在窗口中先人工点选少量明显的灯谱线（点灯谱峰，再点对应参考谱线；右键删除），再点 **Auto identify** —— 人工配对会作为初值。保持 **Show reference** 勾选即可：拟合前参考谱显示在自己的波长轴上，有 ≥2 条锚点建立模型后自动映射到 pixel 轴；需要左右翻转时勾 **Flip x**。保存后 `[slit]` 中的 `identification_file` 指向该 CSV（`peak_snr` 设为同一阈值），流水线再用 `full_linelist` 匹配并重新拟合最终多项式。线表的波长类型（空气或真空）须与所标注的波长一致。（未安装时可用 `PYTHONPATH=src python -m specproc.utils.identify`。）

## 望远镜校准文件

望远镜特定的校准文件存储在 `telescopes/<望远镜名称>/` 目录中。

**文件命名规则**：`<光谱仪>_wlcalib_<YYYYMMDDHH>_<后缀>.fits`

其中：
- `<光谱仪>`：光谱仪名称（例如 hrs 表示高分辨率光谱仪）
- `YYYYMMDDHH`：校准日期和小时
- `<后缀>`：可选后缀（例如 A、B、C 表示不同版本的校准）

**兴隆216 HRS 可用的校准文件**：
- `wlcalib_20141103049.fits` - 2014-11-03 04:50
- `wlcalib_20171202012.fits` - 2017-12-02 01:20
- `wlcalib_20190905028_A.fits` - 2019-09-05 02:50（版本 A）
- `wlcalib_20211123011_A.fits` - 2021-11-23 01:10（版本 A）- **最新**

这些文件包含：
- 2D 波长多项式解：λ(x,y) = Σ p_ij·x^i·y^j
- 已识别的发射谱线及其波长
- 波长拟合的 RMS 值

## 配置

望远镜和光谱仪的配置在 `default_config.cfg` 文件中设置：

```ini
[telescope]
# 望远镜名称（用于查找校准数据）
name = xinglong216hrs

# 光谱仪名称
instrument = hrs

[telescope.linelist]
# 使用的灯类型
linelist_type = ThAr

# 灯谱线文件路径
linelist_path = calib_data/linelists/

# 使用的具体灯谱线文件（可选）
# 兴隆216 HRS 推荐使用：thar-noao.dat
linelist_file = thar-noao.dat

# 使用预先识别的校准文件（可选）
use_precomputed_calibration = no
calibration_path = calib_data/telescopes/xinglong216hrs/

# 使用的具体校准文件（可选）
# 使用最新的：wlcalib_20211123011_A.fits
calibration_file = wlcalib_20211123011_A.fits
```

## 添加自定义校准数据

### 添加新的灯谱线文件：

1. 在 `linelists/` 目录中创建合适的文件
2. 参考同类 `.dat` 文件的格式；不同谱线表的列内容可能不同
3. 在 `default_config.cfg` 中配置：`linelist_file = <文件名>`

### 添加新的望远镜配置：

1. 创建目录：`calib_data/telescopes/<望远镜名称>/`
2. 将校准文件放入该目录，遵循命名规则
3. 在 `default_config.cfg` 中配置：
   ```ini
   [telescope]
   name = <望远镜名称>
   instrument = <光谱仪名称>
   ```

### 添加预计算的波长校准：

1. 将 FITS 校准文件复制到对应的望远镜目录
2. 确保文件遵循命名规则：`<光谱仪>_wlcalib_<YYYYMMDDHH>_<后缀>.fits`
3. 在 `default_config.cfg` 中启用：`use_precomputed_calibration = yes`
4. 指定使用的校准文件：`calibration_file = <文件名>`

## 使用说明

### 处理兴隆216 HRS 数据时的配置：

对于使用 ThAr 灯的兴隆216 HRS 数据，推荐配置：

```ini
[telescope]
name = xinglong216hrs
instrument = hrs

[telescope.linelist]
linelist_type = ThAr
linelist_file = thar-noao.dat
linelist_path = calib_data/linelists/

use_precomputed_calibration = no
calibration_path = calib_data/telescopes/xinglong216hrs/
calibration_file = wlcalib_20211123011_A.fits
```

### 使用最新的校准文件：

如果需要使用最新的校准文件（2021年11月），设置：

```ini
use_precomputed_calibration = yes
calibration_file = wlcalib_20211123011_A.fits
```

### 重新拟合波长校准：

如果需要重新拟合波长校准（使用自己的 ThAr 灯谱图像），设置：

```ini
use_precomputed_calibration = no
```

系统将使用指定的灯谱线文件（`thar-noao.dat`）从 ThAr 灯谱图像中识别发射谱线并拟合波长多项式。

## 注意事项

1. **灯谱线选择**：
   - 兴隆216 HRS 推荐使用 `thar-noao.dat`，这是针对该光谱仪优化的灯谱线
   - 其他灯谱线文件可能不包含所需的特定波长范围

2. **校准文件选择**：
   - 尽量使用与观测时间接近的校准文件
   - 最新的校准文件是 `wlcalib_20211123011_A.fits`（2021年11月）
   - 如果观测时间与校准时间相差较大，考虑重新拟合波长校准

3. **文件格式兼容性**：
   - 谱线目录当前提供纯文本 `.dat` 文件；各文件的列格式可能不同
   - 校准文件必须是 FITS 格式（`.fits`）
