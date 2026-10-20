# 校准数据目录说明

本目录包含 SpecProc 所需的所有校准数据文件。

## 目录结构

```
calib_data/
├── linelists/              # 灯谱发射线目录
│   ├── thar-noao.dat      # ThAr 灯谱线（兴隆216 HRS 推荐）
│   ├── thar.dat           # 标准 ThAr 灯谱线
│   ├── FeAr.dat           # FeAr 灯谱线
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
- `FeAr.dat` - Fe-Ar 空心阴极灯谱线（3000-10,000 Å）。
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
specproc-air2vac input.dat output.dat
specproc-vac2air input.dat output.dat
```

安装或更新 SpecProc 后，这两个命令会添加到当前 Python 环境中（在仓库本地以可编辑模式安装可运行 `python -m pip install -e .`）。默认环境参数为温度 `-t 15` ℃、气压 `-p 760` mmHg、水蒸气分压 `-f 4` mmHg，也可以通过这些参数覆盖。工具只转换每行第一个数值列，并保留注释和后续文本。请使用不同的输出文件名，以免覆盖原始线表。

如需不安装、直接从仓库运行，可使用 `PYTHONPATH=src python -m specproc.utils.airvac input.dat output.dat --type air2vac` 或 `--type vac2air`。

### `thar-noao.dat` 的波长类型

`thar-noao.dat` 文件没有注明波长为空气波长还是真空波长。与 `thar.dat` 比较后，证据支持它使用空气波长：在 `thar.dat` 的波长范围（3184.9486-9664.6983 Å）内，有 2,948 个波长数值与 `thar.dat` 完全一致。IRAF 线表说明指出 `vacthar.dat` 是由 `thar.dat` 转换得到的真空波长表，因此这些完全相同的值与两份原始线表均使用空气波长的判断相符。这是有较强依据的推断，并非 `thar-noao.dat` 文件中的明确声明；如能取得原始来源元数据，应以其为准。

### 手动证认谱线（狭缝光谱）

慕士塔格 1.93 米狭缝光谱仪没有预先计算的定标文件，需要先手动标记少量明显的灯谱线：

```bash
specproc-identify --lamp output/step5_extraction/<灯谱文件>_spec.fits --linelist calib_data/linelists/fear.dat
```

（未安装时可用 `python -m specproc.gui.line_identifier`；`--lamp` 和 `--linelist` 可省略，在窗口中载入。）窗口上方显示抽取的一维灯谱（横轴为像素），下方把线表画成零基线上的 δ 函数（横轴为波长）。可人工点击一条灯谱峰，再点击对应参考谱线添加锚定线；右键可删除。设置 **low/upper** 波长范围和**信噪比阈值**后，即使没有人工配对，也可直接点击 **Auto identify**。程序根据谱线间隔的形态（最强灯谱峰的三线间距之比与线表中三线间距之比相比对，正反色散方向均尝试）寻找初始解，并剔除波长覆盖超出所设范围的解，再迭代匹配高于阈值的灯谱峰。人工配对（如有）会作为锚点保留。也可勾选 **Overlay reference lines**：用当前模型（或 low/upper 范围，波长随像素递减时勾选 **Flip**）把线表以竖线叠加在灯谱上，用 **Shift**/**Stretch**（或方向键：←/→ 平移，↑/↓ 拉伸）手动微调对齐，再点击 **Fit overlay** 匹配谱峰并拟合。保存后将 `[slit]` 中的 `identification_file` 设为该文件（`peak_snr` 设为同一阈值）：流水线据此拟合，再用 `full_linelist` 匹配并重新拟合最终多项式。线表的波长类型（空气或真空）须与所标注的波长一致。

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
