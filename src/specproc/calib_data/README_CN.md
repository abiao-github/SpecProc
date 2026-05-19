# 校准数据目录说明

本目录包含 SpecProc 所需的所有校准数据文件。

## 目录结构

```
calib_data/
├── linelists/              # 灯谱发射线目录
│   ├── thar-noao.dat      # ThAr 灯谱线（兴隆216 HRS 推荐）
│   ├── thar.dat           # 标准 ThAr 灯谱线
│   ├── FeAr.dat           # FeAr 灯谱线
│   ├── thar_lines_3500_4500.csv
│   ├── ar_lines_4000_7000.csv
│   ├── ne_lines_5800_6700.csv
│   ├── he_lines_3800_7300.csv
│   └── fe_lines_4000_4500.csv
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
- `thar-noao.dat` - ThAr 灯谱线（兴隆216 HRS 推荐使用）
- `thar.dat` - 标准 ThAr 灯谱线
- `FeAr.dat` - FeAr 灯谱线
- `thar_lines_3500_4500.csv` - ThAr 灯谱线（CSV格式，3500-4500 Å）
- `ar_lines_4000_7000.csv` - Ar 灯谱线（CSV格式，4000-7000 Å）
- `ne_lines_5800_6700.csv` - Ne 灯谱线（CSV格式，5800-6700 Å）
- `he_lines_3800_7300.csv` - He 灯谱线（CSV格式，3800-7300 Å）
- `fe_lines_4000_4500.csv` - Fe 灯谱线（CSV格式，4000-4500 Å）

**支持的灯类型**：
- `ThAr` - 氩钍灯（阶梯光谱仪最常用）
- `FeAr` - 铁氩灯
- `Ar` - 氩灯
- `Ne` - 氖灯
- `He` - 氦灯
- `Fe` - 铁灯

**文件格式**：
- `.dat` 格式（标准格式，包含波长、强度、注释）
- `.csv` 格式（波长、强度、注释）

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
2. 遵循文件格式（.dat 或 .csv）
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
   - 灯谱线文件支持 `.dat` 和 `.csv` 两种格式
   - 校准文件必须是 FITS 格式（`.fits`）
