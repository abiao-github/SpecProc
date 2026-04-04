# SpecProc 配置文件说明

## 配置文件类型

SpecProc 支持两种配置文件：

### 1. 默认配置文件

**位置**：`SpecProc/default_config.cfg`
**用途**：提供默认参数值
**修改**：不建议直接修改

### 2. 用户配置文件

**位置**：在工作目录中创建 `specproc.cfg`
**用途**：覆盖默认配置，自定义参数
**优先级**：用户配置 > 默认配置

## 配置文件示例

### 完整的用户配置文件

```ini
# SpecProc User Configuration
# Copy this file to specproc.cfg in your working directory

# ============================================================================
# 数据输入配置
# ============================================================================
[data]
# 原始数据目录（相对于当前工作目录）
# 例如：在 /home/user/obs1/ 目录运行 specproc
# 如果 rawpath=20241102_hrs，则数据从 /home/user/obs1/20241102_hrs/ 加载
rawpath = 20241102_hrs

# FITS 头部关键字（通常不需要修改）
statime_key = DATE-OBS
exptime_key = EXPTIME

# 色散方向
# 选项：xr-（x 正方向为色散），xl-（x 负方向为色散）
direction = xr-

# 过扫描区域配置
overscan_start_column = 4097
overscan_method = median

# ============================================================================
# 通用处理选项
# ============================================================================
[reduce]
# 输出目录（相对于当前工作目录）
# 例如：在 /home/user/obs1/ 目录运行 specproc
# 如果 out_path=output，则结果保存到 /home/user/obs1/output/
#
# 输出目录结构：
# output/
#   ├── overscan_corrected/      # 第0步：过扫描校正结果
#   ├── bias_corrected/          # 第1步：偏置减除结果
#   ├── flat_corrected/           # 第2步：平场改正结果
#   ├── background_corrected/     # 第4步：背景扣除结果
#   ├── cosmic_corrected/        # 第5步：宇宙线去除结果（仅科学图像）
#   ├── spectra/                 # 最终一维光谱
#   └── figures/                 # 诊断图像
out_path = output

# 处理模式：normal 或 debug
mode = normal

# 输出图像格式：png、pdf、jpg
fig_format = png

# 最终光谱文件名后缀
oned_suffix = _ods

# 处理核心数
ncores = max

# 宇宙线校正参数
cosmic_enabled = yes
cosmic_sigma = 5.0
cosmic_window = 7

# ============================================================================
# 中间结果保存选项
# ============================================================================
[reduce.save_intermediate]
# 每一步是否保存中间结果
# 可以单独控制每一步是否保存
# 在 GUI 中也应该有对应的勾选选项
save_overscan = yes        # 第0步：过扫描校正
save_bias = yes             # 第1步：偏置减除
save_flat = yes              # 第2步：平场改正
save_background = yes        # 第4步：背景扣除
save_cosmic = yes           # 第5步：宇宙线去除
save_extraction = yes       # 第6步：一维谱提取
save_wlcalib = yes          # 第7步：波长定标
save_deblaze = yes          # 第8步：Blaze 函数改正

# ============================================================================
# 望远镜/光谱仪配置
# ============================================================================
[telescope]
# 望远镜名称（用于查找校准数据）
name = xinglong216hrs

# 光谱仪名称（用于校准文件查找）
instrument = hrs

[telescope.linelist]
# 灯谱线类型
linelist_type = ThAr

# 灯谱线文件路径
linelist_path = calib_data/linelists/

# 使用的具体灯谱线文件（可选）
# 对于兴隆216 HRS：thar-noao.dat 是推荐的
linelist_file = thar-noao.dat

# 是否使用预先识别的校准文件（可选）
use_precomputed_calibration = yes
calibration_path = calib_data/telescopes/xinglong216hrs/

# 使用的具体校准文件（可选）
# 使用最新的：wlcalib_20211123011_A.fits
calibration_file = wlcalib_20211123011_A.fits
```

## 配置参数说明

### 路径参数

| 参数 | 说明 | 示例 |
|------|------|------|
| `rawpath` | 原始数据目录（相对路径） | `20241102_hrs` |
| `out_path` | 输出目录（相对路径） | `output` |

### 路径相对性说明

所有路径都是相对于**当前工作目录**的：

```bash
# 假设工作目录是 /home/user/obs1/2024/11/02
cd /home/user/obs1/2024/11/02

# 配置文件：
[data]
rawpath = 20241102_hrs
[reduce]
out_path = output

# 实际使用的路径：
# 输入：/home/user/obs1/2024/11/02/20241102_hrs/
# 输出：/home/user/obs1/2024/11/02/output/
```

### 中间结果保存选项

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `save_overscan` | 保存过扫描校正结果 | yes |
| `save_bias` | 保存偏置减除结果 | yes |
| `save_flat` | 保存平场改正结果 | yes |
| `save_background` | 保存背景扣除结果 | yes |
| `save_cosmic` | 保存宇宙线去除结果 | yes |
| `save_extraction` | 保存一维谱提取结果 | yes |
| `save_wlcalib` | 保存波长定标结果 | yes |
| `save_deblaze` | 保存 Blaze 函数改正结果 | yes |

如果某一步设置为 `no`，则：
- ❌ 不会在 `output/` 下创建对应的子目录
- ❌ 不会保存该步的中间结果和图像
- ✅ 在 GUI 中应该有对应的取消勾选选项

## 使用建议

### 推荐配置（兴隆216 HRS）

```ini
[data]
rawpath = 20241102_hrs

[reduce]
out_path = output

[reduce.save_intermediate]
save_overscan = yes
save_bias = yes
save_flat = yes
save_background = yes
save_cosmic = yes
save_extraction = yes
save_wlcalib = yes
save_deblaze = yes

[telescope]
name = xinglong216hrs
instrument = hrs

[telescope.linelist]
linelist_type = ThAr
linelist_file = thar-noao.dat
use_precomputed_calibration = yes
calibration_file = wlcalib_20211123011_A.fits
```

### 节省空间配置

如果不需要保存所有中间结果：

```ini
[reduce.save_intermediate]
save_overscan = no     # 不保存过扫描结果
save_bias = no           # 不保存偏置结果
save_flat = no            # 不保存平场结果
save_background = no      # 不保存背景结果
save_cosmic = no         # 不保存宇宙线结果
save_extraction = no      # 不保存提取结果
save_wlcalib = no         # 不保存波长定标结果
save_deblaze = no         # 不保存 Blaze 改正结果

# 这样只会保存最终的光谱文件到 output/spectra/
```

## GUI 配置

### 配置文件对应

GUI 中应该有对应的设置界面：

1. **路径设置**
   - 原始数据目录
   - 输出目录

2. **保存选项**
   - 每一步的"保存中间结果"复选框
   - 对应配置文件中的 `save_*` 参数

3. **校准设置**
   - 望远镜/光谱仪选择
   - 灯谱线选择
   - 是否使用预计算校准

### 配置优先级

1. 用户配置文件（`specproc.cfg`）- 最高优先级
2. 默认配置文件（`default_config.cfg`）- 较低优先级

## 常见问题

### Q: 如何使用绝对路径？

A: 在配置文件中直接使用绝对路径：
```ini
[data]
rawpath = /absolute/path/to/rawdata
[reduce]
out_path = /absolute/path/to/output
```

### Q: 可以混合使用相对路径和绝对路径吗？

A: 可以！例如：
```ini
[data]
rawpath = 20241102_hrs              # 相对路径

[reduce]
out_path = /absolute/path/to/output  # 绝对路径
```

### Q: 配置文件放在哪里？

A: 配置文件放在你的**工作目录**中，而不是 SpecProc 源代码目录中：
```bash
# ❌ 错误
/path/to/SpecProc/specproc.cfg

# ✅ 正确
/home/user/obs1/specproc.cfg
```

### Q: 如何控制某一步不保存？

A: 在配置文件中设置对应的 `save_*` 参数为 `no`：
```ini
[reduce.save_intermediate]
save_bias = no  # 不保存偏置减除结果
```

或在 GUI 中取消勾选对应的"保存中间结果"选项。
