# SpecProc Quick Start Guide

## 快速开始

### 安装 SpecProc（5分钟）

```bash
# 1. 进入 SpecProc 目录
cd /Users/abiao/Documents/资料存档/软件工具/gitrepo/SpecProc

# 2. 运行安装脚本
./install.sh

# 3. 激活环境
conda activate specproc

# 4. 启动软件
python run.py
```

### 命令行快速执行

```bash
# 激活环境
conda activate specproc

# GUI 模式（默认，可在任意目录运行）
specproc

# CLI 模式（可在任意目录运行）
specproc --mode cli

# 指定配置文件
specproc --config /path/to/config.cfg

# 查看帮助
specproc --help
```

**注意**：`specproc` 命令可以在任意目录下使用，无需进入 SpecProc 目录。

## 基本使用流程

### GUI 模式

1. **启动软件**
   ```bash
   conda activate specproc
   specproc
   ```

   或显式指定 GUI 模式：
   ```bash
   specproc --mode gui
   ```

2. **选择文件**
   - 点击 "Select Bias Files" 选择 bias 文件
   - 点击 "Select Flat Files" 选择 flat 文件
   - 点击 "Select Calibration Files" 选择 ThAr 灯谱文件
   - 点击 "Select Science Files" 选择科学图像文件

3. **运行处理**
   - 点击 "Run Full Pipeline" 运行完整流程
   - 或逐个点击步骤按钮分步执行

4. **查看结果**
   - 结果保存在 `output/` 目录
   - 中间结果保存在 `output/midpath/`

### CLI 模式

```bash
# 运行 CLI（可在任意目录）
specproc --mode cli

# 或指定配置文件
specproc --mode cli --config /path/to/config.cfg

# 按提示选择文件：
# - 选择 Bias 文件（例如：1-3）
# - 选择 Flat 文件（例如：4-6）
# - 选择 Science 文件（例如：8-9）
# - 选择 Calibration 文件（例如：7）
# - 选择处理步骤（回车执行全部）
```

## 数据准备

### 工作目录设置

**重要**：在你的工作目录中运行 SpecProc，而不是在源代码目录中运行。

### 正确的使用流程

```bash
# 1. 创建工作目录（用于观测项目）
mkdir -p ~/projects/2024/obs1
cd ~/projects/2024/obs1

# 2. 创建数据处理所需的子目录
mkdir -p rawdata output/midpath output/spectra output/figures

# 3. 将 FITS 数据文件放到 rawdata 目录
cp /somewhere/bias_*.fits ./rawdata/
cp /somewhere/flat_*.fits ./rawdata/
cp /somewhere/thar_*.fits ./rawdata/
cp /somewhere/science_*.fits ./rawdata/

# 4. 创建用户配置文件（可选）
cp /path/to/SpecProc/default_config.cfg ./specproc.cfg
# 编辑配置：nano ./specproc.cfg

# 5. 在工作目录中运行 SpecProc
specproc --config ./specproc.cfg
```

### 目录结构

**工作目录**（你处理数据的地方）：
```
~/projects/2024/obs1/       # 你的工作目录
├── rawdata/                  # 原始 FITS 数据
│   ├── bias_*.fits
│   ├── flat_*.fits
│   ├── thar_*.fits      # ThAr 灯谱
│   └── science_*.fits   # 科学图像
├── output/                   # 输出结果（自动生成）
│   ├── midpath/             # 中间结果
│   ├── spectra/             # 一维光谱
│   └── figures/             # 诊断图
├── specproc.cfg              # 用户配置文件（可选）
└── ...
```

**注意**：
- ❌ 不要在 SpecProc 源代码目录（`/path/to/SpecProc`）中运行
- ✅ 在你的工作目录中运行 `specproc`
- ✅ `rawdata` 和 `output` 会在你的工作目录中自动创建

## 配置文件

### 使用默认配置

直接使用 `default_config.cfg`，无需修改。

### 自定义配置

```bash
# 复制默认配置
cp default_config.cfg specproc.cfg

# 编辑配置
nano specproc.cfg
```

### 常用配置项

```ini
[data]
rawpath = ./rawdata           # 原始数据路径

[telescope]
name = xinglong216hrs         # 望远镜名称
instrument = hrs               # 光谱仪名称

[telescope.linelist]
linelist_type = ThAr           # 灯类型
linelist_file = thar-noao.dat  # 灯谱线文件

use_precomputed_calibration = yes  # 使用预计算校准
calibration_file = wlcalib_20211123011_A.fits
```

## 处理步骤说明

### 8阶段处理流程

1. **Overscan correction** - 过扫描校正（所有图像）
2. **Bias subtraction** - 偏置减除
3. **Flat fielding & order tracing** - 平场改正与阶序追踪
4. **Background subtraction** - 背景扣除
5. **Cosmic ray correction** - 宇宙线去除（仅科学图像）
6. **1D spectrum extraction** - 一维谱提取
7. **Wavelength calibration** - 波长定标
8. **De-blazing** - Blaze 函数改正

### CLI 模式步骤选择

```bash
# 执行所有步骤（推荐）
[按回车]

# 只执行前3步
0-2

# 执行特定步骤
0,1,2,4,5
```

## 输出结果

### 输出文件

```
output/
├── midpath/
│   ├── bias_master.fits          # 主 bias
│   ├── flat_master.fits          # 主 flat
│   ├── wlcalib.fits              # 波长校准
│   └── apertures.fits            # 阶序追踪结果
├── spectra/
│   └── science_001_ods.fits      # 一维光谱
└── figures/
    └── *.png                     # 诊断图
```

### 查看结果

```bash
# 查看一维光谱
ls output/spectra/

# 使用 Python/Astropy 查看光谱
python -c "from astropy.io import fits; hdul = fits.open('output/spectra/science_001_ods.fits'); print(hdul.info())"
```

## 常见问题

### 1. ImportError: No module named 'PyQt5'

```bash
conda activate specproc
pip install PyQt5
```

### 2. 配置文件找不到

```bash
cd SpecProc
ls -la default_config.cfg
```

### 3. 原始数据路径错误

编辑 `specproc.cfg`:
```ini
[data]
rawpath = /absolute/path/to/rawdata
```

### 4. 环境未激活

```bash
conda activate specproc
```

## 获取帮助

- **完整安装指南**: `INSTALLATION_GUIDE_CN.md` (中文)
- **英文安装指南**: `INSTALLATION_GUIDE.md`
- **项目文档**: `README.md`
- **处理流程图**: `PIPELINE_FLOWCHART.md`
- **校准数据**: `calib_data/README_CN.md`

## 下一步

1. 准备你的 FITS 数据文件
2. 配置文件设置（如需要）
3. 运行处理流程
4. 检查输出结果

祝你使用愉快！
