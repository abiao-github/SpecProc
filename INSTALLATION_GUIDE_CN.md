# SpecProc 安装和使用指南

## 系统要求

- Python 3.8 或更高版本
- Anaconda 或 Miniconda（推荐）

## 安装步骤

### 方法一：使用 Anaconda 环境（推荐）

#### 1. 创建并激活 Anaconda 环境

```bash
# 创建新的 conda 环境（Python 3.8）
conda create -n specproc python=3.8

# 激活环境
conda activate specproc
```

#### 2. 安装依赖包

```bash
# 安装 PyQt5（如果系统是 macOS，可能需要先安装 Qt）
conda install -c conda-forge pyqt5

# 安装科学计算包
conda install numpy scipy astropy matplotlib

# 或者使用 pip 安装
pip install numpy scipy astropy matplotlib PyQt5 configparser
```

#### 3. 安装 SpecProc

```bash
# 进入 SpecProc 目录
cd /Users/abiao/Documents/资料存档/软件工具/gitrepo/SpecProc

# 开发模式安装（推荐，这样修改代码后无需重新安装）
pip install -e .
```

#### 4. 验证安装

```bash
# 检查命令是否可用
specproc --help

# 查看版本
specproc --version

# 测试在任意目录下运行
cd /tmp && specproc --help
```

如果命令不可用，请确保：
1. 已成功安装 SpecProc：`pip install -e .`
2. 在正确的 conda 环境中：`conda activate specproc`
3. conda 环境 bin 目录在 PATH 中

### 方法二：直接使用 pip 安装

如果不想使用 conda 环境：

```bash
# 进入 SpecProc 目录
cd /Users/abiao/Documents/资料存档/软件工具/gitrepo/SpecProc

# 安装依赖
pip install -r requirements.txt  # 如果有 requirements.txt
# 或者手动安装
pip install numpy scipy astropy matplotlib PyQt5 configparser

# 安装 SpecProc
pip install -e .
```

## 使用方式

安装完成后，可以在任意目录下使用 `specproc` 命令。

### 方式一：GUI 模式（图形界面）- 默认

```bash
# 启动 GUI（默认模式）
specproc

# 指定配置文件
specproc --config /path/to/config.cfg

# 显式指定 GUI 模式
specproc --mode gui
```

**GUI 功能**：
- 选择 bias、flat、波长校准（ThAr灯谱）、科学图像文件
- 实时进度显示
- 处理日志输出
- 一步执行或分步执行处理流程

### 方式二：CLI 模式（命令行）

```bash
# 启动 CLI 模式
specproc --mode cli

# 指定配置文件
specproc --mode cli --config /path/to/config.cfg
```

### 方式三：使用 run.py（开发调试）

```bash
# 在 SpecProc 目录下使用 run.py
cd /Users/abiao/Documents/资料存档/软件工具/gitrepo/SpecProc
python run.py              # GUI 模式
python run.py --mode cli   # CLI 模式
```

**注意**：推荐使用 `specproc` 命令，可以在任意目录下运行。`python run.py` 主要用于开发调试。

**CLI 交互流程**：
1. 选择原始数据路径（默认 `./rawdata`）
2. 选择 Bias 文件
3. 选择 Flat 文件
4. 选择 Science 文件
5. 选择 Wavelength Calibration 文件（ThAr 灯谱）
6. 选择处理步骤（0-5，或回车执行全部）

**处理步骤说明**：
- 0: Overscan correction（过扫描校正）
- 1: Bias correction（偏置减除）
- 2: Flat fielding & order tracing（平场改正与阶序追踪）
- 3: Wavelength calibration（波长定标）
- 4: Background subtraction（背景扣除）
- 5: Spectrum extraction（一维谱提取）

## 配置文件

### 默认配置

默认配置文件位于：`SpecProc/default_config.cfg`

### 用户配置

可以在项目根目录创建 `specproc.cfg` 来覆盖默认配置：

```bash
# 在 SpecProc 目录下创建配置文件
cp default_config.cfg specproc.cfg

# 编辑配置文件
nano specproc.cfg  # 或使用其他编辑器
```

### 重要配置项

```ini
[data]
# 原始数据路径
rawpath = ./rawdata

[telescope]
# 望远镜名称
name = xinglong216hrs

# 光谱仪名称
instrument = hrs

[telescope.linelist]
# 灯谱线类型
linelist_type = ThAr

# 灯谱线文件路径
linelist_path = calib_data/linelists/
linelist_file = thar-noao.dat

# 是否使用预计算的校准文件
use_precomputed_calibration = yes
calibration_path = calib_data/telescopes/xinglong216hrs/
calibration_file = wlcalib_20211123011_A.fits
```

## 数据准备

### 工作目录设置

**重要**：SpecProc 应该在你的工作目录中运行，而不是在源代码目录中运行。

### 正确的工作流程

```bash
# 1. 创建你的工作目录（例如用于某个观测项目）
mkdir -p /my/project/2024/obs1
cd /my/project/2024/obs1

# 2. 创建数据处理所需的子目录
mkdir -p rawdata output/midpath output/spectra output/figures

# 3. 将 FITS 数据文件放到 rawdata 目录
# （从其他位置复制或移动）
cp /somewhere/bias_*.fits ./rawdata/
cp /somewhere/flat_*.fits ./rawdata/
cp /somewhere/science_*.fits ./rawdata/

# 4. 创建用户配置文件（可选）
# 从 SpecProc 源码目录复制默认配置
cp /path/to/SpecProc/default_config.cfg ./specproc.cfg

# 5. 在工作目录中运行 SpecProc
specproc --config ./specproc.cfg
```

### 目录结构

**工作目录结构**（你处理数据的地方）：
```
/my/project/2024/obs1/          # 你的工作目录
├── rawdata/                      # 原始 FITS 数据
│   ├── bias_*.fits
│   ├── flat_*.fits
│   ├── thar_*.fits
│   └── science_*.fits
├── output/                       # 处理结果（自动生成）
│   ├── midpath/                 # 中间处理结果
│   ├── spectra/                 # 最终一维光谱
│   └── figures/                 # 诊断图像
├── specproc.cfg                 # 用户配置文件（可选）
└── ...
```

**注意**：
- ❌ **不要**在 SpecProc 源代码目录（`/path/to/SpecProc`）中运行 `specproc`
- ✅ **应该**在你的工作目录中运行 `specproc`
- ✅ `rawdata` 和 `output` 会创建在你的工作目录中

### 放置数据文件

将 FITS 数据文件按类型分类到 `rawdata` 目录：

```
rawdata/
├── bias_001.fits
├── bias_002.fits
├── bias_003.fits
├── flat_001.fits
├── flat_002.fits
├── flat_003.fits
├── thar_001.fits        # ThAr 灯谱文件
├── science_001.fits     # 科学图像
├── science_002.fits
└── ...
```

## 运行示例

### 示例1：GUI 模式完整流程

```bash
# 1. 激活 conda 环境
conda activate specproc

# 2. 进入 SpecProc 目录
cd /Users/abiao/Documents/资料存档/软件工具/gitrepo/SpecProc

# 3. 启动 GUI
python run.py

# 4. 在 GUI 中：
#    - 点击 "Select Bias Files" 选择 bias 文件
#    - 点击 "Select Flat Files" 选择 flat 文件
#    - 点击 "Select Calibration Files" 选择 ThAr 灯谱文件
#    - 点击 "Select Science Files" 选择科学图像文件
#    - 点击 "Run Full Pipeline" 运行完整处理流程
```

### 示例2：CLI 模式完整流程

```bash
# 1. 激活 conda 环境
conda activate specproc

# 2. 进入 SpecProc 目录
cd /Users/abiao/Documents/资料存档/软件工具/gitrepo/SpecProc

# 3. 运行 CLI
python run.py --mode cli

# 4. 按提示操作：
#    - 输入 rawdata 路径（或回车使用默认）
#    - 选择 bias 文件（例如：1-3）
#    - 选择 flat 文件（例如：4-6）
#    - 选择 science 文件（例如：8-9）
#    - 选择 calibration 文件（例如：7）
#    - 选择处理步骤（回车执行全部）
```

## 输出结果

处理完成后，结果保存在配置文件中 `out_path` 指定的目录（默认为 `output/`）：

### 输出目录结构

```
output/
├── overscan_corrected/      # 第0步：过扫描校正结果
│   ├── 20241102001.fits
│   ├── 20241102002.fits
│   └── ...
├── bias_corrected/          # 第1步：偏置减除结果
│   ├── bias_master.fits
│   └── ...
├── flat_corrected/           # 第2步：平场改正结果
│   ├── flat_master.fits
│   ├── apertures.fits
│   └── ...
├── background_corrected/     # 第4步：背景扣除结果
│   └── ...
├── cosmic_corrected/        # 第5步：宇宙线去除结果（仅科学图像）
│   └── ...
├── spectra/                 # 最终一维光谱
│   ├── science_001_ods.fits
│   ├── science_002_ods.fits
│   └── ...
└── figures/                 # 诊断图
    ├── *.png
    ├── *.pdf
    └── ...
```

### 控制每一步是否保存

在配置文件中可以控制每一步是否保存中间结果：

```ini
[reduce.save_intermediate]
# 每一步是否保存中间结果（yes 或 no）
save_overscan = yes        # 第0步：过扫描校正
save_bias = yes             # 第1步：偏置减除
save_flat = yes              # 第2步：平场改正
save_background = yes        # 第4步：背景扣除
save_cosmic = yes           # 第5步：宇宙线去除
save_extraction = yes       # 第6步：一维谱提取
save_wlcalib = yes          # 第7步：波长定标
save_deblaze = yes          # 第8步：Blaze 函数改正
```

**注意**：
- 如果某一步设置为 `no`，则不会在 `output/` 下创建对应的子目录
- GUI 中也应该有对应的勾选选项
- 默认所有步骤都保存

## 常见问题

### 1. ImportError: No module named 'PyQt5'

```bash
# 安装 PyQt5
conda install -c conda-forge pyqt5
# 或
pip install PyQt5
```

### 2. 权限错误（macOS）

如果遇到 Qt 相关的权限错误：

```bash
# 重新安装 PyQt5
pip uninstall PyQt5
pip install PyQt5
```

### 3. 配置文件找不到

确保配置文件在 SpecProc 根目录下：

```bash
cd /Users/abiao/Documents/资料存档/软件工具/gitrepo/SpecProc
ls -la default_config.cfg specproc.cfg
```

### 4. 原始数据路径错误

在配置文件中设置正确的路径：

```ini
[data]
rawpath = /path/to/your/rawdata
```

## 卸载

```bash
# 卸载 SpecProc
pip uninstall SpecProc

# 或删除 conda 环境
conda deactivate
conda remove -n specproc --all
```

## 下一步

- 查看 `README.md` 了解更多功能
- 查看 `PIPELINE_FLOWCHART.md` 了解处理流程
- 查看 `calib_data/README_CN.md` 了解校准数据配置

## 技术支持

如有问题，请检查：
1. Python 版本是否 >= 3.8
2. 所有依赖包是否正确安装
3. 配置文件路径是否正确
4. 原始数据文件是否存在
