# SpecProc: 阶梯光栅光谱 FITS 数据处理图形界面管线

![SpecProc GUI](docs/SpecProc.png)

一个功能完整、基于 PyQt 的图形化交互式管线，专门用于处理阶梯光栅光谱仪（Echelle Spectrograph）的 FITS 数据。

## 目录

- [功能特性](#功能特性)
- [安装指南](#安装指南)
- [快速开始](#快速开始)
- [配置说明](#配置说明)
- [数据处理流程](#数据处理流程)
- [使用方法](#使用方法)
- [校准数据](#校准数据)
- [故障排除](#故障排除)
- [文档指南](#文档指南)
- [项目结构](#项目结构)
- [开源协议](#开源协议)
- [参与贡献](#参与贡献)
- [致谢](#致谢)

## 功能特性

### 完整的处理管线

**8 步全自动化光谱处理流程**：

1. **基础预处理 (Basic Pre-processing)** - 包含过扫区(Overscan)扣除、本底(Bias)扣除以及宇宙线去除。
2. **级次寻迹 (Orders Tracing)** - 生成主平场并追踪阶梯光栅的各个衍射级次。
3. **散射光扣除 (Scattered Light Subtraction)** - 使用 2D 卷积或样条函数建立并扣除级次间的杂散光背景。
4. **二维平场校正 (2D Flat-Field Correction)** - 高精度逐像素二维平场校正。
5. **一维光谱抽取 (1D Spectrum Extraction)** - 采用简单求和(Sum)或最优提取(Optimal)算法提取 1D 谱。
6. **去闪耀 (De-blazing)** - 闪耀函数改正，拉平光谱能量分布。
7. **波长定标 (Wavelength Calibration)** - 基于特征峰的模式盲配(Pattern Matching)与二维多项式曲面拟合定标

### 交互式图形界面

基于 PyQt5 的用户界面提供：
- 便捷的 Bias、Flat 和 Science FITS 文件管理
- 实时处理进度跟踪
- 详尽的处理日志与报错诊断
- 支持一键全自动执行或单步选中执行

## 安装指南

### 推荐方式：使用 Conda (本地源码安装)

```bash
# 创建全新的 conda 环境
conda create -n specproc python=3.8
conda activate specproc

# 安装科学计算与 GUI 依赖
conda install -c conda-forge pyqt5 numpy scipy astropy matplotlib -y

# 进入源码目录并以开发模式安装
cd /path/to/SpecProc
pip install -e .

# 启动程序
specproc
```

## 快速开始

### 工作目录设置

**极其重要**：`specproc` 必须在您的观测数据工作目录下运行，**绝不能**在 SpecProc 源码目录下运行！

### 标准工作流

```bash
# 1. 创建您的数据处理工作目录 (例如针对某次观测)
mkdir -p /myworkspace
cd /myworkspace

# 2. 创建一个目录用于存放 FITS 数据
mkdir -p 20241102_hrs output

# 3. 将您的 FITS 文件复制到该目录中
cp /somewhere/bias_*.fits ./20241102_hrs/
cp /somewhere/flat_*.fits ./20241102_hrs/
cp /somewhere/thar_*.fits ./20241102_hrs/
cp /somewhere/science_*.fits ./20241102_hrs/

# 4. 拷贝默认配置文件到当前目录并重命名 (可选)
cp /path/to/SpecProc/default_config.cfg ./specproc.cfg

# 5. 在您的工作目录启动管线
specproc --config ./specproc.cfg
```

### 生成的目录结构

一旦运行成功，您的工作目录将呈现如下结构：
```text
/myworkspace/                  # 您的当前工作目录
├── 20241102_hrs/              # 输入的 FITS 文件
├── output/                    # 处理结果 (自动生成)
│   ├── step1_basic/            # 步骤 1: 基础预处理 (脱偏置、宇宙线)
│   ├── step2_flat/             # 步骤 2: 主平场与寻迹边界、诊断图
│   ├── step3_scatterlight/     # 步骤 3: 杂散光背景模型
│   ├── step4_flat_corrected/   # 步骤 4: 2D 平场校正后的图像、Blaze模型
│   ├── step5_extraction/       # 步骤 5: 抽取出的 1D 光谱 (未波长定标)
│   ├── step6_deblazing/        # 步骤 6: 去闪耀后的 1D 光谱
│   ├── step7_wavelength/       # 步骤 7: 波长定标后的 1D 光谱及诊断PDF
│   └── step8_stitching/        # 步骤 8: 最终拼接完成的连续 1D 光谱
└── specproc.cfg              # 您的自定义配置文件
```

## 配置说明

您可以通过修改工作目录下的 `specproc.cfg` 来调整管线行为。

### 常见参数解析

```ini
[reduce]
# 总体输出目录
output_path = output

[data]
# 是否开启过扫区扣除 (设为 -1 禁用)
overscan_start_column = -1

# 探测器上下两半的分割行数 (Xinglong 2.16m 为 2068)
detector_split_row = 2068

[reduce.flat]
# Blaze 函数的 B-Spline 平滑度因子 (越大越平滑刚硬，越小越贴合局部特征)
blaze_smooth_factor = 1.0
# 横截面 Profile 沿 X 轴方向的切片数
n_profile_segments = 100
# 自动识别并使用刚性多项式对抗干涉条纹的红端级次数
fringe_orders = 20

[reduce.background]
# 杂散光拟合方法: convolution 或 column_spline
method = convolution
# 高斯核大小 (控制横向杂散光波纹的贴合度)
kernel_sigma_x = 13.0
kernel_sigma_y = 13.0

[reduce.wlcalib]
# 盲配的 RMS 误差容忍度 (埃)
rms_threshold = 0.5
```

## 数据处理流程

#### STEP 1: 基础预处理 (Basic Pre-processing)
- **输入**: 原始 FITS 文件 (bias, flat, ThAr, science)
- **处理**:
  - 提取过扫区 (读出偏置区域)
  - 计算中值或多项式拟合
  - 从图像中扣除过扫区偏置
  - 合并多个 bias 图像 (均值/中值)
  - 生成 master bias
  - 从 science/flat/ThAr 图像中扣除 master bias
  - 使用 L.A.Cosmic 识别并去除宇宙线 (仅限科学图像)
- **输出**: 预处理后的图像 (已做过扫区、本底、宇宙线校正)
- **说明**: 基础校正，为寻迹和提取准备数据。

#### STEP 2: 级次寻迹 (Orders Tracing)
- **输入**: 预处理后的 flat 图像
- **处理**:
  - 合并 flat 图像
  - 生成 master flat
  - 识别阶梯光栅衍射级次
  - 为每个级次拟合多项式迹线
  - 提取闪耀分布 (blaze profiles)
- **输出**: Master flat, 孔径边界 (apertures) 和 闪耀曲线
- **说明**: 为后续步骤提供级次边界和分布特征。

#### STEP 3: 散射光扣除 (Scattered Light Subtraction)
- **输入**: 预处理后的 science 图像
- **处理**:
  - 使用 2D 卷积或样条函数评估背景散射光
  - 从科学图像中扣除背景模型
- **输出**: 扣除背景后的图像
- **说明**: 消除级次间的杂散光。

#### STEP 4: 二维平场校正 (2D Flat-Field Correction)
- **输入**: 扣除背景的 science 图像
- **处理**:
  - 生成 2D 像素平场校正映射 (pixel-to-pixel flat)
  - 将 2D 平场校正应用于科学图像
- **输出**: 已做 2D 平场校正的图像
- **说明**: 纠正像素级的灵敏度差异。

#### STEP 5: 一维光谱抽取 (1D Spectrum Extraction)
- **输入**: 2D 平场校正图像
- **处理**:
  - 抽取每个级次的 1D 光谱
  - 方法: 简单求和 (Sum) 或 最优提取 (Optimal, Horne 1986)
  - 计算提取误差
- **输出**: SpectraSet (像素坐标系)
- **说明**: 将二维弯曲迹线坍缩为一维像素光谱。

#### STEP 6: 去闪耀 (De-blazing)
- **输入**: 提取的 1D 光谱 (像素坐标系)
- **处理**:
  - 读取平场的闪耀函数
  - 匹配对应级次
  - 除以闪耀函数: F_corrected(λ) = F_observed(λ) / B(λ)
  - 归一化为单位连续谱
- **输出**: 去闪耀的光谱
- **说明**: 校正光栅闪耀函数的包络效应。

#### STEP 7: 波长定标 (Wavelength Calibration)
- **输入**: 去闪耀后的 1D 光谱 (像素坐标系)
- **处理**:
  - 步骤 1: 标定 ThAr 灯谱
    - 提取 1D 光谱
    - 证认发射线
    - 拟合 2D 波长多项式: λ(x,y) = Σ p_ij·x^i·y^j
  - 步骤 2: 应用于科学光谱
    - 将像素坐标转换为波长单位
- **输出**: 完成波长定标的 1D 光谱
- **说明**: 建立各级次的物理波长标尺。

#### STEP 8: 级次拼接 (Order Stitching)
- **输入**: 已定标的一维多级次光谱
- **处理**: 将多个重叠的级次在重合波段依据信噪比进行加权融合，输出一条便于后续进行吸收线分析的连续一维物理光谱。
- **输出**: 最终的连续 1D 光谱
- **说明**: 生成可用于科学分析的最终数据产品。

## 使用方法

### GUI 模式 (默认)

```bash
# 启动 GUI (默认模式)
specproc

# 或明确指定 GUI 模式
specproc --mode gui

# 使用自定义配置文件
specproc --config /path/to/config.cfg
```
**操作流程**:
1. 在左侧栏分别添加对应的 Bias、Flat、Calibration(ThAr) 和 Science 图像。
2. （可选）点击右上角 "Settings" 确认望远镜台址和参数。
3. 勾选需要执行的步骤，点击 **"Run All Steps"** 一键到底，或点击 "Run Selected Steps" 单步调试。

### 命令行批处理模式 (CLI)
```bash
specproc --mode cli --config ./specproc.cfg
```
无需图形界面，适合部署在远端服务器批量执行常规管线作业。