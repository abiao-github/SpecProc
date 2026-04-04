# SpecProc: 光谱数据处理图形界面工具

一个完整的基于 PyQt 的图形界面工具，用于处理阶梯光谱仪的 FITS 数据。受 [gamse](https://github.com/wangleon/gamse) 包启发。

## 目录

- [功能特性](#功能特性)
- [安装](#安装)
- [快速开始](#快速开始)
- [配置](#配置)
- [处理流程](#处理流程)
- [使用方法](#使用方法)
- [校准数据](#校准数据)
- [常见问题](#常见问题)
- [文档](#文档)

## 功能特性

### 完整的处理流程

**8 阶段自动化光谱处理**：

1. **过扫描校正** - 过扫描改正（所有图像）
2. **偏置减除** - 本底减除（使用均值/中值合并）
3. **平场改正与阶序追踪** - 平场改正与阶序追踪
4. **背景扣除** - 背景扣除
5. **宇宙线去除** - 宇宙线去除（仅科学图像）
6. **一维谱提取** - 一维谱提取
7. **波长定标** - 波长定标（对提取的 1D 光谱应用）
8. **Blaze 函数改正** - Blaze 函数改正

### 图形界面

基于 PyQt5 的用户界面，支持：
- bias、flat、science 图像文件管理
- 实时进度跟踪
- 处理日志和诊断
- 一键执行完整流程或分步执行

### 主要特性

- **配置驱动**：使用 INI 格式配置文件，易于参数调整
- **通用光谱仪支持**：可配置不同的阶梯光谱仪
- **灵活输出**：控制是否保存每一步的中间结果
- **命令行界面**：除 GUI 外，还支持 CLI 批处理

## 安装

选择以下任一安装方法：

### 方法 1：使用 pip

从 PyPI 直接安装 SpecProc。

```bash
# 安装 SpecProc
pip install specproc

# 启动应用程序
specproc
```

**说明：**
- 需要 Python 3.7+
- 依赖包会自动从 PyPI 安装
- 卸载使用 `pip uninstall specproc`

### 方法 2：使用 conda

Conda 提供完整的环境和所有依赖。有两种安装选项：

#### 选项 2.1：在已有的 conda 环境中安装

```bash
# 激活你的 conda 环境
conda activate your_environment

# 从 conda-forge 安装 SpecProc
conda install -c conda-forge specproc

# 启动应用程序
specproc
```

#### 选项 2.2：创建新的 conda 环境

```bash
# 为 SpecProc 创建新的 conda 环境
conda create -n specproc python=3.8
conda activate specproc

# 安装 SpecProc 和所有依赖
conda install -c conda-forge specproc

# 启动应用程序
specproc
```

**说明：**
- 推荐给希望使用独立环境的用户
- 所有依赖由 conda 管理
- 支持 Python 3.7-3.11

### 方法 3：从源代码安装

运行安装脚本从本地源代码安装 SpecProc 和所有依赖。

```bash
# 进入 SpecProc 目录
cd /path/to/SpecProc

# 添加可执行权限（如需要）
chmod +x install.sh

# 运行安装脚本
./install.sh

# 启动应用程序
specproc
```

**说明：**
- 安装脚本自动处理所有依赖
- 自动检测可用的包管理器（pip 或 conda）
- 将 SpecProc 安装到你的系统
- 适用于从本地源代码进行自动化安装

## 快速开始

### 工作目录设置

**重要**：SpecProc 应该在你的工作目录中运行，而不是在源代码目录中运行。

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
├── output/                    # 处理结果（自动生成）
│   ├── overscan_corrected/      # 第0步：过扫描校正结果
│   ├── bias_corrected/          # 第1步：偏置减除结果
│   ├── flat_corrected/           # 第2步：平场改正结果
│   ├── background_corrected/     # 第4步：背景扣除结果
│   ├── cosmic_corrected/        # 第5步：宇宙线去除结果（仅科学图像）
│   ├── spectra/                 # 最终一维光谱
│   └── figures/                 # 诊断图像
├── specproc.cfg              # 用户配置文件（可选）
└── ...
```

**注意**：
- ❌ 不要在 SpecProc 源代码目录（`/path/to/SpecProc`）中运行
- ✅ 在你的工作目录中运行 `specproc`
- ✅ `rawdata` 和 `output` 会在你的工作目录中自动创建

## 配置

### 配置文件类型

#### 默认配置文件

**位置**：`SpecProc/default_config.cfg`
**用途**：提供默认参数值
**修改**：不建议直接修改

#### 用户配置文件

**位置**：工作目录中的 `specproc.cfg`
**用途**：覆盖默认配置，自定义参数
**优先级**：用户配置 > 默认配置

### 路径配置

#### 数据路径

```ini
[data]
# 原始数据目录（相对于当前工作目录）
# 例如：在 /home/user/obs1 目录运行 specproc
# 如果 rawpath=20241102_hrs，则数据从 /home/user/obs1/20241102_hrs/ 加载
rawpath = 20241102_hrs
```

#### 输出路径

```ini
[reduce]
# 输出目录（相对于当前工作目录）
# 例如：在 /home/user/obs1 目录运行 specproc
# 如果 out_path=output，则结果保存到 /home/user/obs1/output/
#
# 输出目录结构：
# output/
#   ├── overscan_corrected/      # 第0步结果
#   ├── bias_corrected/          # 第1步结果
#   ├── flat_corrected/           # 第2步结果
#   ├── background_corrected/     # 第4步结果
#   ├── cosmic_corrected/        # 第5步结果（仅科学图像）
#   ├── spectra/                 # 最终一维光谱
#   └── figures/                 # 诊断图像
out_path = output
```

#### 路径相对性说明

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

### 中间结果保存

控制每一步是否保存中间结果：

```ini
[reduce.save_intermediate]
# 每一步是否保存中间结果
# 可以单独控制每一步是否保存
# 默认所有步骤都保存
save_overscan = yes        # 第0步：过扫描校正
save_bias = yes             # 第1步：偏置减除
save_flat = yes              # 第2步：平场改正
save_background = yes        # 第4步：背景扣除
save_cosmic = yes           # 第5步：宇宙线去除
save_extraction = yes       # 第6步：一维谱提取
save_wlcalib = yes          # 第7步：波长定标
save_deblaze = yes          # 第8步：Blaze 函数改正
```

**效果**：
- 如果某一步设置为 `no`，则不会在 `output/` 下创建对应的子目录
- 在 GUI 中应该有对应的取消勾选选项
- 默认所有步骤都保存

### 望远镜和校准配置

```ini
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

## 处理流程

### 完整的 8 阶段处理流程

```mermaid
flowchart TD
    Start([开始]) --> Stage0
    subgraph Stage0 [第0步: 过扫描校正]
        A0[读取原始 FITS 文件] --> A1[提取 overscan 区域]
        A1 --> A2[计算中值或多项式拟合]
        A2 --> A3[从图像减去 overscan 偏置]
        A3 --> End0[输出: overscan 校正的图像]
    end
    End0 --> Stage1

    subgraph Stage1 [第1步: 偏置减除]
        B0[读取 overscan 校正的图像] --> B1[合并 bias 帧]
        B1 --> B2[生成 master bias]
        B2 --> B3[从图像减去 master bias]
        B3 --> End1[输出: 偏置校正的图像]
    end
    End1 --> Stage2

    subgraph Stage2 [第2步: 平场改正与阶序追踪]
        C0[读取偏置校正的平场图像] --> C1[合并平场图像]
        C1 --> C2[生成 master flat]
        C2 --> C3[检测阶梯阶序]
        C3 --> C4[为每个阶序拟合多项式轨迹]
        C4 --> C5[提取 blaze 函数]
        C5 --> End2[输出: master flat 和阶序]
    end
    End2 --> Stage3

    subgraph Stage3 [第3步: 背景扣除]
        D0[读取偏置校正的图像] --> D1[估计背景]
        D1 --> D2{方法?}
        D2 -- 2D 多项式 --> D3a[拟合 2D 多项式]
        D2 -- 中值滤波 --> D3b[应用中值滤波]
        D3a --> D4[减去背景]
        D3b --> D4
        D4 --> End3a[输出: 背景扣除的图像]
        D4 --> End3b[继续到第4步]
    end
    End3a --> Stage4a
    End3b --> Stage4b

    subgraph Stage4 [第4步: 宇宙线去除]
        E0[读取背景扣除的图像] --> E1[检测宇宙线]
        E1 --> E2[均值 + σ×std 阈值]
        E2 --> E3[识别宇宙线像素]
        E3 --> E4[用中值滤波插值]
        E4 --> End4[输出: 宇宙线校正的图像]
    end
    End4 --> Stage5

    subgraph Stage5 [第5步: 一维谱提取]
        F0[读取宇宙线校正的图像] --> F1{提取方法?}
        F1 -- 求和提取 --> F2a[简单孔径求和]
        F1 -- 最优提取 --> F2b[最优提取 Horne 1986]
        F2a --> F3[提取每个阶序的 1D 光谱]
        F2b --> F3
        F3 --> F4[计算提取误差]
        F4 --> End5[输出: SpectraSet (像素空间)]
    end
    End5 --> Stage6

    subgraph Stage6 [第6步: 波长定标]
        G0[第一步: ThAr 灯谱定标] --> G1[提取 ThAr 1D 光谱]
        G1 --> G2[识别发射谱线]
        G2 --> G3[拟合 2D 波长多项式 λ x,y = Σ p_ij·x^i·y^j]
        G3 --> G4[建立像素→波长映射]
        G4 --> G5[第二步: 应用到 science 光谱]
        G5 --> G6[转换像素坐标为波长单位]
        G6 --> End6[输出: 波长定标的 1D 光谱]
    end
    End6 --> Stage7

    subgraph Stage7 [第7步: Blaze 函数改正]
        H0[读取波长定标的光谱] --> H1[读取 flat 光谱的 Blaze 函数]
        H1 --> H2{阶序匹配}
        H2 --> H3[匹配对应阶序的 B λ]
        H3 --> H4[除以 Blaze 函数 F_corrected λ = F_observed λ / B λ]
        H4 --> H5[归一化到单位连续谱]
        H5 --> End7[输出: 最终校准光谱]
    end
    End7 --> Final

    Final([结束]) --> Output[输出文件 output/spectra/*.fits output/midpath/ output/figures/*.png]

    style Stage0 fill:#e1f5ff
    style Stage1 fill:#fff4e1
    style Stage2 fill:#e8f5e9
    style Stage3 fill:#fce4ec
    style Stage4 fill:#f3e5f5
    style Stage5 fill:#fff9c4
    style Stage6 fill:#ffccbc
    style Stage7 fill:#d1c4e9
```

### 阶段说明

#### 第0步：过扫描校正
- **输入**：原始 FITS 文件（bias, flat, ThAr, science）
- **处理**：
  - 提取 overscan 区域（读出偏置区域）
  - 计算中值或多项式拟合
  - 从图像减去 overscan 偏置
- **输出**：Overscan 校正的图像
- **说明**：必须是第一步，应用于所有图像类型

#### 第1步：偏置减除
- **输入**：Overscan 校正的图像
- **处理**：
  - 合并多个 bias 帧（均值/中值）
  - 生成 master bias
  - 从 science/flat/ThAr 减去 master bias
- **输出**：偏置校正的图像
- **说明**：Bias 是 0 秒曝光，不需要宇宙线去除

#### 第2步：平场改正与阶序追踪
- **输入**：偏置校正的平场图像
- **处理**：
  - 合并平场图像
  - 生成 master flat
  - 检测阶梯阶序
  - 为每个阶序拟合多项式轨迹
  - 提取 blaze 函数
- **输出**：Master flat、阶序和 blaze 函数
- **说明**：为后续步骤提供阶序和 blaze 函数

#### 第3步：背景扣除
- **输入**：偏置校正的图像
- **处理**：
  - 使用 2D 多项式或中值滤波估计背景
  - 从图像减去背景
- **输出**：背景扣除的图像
- **说明**：应用于科学图像，在宇宙线去除后

#### 第4步：宇宙线去除
- **输入**：背景扣除的图像（仅科学图像）
- **处理**：
  - 使用 sigma 阈值检测宇宙线
  - 用中值滤波插值
- **输出**：宇宙线校正的图像
- **说明**：仅应用于科学图像（长曝光）

#### 第5步：一维谱提取
- **输入**：宇宙线校正的图像
- **处理**：
  - 为每个阶梯阶序提取 1D 光谱
  - 方法：求和提取或最优提取（Horne 1986）
  - 计算提取误差
- **输出**：像素空间的 1D 光谱
- **说明**：依赖于第2步的阶序

#### 第6步：波长定标
- **输入**：提取的 1D 光谱（像素空间）
- **处理**：
  - 第一步：定标 ThAr 灯谱
    - 提取 1D 光谱
    - 识别发射谱线
    - 拟合 2D 波长多项式 λ(x,y) = Σ p_ij·x^i·y^j
  - 第二步：应用到 science 光谱
    - 将像素坐标转换为波长单位
- **输出**：波长定标的 1D 光谱
- **说明**：必须在一维谱提取之后

#### 第7步：Blaze 函数改正
- **输入**：波长定标的 1D 光谱
- **处理**：
  - 从平场读取 blaze 函数（第2步）
  - 匹配阶序
  - 除以 blaze 函数：F_corrected(λ) = F_observed(λ) / B(λ)
  - 归一化到单位连续谱
- **输出**：最终校准的光谱
- **说明**：必须在波长定标之后

## 使用方法

### GUI 模式（默认）

```bash
# 启动 GUI（默认模式）
specproc

# 或显式指定 GUI 模式
specproc --mode gui

# 使用自定义配置文件
specproc --config /path/to/config.cfg
```

**GUI 工作流程**：
1. 选择 bias 文件
2. 选择 flat 文件
3. 选择校准文件（ThAr 灯谱）
4. 选择 science 文件
5. 点击"运行完整流程"或分步执行
6. 实时查看进度
7. 在 output 目录查看结果

### CLI 模式（命令行）

```bash
# 运行 CLI 模式
specproc --mode cli

# 使用自定义配置
specproc --mode cli --config /path/to/config.cfg
```

**CLI 工作流程**：
1. 按提示选择文件
2. 选择处理步骤（0-7，或回车执行全部）
3. 监控控制台进度
4. 在 output 目录查看结果

## 校准数据

### 目录结构

```
calib_data/
├── linelists/              # 灯谱发射线目录
│   ├── thar-noao.dat      # ThAr 灯谱线（兴隆216 HRS 推荐）
│   ├── thar.dat           # 标准 ThAr 灯谱线
│   └── FeAr.dat           # FeAr 灯谱线
└── telescopes/             # 望远镜特定校准文件
    ├── generic/           # 通用配置模板
    └── xinglong216hrs/    # 兴隆216望远镜
        ├── wlcalib_20141103049.fits
        ├── wlcalib_20171202012.fits
        ├── wlcalib_20190905028_A.fits
        └── wlcalib_20211123011_A.fits
```

### 灯谱线文件

**可用的灯谱线文件**：
- `thar-noao.dat` - ThAr 灯谱线（兴隆216 HRS 推荐）
- `thar.dat` - 标准 ThAr 灯谱线
- `FeAr.dat` - FeAr 灯谱线

**支持的灯类型**：
- `ThAr` - 氩钍灯（阶梯光谱仪最常用）
- `FeAr` - 铁氩灯
- `Ar` - 氩灯
- `Ne` - 氖灯
- `He` - 氦灯
- `Fe` - 铁灯

### 望远镜校准文件

**兴隆216 HRS 可用的校准文件**：
- `wlcalib_20141103049.fits` - 2014-11-03 04:50
- `wlcalib_20171202012.fits` - 2017-12-02 01:20
- `wlcalib_20190905028_A.fits` - 2019-09-05 02:50（版本 A）
- `wlcalib_20211123011_A.fits` - 2021-11-23 01:10（版本 A）- **最新**

### 使用选项

#### 使用预计算的校准文件（推荐）

```ini
use_precomputed_calibration = yes
calibration_file = wlcalib_20211123011_A.fits
```

#### 重新拟合波长定标

```ini
use_precomputed_calibration = no
linelist_file = thar-noao.dat
```

### 添加自定义校准数据

#### 添加新的灯谱线文件：

1. 在 `linelists/` 目录创建文件
2. 遵循文件格式（波长、强度、注释）
3. 在配置文件中配置：`linelist_file = <文件名>`

#### 添加新的望远镜：

1. 创建目录：`calib_data/telescopes/<望远镜名称>/`
2. 放置校准文件，使用正确的命名规则
3. 在配置文件中配置：
   ```ini
   [telescope]
   name = <望远镜名称>
   instrument = <光谱仪名称>
   ```

## 常见问题

### ImportError: No module named 'PyQt5'

```bash
conda activate specproc
pip install PyQt5
```

### specproc 命令找不到

```bash
conda activate specproc
pip install -e .
```

### 配置文件找不到

```bash
# 复制默认配置
cp /path/to/SpecProc/default_config.cfg ./specproc.cfg
```

### GitHub 大文件错误

**错误**：`File exceeds GitHub's file size limit of 100.00 MB`

**解决方案**：大 FITS 文件不应该提交。使用 `.gitignore` 排除它们。

**防止将来添加**：
- 在 `.gitignore` 中添加输出目录
- 在单独的工作目录运行 SpecProc，而不是在源代码目录

### 处理错误

1. **缺少 bias 文件**：Bias 校正是可选的，但推荐使用
2. **缺少 flat 文件**：阶序追踪需要
3. **缺少校准文件**：波长定标需要

## 文档

- 查看 [calib_data/README_CN.md](calib_data/README_CN.md) 了解校准数据配置
- 查看 [CONFIGURATION_GUIDE_CN.md](CONFIGURATION_GUIDE_CN.md) 了解详细配置指南
- 查看 [PIPELINE_FLOWCHART.md](PIPELINE_FLOWCHART.md) 了解详细处理流程

## 项目结构

```
SpecProc/
├── README.md                    # 主文档（英文）
├── README_CN.md                # 主文档（中文）
├── DOCUMENTATION.md             # 文档索引
├── INSTALLATION_GUIDE.md       # 安装指南（英文）
├── INSTALLATION_GUIDE_CN.md    # 安装指南（中文）
├── QUICK_START.md               # 快速开始（中文）
├── QUICK_REFERENCE.md           # 快速参考（英文）
├── CONFIGURATION_GUIDE_CN.md   # 配置指南（中文）
├── default_config.cfg           # 默认配置
├── specproc.cfg.example         # 用户配置示例
├── calib_data/
│   ├── README.md                # 校准数据指南（英文）
│   ├── README_CN.md            # 校准数据指南（中文）
│   ├── linelists/               # 灯谱线文件
│   └── telescopes/              # 望远镜校准文件
├── src/                         # 源代码
│   ├── gui/                     # GUI 模块
│   ├── core/                    # 核心处理
│   ├── config/                  # 配置管理
│   ├── utils/                   # 工具函数
│   └── plotting/                # 绘图功能
├── install.sh                   # 安装脚本
├── requirements.txt              # Python 依赖
├── setup.py                     # 安装配置
├── run.py                      # 主入口点
└── test_*.py                    # 测试文件
```

## 许可证

查看 LICENSE 文件了解详情。

## 贡献

欢迎贡献！请随时提交 Pull Request。

## 致谢

- 受 [gamse](https://github.com/wangleon/gamse) 包启发
- 使用 PyQt5、NumPy、SciPy 和 Astropy 构建

## 支持

如有问题和疑问，请在 GitHub 上提交 Issue。
