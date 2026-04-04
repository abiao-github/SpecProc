# SpecProc Documentation

## Project Documentation

SpecProc maintains essential documentation in both English and Chinese.

## Core Documents

### User Guides (Essential)

1. **[README.md](README.md)** - Main Documentation
   - Project overview and features
   - 8-stage processing pipeline
   - Project structure
   - Installation and usage

2. **[INSTALLATION_GUIDE.md](INSTALLATION_GUIDE.md)** - Installation Guide (English)
   - Detailed installation steps
   - Dependencies
   - Configuration
   - Troubleshooting

3. **[INSTALLATION_GUIDE_CN.md](INSTALLATION_GUIDE_CN.md)** - 安装指南（中文）
   - 详细的安装步骤
   - 依赖包说明
   - 配置文件说明
   - 常见问题解决

4. **[QUICK_START.md](QUICK_START.md)** - Quick Start Guide (Chinese)
   - 5-minute quick installation
   - Basic usage workflow
   - Data preparation
   - Common issues

5. **[QUICK_REFERENCE.md](QUICK_REFERENCE.md)** - Quick Reference (English)
   - Command cheat sheet
   - Configuration reference
   - Common issues
   - File naming conventions

6. **[CONFIGURATION_GUIDE_CN.md](CONFIGURATION_GUIDE_CN.md)** - 配置文件说明（中文）
   - 配置文件类型和优先级
   - 完整的配置文件示例
   - 所有参数说明
   - 路径相对性说明
   - 中间结果保存选项
   - GUI 配置说明

### Technical Reference

6. **[PIPELINE_FLOWCHART.md](PIPELINE_FLOWCHART.md)** - Processing Flowchart
   - Complete 8-stage pipeline
   - Mermaid flowchart
   - Data flow

### Calibration Data

7. **[calib_data/README.md](calib_data/README.md)** - Calibration Data Guide (English)
   - Calibration data directory structure
   - Lamp line lists
   - Telescope calibration files
   - Configuration methods

8. **[calib_data/README_CN.md](calib_data/README_CN.md)** - 校准数据说明（中文）
   - 校准数据目录结构
   - 灯谱线文件说明
   - 望远镜校准文件说明
   - 配置方法

## 配置和脚本

### 配置文件

- **default_config.cfg** - 默认配置文件
- **specproc.cfg** - 用户配置文件（可选）

### 安装脚本

- **install.sh** - 自动安装脚本

### 依赖文件

- **requirements.txt** - Python 依赖列表
- **setup.py** - 安装配置

## 运行脚本

- **run.py** - 主运行脚本

## 测试文件

- **test_core.py** - 核心功能测试
- **test_gui_init.py** - GUI 初始化测试
- **test_gui_steps.py** - GUI 步骤测试
- **test_overscan.py** - 过扫描校正测试
- **test_stages.py** - 处理阶段测试

## 其他文件

- **LICENSE** - 许可证
- **.gitignore** - Git 忽略规则
- **SpecProc.code-workspace** - VS Code 工作区配置（在父目录）

## 文档使用指南

### 新用户

1. 先阅读 **README.md** 了解项目
2. 按照 **INSTALLATION_GUIDE_CN.md** 安装
3. 使用 **QUICK_START.md** 快速开始

### 高级用户

1. 查看 **PIPELINE_FLOWCHART.md** 了解处理流程
2. 参考 **calib_data/README_CN.md** 配置校准数据
3. 编辑 **default_config.cfg** 或 **specproc.cfg** 自定义参数

### 开发者

1. 阅读 **README.md** 中的项目结构
2. 查看 **setup.py** 了解安装配置
3. 运行测试文件验证功能

## 项目结构

```
SpecProc/
├── README.md                          # 主文档
├── INSTALLATION_GUIDE_CN.md          # 安装指南（中文）
├── QUICK_START.md                    # 快速开始（中文）
├── PIPELINE_FLOWCHART.md             # 处理流程图
├── default_config.cfg                # 默认配置
├── specproc.cfg                     # 用户配置
├── install.sh                      # 安装脚本
├── requirements.txt                 # 依赖列表
├── setup.py                       # 安装配置
├── run.py                         # 主运行脚本
├── LICENSE                        # 许可证
├── .gitignore                    # Git 忽略
├── calib_data/                     # 校准数据
│   └── README_CN.md                # 校准数据说明
├── src/                           # 源代码
│   ├── gui/                       # GUI 模块
│   ├── core/                      # 核心处理
│   ├── config/                    # 配置管理
│   ├── utils/                     # 工具函数
│   └── plotting/                  # 绘图功能
└── test_*.py                     # 测试文件
```

## 快速链接

### 安装和运行
```bash
# 安装
./install.sh

# 运行 GUI
specproc

# 运行 CLI
specproc --mode cli
```

### 文档阅读
- [README.md](README.md) - 从这里开始
- [INSTALLATION_GUIDE_CN.md](INSTALLATION_GUIDE_CN.md) - 详细安装
- [QUICK_START.md](QUICK_START.md) - 快速上手
- [PIPELINE_FLOWCHART.md](PIPELINE_FLOWCHART.md) - 了解流程
- [calib_data/README_CN.md](calib_data/README_CN.md) - 配置校准

## 注意

所有文档都是核心必要文档，不包含冗余或历史记录。如有需要，请查看 Git 历史记录。
