# SpecProc 数据处理流程图

## 完整处理流程

```mermaid
flowchart TD
    Start([开始]) --> Stage0
    
    subgraph Stage0 [STAGE 0: 过扫描校正]
        A0[读取原始 FITS 文件<br/>bias, flat, ThAr, science] --> A1[提取 overscan 区域]
        A1 --> A2[计算 overscan 中值/多项式拟合]
        A2 --> A3[从图像减去 overscan 偏置]
        A3 --> End0[输出: overscan 校正的图像]
    end
    
    End0 --> Stage1
    
    subgraph Stage1 [STAGE 1: 偏置减除]
        B0[读取 bias 图像] --> B1[合并多个 bias 帧<br/>mean/median combine]
        B1 --> B2[生成 master bias]
        B2 --> B3[从 science/flat/ThAr 减去 master bias]
        B3 --> End1[输出: 偏置校正的图像]
    end
    
    End1 --> Stage2
    
    subgraph Stage2 [STAGE 2: 平场改正与阶序追踪]
        C0[读取 flat 图像] --> C1[合并多个 flat 帧]
        C1 --> C2[归一化 master flat]
        C2 --> C3[检测 echelle 阶序<br/>cross-correlation]
        C3 --> C4[追踪每个阶序的位置<br/>polynomial fitting]
        C4 --> C5[提取灵敏度图和 Blaze 函数]
        C5 --> End2[输出: FlatField + ApertureSet]
    end
    
    End2 --> Stage3
    
    subgraph Stage3 [STAGE 3: 背景扣除]
        D0[读取 bias 校正的 science 图像] --> D1[估计级间背景]
        D1 --> D2[2D 多项式拟合或中值滤波]
        D2 --> D3[减去背景模型]
        D3 --> End3[输出: 背景校正的图像]
    end
    
    End3 --> Stage4
    
    subgraph Stage4 [STAGE 4: 宇宙线去除<br/>仅科学图像]
        E0[读取背景校正的 science 图像] --> E1{启用宇宙线校正?}
        E1 -- 是 --> E2[计算阈值: mean + σ×std]
        E2 --> E3[检测宇宙线像素]
        E3 --> E4[中值滤波替换]
        E4 --> End4a[输出: 宇宙线校正的图像]
        E1 -- 否 --> End4b[跳过校正]
    end
    
    End4a --> Stage5
    End4b --> Stage5
    
    subgraph Stage5 [STAGE 5: 一维谱提取]
        F0[读取 cosmic 校正的图像] --> F1{提取方法?}
        F1 -- Sum extraction --> F2a[简单孔径求和]
        F1 -- Optimal extraction --> F2b[最优提取<br/>Horne 1986]
        F2a --> F3[提取每个阶序的 1D 光谱]
        F2b --> F3
        F3 --> F4[计算提取误差]
        F4 --> End5[输出: SpectraSet<br/>(像素空间)]
    end
    
    End5 --> Stage6
    
    subgraph Stage6 [STAGE 6: 波长定标]
        G0[第一步: ThAr 灯谱定标] --> G1[提取 ThAr 1D 光谱]
        G1 --> G2[识别发射谱线]
        G2 --> G3[拟合 2D 波长多项式<br/>λ(x,y) = Σ p_ij·x^i·y^j]
        G3 --> G4[建立像素→波长映射]
        G4 --> G5[第二步: 应用到 science 光谱]
        G5 --> G6[转换像素坐标为波长单位]
        G6 --> End6[输出: 波长定标的 1D 光谱]
    end
    
    End6 --> Stage7
    
    subgraph Stage7 [STAGE 7: Blaze 函数改正]
        H0[读取波长定标的光谱] --> H1[读取 flat 光谱的 Blaze 函数]
        H1 --> H2{阶序匹配}
        H2 --> H3[匹配对应阶序的 B λ]
        H3 --> H4[除以 Blaze 函数<br/>F_corrected λ = F_observed λ / B λ]
        H4 --> H5[归一化到单位连续谱]
        H5 --> End7[输出: 最终校准光谱]
    end
    
    End7 --> Final
    
    Final([结束]) --> Output[输出文件<br/>output/spectra/*.fits<br/>output/midpath/<br/>output/figures/*.png]
    
    style Stage0 fill:#e1f5ff
    style Stage1 fill:#fff4e1
    style Stage2 fill:#e8f5e9
    style Stage3 fill:#fce4ec
    style Stage4 fill:#f3e5f5
    style Stage5 fill:#fff9c4
    style Stage6 fill:#ffccbc
    style Stage7 fill:#d1c4e9
```

## 处理步骤详细说明

### STAGE 0: 过扫描校正 (Overscan Correction)
- **输入**: 原始 FITS 文件 (bias, flat, ThAr, science)
- **处理**:
  - 提取 overscan 区域 (读出偏置区域)
  - 计算中值或多项式拟合
  - 从图像减去 overscan 偏置
- **输出**: Overscan 校正的图像
- **说明**: 必须是第一步,应用于所有图像类型

### STAGE 1: 偏置减除 (Bias Correction)
- **输入**: Overscan 校正的图像
- **处理**:
  - 合并多个 bias 帧 (mean/median)
  - 生成 master bias
  - 从 science/flat/ThAr 减去 master bias
- **输出**: 偏置校正的图像
- **说明**: Bias 是 0 秒曝光,不需要宇宙线去除

### STAGE 2: 平场改正与阶序追踪 (Flat Fielding & Order Tracing)
- **输入**: 偏置校正的 flat 图像
- **处理**:
  - 合并多个 flat 帧
  - 归一化 master flat
  - 检测 echelle 阶序位置
  - 追踪每个阶序 (polynomial fitting)
  - 提取灵敏度图和 Blaze 函数
- **输出**: FlatField, ApertureSet
- **说明**: 合并时自动去除宇宙线

### STAGE 3: 背景扣除 (Background Subtraction)
- **输入**: 偏置校正的 science 图像
- **处理**:
  - 估计级间背景
  - 2D 多项式拟合或中值滤波
  - 减去背景模型
- **输出**: 背景校正的图像
- **说明**: 在宇宙线去除之前进行

### STAGE 4: 宇宙线去除 (Cosmic Ray Correction)
- **输入**: 背景校正的 science 图像
- **处理**:
  - 计算阈值 (mean + σ×std)
  - 检测宇宙线像素
  - 中值滤波替换
- **输出**: 宇宙线校正的图像
- **说明**: 仅对 science 图像进行

### STAGE 5: 一维谱提取 (Spectrum Extraction)
- **输入**: 宇宙线校正的图像
- **处理**:
  - Sum extraction: 简单孔径求和
  - Optimal extraction: 最优提取 (Horne 1986)
  - 提取每个阶序的 1D 光谱
  - 计算提取误差
- **输出**: SpectraSet (像素空间)
- **说明**: 在波长定标之前进行

### STAGE 6: 波长定标 (Wavelength Calibration)
- **输入**: 提取的 1D 光谱 (像素空间)
- **处理**:
  - **第一步**: ThAr 灯谱定标
    - 提取 ThAr 1D 光谱
    - 识别发射谱线
    - 拟合 2D 波长多项式 λ(x,y) = Σ p_ij·x^i·y^j
    - 建立像素→波长映射
  - **第二步**: 应用到 science 光谱
    - 转换像素坐标为波长单位
- **输出**: 波长定标的 1D 光谱
- **说明**: 必须在光谱提取之后进行

### STAGE 7: Blaze 函数改正 (De-blazing)
- **输入**: 波长定标的光谱
- **处理**:
  - 匹配对应阶序的 Blaze 函数 B(λ)
  - 除以 Blaze 函数: F_corrected(λ) = F_observed(λ) / B(λ)
  - 归一化到单位连续谱
- **输出**: 最终校准光谱
- **说明**: 必须在波长定标之后进行,因为 Blaze 函数是波长的函数

## 数据流向

```
原始图像 (bias, flat, ThAr, science)
    ↓
[STAGE 0] Overscan 校正
    ↓
[STAGE 1] Bias 减除
    ↓
         ├→ Flat → [STAGE 2] 平场改正与阶序追踪 → FlatField + ApertureSet
         │
         └→ Science → [STAGE 3] 背景扣除
                          ↓
                     [STAGE 4] 宇宙线去除
                          ↓
                     [STAGE 5] 一维谱提取 (像素空间)
                          ↓
                     [STAGE 6] 波长定标 (波长空间)
                          ↓
                     [STAGE 7] Blaze 函数改正
                          ↓
                     最终校准光谱
```

## 关键依赖关系

| 阶段 | 依赖 | 输出 |
|------|------|------|
| STAGE 0 | 无 (第一步) | Overscan 校正的图像 |
| STAGE 1 | STAGE 0 | 偏置校正的图像 |
| STAGE 2 | STAGE 1 (flat) | FlatField, ApertureSet |
| STAGE 3 | STAGE 1 (science) | 背景校正的图像 |
| STAGE 4 | STAGE 3 | 宇宙线校正的图像 |
| STAGE 5 | STAGE 4, STAGE 2 | SpectraSet (像素空间) |
| STAGE 6 | STAGE 5 | 波长定标的光谱 |
| STAGE 7 | STAGE 6, STAGE 2 | 最终校准光谱 |

## 配置参数

```ini
[reduce]
cosmic_enabled = yes              # 宇宙线校正开关
cosmic_sigma = 5.0                # 宇宙线检测 sigma 阈值
cosmic_window = 7                 # 中值滤波窗口大小

[reduce.bias]
combine_method = median           # Bias 合并方法

[reduce.flat]
combine_method = median           # Flat 合并方法

[reduce.trace]
degree = 3                       # 阶序追踪多项式阶数

[reduce.extract]
method = sum                      # 提取方法: sum/optimal
lower_limit = -5.0                # 提取下限 (像素)
upper_limit = 5.0                 # 提取上限 (像素)

[reduce.wlcalib]
linelist = ThAr                   # 校准灯类型
xorder = 4                        # 色散方向多项式阶数
yorder = 4                        # 空间方向多项式阶数
rms_threshold = 0.1                # RMS 阈值 (Å)

[reduce.background]
method = 2d_poly                  # 背景估计方法
poly_order = 2                    # 多项式阶数
```

## 输出文件

```
output/
├── midpath/                      # 中间结果
│   ├── bias_master.fits          # Master bias
│   ├── flat_master.fits          # Master flat
│   ├── apertures.fits            # 阶序定义
│   ├── wlcalib.fits              # 波长定标系数
│   └── cosmic_corrected.fits     # 宇宙线校正图像
├── spectra/                      # 最终光谱
│   └── science_filename.fits     # 校准后的 1D 光谱
└── figures/                      # 诊断图
    ├── bias_master.png
    ├── flat_master.png
    ├── orders.png
    ├── wlcalib.png
    └── spectra.png
```

## 注意事项

1. **Overscan correction 必须是第一步**: 应用于所有图像类型
2. **Bias 不需要宇宙线校正**: 0 秒曝光,几乎没有 cosmic ray
3. **Flat 不需要宇宙线校正**: 合并时自动去除
4. **仅 Science 需要宇宙线校正**: 长曝光单帧,必须去除
5. **波长定标必须在提取之后**: 只能对 1D 光谱进行
6. **Blaze 改正必须在波长定标之后**: Blaze 函数是波长函数

## 参考资料

- gamse: https://github.com/wangleon/gamse
- Horne 1986: PASP 98, 609 (最优提取)
- ESO Pipelines: 标准处理流程
- IRAF echelle: 传统处理流程
