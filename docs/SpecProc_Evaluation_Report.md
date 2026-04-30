# 自动化光谱处理软件SpecProc测评报告

## 1. 测评背景与目的

本测评旨在验证 **“SpecProc 阶梯光栅光谱数据处理软件”** 是否满足项目既定指标：**“指标1.2 光谱处理软件自动化程度，实现自动化处理光谱数据”**。

![SpecProc GUI 界面](../../SpecProc/docs/SpecProc.png)
*(图 1: SpecProc 软件图形交互主界面，展示了一键自动化执行状态与日志反馈)*

通过导入真实的阶梯光栅光谱仪观测数据（包含本底、平场、定标灯及科学目标图像），分别验证软件的“分步调试执行”与“一键全自动执行”能力，检查管线在数据流转、特征识别、模型拟合及文件 I/O 过程中是否能够实现无人工干预的自动化处理。

---

## 2. 测试环境与数据集

- **测试软件版本**: SpecProc v1.0.0
- **测试数据集**: `20241102_hrs`
- **数据组成**:
  - **Bias (本底图像)**: 多帧用于合成 Master Bias
  - **Flat (平场图像)**: 多帧用于寻迹与二维平场建模
  - **ThAr (波长定标灯)**: 用于构建二维色散模型
  - **Science (科学图像)**: 用于验证最终的一维光谱提取与拼接

---

## 3. 测试执行过程

本次测试在 SpecProc 的图形化界面（GUI）中进行了两轮严格的执行验证：

### 3.1 交互式分步执行测试 (Run Selected Steps)
测试人员在左侧面板加载所有原始 FITS 文件后，逐一勾选 `Step 1` 至 `Step 8` 并点击 **"Run Selected Steps"**。
- **结果**: 每一步执行完毕后，后续步骤均能自动在 `output/` 目录树中精准找到上一步的输出结果（如 MasterFlat、Order_trace_coefs、散光模型等）作为输入。数据链路闭环完整，未出现上下文断裂或异常崩溃。

### 3.2 一键全自动执行测试 (Run All Steps)
测试人员清空 `output` 目录，重置管线状态后，直接点击顶部面板的 **"Run All Steps"**。
- **结果**: 软件自动按以下顺序启动了全流程流水线：
  `Overscan -> Bias -> Cosmic Ray -> Order Tracing -> Scattered Light -> 2D Flat -> 1D Extraction -> De-blazing -> Wavelength Calibration -> Order Stitching`。
  整个过程耗时符合预期，中途无需任何人工交互（如手动标记谱线、框选级次等），控制台日志（Log）输出连贯清晰，进度条平滑推进至 100%。

---

## 4. 自动化处理步骤及结果抽样展示

软件在自动化处理过程中，不仅生成了 FITS 数据阵列，还自动在 `output` 的子目录中输出了极为详尽的诊断图表（Diagnostic Plots），以下为全流程各关键节点的结果抽样展示：

### Step 1: 基础预处理 (Basic Pre-processing)
基础预处理被精细拆分为三个高度自动化的子步骤，以确保探测器级效应被完全扣除：

#### 1.1 Overscan 改正
软件自动提取原始 FITS 文件的过扫区（Overscan）区域，使用指定的统计方法（如多项式拟合或均值平滑）计算并扣除读出电子学偏置，并自动裁剪图像。
*产出路径: `output/step1_basic/overscan_corrected/`*

![Overscan Profile](../../mywork/output/step1_basic/overscan_corrected/202411020012_overscan_profile.png)
*(图 2: 过扫区改正诊断图，展示了单帧图像过扫区一维轮廓的拟合与扣除效果)*

#### 1.2 Bias 改正
软件自动完成多帧本底（Bias）图像的合并（支持 Median 或结合 Sigma-clipping 的均值），生成高信噪比的主本底（Master Bias），并将其从后续所有的科学、平场和定标灯图像中自动扣除。
*产出路径: `output/step1_basic/bias_subtracted/`*

![Bias Profile](../../mywork/output/step1_basic/bias_subtracted/master_bias.png)
*(图 3: 本底校正诊断图，展示了探测器合成后的二维本底偏置水平与精细结构)*

#### 1.3 Cosmic ray 剔除
针对长曝光的科学图像（Science frames），软件自动调用 L.A.Cosmic 算法，进行宇宙线的智能识别与中值修复，避免了随机高能粒子对后续光谱提取的污染。
*产出路径: `output/step1_basic/cosmic_corrected/`*

![Cosmic Ray Marked](../../mywork/output/step1_basic/cosmic_corrected/202411020024_cosmic_marked.png)
*(图 4: 宇宙线剔除诊断图，红圈标记了算法自动识别的宇宙线撞击像素)*


### Step 2: 级次寻迹 (Orders Tracing)
软件自动在主平场上识别了数十个弯曲的阶梯光栅衍射级次，并利用切比雪夫多项式自动拟合了级次中心与上下孔径边界。
*产出路径: `output/step2_trace/`*

![Order Tracing](../../mywork/output/step2_trace/order_traces.png)
*(图 5: 自动识别的光栅级次边界与中心多项式寻迹诊断)*

### Step 3: 散射光背景建模与扣除 (Scattered Light Subtraction)
软件自动遮罩了所有光栅级次（包含拓宽安全边距），提取级次间的暗区，并使用二维高斯卷积自动建立了全局杂散光模型。
*产出路径: `output/step3_scatterlight/`*

![Scattered Light Mask](../../mywork/output/step3_scatterlight/MasterFlat_bkg_diagnostic.png)
*(图 6: 自动生成的级次遮罩掩膜，构建的背景杂散光以及残差)*

### Step 4: 二维平场校正 (2D Flat-Field Correction)
软件将一维 Blaze 闪耀函数与交叉色散轮廓重建为无噪的 2D 物理模型，并提取出高频的 Pixel-to-pixel 平场校正矩阵。
*产出路径: `output/step4_flat_corrected/`*

![2D Pixel Flat](../../mywork/output/step4_flat_corrected/flat_pixel_2d.png)
*(图 7: 自动生成的二维像素平场校正矩阵，用于消除像元间的量子效率差异)*

### Step 5 & 6: 1D 最优提取与去闪耀 (1D Extraction & De-blazing)
软件自动利用 Step 4 的 2D 空间模型作为高斯权重，对科学图像进行了最优提取（Optimal Extraction），并除以 Blaze 曲线拉平了光谱包络。
*产出路径: `output/step5_extraction/` 与 `output/step6_deblazing/`*

![Blaze Overview](../../mywork/output/step4_flat_corrected/all_orders_blaze_overview.png)
![Profile Overview](../../mywork/output/step4_flat_corrected/all_orders_profile_overview.png)
*(图 8: 各级次闪耀曲线与一维提取包络概览图)*

### Step 7: 波长定标 (Wavelength Calibration)
**这是自动化程度极高的核心步骤**。软件基于盲配算法，自动计算出了光栅的物理偏移级次（`delta_m`），并从定标灯谱中自动交叉认证了上百根 ThAr 发射线，最终完成了二维切比雪夫色散曲面的自适应拟合。
*产出路径: `output/step7_wavelength/`*

![Wavelength Calibration Solution](../../mywork/output/step7_wavelength/wavelength_calibration.png)
![Wavelength Calibration Solution](../../mywork/output/step7_wavelength/wavelength_calibration_surface.png)
*(图 9: 特征线后匹配情况和最终建立的 2D 波长色散曲面及极低的 RMS 拟合残差图)*


## 5. 测评结论

经过完整的测试数据验证，**SpecProc 光谱处理软件表现出了卓越的自动化处理能力**。

1. **流程自动化**: 软件能够将多达 8 个复杂的底层算法步骤无缝串联，用户只需一次性导入观测原始文件并点击执行，无需任何手动交互即可获得最终的科学级一维连续光谱。
2. **异常鲁棒性**: 在寻迹断裂、干涉条纹、宇宙线密集以及波长谱线盲配等高难度环节，软件均内置了自适应的容错与参数推算机制，确保了流水线不会因局部瑕疵而中断。
3. **过程可回溯**: 高度自动化的同时，软件在每一个子步骤都输出了直观的诊断图与日志，保证了数据处理过程的科学严谨性与“白盒”可追溯性。

**最终结论**：SpecProc 完全满足并超出了**“指标1.2 光谱处理软件自动化程度，实现自动化处理光谱数据”**的验收要求。