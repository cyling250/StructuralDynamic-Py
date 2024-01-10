# Python在结构动力计算中的应用（StructuralDynamic-Py）

## 1 说明

本项目为建筑工业出版社《Python在结构动力计算中的应用》书本的配套代码。

本书中的所有代码都省略了库导入部分，所有代码的库导入部分都仅在代码文件中写出，书中的所有出现的代码都在本项目中有对应文件，并且本项目还自行编写了一些能够帮助读者快速出图的绘图函数和一些额外小demo。

本书秉承书本重原理，代码重实践的思路，请各位读者将书本与代码联系起来阅读，而非一味看重书本或者一味只学代码。

<!--TODO 添加一部分的前言部分-->

## 2 文件架构

### 2.1 libs文件夹

[libs文件夹](libs)主要包含可调用的库文件，对应书中的各种算法。在执行[parts文件](parts)中的各demo时，
需要调用这些库。

libs文件夹中一共包含10个.py文件，这些文件的主要内容介绍如下（按照字母排序）：

[abaqus.py](libs/abaqus.py):主要存放第七章中Python在ABAQUS中的应用相关内容的函数。

[beam6.py](libs/beam6.py):主要存放等截面六自由度梁单元模型。

[damping.py](libs/damping.py):主要存放阻尼矩阵构建的相关程序。

[figtools.py](libs/figtools.py):主要存放一些自定义的绘图相关函数。

[layer_shear.py](libs/layer_shear.py):主要存放层剪切模型的刚度矩阵、质量矩阵构建函数。

[model_analysis.py](libs/modal_analysis.py):主要存放模态分析方法的相关程序，如振型叠加法，反应谱计算等。

[read_wave.py](libs/read_wave.py):主要存放地震波处理的相关程序。

[sbs_integration_linear](libs/sbs_integration_linear.py):主要存放线性状态下逐步积分法的相关程序。

[sbs_integration_nonlinear.py](libs/sbs_integration_nonlinear.py):主要存放非线性状态下逐步积分法的相关程序。

[vibration_mode.py](libs/vibration_mode.py):主要存放模态分析的相关程序。

### 2.2 parts文件夹

[parts文件夹](parts)
主要包含书中的各种例子，其下的子目录命名格式为part+章节数字，例如第三章的相关代码就存放在[parts/part3](parts/part3)
目录下。

parts文件夹下的文件命名规则为：按照书中出现的例子进行编号，例如在书中
<!--TODO 等书本例子顺序敲定后重新对文件进行命名-->

### 2.3 res文件夹

[res文件夹](res)
主要存放的是一些在编程中用到的资源文件和程序的输出文件。例如本书中经常使用到的地震波ELCENTRO就存放在[res/ELCENTRO.DAT](res/ELCENTRO.DAT)
文件下。

## 3 各函数的详细对应关系(按照书中章节排序)

### 3.1 第三章

<div class="center">

|  对应章节   |     描述     |                          对应文件或函数                          |
|:-------:|:----------:|:---------------------------------------------------------:|
| 3.1.1_1 |   层剪切模型    |           [layer_shear.py](libs/layer_shear.py)           |
| 3.1.1_1 |    杆系模型    |                 [beam6.py](libs/beam6.py)                 |
| 3.3.2_1 | Rayleigh阻尼 |       [damping.py:<br/>rayleigh()](libs/damping.py)       |
| 3.3.2_1 | Caughey阻尼  |       [damping.py:<br/>caughey()](libs/damping.py)        |
| 3.3.2_2 |   非比例阻尼    | [damping.py:<br/>damping_nonclassical()](libs/damping.py) |
|   3.4   |  地震波读取函数   | [read_wave.py:<br/>read_quake_wave()](libs/read_wave.py)  |

</div>

### 3.2 第四章

<div class="center">

| 对应章节  |       描述       |                                   对应文件或函数                                    |
|:-----:|:--------------:|:----------------------------------------------------------------------------:|
|  4.1  |   Rayleigh法    |       [vibration_mode.py:<br/>rayleigh_psi()](libs/vibration_mode.py)        |
|  4.2  | Rayleigh-Ritz法 |       [vibration_mode.py:<br/>rayleigh_ritz()](libs/vibration_mode.py)       |
|  4.3  |   荷载相关Ritz向量   | [vibration_mode.py:<br/>load_depended_ritz_vector()](libs/vibration_mode.py) |
| 4.4.1 |   矩阵迭代法：基本模态   |     [vibration_mode.py:<br/>mat_iterate_base()](libs/vibration_mode.py)      |
| 4.4.2 |   矩阵迭代法：高阶模态   |     [vibration_mode.py:<br/>mat_iterate_high()](libs/vibration_mode.py)      |
| 4.4.3 |  矩阵迭代法：最高阶模态   |    [vibration_mode.py:<br/>mat_iterate_highest()](libs/vibration_mode.py)    |
|  4.5  |     子空间迭代法     |    [vibration_mode.py:<br/>subspace_iteration()](libs/vibration_mode.py)     |
|  4.6  |   Lanczos方法    |          [vibration_mode.py:<br/>lanczos()](libs/vibration_mode.py)          |
|  4.7  |  Dunkerley方法   |         [vibration_mode.py:<br/>dunkerley()](libs/vibration_mode.py)         |
|  4.8  |   Jacobi迭代法    |          [vibration_mode.py:<br/>jacobi()](libs/vibration_mode.py)           |
|   -   |    振型归一化函数     |         [vibration_mode.py:<br/>normalize()](libs/vibration_mode.py)         |

</div>

### 3.3 第五章

<div class="center">

| 对应章节  |        描述         |                                    对应文件或函数                                     |
|:-----:|:-----------------:|:------------------------------------------------------------------------------:|
| 5.1.2 |       快速FFT       |           [modal_analysis.py:<br/>fourier()](libs/modal_analysis.py)           |
| 5.1.3 | 多自由度实模态<br/>振型分解法 |     [modal_analysis.py:<br/>modal_superposition()](libs/modal_analysis.py)     |
| 5.2.2 | 多自由度复模态<br/>振型分解法 | [modal_analysis.py:<br/>complex_modal_superposition()](libs/modal_analysis.py) |
| 5.3.2 |     伪加速度地震反应谱     |        [modal_analysis.py:<br/>pse_spectrum()](libs/modal_analysis.py)         |
| 5.3.2 |    绝对加速度地震反应谱     |        [modal_analysis.py:<br/>abs_spectrum()](libs/modal_analysis.py)         |
| 5.4.3 |       规范设计谱       |       [modal_analysis.py:<br/>design_spectrum()](libs/modal_analysis.py)       |
| 5.4.4 |       SRSS法       |            [modal_analysis.py:<br/>srss()](libs/modal_analysis.py)             |
| 5.4.4 |       CQC法        |             [modal_analysis.py:<br/>cqc()](libs/modal_analysis.py)             |
| 5.4.5 |      振型参与质量       |         [modal_analysis.py:<br/>modal_mass()](libs/modal_analysis.py)          |
| 5.4.5 |     振型分解反应谱法      |   [modal_analysis.py:<br/>modal_response_spectrum()](libs/modal_analysis.py)   |

</div>

### 3.4 第六章

<div class="center">

|  对应章节   |                描述                 |                                             对应文件或函数                                              |
|:-------:|:---------------------------------:|:------------------------------------------------------------------------------------------------:|
|  6.1.1  |           Duhamel积分解析式型           |         [sbs_integration_linear.py:<br/>duhamel_parse()](libs/sbs_integration_linear.py)         |
|  6.1.1  |           Duhamel积分数值型            |       [sbs_integration_linear.py:<br/>duhamel_numerical()](libs/sbs_integration_linear.py)       |
|    -    |      多线程Duhamel积分<br/>线程调用函数      |             [sbs_integration_linear.py:<br/>func()](libs/sbs_integration_linear.py)              |
|    -    |         多线程的Duhamel积分数值型          | [sbs_integration_linear.py:<br/>duhamel_numerical_multprocess()](libs/sbs_integration_linear.py) |
|  6.1.2  |               分段解析法               |       [sbs_integration_linear.py:<br/>segmented_parsing()](libs/sbs_integration_linear.py)       |
|    -    |             单自由度中心差分法             |   [sbs_integration_linear.py:<br/>center_difference_single()](libs/sbs_integration_linear.py)    |
|  6.1.3  |             多自由度中心差分法             |  [sbs_integration_linear.py:<br/>center_difference_multiple()](libs/sbs_integration_linear.py)   |
|    -    |         单自由度Newmark-beta法         |      [sbs_integration_linear.py:<br/>newmark_beta_single()](libs/sbs_integration_linear.py)      |
|  6.1.4  |         多自由度Newmark-beta法         |     [sbs_integration_linear.py:<br/>newmark_beta_multiple()](libs/sbs_integration_linear.py)     |
|  6.1.5  |         多自由度Wilson-theta法         |     [sbs_integration_linear.py:<br/>wilson_theta_multiple()](libs/sbs_integration_linear.py)     |
| 6.2.2_1 |             刚度退化双折线模型             |        [sbs_integration_nonlinear.py:<br/>Bilinear2()](libs/sbs_integration_nonlinear.py)        |
| 6.2.2_2 |             刚度退化三折线模型             |        [sbs_integration_nonlinear.py:<br/>Bilinear3()](libs/sbs_integration_nonlinear.py)        |
| 6.2.2_4 |            Bouc-wen模型             |         [sbs_integration_nonlinear.py:<br/>Boucwen()](libs/sbs_integration_nonlinear.py)         |
| 6.2.2_5 |     多自由度处理方法<br/>(绝对坐标转相对坐标)      |  [sbs_integration_nonlinear.py:<br/>bili_mdof();bouc_mdof()](libs/sbs_integration_nonlinear.py)  |
| 6.2.2_6 |   非线性多自由度双折线<br/>Newmark-beta法    |      [sbs_integration_nonlinear.py:<br/>newmark_bili()](libs/sbs_integration_nonlinear.py)       |
|    -    | 非线性多自由度Bouc-wen<br/>Newmark-beta法 |      [sbs_integration_nonlinear.py:<br/>newmark_bouc()](libs/sbs_integration_nonlinear.py)       |
|    -    |       非线性多自由度双折线<br/>中心差分法        |       [sbs_integration_nonlinear.py:<br/>center_bili()](libs/sbs_integration_nonlinear.py)       |
|    -    |     非线性多自由度Bouc-wen<br/>中心差分法     |       [sbs_integration_nonlinear.py:<br/>center_bouc()](libs/sbs_integration_nonlinear.py)       |

</div>

### 3.5 第七章

<div class="center">

| 对应章节  |         描述         |           对应文件或函数           |
|:-----:|:------------------:|:---------------------------:|
| 7.3.3 | 自定义ABAQUS-Python函数 | [absqus.py](libs/abaqus.py) |
</div>

### 3.6 第八章

<div class="center">

| 对应章节  |                 描述                 |                               对应文件或函数                               |
|:-----:|:----------------------------------:|:-------------------------------------------------------------------:|
| 8.1.1 |         支持向量机(SVM)鸢尾花分类任务          |         [SVM_iris:</br>8.1.1_1.py](parts/part8/8.1.1_1.py)          |
| 8.1.1 |          随机森林(RF)鸢尾花分类任务           |          [RF_iris:</br>8.1.1_2.py](parts/part8/8.1.1_2.py)          |
| 8.1.1 |          XGBoost波士顿房价回归任务          |      [XGBoost_boston:</br>8.1.1_3.py](parts/part8/8.1.1_3.py)       |
| 8.1.2 |            ANN波士顿房价回归任务            |        [ANN_boston:</br>8.1.2_1.py](parts/part8/8.1.2_1.py)         |
| 8.1.3 |      卷积神经网络(CNN)模型示例-Pytorch       |         [CNN_model:</br>8.1.3_1.py](parts/part8/8.1.3_1.py)         |
| 8.1.3 |    长短期记忆神经网络(LSTM)模型示例-Pytorch     |        [LSTM_model:</br>8.1.3_2.py](parts/part8/8.1.3_2.py)         |
| 8.2.3 |             数据清洗-重复值处理             |   [Data_clean_duplicates:</br>8.2.3_1.py](parts/part8/8.2.3_1.py)   |
| 8.2.3 |             数据清洗-缺失值处理             |     [Data_clean_fillna:</br>8.2.3_2.py](parts/part8/8.2.3_2.py)     |
| 8.2.3 |             数据清洗-异常值处理             |     [Data_clean_dropna:</br>8.2.3_3.py](parts/part8/8.2.3_3.py)     |
| 8.2.4 | 10折交叉验证(10-folds cross-validation) | [10-folds cross-validation:</br>8.2.4_1.py](parts/part8/8.2.4_1.py) |
| 8.3.7 |         SVM输电线路脱冰跳跃高度预测任务          |   [SVM_transmission line:</br>8.3.7_1.py](parts/part8/8.3.7_1.py)   |
| 8.3.8 |         随机森林输电线路脱冰跳跃高度预测任务         |   [RF_transmission line:</br>8.3.7_2.py](parts/part8/8.3.7_2.py)    |
| 8.3.9 |       XGBoost输电线路脱冰跳跃高度预测任务        | [XGBoost_transmission line:</br>8.3.7_3.py](parts/part8/8.3.7_3.py) |




</div>
