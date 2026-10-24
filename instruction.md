# 漂移 Maxwell 群的有限频窗动理学谱

为基于 PlasmaPy 的一维静电等离子体分析实现 `kinetic_spectrum`，给出指定复频率窗口内的全部离散模、残数，以及这些极点对时域响应的贡献。输入群分布及工况为合成数据。

在 `/app/src/kinetic_modes` 实现以下接口，可以在 `/app/src` 增加模块。容器已安装 PlasmaPy 2026.2.0、NumPy、SciPy 等依赖，无需外部数据。不得修改调用者的输入。

```python
kinetic_spectrum(fractions, thermal_widths, drift_frequencies,
                 window, sample_times)
```

## 科学定义

考虑均匀、无磁、无碰撞电子等离子体，正电背景固定。电子速度分布为若干一维 Maxwell 群的和。总电子等离子体频率为 ωp，实波数 k>0。第 j 群的数密度份额为 fj，速度标准差为 σj，平均漂移速度为 uj。输入均已归一化：

```text
z = ω/ωp
aj = k*σj/ωp        # thermal_widths
dj = k*uj/ωp        # drift_frequencies
τ = ωp*t           # sample_times
Σj fj = 1
```

时间因子采用 `exp(-iωt)`。Landau 延拓的等离子体色散函数及介电函数定义为

```text
Z(ξ) = i*sqrt(π)*exp(-ξ²)*erfc(-iξ)
ξj(z) = (z-dj)/(sqrt(2)*aj)
ε(z) = 1 + Σj fj/aj² * [1 + ξj(z)*Z(ξj(z))]
```

erfc 为整函数的解析延拓。所有频率和残数都对归一化变量 z 定义，不能换成对 ξ 或物理频率 ω 的导数。

`window=[left,right,bottom,top]` 定义开矩形 `left<Re(z)<right`、`bottom<Im(z)<top`。求其中 ε 的全部零点 zn；本题支持域内它们均为简单零点。每个模的残数为

```text
rn = 1 / ε'(zn)
```

本接口的响应仅为该有限窗口内的极点和：

```text
R(τ)  = Σn rn*exp(-i*zn*τ)
R'(τ) = Σn (-i*zn)*rn*exp(-i*zn*τ)
```

它不包含窗口外的极点或其他反演贡献。增长模和阻尼模均计入，只要位于窗口内。

## 输入、输出与误差

前三个参数为等长、非空的一维有限实数数值序列。`fractions` 严格为正，和与 1 的偏差不超过 `1e-12`；`thermal_widths` 严格为正。`window` 为四个有限实数，左右、上下端点分别严格递增。`sample_times` 为非空的一维非负有限实数序列，允许重复或不按顺序排列。

返回 dict，恰含下列字段；数组也可用同形状数值序列表示。M 为窗口内模数，K 为请求时间数。

|字段|形状|含义|
|---|---|---|
|`roots`|`(M,)`|全部 zn，不重复，顺序不限|
|`residues`|`(M,)`|与 roots 同序的 rn|
|`response`|`(K,)`|按 sample_times 原顺序的 R|
|`response_derivative`|`(K,)`|同顺序的 R'|

所有结果必须有限。没有模时，roots、residues 保留空轴，两个响应为零。不得只返回增长最快的模、每群一个模或实轴近似模。每个返回零点与期望零点按一一对应比较，对应残数沿用同一配对；不强制排序。同一输入重复调用、群顺序改变或将完全相同的群拆成数密度份额之和相同的群，物理结果应一致。

复数误差取模。roots 使用 `atol=2e-8, rtol=2e-8`，residues、response 和 response_derivative 使用 `atol=2e-8, rtol=2e-7`，均按 `|实际-期望|≤atol+rtol*|期望|` 判断。模式数及返回形状精确匹配。

## 支持域与运行条件

群数为 1 至 4，`0.035≤aj≤1.2`、`|dj|≤2`，K≤16、`0≤τ≤10`。矩形端点满足 `-4≤left<right≤4`、`-2≤bottom<top≤1`。窗口内至多 96 个模，任意两模距离至少 `1e-5`，残数模不超过 1000；介电零点与窗口边界的距离至少 `1e-6`。调用方保证这些正常分支条件，不要求诊断所有支持域外的退化情形。

空群、形状不匹配、非有限或非实输入、非正份额/热宽、份额未归一化、非法矩形，以及空、负或非有限时间，应抛出 ValueError，且不修改输入。

每次调用限 180 秒、2 CPU、4 GiB，可使用临时目录。只有 `/app/src` 源码进入独立验收环境。
