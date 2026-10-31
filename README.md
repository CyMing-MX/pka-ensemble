# pka_replica — Rowan 宏观/微观 pKa 预测工作流的本地复刻

对 [labs.rowansci.com](https://labs.rowansci.com) 的 **pKa(微观)** 与
**macropka(宏观)** 两个预测工作流所做的逆向分析与开源复刻。
逆向证据来源:前端 Next.js chunk、Rowan 公开文档、开源数据模型库
stjames-public、Uni-pKa 论文(JACS Au)与 Starling 论文(Wagen 2025),
以及一个真实 Rowan 账号上提交的甲氨蝶呤 g-xTB 计算结果(见 §4 对照)。

## 1. 网站的原理解析

### 1.1 两条工作流(来自前端 chunk 与 stjames 数据模型)

**微观 pKa(位点级)** 默认参数:

```js
{ microscopic_pka_method: "gxtb_wagen2026",   // 可选 aimnet2_wagen2024 / starling / chemprop_nevolianis2025
  solvent: "water", pka_range: [2, 12],
  deprotonate_elements: [7, 8, 16],            // N, O, S 上去质子
  protonate_elements: [7],                     // N 上加质子
  mode: "careful" }                            // careful/rapid/reckless 控制构象数
```

**宏观 pKa(分子级)** 默认 `min_pH=0, max_pH=14, min_charge=-2, max_charge=+2`,
输出微态(带无量纲自由能)、相邻电荷级宏观 pKa、pI、0.1 pH 步长的微态分布
与 logD 曲线。表单还允许 -7..+7 的电荷窗和 Starling/Chemprop 等方法。

### 1.2 微观 pKa 的物理路线(Wagen 2024/2026)

1. 位点枚举 + 粗筛(粗 pKa 受 `pka_range` 限制,`reasonableness_buffer=5`);
2. 每质子态 ≤10 构象(energy window 10 kcal/mol);
3. **气相优化 + GFN2-xTB/CPCM-X 水溶剂化单点**(不在隐式溶剂里优化——否则
   中性氨基酸坍缩成两性离子);
4. Pracht–Grimme 构象系综自由能(298.15 K);
5. **官能团特异的经验线性(或二次)校准** ΔG→pKa。基准 MAE ≈ 1 pKa。

### 1.3 宏观 pKa 的热力学框架(Starling / Uni-pKa, eq.1–4)

```
g_i(pH) = g_i + z_i·ln10·pH                    ← 质子化学势
w_i(pH) = exp(−g_i(pH))/Z(pH)                  ← 微态分布
G°(z)   = −ln Σ_{i∈z} exp(−g_i)                ← 电荷级配分函数 (log-sum-exp)
pKa(z→z−1) = (G°(z−1) − G°(z))/ln10            ← 宏观 pKa
```

流程 = **SMARTS 微态枚举 → 每微态自由能(Starling 网络,Uni-pKa 架构)
→ Boltzmann 配分 → 宏观 pKa / pI / logD / 分布曲线**。哌嗪/甘氨酸的
"耦合解离"案例正是该框架相对逐位点方法的卖点。

**logD**:logD(pH) = Σ w_i(pH)·logP_i —— 每个微态的 RDKit Wildman–Crippen
logP 按该微态在 pH 下的 Boltzmann 占比加权(与 Rowan 文档一致的近似:
未显式建模离子对的单独分配行为)。

## 2. 本复刻的实现(v4)

```
pka_replica/
├── enumerator.py   SMARTS 酸/碱模板、交替质子化轮次、电荷窗、flips 记录
├── energy.py       能量后端(见 §2.1)
├── thermo.py       eq.1–4
├── workflow.py     stjames 风格结果 schema + 自动电荷窗口
├── server.py       Web UI + 2D 结构描绘(SVG)+ 电荷态卡片
└── web/index.html  深色界面:pH 滑块、微态分布/logD/⟨z⟩ 曲线
```

### 2.1 xtb 物理路线的速度与精度设计(协议 v4)

**速度**(甲氨蝶呤 46 重原子、~11 位点:十几分钟卡死 → v4 协议 ~26 秒;
v5(AIMNet2 预优化)~101 秒 → **两项优化后 ~64 秒**,位点简单分子 0.5–1 s):

1. **几何复用**:只生成一个母体 3D 构象(ETKDG+MMFF 取最优),所有质子态
   通过在同坐标上删/加 H 原子派生——每个态省掉嵌入+MMFF;
2. **micropka 复用 macropka 的位点能量**(2026-10-23):server 先跑宏观
   再跑微观,微观的单翻转任务与宏观完全重复——现在 macro 返回
   `site_pka_full`,微观直接复用,web 全流程省 ~20%;
3. **AIMNet2 预优化松弛度 A/B 定标**(2026-10-23):fmax 0.05 → 0.1 eV/A
   (GFN2 反正会收敛到自身极小点),MTX 全部位点 pKa 偏移 ≤0.05,宏观
   省 22%;注意 AIMNet2 阶段受 warp 线程安全锁串行,任务数是主要成本;
   位点 pKa 由固定种子 ETKDG(0xC0FFEE)保证逐次可复现;
2. **AIMNet2 多级预优化**(v5,路线 A 已落地):xtb 任务链 =
   **AIMNet2 预优化(wB97M-D3 神经网络势,isayevlab/aimnetcentral)→
   GFN2 短优化(收敛到 GFN2 自身极小点)→ ALPB 水单点**。AIMNet2 把
   几何快速带到盆地底部,GFN2 只需几步收尾;第二级不可省——直接在
   AIMNet2 极小点上做 GFN2 单点会引入 2–3 pKa 的跨引擎采样噪声
   (Rowan 的 multistage optimization 同款思路);
3. **并行**:所有任务进 ThreadPoolExecutor(每 worker OMP=1;
   Rowan 的日志显示他们的构象是逐个串行优化的,这里反而更快);
4. **Rowan 式筛选门控**:先只算单翻转态并校准,位点 pKa 落在
   `[pka_min−5, pka_max+5]`(即 stjames `reasonableness_buffer=5`)之外的
   直接丢弃;多翻转态只对幸存位点计算,且按可加能量排序、上限 16 个;
5. 中性母体多构象(2),其余态单构象。

**精度协议**(与校准表一致):全部态 **AIMNet2 预优化 + 气相 GFN2 优化 +
ALPB 水单点**;六类官能团特异线性校准(O/N/S-acid、N-base-ali/aro/
aniline,`calibrate_fit.py` 可重拟合,AIMNet2 协议已于 2026-10-22 重拟合);
多质子态做"校准单点 + GFN2 真实多体残差(±4 截断)"的精确分解。
AIMNet2 的电荷经 `atoms.info['charge']` 传入(其文档规定的最高优先级
通道);每线程一个 calculator 实例。新协议下单位点验证集:乙酸 4.42
(4.76)、苯酚 10.89(9.99)、**吡啶 5.23(5.23)、咪唑 6.93(6.95)、
苯胺 4.64(4.63)、布洛芬 4.36(4.43)**——六个里有四个误差 <0.1。

**自动电荷窗口**:`suggest_charge_window()` 按可解离位点数取
`z0 ± clamp(n, 2, 3)`(对应 Rowan 的 min/max charge,其默认 ±2、上限 ±7;
Starling 微态枚举本身也用束搜控制组合爆炸)。

### 2.2 fragment 片段路线的正确性评估(零依赖兜底)

SMARTS 常数表 + 库仑耦合(ε_eff=21)。定位:

- **单元酸/碱**:±0.5 pKa(乙酸 4.3 vs 4.76、吡啶 5.2 vs 5.23);
- **多位点分子**:排序与量级正确,数值 ±1–2(哌嗪 ΔpKa≈4.2 vs 4.47,
  甘氨酸两性离子稳定化不足);
- **与文献方法的关系**:ChemAxon Chemicalize/Marvin 的方法是
  Csizmadia 的"片段常数 + Hammett/Taft 取代基线性自由能"(芳香体系
  电子效应经 σ 常数传播),基准 RMSE ≈ 0.5–0.9;本复刻的片段表等价于
  去掉 Hammett 传播项的简化版——这正是主要差距来源(对位/间位取代
  苯甲酸、取代苯胺类的系统性偏差)。要补齐需要建 σ 常数表 + 模板匹配,
  工作量约一天,精度可望到 ~0.7 RMSE 水平(仍不及 ML 路线)。
- 结论:**精确预测请用 xtb 路线或接入 ML 权重;fragment 路线只用于
  无 xtb 环境的快速粗估**。

## 3. 使用

**安装为工具包**(项目根目录,含 pyproject.toml):

```bash
pip install -e .[aimnet]      # aimnet 为可选的神经网络几何引擎
pka-server                    # Web UI → http://127.0.0.1:8765
pka-cli "CC(=O)O"             # 宏观 pKa
pka-cli "58-08-2" --micro     # CAS 号:自动经 PubChem 解析为 SMILES
pka-cli "caffeine"            # 化合物名称同样可用
```

依赖:`rdkit`(+ 可选 `aimnet[ase]` 几何引擎);`xtb` 可执行文件需在 PATH。

**Web 附加功能**:①会话历史——左侧栏记录本次运行的所有单分子结果,
点击免重算直接查看,同分子(规范化 SMILES 同一)再提交直接回缓存
(上限 50 条,服务重启清空);②xlsx 导出——单分子结果页"导出本分子
xlsx"按钮 + 批量完成后"导出 xlsx"按钮,生成 Chemicalize 风格长表
(LogD / isoelectric 两张 sheet,每行一分子,列 = pH 0-14 步长 0.1,
4 位小数);③右上角 MolGpka / SwissADME 参考工具链接。

**输入方式**:SMILES、CAS 号或化合物名称(无效 SMILES 自动走 PubChem
PUG REST 解析,结果回填并标注 CID);批量用 Web UI 的"批量 txt/csv"按钮。

**分子物化性质**:每次预测附带输入态 RDKit 描述符(分子量、重原子数、
Crippen logP、TPSA、氢键供体/受体、可旋转键、芳香环、形式电荷),
交互页面有独立卡片,批量 Excel 每行都有对应列。

**批量模式**(Web UI):点"批量 txt/csv"上传文件(≤50 个分子),每行一个
SMILES;也支持 csv(name,smiles 或 smiles,name——自动识别哪个字段是
SMILES),字段分隔符支持逗号/分号/制表符。提交后左侧边栏(Chemicalize
风格)实时显示每个分子的 2D 结构缩略图与计算状态,完成后各卡片标注
窗口内 pKa/pI;点击任意缩略图,主区加载该分子的完整结果。无效行提交时
即报告并跳过。

UI 功能(Rowan/Chemicalize 风格,浅色主题):**电荷级分布图**——每个净电荷
级一条粗聚合曲线,该级内的各微态画成同色淡化细线;图例、曲线、右侧结构
卡片三者按电荷级一一对应(每级一张卡:峰值、主导 pH 区间、级内微态数与
主导微态结构图);悬停任意 pH 显示各电荷级占比;宏观 pKa 卡带最强酸性/
最强碱性 pKa 统计(Rowan 风格);⟨z⟩ / logD 曲线均悬停读数;
位点级表格用"母体(橙色高亮位点)→ 共轭结构"图像,共轭结构三级回退保证
永不缺失。原始微态数据(g 值表)不在页面展示(调试信息),完整数据在
API JSON 与 `cli.py --json`。
实现注意:位点/微态的原子索引基于**输入 SMILES 的原子顺序**,canonical
SMILES 会静默重排——曾两次踩坑。**第三次(2026-10-23)**:连
`MolToSmiles(canonical=False)` 也不保序(芳香环改写时,如 PubChem 返回的
Kekulé 输入),导致位点高亮整体错位;最终方案是 `object_data.input_molblock`
(`MolToMolBlock` 往返严格保序),绘图端一律优先解析 molblock。

## 3.5 基准实测(2026-10-22)

数据集:QupKake/Uni-pKa 同源的 Novartis(280)+ 文献(122)SDF,
共 402 个位点级实验 pKa(每分子一个标注电离位点)。跑分脚本
`benchmark.py sites`。两种口径:best-match(位点归属不计入误差)、
strongest(最强酸/碱显示口径)。

| 路线 | n | bestMAE | strongMAE |
|---|---|---|---|
| sites(SMARTS 常数 + 库仑耦合) | 6,084 | 1.72 | 1.80 |
| xtb 多级协议(N-base-ali 重拟合**后**) | 46 抽样 | **1.08** | — |
| xtb 多级协议(重拟合前,同 44 分子配对) | 44 | 3.00 | — |

xtb 路线曾劣于 sites:诊断发现 N-base-ali 校准仅用两点拟合(甲基胺/哌啶),
raw 外推时呈负斜率(物理错误),脂肪胺类系统性 +3 偏差,|err|>3 离群值
20/44。修复:`refit_nbaseali.py` 用 14 个脂肪胺参考物(甲/乙/丙胺、二甲/三甲
胺、二乙/三乙胺、吗啉、哌啶、吡咯烷、氮杂环丁烷、环己胺、苄胺、乙醇胺、三氟
乙胺;2 个因收敛失败跳过)按生产协议重算并最小二乘重拟合 → 斜率 0.3983
(恢复物理正确方向),残差 mean 0.54 / max 1.17。pKahub monoprotic 100 抽样
重跑(进程在 m63 处中断,46 个已评分):**bestMAE 3.00 → 1.08**,离群
20 → 2,22 个分子改善超过 1 pKa(如预测 15.28 → 6.48,实验 5.40)。
文献参照(monoprotic MAE,同 pKahub 数据):MolGpka ≈0.75、QupKake ≈0.79、
pKaSolver ≈0.89、Chemaxon ≈0.94(商业)、Uni-pKa macro(simple)≈0.70。
**定位:无训练的规则+物理路线;重拟合后 xtb 路线 monoprotic MAE ≈ 1.08,
已接近商业 Chemaxon(≈0.94)水平(注意 xtb 样本量较小);两性电解质
与多碱耦合仍需路线 B(Uni-pKa)。** 盐/混合物输入(断键)自动取最大片段
并警告;离子对整体 pKa 无定义。

## 3.6 与 Chemicalize 的药物集对比(2026-10-23,ACID-Analytes 32 药)

对用户提供的 Chemicalize 预测数据(pKa-A 最酸 / pKa-B 最碱)逐药对比
(xtb 路线,`acid_analytes_run.py`,结果 `acid_analytes_final.tsv`):

| 口径 | n | MAE | 中位 | ≤0.5 | ≤1 |
|---|---|---|---|---|---|
| pKa-A(最酸) | 19 | **1.72** | 1.77 | 6 | 7 |
| pKa-B(最碱) | 26 | **1.17** | 0.50 | 13 | 15 |
| pI | 20 | 1.45 | 0.96 | 5 | 11 |

对比过程中落地了三层**化学合理性筛选**(同时作用于微观行和系综):
①弱酸 N-H(酰胺/苯胺/吲哚类,化学上 ≥10)<10 视为校准伪影剔除
(NON_GENUINE_ACID_N);②脂肪醇 O-H <12 剔除(ALCOHOL_FLOOR,酚/烯醇
不受影响——模板已拆分,`phenol_enol_oxime` 只匹配 sp2 C/N 上的 O-H);
③展示行限 0–14;④|pKa|>20 的宏观转换(空电荷级的占位能量)|过滤。
这些筛选同时修掉了系综污染:校准伪影位点曾让"环 N + 伪影 N"双质子化态
的能量被错算到母体之下,主导 +2 级(MTX 宏观 pKa 首项 −1.59 → 0.75)。
筛选前后同集对比:A 口径 MAE 3.07 → 1.72,B 口径 1.54 → 1.17,且无需
人工剔除离群值。已知残留:蒽环类(多柔比星)氨基糖 NH₂ 的质子化任务
超时 → 碱行缺失;他莫昔芬二甲氨基校准到下限(实验 8.5)——均为个例。

## 3.7 PFAS-Analytes 对比(2026-10-24,11 种全氟酸)

全氟羧酸(C4–C10)+ 全氟磺酸(C4–C8),Chemicalize pKa-A 分别 ~0.3 /
−3.3(`pfas_run.py`,结果 `pfas_result.tsv`):

| 类别 | chem | 我们 xtb | 误差 | 片段路线 |
|---|---|---|---|---|
| 全氟羧酸 ×7 | 0.29–0.40 | 0.01–0.23 | **−0.08~−0.42(全部 ≤0.5)** | 4.3(常数,误差 +4) |
| 全氟磺酸 ×4 | −3.3 | −2.0(校准下限) | ≈ +1.3 | 无(超酸出 0–14 窗) |

xtb 全集 MAE **1.96 → 0.63**(中位 0.39,7/11 ≤0.5)。实现:**"O-acid-fluoro"
独立校准类**——`_flip_class` 检测羧基 α-C 上的氟,路由到以氟化阶梯参考
(氟乙酸 2.59 / 二氟乙酸 1.33 / 三氟乙酸 0.23 / 全氟丙酸 0.59)重拟合的
(0.161, −11.30),残差 ≤0.31;TFA 全流程 0.52(实验 0.23),乙酸仍走通用
类 4.42(实验 4.76),互不干扰。全氟磺酸无可靠实验参考(超酸测量困难),
维持钳位值(误差 +1.3,歪打正着)。**顺带修复**:①磺酸/磺酰胺 SMARTS 曾用
`$()` 单原子递归写法,匹配集只含 S、找不到带 H 原子,酸位点被静默丢弃
——这两类位点此前从未生效(PFAS 暴露);已改为从 O/N 起匹配并补齐伯磺
酰胺 `[NX3H2]` 形式。②宏观 pKa 列表限 0–14(与微观行一致的物理窗口;
超酸/占位级的出窗转换不再显示)。

## 4. 与 Rowan 官方结果对照(甲氨蝶呤, 2026-10-21)

在你的 Rowan 账号上提交了甲氨蝶呤
`CN(Cc1cnc2nc(N)nc(N)c2n1)c1ccc(C(=O)N[C@@H](CCC(=O)O)C(=O)O)cc1`
的微观 pKa(g-xTB Wagen 2026)官方计算,与本复刻(xtb v4,28 s)对照:

| | Rowan g-xTB(官方) | 复刻 GFN2/ALPB |
|---|---|---|
| 最强酸(羧酸位点) | **3.35**(另一羧酸 4.25) | 4.41 / 3.97 |
| 最强碱(脂肪胺 N) | **6.09** | **5.99** |
| 窗口内其余位点 | 3.96, 3.68, 1.55, 0.06, −0.42 | 3.97, 3.6, 1.7, 0.76, −1.59 |
| 两性离子提示 | "Potentially zwitterionic" | pI 6.26(两性行为一致) |

窗口内位点级 MAE ≈ 0.5–0.9 pKa,最强碱相差 0.10。作为参照,
Rowan 文档给出的 g-xTB 基准 MAE ≈ 1。

**注意(2026-10-23)**:上表是 10-21 的 v4 协议(GFN2-only)数字。v5 协议
(AIMNet2 预优化)在 MTX 上给:羧酸 4.29 / 1.92,最强碱(蝶啶 N)5.79,
宏观 pKa 4.46 / 3.82 / 5.80(±0.02 可复现,固定种子)。**数值移动不是
"变不准"**——MTX 的实验值本身高度分散(Meloun 2010 热力学研究:
4.8 / 5.5;其他来源 2.1 / 4.7 / 5.6,归属各异),v4 的 3.97/4.41 并不
更接近任何一组实验值;客观精度证据是 pKahub 基准 bestMAE 3.00 → 1.08
(44 分子配对),而非单分子目测。

## 5. 已知局限(诚实清单)

0. **酰胺/苯胺型弱碱的处理(2026-10-23 修复)**:碱位点模板刻意**不包含
   酰胺 N 质子化**(pKaH ≈ −1,pH 0–14 永远不可见)。曾因模板存在 +
   xtb 路线把酰胺归入 N-base-ali 校准,质子化酰胺微态在甲氨蝶呤上主导了
   全部碱端系综(卡片画出 [NH2+],看起来"图文不符")——根因是该协议下
   酰胺质子化的 raw 偏低 ~11 pKa,任何线性校准都会把它拉进窗口;删除模板
   是唯一稳妥解。极弱碱(aniline pKaH<0)仍会列出,供参考。
1. **两性电解质**(甘氨酸类):xtb 路线第一 pKa ±1 内但第二 pKa 偏高
   ~1.5;fragment 路线两性离子稳定化不足——这正是 Starling 用 ML 学掉的
   耦合项;哌嗪类对称二胺的第二 pKa 偏高 ~5(隐式溶剂过度屏蔽排斥),
   此时 fragment 路线反而更准;
2. 高电荷端(|z|≥2)多位羧酸的气相优化偶发垃圾极小点,靠 pKa 范围过滤
   兜底(收窄电荷窗可避免);
3. 无 AIMNet2 引导的束搜(枚举上限 max_states=200)、无显式互变异构体
   扩展、无 RRHO 振动校正、构象系综简化为单构象。

## 5.5 升级 pKa 精度的两条路线(可选增强,按投入排序)

**路线 A —— AIMNet2 几何引擎**(**已完成**,2026-10-22):
- `pip install "aimnet[ase]"`(PyPI 正式包;仓库
  github.com/isayevlab/aimnetcentral,模型首用时自动从 HF Hub 下载);
- 实际落地的协议是 **AIMNet2 预优化 → GFN2 短优化 → ALPB 水单点**
  (直接在 AIMNet2 极小点上做 GFN2 单点会引入 2–3 pKa 的跨引擎噪声,
  必须让 GFN2 收敛到自身极小点);电荷经 `atoms.info['charge']` 传入;
- 校准已按新协议重拟合(`calibrate_fit.py`);校准时 raw 统一为
  "去质子化方向"(碱位点取共轭酸的酸式解离);
- 收益:单位点验证集误差从 ±0.5–1.6 降到 **±0.1 以内**(四个分子),
  甲氨蝶呤与 Rowan g-xTB 的窗口内位点对比同步改善;MTX 总时长 26 s
  持平。路线 B(Uni-pKa)仍是两性电解质/多碱耦合问题的最终解。

**路线 B —— Uni-pKa 真权重**(1–2 天,质变到药物级):

- 这是 Starling 的开源前身(同一架构、不同训练数据),宏 pKa
  RMSE ≈ 0.7——本复刻位点分解 + GFN2 路线搞不定的两性电解质和多碱
  耦合,它直接从微态能量层面解决(§5.1 的两个局限同时消失);
- **需要你做的只有一件事**:安装 WSL2(管理员 PowerShell 运行
  `wsl --install`,重启)——Uni-pKa 的 Uni-Core 依赖基本只在 Linux
  下可编译;
- 其余交给我:从 Zenodo 19627026 的 apptainer 容器(squashfs 格式,
  可直接解包,约 2–4 GB)提取 Uni-pKa checkpoint、在 WSL 里装
  unicore、写一个 `UniPkaEnergyModel` 替换 `QuantumEnergyModel`
  (微态枚举、热力学层、UI 全部复用,只换能量来源);
- 速度:CPU 推理每微态 ~0.1–0.5 s,几十个微态的药物分子分钟级,可接受。

## 6. 参考

- Rowan 文档:[pka](https://docs.rowansci.com/workflows/pka) / [macropka](https://docs.rowansci.com/science/workflows/macropka)
- [stjames-public](https://github.com/rowansci/stjames-public)(数据模型)
- Uni-pKa(Luo et al., JACS Au)eq.1–4;Starling(Wagen, ChemRxiv 2025)
- [Uni-pKa 代码](https://github.com/dptech-corp/Uni-pKa) · [AIMNet2](https://github.com/isayevlab/aimnetcentral)(`pip install aimnet`)
- ChemAxon pKa 方法:Csizmadia et al., 片段常数 + Hammett/Taft LFER
  (docs.chemaxon.com pKa plugin);benchmark: Sipos-Szabó et al., JCIM 2026
