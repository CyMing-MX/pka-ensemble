# pka_replica — 本地有机分子 pKa 预测工作流

开源的水相有机分子 pKa 预测工具:**位点级(微观)与分子级(宏观)两条计算
路线、微态系综热力学、Web 界面、批量计算与 xlsx 报表**,并提供免安装的
Windows 便携版(见 Releases)。

## 1. 它怎么算

### 1.1 微态系综热力学(宏观 pKa / pI / logD)

遵循 Uni-pKa(JACS Au)与 Starling(Wagen 2025)的热力学框架(eq. 1–4):

```
g_i(pH) = g_i + z_i·ln10·pH                    ← 质子化学势
w_i(pH) = exp(−g_i(pH))/Z(pH)                  ← 微态分布
G°(z)   = −ln Σ_{i∈z} exp(−g_i)                ← 电荷级配分函数 (log-sum-exp)
pKa(z→z−1) = (G°(z−1) − G°(z))/ln10            ← 宏观 pKa
```

流程:SMARTS 模板微态枚举 → 每微态自由能 → Boltzmann 配分 → 宏观 pKa、
pI、微态分布与 logD 曲线。logD(pH) = Σ w_i(pH)·logP_i(每微态的
Wildman–Crippen logP 按占比加权)。

### 1.2 两条能量路线

**sites 路线(零依赖,毫秒级)**:SMARTS 官能团常数表 + 库仑耦合
(ε_eff=21)。方法学上属于"片段常数"家族(与 ChemAxon pKa 插件的
Csizmadia 片段常数 + Hammett/Taft LFER 同类,本实现是其不含 σ 传播项的
简化版)。适合常见官能团的快速粗估。

**xtb 路线(物理,推荐)**:逐微态热力学循环:

1. ETKDG+MMFF 母体构象,全部质子态由同一坐标删/加 H 派生(几何复用);
2. (可选)AIMNet2 神经网络势(wB97M-D3)气相预优化,把几何快速带到
   盆地底部;
3. GFN2-xTB 气相优化——必须收敛到 GFN2 自身极小点(跨引擎极小点采样会
   引入 2–3 pKa 噪声);
4. GFN2-xTB + ALPB 水单点;
5. 六类官能团特异线性校准(O/N/S-acid、N-base-ali/aro/aniline、
   O-acid-fluoro)ΔG→pKa;多质子态做"校准加和 + 真实多体残差(±4 截断)"
   分解;
6. 筛选门控:单翻转态先算,位点 pKa 落在窗口 ±5 之外直接丢弃;化学合理
   性过滤(弱酸 N-H/脂肪醇的校准伪影、pH 0–14 窗口)。

校准参数可用 `calibrate_fit.py` 按上述协议重拟合。

## 2. 精度

所有基准为**实验值对照**(本工具不含任何训练/拟合数据泄漏;校准仅用
少量小分子参考物,详见仓库内 `refit_*.py` / `calibrate_fit.py`)。

**pKahub 基准**(Sipos-Szabó et al., JCIM 2026,monoprotic 子集):

| 路线 | n | best-match MAE |
|---|---|---|
| xtb 多级协议 | 46 抽样 | **1.08** |
| sites | 6,084 | 1.72 |

文献参照(同数据):MolGpka ≈0.75、QupKake ≈0.79、pKaSolver ≈0.89、
Chemaxon ≈0.94(商业)、Uni-pKa macro ≈0.70。定位:无训练的规则+物理
路线,接近商业工具水平(注意 xtb 样本量较小)。

**药物集对比**(32 个肿瘤药物,vs Chemicalize 预测值):

| 口径 | n | MAE | 中位 |
|---|---|---|---|
| pKa-A(最酸) | 19 | 1.72 | 1.77 |
| pKa-B(最碱) | 26 | 1.17 | 0.50 |
| pI | 20 | 1.45 | 0.96 |

**PFAS 对比**(11 种全氟羧酸/磺酸,vs Chemicalize):全氟羧酸 ×7 误差全部
≤0.5,全集 MAE 0.63(独立的 O-acid-fluoro 校准类)。

单位点抽检:乙酸 4.42(4.76)、苯酚 10.89(9.99)、吡啶 5.23(5.23)、
咪唑 6.93(6.95)、苯胺 4.64(4.63)、布洛芬 4.36(4.43)。跑分脚本:
`benchmark.py` / `benchmark_pkahub.py` / `acid_analytes_run.py` / `pfas_run.py`。

## 3. 使用

**便携版(推荐,Releases 页下载)**:解压 → 双击 `pka-server.exe` →
浏览器自动打开 `http://127.0.0.1:8765`。命令行:`pka-cli.exe "58-08-2"`。
完整说明见包内 `操作手册.md`。

**源码安装**:

```bash
pip install -e .[aimnet]      # aimnet 为可选的神经网络几何引擎
pka-server                    # Web UI → http://127.0.0.1:8765
pka-cli "CC(=O)O"             # 宏观 pKa
pka-cli "58-08-2" --micro     # CAS 号:自动经 PubChem 解析为 SMILES
```

依赖:`rdkit`(+ 可选 `aimnet[ase]`);`xtb` 可执行文件需在 PATH
(便携版已捆绑)。

**输入**:SMILES、CAS 号或化合物名称(无效 SMILES 自动经 PubChem 解析,
需联网);批量支持 txt/csv 上传(≤50 分子,逗号/分号/制表符分隔)。

**功能**:电荷级分布卡(每净电荷级一条曲线 + 主导微态结构图,悬停读数)、
位点级表格(母体→共轭结构图,高亮可解离原子)、⟨z⟩/logD 悬停曲线、
会话历史(同分子免重算)、PubChem+RDKit 双源物化性质卡、Chemicalize 风格
xlsx 长表导出(LogD / isoelectric,pH 0–14 步长 0.1)、MolGpka /
SwissADME 交叉验证友链。

## 4. 已知局限

1. **两性电解质与多碱耦合**:第二解离常受隐式溶剂的耦合误差影响
   (哌嗪类对称二胺第二 pKa 偏高 ~5,此时 sites 路线反而更准);微态系综
   的多体质子耦合按可加 + 有界残差近似;
2. 高电荷端(|z|≥2)多位羧酸偶发垃圾极小点,靠 pKa 窗口过滤兜底;
3. 枚举上限 200 微态、无显式互变异构体扩展、无 RRHO 振动校正、
   构象系综简化为单构象;
4. 校准为经验线性外推——超出参考物化学空间(如全新骨架)时精度不可靠;
   超酸/超碱(pKa 出 0–14 窗)仅报告钳位值。

## 5. 许可

本项目代码 MIT;分发二进制中捆绑的第三方组件(xtb LGPL-3.0、RDKit
BSD-3-Clause、PyTorch BSD-3-Clause、warp、AIMNet2 等)依其自身许可,
见包内 `LICENSES/` 目录。

## 6. 参考

- Uni-pKa:Luo et al., _JACS Au_(微态系综热力学 eq.1–4);
  [代码](https://github.com/dptech-corp/Uni-pKa)
- Starling:Wagen et al., ChemRxiv 2025(多级优化与校准思路)
- AIMNet2 neural potential:[isayevlab/aimnetcentral](https://github.com/isayevlab/aimnetcentral)
- xTB:GFN2-xTB 与 ALPB 溶剂模型,Bannwarth et al.
- pKa 基准:Sipos-Szabó et al., _JCIM_ 2026(pKahub)
- 片段法方法学:Csizmadia et al.(ChemAxon pKa 插件);Wildman–Crippen logP
