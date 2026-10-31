# pka_replica 便携版操作手册

**版本 1.0 · Windows 10/11 · 无需安装、无需 Python**

pka_replica 是 Rowan Scientific 宏观/微观 pKa 预测工作流的本地开源复刻:
SMARTS 片段路线 + GFN2-xTB 量子物理路线(+ AIMNet2 神经网络势预优化,
取决于你下载的包),支持单分子预测、批量 txt/csv、PubChem 名称/CAS 解析、
物化性质卡、logD/等电点曲线与 xlsx 长表导出。

---

## 1. 选哪个包?

| 包 | 大小(解压后) | 计算路线 | 适合 |
|---|---|---|---|
| `pka_replica_v1.0_lite.zip` | ~300 MB | sites + xtb | 大多数场景;磁盘紧张 |
| `pka_replica_v1.0_full.zip` | ~1.5 GB | + AIMNet2 预优化(与开发者本机完全一致) | 想要最高精度 |
| `pka_replica_v1.0_ai-addon.zip` | ~1 GB | —— 仅含 `_ai` 文件夹 | **已下轻量版,后想升级** |

升级方法:把 `ai-addon.zip` 里的 `_ai` 文件夹整个解压到 pka_replica 根目录
(`pka-server.exe` 旁边),重新启动即变为完整版。删掉 `_ai` 文件夹即还原轻量版。

## 2. 快速开始

1. 解压 zip 到任意**本地磁盘**目录(路径建议不含中文和空格;不要解压到网络盘);
2. 双击 **`pka-server.exe`** —— 会弹出一个黑色控制台窗口并自动打开浏览器
   `http://127.0.0.1:8765`(端口被占用时自动顺延,以控制台显示的地址为准);
3. 在输入框粘贴 SMILES / CAS 号 / 化合物名称(后两者经 PubChem 解析,需联网),
   点"预测 pKa"。快捷标签(乙酸、苯酚……)可一键试算;
4. **关闭程序 = 关掉那个黑色窗口**。历史记录保存在内存,关闭即清空。

> 安全提示:程序只监听本机 127.0.0.1,不对外网开放;如 Windows 防火墙首次
> 弹窗,选择"取消/仅专用网络"即可正常使用本机访问。
>
> 高级:设置环境变量 `PKA_PORT` 可改起始端口(默认 8765,被占自动顺延)。

## 3. 命令行用法

在解压目录打开 CMD 或 PowerShell:

```
pka-cli.exe "CC(=O)O"                      :: 宏观 pKa(乙酸)
pka-cli.exe "58-08-2"                      :: CAS 号 → PubChem → pKa
pka-cli.exe "caffeine" --micro             :: 化合物名,位点级 pKa
pka-cli.exe "NCC(=O)O" --backend sites     :: 快速片段路线(毫秒级)
pka-cli.exe "CC(=O)O" --json out.json      :: 完整结果存 JSON
```

常用参数:`--backend auto|xtb|sites`、`--micro`(只算位点)、
`--min-charge/--max-charge`(默认 auto)、`--json 文件名`。

## 4. 计算路线怎么选

界面 backend 下拉框:

- **auto**(默认):自动探测——有 xtb 用物理路线,没有则用 sites;
- **xtb**:GFN2-xTB/ALPB 逐微态量子计算,慢(小分子数秒~大分子 1-2 分钟)
  但普适;若 `_ai` 文件夹存在,自动叠加 AIMNet2 预优化;
- **sites**:SMARTS 片段库,毫秒级,常见官能团够用,新化学空间会失准。

**首次用 xtb 路线的提示**:首次运行 AI 增强会进行内核编译,头 10–30 秒
可能无输出,属正常现象;之后每次启动都很快。

## 5. 功能速览

- **历史记录**(左侧栏):本次运行内的所有单分子结果,点击免重算直接查看;
  相同分子再次提交直接命中缓存(换 backend 会真正重算);上限 50 条,关程序清空;
- **批量计算**:"批量 txt/csv" 上传文件,每行一个 SMILES/CAS/名称(支持
  逗号/分号/制表符分隔的 `名称,SMILES` 格式),完成后左侧栏可逐个查看,
  标题旁"导出 xlsx"一次导出全部;
- **导出 xlsx**:单分子结果页"导出本分子 xlsx";格式为 Chemicalize 风格
  长表(LogD / isoelectric 两张 sheet,pH 0–14 步长 0.1);
- **物化性质卡**:PubChem 在线数据(XLogP/TPSA/电荷/HBD/HBA 等)与
  本地 RDKit 计算(logP Crippen 等)分区展示;断网时自动只显示本地段;
- **参考链接**:右上角 MolGpka / SwissADME 交叉验证入口。

## 6. 常见问题

**Q:双击 exe 后黑色窗口一闪而过?**
看窗口内报错文字。最常见:解压不完整(重新解压)、杀毒软件拦截(加入
白名单后重试)。

**Q:提示"xtb executable not found"?**
确认 `bin\xtb.exe` 存在于解压目录;不要把 exe 单独拖出目录使用。

**Q:PubChem 解析一直失败?**
需要联网;公司内网/代理环境可能拦截,可改用直接粘贴 SMILES。

**Q:计算很慢?**
xtb 物理路线对大分子(50+ 重原子)本就需要 1–2 分钟;批量任务逐个排队。
赶时间用 `--backend sites` 或界面 sites 选项。

**Q:两性离子/盐类输入结果为空或奇怪?**
直接贴中性形式(如 `NCC(=O)O` 而非 `[NH3+]CC(=O)[O-]`);两性离子输入
现已支持,但极端电荷态(±3 以上)仍建议手动收窄电荷窗。

**Q:如何彻底卸载?**
删除解压目录即可;程序不写注册表,仅在首次运行 AI 增强时于
`%LOCALAPPDATA%\NVIDIA\warp` 留有约 1 MB 内核缓存,可一并删除。

## 7. 目录结构

```
pka-server.exe        双击启动 Web 界面
pka-cli.exe           命令行工具
_internal\            运行时库(勿动)
bin\                  xtb 引擎(勿动)
_ai\                  [仅完整版] AI 几何引擎(勿动;删除即变轻量版)
LICENSES\             第三方组件许可
```

## 8. 已知局限(诚实清单)

- 两性电解质第二 pKa、哌嗪类对称二胺偏差较大(隐式溶剂的固有局限);
- 全氟磺酸等超酸超出 0–14 展示窗口,按钳位值报告;
- 构象系综为单构象,无 RRHO 振动校正;
- 精度基准与出处见项目仓库 README(32 药 vs Chemicalize:最碱 MAE 1.17 /
  最酸 1.72;PFAS MAE 0.63)。
