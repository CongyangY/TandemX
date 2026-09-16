# TandemX：面向 Genome Research Methods 水平的第一轮差距审计与实验方案

**日期：**2026-09-16。**性质：**只读代码、数据表和文献审计后的方案；尚未实现新算法、运行新基准、打开新留出集或修改稿件。**仓库起点：**`2676aad`。本文件提出未来可预注册的实验与门槛，所有预测的收益均不是现有结果。最终是否达到目标水平只能由冻结后的独立结果判断。

## 0. 项目目标与当前判决

要回答的问题是：在复杂植物基因组中，TandemX 能否从原始 HiFi/ONT 长读段形成可解释的串联重复家族、可校准的丰度与组装表示不足判断，并在至少一个共同科学终点上取得相对于当前强基线的、限定输入和真值的优势。**当前答案：尚不能。** 已有可运行的家族级审计工作流、合格的来源记录、两种材料的回顾性方向证据；尚缺严格物理/编辑真值、真实材料同终点强基线、校准的不确定度、经验证的 HOR 顺序、ONT 去新主链与完整大基因组工程验收。现有 871 项 pytest 通过只说明代码回归，不是科学验证。

当前最合理的论文定位是“有来源和失败状态的家族级审计研究系统”。Methods 论文所需的核心算法贡献尚未建立。新阶段允许在**分离的研究原型**中检验方法假说；旧的生产估计器、旧实验、失败记录、v5 稿均保留原样。暂停 v5 修辞润色、Bioconda、Zenodo 和正式 tag，不把审计自动解释为对这些发布动作的批准。

### 证据基线及其反例

| 终点 | 已有证据 | 限制或反例 |
| --- | --- | --- |
| 发现/模拟定位 | 三个独立 10-Mb 模拟基因组的 55 个种植家族均在 5× 找到，20× 达到预设目录稳定；模拟定位平均 recall 0.981、precision 0.9996。 | 该模拟不覆盖真实复杂 HOR、多位点同源与植株/单倍型异质性；20× 不是植物通用深度。 |
| 冻结丰度与二元分类 | 模拟生产估计器 MARE 0.363；另一个留出二元测试 sensitivity 0.654、FPR 0.016、precision 0.983。 | `copy_number_interval_low/high` 是诊断词 10–90% 离散，不是校准 CI；词“唯一”只在当前目录内检验。 |
| 真实历史组装代理 | Ey15-2 15-kb 合格家族 8 TP/0 FN/0 FP/11 TN；Macadamia 2/0/0/41。 | 分母分别只含 19/43 个新组装已表示家族，排除新组装完全缺失者；新组装与 HiFi 有共享证据；5-kb 阈值出现 FN/FP；Macadamia 阳性仅 2 个。 |
| 差额幅度/独立读段 | Macadamia 读段–旧组装差额 2.119 Mb，对新旧组装增量 0.159 Mb；六个预选新组装差额候选中 3 个获独立方向支持、3 个未解决。 | 增量不是物理真实缺失量；跨 k/平台冲突的 `TXF000695` 不能认定塌缩。无个体重复或物理拷贝真值。 |
| 正式模拟竞争者 | 72/72 格完成；TandemX 54/54 planted family conditions，SRF k151 45/54、k101 53/54。 | 普通竞争性映射在六条件中正家族 MARE 均较低；共享片段阴性归属 TandemX 平均 30,887 bp，SRF k151 0、k101 227。其映射输入用了 TandemX 目录，去新召回不适用。 |
| 真实历史竞争者 | SRF/普通映射六格方案已冻结。 | 首格 SRF k151 原生组装到 7,200 s 技术超时，其余五格未运行；无真实历史同终点优劣结果。 |
| 性能 | Morex 1.170 Gb 输入，TandemX 1,909.525 s/913.891 MiB；TRF 2,409.145 s/286.172 MiB；TideHunter 566.196 s/580.375 MiB。 | 同机有并发，非多重复隔离排名；TandemX 不在速度或内存上统一领先；完整审计量化需数小时。现有 `--chunk-size/--chunk-bases` 是投喂上限而非阶段内检查点。 |

来源：[现状总账](current_status.md)、[专家状态及资源表入口](tandemx_expert_status_20260916_zh.md)、[正式 SRF 结果](srf_formal_results_20260910.md)、[投稿方法缺口](submission_method_gaps_20260910.md)、[v5 保留稿](../paper/0910/manuscript_v5.md)。`README.md` 部分 MVP 和稿件版本语句落后于后续真实实验；发布前应修正，而本轮不以修改它来掩盖证据差距。

## 1. 方法新意与 2026 竞争格局的修正

**当前真正可辨识的差异**是把一个操作性家族的读段发现、估计丰度、组装坐标、读段—组装差额与不确定/失败状态串到可审计记录。该整合本身有使用价值，但周期种子、序列聚类、k-mer 深度中位数、普通映射、标签序列图和阈值分类不能仅因被放在同一流水线而视作新算法。尤其不能写“首次连接 reads 与 assembly”：SRF 已从准确读段/组装重建卫星重复并比较表示，HiCAT-human 有 reads 与 assembly_match 路径，Merqury 有 read–assembly k-mer copy 谱。

**现行竞品需要重估。** [CENdetectHOR 官方 pipeline](https://github.com/CENdetectHOR/CENdetectHOR_pipeline) 接受单染色体/着丝粒 assembly FASTA，卫星 consensus 可选，执行高重复区定位、单体长度推断、单体提取、聚类与 HOR/家族分析，并配合 PhyloTreeGUI。其 [2025 年方法预印本](https://doi.org/10.1101/2025.01.07.631657) 含 human 与 Arabidopsis；[2026-04 的软件库 0.0.6](https://pypi.org/project/cen-detect-hor/) 是当前软件版本线索，**不是 2026 年方法论文**；pipeline 的运行 commit 仍须另行固定。它的官方输入路径不提供原始 FASTQ 丰度或 read–assembly 差额，故这两项应记为任务范围外，不能记低准确率。现行 [TideCluster](https://github.com/kavonrtep/TideCluster) 包括逐阵列单体/HOR 线索、TRC/superfamily、跨样本 assembly 比较和报告；不能再写其“没有 HOR、跨样本、不确定标签或大输入能力”。[HiCAT-human](https://github.com/865699871/HiCAT-human) 已有 HiFi reads 加组装路径，但偏人类着丝粒，不能当成植物无先验定量性能已经验证。[HiCAT 原版](https://pmc.ncbi.nlm.nih.gov/articles/PMC10053651/) 也在 Col-CEN 检查过 HOR；植物结构注释不是空白。

### 按任务决定可比性（✓ 原生；△ 有配套步骤或限定模式；— 设计范围外；? 当前能力尚待正式验证）

| 工具 | 输入及读段去新 | 单体/家族/HOR | 家族读段丰度 | assembly 阵列 | read–assembly 缺额 | 跨样本/不确定性 |
| --- | --- | --- | --- | --- | --- | --- |
| TandemX 当前 | HiFi reads ✓；ONT 去新 ? | 单体/操作家族 ✓；HOR 顺序 — | 诊断 k-mer 点估计 ✓；校准 CI — | ✓ | 家族差额/状态 ✓；物理塌缩真值 — | 基本 cohort △；校准联合 CI — |
| [TRF](https://github.com/Benson-Genomics-Lab/TRF) | FASTA 逐序列阵列/周期 ✓ | 跨读段家族/HOR — | — | 可逐 assembly 运行 ✓ | — | — |
| [TideHunter](https://github.com/yangao07/TideHunter) | 长读段逐条阵列/共识 ✓ | 原生跨读段家族/HOR — | 读内 copyNum ≠ 基因组家族 CN | 可逐 assembly 运行 ✓ | — | — |
| [SRF](https://github.com/lh3/srf) | 高精度读段/短读/组装去新 ✓ | 卫星单位与可能 HOR ✓ | 完整映射流水线的 mapped-bp 丰度 △ | 映射至组装 △ | 文献已有差额观察，内建统一判定待核 △ | k/count 敏感，失败/空目录需区分 |
| [TRASH/TRASH2](https://github.com/vlothec/TRASH) | assembly 工具，raw-read 去新 — | 阵列/单体/结构 ✓；TRASH2 的 HOR 脚本另定版本 | — | ✓ | — | 结构/跨区域有能力，跨样本终点需细分 |
| [TideCluster](https://github.com/kavonrtep/TideCluster) | assembly ✓，raw-read 去新 — | TRC/superfamily、逐阵列 KITE/HOR 线索 ✓ | 组装覆盖量 ≠ read CN | ✓ | — | 比较脚本与 HTML ✓；结构置信标签是操作规则 |
| [CENdetectHOR](https://github.com/CENdetectHOR/CENdetectHOR_pipeline) | assembly FASTA ✓，raw FASTQ — | 单体、家族、HOR、变体 ✓ | — | ✓ | — | GUI/图谱 ✓；手工筛选成本须记录 |
| [HiCAT-human](https://github.com/865699871/HiCAT-human) | 人类 HiFi reads + assembly ✓ | 人类 HOR ✓ | 不能预设为校准植物 CN | ✓ | reads/assembly 注释相连 △ | 多人样本 ✓；植物迁移 ? |
| [Merqury](https://github.com/marbl/merqury) | k-mer 光谱，非单体去新 | — | 全局 k-mer multiplicity ✓ | 全局/局部 copy 谱 △ | 非自动家族标签/精确塌缩分类 | 图谱 ✓ |

**不可强行合并的评分：**assembly-only 工具在 raw-read 去新记 `N/A`；普通映射拿到 TandemX 目录后发现召回记 `N/A`；TideHunter 读内 copyNum 不当作基因组 CN；SRF 的 mapped bp、TandemX 的 diagnostic bp 与物理缺失 bp 分列；CENdetectHOR 的 HOR 真值评分不与现有 `family_hierarchy.tsv` 的候选整数倍边比较；超时、输入不合格、存储失败、合法空输出分列状态。所有工具保留默认与同等预算调参轨道，并计完整依赖流水线、人工选择、准备和临时盘成本。

为避免上表合并列遮蔽任务差别，正式比较协议逐项登记以下 **12 个任务维度**。其中 `N/A` 必须在运行前依据官方输入/输出定义判定，不能事后因结果不好而宣布不适用。

| 维度 | 同终点可比的主要路径 | 单独报告的界限 |
| --- | --- | --- |
| 读段 discovery | TandemX、SRF；TRF/TideHunter 的逐读段原生单位发现 | TRF/TideHunter 如需额外跨读段合并，完整流水线才可评家族召回；assembly-only 为 N/A |
| monomer identification | TRF、TideHunter、SRF；assembly 轨加 CENdetectHOR、TideCluster、TRASH、HiCAT | 真实单位长度及谐波需要等价规则 |
| family clustering | TandemX 与预注册的 SRF 重复单位目录→跨序列家族合并适配器；assembly 轨单列 | [SRF 官方说明](https://github.com/lh3/srf)提醒同类可输出多个相近序列；其原生单位目录与额外聚类结果分开评分并计适配器成本 |
| HOR detection | CENdetectHOR、TideCluster/KITE、HiCAT、TRASH 相关版本；TandemX 只有候选层级 | 当前 TandemX 无逐拷贝 HOR 分数 |
| abundance estimation | TandemX、普通竞争映射、SRF mapped-bp 完整流程 | 分别报告 diagnostic bp、assigned bp、物理 bp；不可混为一列 |
| assembly localization | TandemX 与 assembly 工具共同 FASTA 或同目录映射轨 | 坐标、边界与阵列身份分开 |
| assembly collapse detection | TandemX 与同目录映射/局部审计适用流水线 | assembly-only HOR 注释在读段差额终点为 N/A |
| copy-number estimation | 以独立真值评家族总 CN；Merqury 为 k-mer 光谱参照 | TideHunter 读内 copyNum、组装拷贝数和物理家族 CN 不可等同 |
| cross-sample comparison | TandemX cohort、TideCluster 比较脚本、适用的 HiCAT-human 轨 | 对比输入、家族身份统一规则和材料配对 |
| large-genome scalability | 同机、同输入、同预算的完整适用流水线 | 文献最大样本规模只证明作者运行过，不能代替本机资源曲线 |
| uncertainty estimation | TandemX 状态/词离散、TideCluster 结构强弱标签及可输出 CI 的新方法 | 规则性置信标签与有覆盖率校准的区间分开 |
| visualization/reproducibility | TandemX 报告、TideCluster HTML、CENdetectHOR GUI 等各自原生输出 | 可编辑图、自动重跑、版本/人工整理成本分别计量 |

**CENdetectHOR 专项共同任务：**在其公开 human 与 Arabidopsis 用例及封存植物组装上，固定 assembly FASTA、可选 consensus 的供给轨、版本/commit、资源预算与同一人工整理规则，分开计 repeat-region detection、monomer size、monomer decomposition、HOR composition、HOR length、variant HOR、wall/CPU/RSS/临时盘。**共同准确率排名只允许两方获得同一 assembly 与同等先验；若 TandemX 须先用原始 reads 生成目录，该轨道标成 `reads+assembly` 的扩展输入示范，不能与纯 assembly CENdetectHOR 构成算法优劣排名。**公开用例是复现/阳性对照，若用于方法选择就不再是独立留出。单染色体/着丝粒裁剪和 TandemX 对应输入准备均计时保存；如需人工 curated HOR，原始输出和人工改动分开评分。读段丰度与差额任务对 CENdetectHOR 列 `N/A`。

## 2. 十二字段的可执行差距登记册

每条均用同一字段，避免只列愿望清单。下述数值门槛是**第一轮建议**，需先在开发集测量基准变异和统计功效、形成冻结协议与源码/输入哈希，**在任何新留出结果揭盲前**固定；若不能先固定，就不能将后续结果叫独立验证。

### G1：校准的家族丰度与差额幅度

1. **Current claim：**诊断 k-mer 点估计与 assembly 表示量能够指出可能家族差额。
2. **Current evidence：**模拟 MARE 0.363；Macadamia 方向二元分类 2/0/0/41，但读段差额 2.119 Mb 对组装增量 0.159 Mb；`TXF000695` 跨 k 剧烈反转。
3. **Weakness：**词仅在当前目录内排他；同源背景、单体变异、错误、生物/技术覆盖偏倚可使中位数失真；10–90% 词深度离散不是置信区间；新组装增量不是物理真值。
4. **Required experiment：**在完整背景、近缘家族、无家族 decoy、不同 GC/覆盖/读长/错误/平台的精确真值中，分别测已知目录和端到端未知目录；独立物理材料只作另一层验证。
5. **Required algorithm change：**研究原型中比较 (i) 有辨识性约束的共享特征混合分解与 (ii) 带 `unknown` 拒绝归属的逐读段竞争性周期占据。允许输出 family group/未解决；95% 区间从独立 read/molecule 分层重抽样得出。只有辨识矩阵有足够秩/独特锚、误差/曝光可估时才输出单家族精值。NNLS、Poisson/NB、EM/Bayes 是模型备选与消融，不自动构成新意。
6. **Comparator：**冻结的 TandemX A0；使用相同家族目录的普通竞争性映射 A1；SRF 本地完整流水线的原生 mapped-bp 单独报告，不偷换测量定义。
7. **Metric：**真实丰度大于零的家族报告配对 `|log2(估计/真实)|`、MARE 与偏差；真实零丰度/decoy 另报假阳性赋值、分配 bp 和拒判，不对零取对数或除零。另报符号反转、负读误归属、未分配率、95% 区间经验覆盖/宽度和计算资源；完整分母及家族、供体、基因组层并列。
8. **Ground truth：**合成来源坐标与实体 bp/完整 copies；真实 T2T 经精确编辑的**编辑量**真值；真实生物绝对量仅在独立可跨阵列分子或 ddPCR/qPCR 等验证足够时作校准真值。
9. **Development or held-out：**先按物种/供体/相关家族组整体分离，旧 Ey15/Macadamia 和已观察种子只作开发/压力测试；全新物种/家族组一次性留出。
10. **Pass/fail criterion：**可辨识性必须在运行前仅依输入设计矩阵/预置序列锚与冻结容差定义，不能看误差后筛选。所有预注册家族进入主分母，主指标同时惩罚错误值和拒判；预定可辨识层另报配对主误差相对**开发集选定且锁定的主基线**下降 ≥15% 这一候选门槛，并报告不可辨识层的组级结果、拒判率及覆盖。多基线若都作正式比较，须预先设同时推断/多重比较规则。95% 区间的经验覆盖、宽度和负读误归属/FPR 预先定非劣界；样本量须先按独立模拟基因组/供体簇与最小有意义差异做功效计算，足够时才用预注册的聚类区间作正式通过判决。若只有 2–4 独立材料，只作描述性 feasibility，绝不凭 family 数量生成窄 CI 或 Methods 优势结论。
11. **Estimated compute/data requirement：**先 10-Mb/100-Mb 梯度，随后至少两个未用材料的 20–30× HiFi 与独立平台；每 1-Gb 基因组约 20–30 Gb/read platform 基础读段量，临时文件/索引需另行预估；T7 稳定性是先决门。
12. **Publication relevance：**若对强映射基线有独立、可解释的校准增益，才可能成为第一项 Methods 核心贡献；仅优化 Macadamia 的代理增量数值不算成功。

### G2：精确受控塌缩真值及完整分母

1. **Current claim：**旧/新组装二元代理评价显示有目标发现能力。
2. **Current evidence：**Ey15-2 8/0/0/11、Macadamia 2/0/0/41，主入选要求新组装至少 15 kb。
3. **Weakness：**完全缺席的新组装家族被排除；新组装真值与 HiFi 来源相关；阈值变为 5 kb 后出现 FN/FP；没有精确减失 bp 真值。
4. **Required experiment：**先核验未编辑 T2T 与原始 reads 的同供体/单倍型和局部一致性，并标记本来已疑似塌缩或不可辨识的区域。保持真实 reads 不变，在符合基线资格的组装副本上做 0/25/50/75/100% 单位收缩/删除；记录编辑前后完整单位、残单位、断点、序列哈希和双向坐标链。设未编辑、非重复等长删除、扩张、多位点同源与错配供体阴性；另做不同覆盖及组装器的自然塌缩实验，分表报告。
5. **Required algorithm change：**本项主工作是独立 benchmark 生成/评分器，**不是** TandemX 新推断算法；须让主评分器包含 100% 删除而不使用事后新组装 15-kb 过滤。工具输出未判定及无映射要有稳定 schema。
6. **Comparator：**TandemX A0/候选推断；同目录 ordinary mapping + localization；适用时 TandemTools/长读局部审计；SRF 以合法本机完整流程配对。assembly-only HOR 工具只参与它们支持的结构/定位端点。
7. **Metric：**按 family、locus、editing event 分开计算 sensitivity/FPR/PPV、AUPRC、缺失 bp MAE/偏差、阵列定位 IoU/边界误差、拒判率、全缺席召回；未运行/失效/不可辨识不填 0。
8. **Ground truth：**精确序列编辑账本仅为**编辑前后差值/注入缺失**真值；即使基线读段一致，未编辑 T2T 仍是参考代理，不自动等于生物绝对 CN。若未编辑阵列与供体读段没有独立锚定长分子或等效证据验收，可评价“能否发现注入编辑”，不能把原始读段与编辑组装的差额误差称为物理缺失量误差。基线合格时对 `read–assembly deficit` 的编辑诱发增量可用配对设计评价；自然组装塌缩需独立 flank-anchored 比对到冻结真值才可判断。
9. **Development or held-out：**按完整材料/单倍型/物种分离，不能同一材料随机分阵列；编辑参数先在开发材料冻结，然后对封存材料一次性实施与评分。
10. **Pass/fail criterion：**必须覆盖 0–100% 梯度及全部预定负例，零/完全缺席家族保留在分母；主分类检验同一完整预定 family/locus 分母上的配对 **ΔAUPRC = TandemX − 开发集锁定的主基线**，FPR 和注入缺失增量 MAE 另立门，不把绝对 AUPRC 的 CI 下界 >0 当成绩。最小有意义 Δ、FPR 容差、独立材料数与按供体/基因组聚类的区间方法均由开发集功效分析在留出前锁定；多基线需同时推断。独立簇不足则只能报告可行性，不能用编辑事件数冒充独立样本作显著性结论。缺失量门未过仅保留二元审计；来源/真值不完整为技术阻断，不能把阳性数量写作零。
11. **Estimated compute/data requirement：**10 Mb 合成加约 100–150 Mb 真 T2T 试验；继而 ≥1 Gb 材料及独立留出。每材料真实 reads 的 20–30× 输入与成对组装、编辑版本、索引、日志；多组装与 1-Gb 级需数十至数百 GB 暂存，开跑前按实际文件/磁盘证明。
12. **Publication relevance：**是识别方法是否真能测塌缩的必要验收基础，可作为共享资源/benchmark 贡献；仅制造真值不等同新的 TandemX 算法。

### G3：逐拷贝层级/HOR 架构与读段互证

1. **Current claim：**`family_hierarchy.tsv` 保存相关家族和可能整数周期倍数。
2. **Current evidence：**代码是家族代表序列成对边，测试含 171/342/684-bp 关系；当前无逐拷贝顺序或跨分子支持。
3. **Weakness：**仅整数长度比/局部同源可由谐波、共享片段或部分单体造成；不能命名 HOR、变体、嵌套结构。
4. **Required experiment：**精确 AB/ABC/ABAC、变体、倒向、嵌套、谐波、同质单体、截断与结构塌缩模拟；真 T2T 单体变体的半合成标签串；human/Arabidopsis 与 CENdetectHOR 作者公开输入上同 assembly 任务；新植物材料读段/组装不一致单列。
5. **Required algorithm change：**研究原型先逐拷贝 indel-aware 分解并保留相位/方向/次优，再将独立分子的有序标签串建稀疏支持图；以最短循环单位和邻接/held-out 标签预测检验局部周期，输出变体、等价模型及拒判状态。不能把旧 `family_hierarchy.tsv` 重新标名为 HOR。
6. **Comparator：**现行 CENdetectHOR 的 assembly FASTA 路径（无先验与同源共识供给两轨）、TideCluster/KITE、HiCAT/适用的 HORmon 或 StringDecomposer 任务轨，只在共同 assembly 输入的 HOR 任务比较。read-supported discordance 的同输入基线需给现有结构工具输出再加常规原始读段比对/局部检查，或选择可接收相同 reads+assembly 的工具；纯 assembly-only 输出在此端点记 N/A。
7. **Metric：**copy 边界/标签聚类、最短循环+RC 单位、HOR/变体/嵌套 precision/recall、邻接边、拒判率、read/assembly 不一致检测、时间/RSS；人工 GUI 与依赖准备时间另计。
8. **Ground truth：**预置 copy/order/坐标的 exact 和半合成；公开 human/Arabidopsis 注释为参考标签与复现实验，不等同完整绝对真值；真实 read-spanning 证据独立标注。
9. **Development or held-out：**同一单体家族的所有谐波/变体、同供体 reads/assembly 必须同 split；新物种/结构类一次性留出。CHM13/Col-CEN 若曾被方法选择看过，作为比较 positive control，不能伪称最终植物留出。
10. **Pass/fail criterion：**assembly-only 共同任务须两方使用相同 assembly 和先验，才能以结构 precision/recall、边界、变体及完好组装回退与 CENdetectHOR 等公平比；如果 TandemX 的目录必须由 reads 生成，该结构比较只能做扩展输入示范，不能形成直接优劣判断。read-supported 错误结构检测只对拥有**相同 reads+assembly 输入**的常规映射审计/合适工具比，并做移除 reads 的自家消融。HOR precision ≥0.95、recall ≥0.80、non-HOR FPR ≤0.05 是待开发集验证的候选目标；功效、独立材料数、最小改进和等价容差均须揭盲前固定，不足时只能报可行性。全部拒判与不可辨识结构保留分母；若无同输入可辨识增益，停止 HOR Methods 主张。
11. **Estimated compute/data requirement：**先 10-Mb 可控结构、100-Mb 半合成；之后至少一个现行 CENdetectHOR 公共 human/Arabidopsis positive control 和 ≥2 个未用植物结构材料；按 array 与分子分块，限制周期枚举及候选图大小。
12. **Publication relevance：**若实现 read 侧对 assembly 错误结构的独立增益，可成为第二项核心贡献；单在 assembly 注释终点与 CENdetectHOR 持平只能称竞争性，不是优势。

### G4：HiFi/ONT 适用域和生物学前瞻证据

1. **Current claim：**当前主发现基于 HiFi，Macadamia ONT 给个别家族方向级交叉检查。
2. **Current evidence：**十完整 HiFi 文库/八报告物种共 262.731 Gb 完整输入 QC；六预选家族的 3 支持/3 未解决；YSD56 785-bp 候选有套嵌深度、组装和 TRF/TideHunter 支持，但非全球新家族/独立重复。
3. **Weakness：**没有 ONT 原始读段去新+定量的一般误差模型与真值测试；十物种可运行不是生物验证；旧候选已被查看，不能回收成新前瞻留出。
4. **Required experiment：**真正前瞻盲测优先新获得、未参与开发的至少 2 物种×2 独立供体植株；每植株逐级记录采样、组织、提取、文库、平台、单倍型和组装版本，并预留独立提取。各供体 HiFi、ONT 超长读段、PCR-free WGS 与锚定分子/定量验证来源必须审计。独立保管者先冻结已知家族、全目录候选、阴性、匹配规则、阈值、统计单位与全部失败状态；TandemX 先以读段冻结预测，再揭盲。已知家族恢复、读段差额方向、物理阵列长度/家族 CN 是三个独立终点，不能互相代替。
5. **Required algorithm change：**仅当 G1/G3 留出通过后，限定加入经验证的平台误差模型/证据融合；不为了个别家族事后调阈。生物实验本身不要求新软件功能。
6. **Comparator：**适用读段去新工具 SRF、TideHunter；同目录普通映射；组装结构工具只对重叠端点比较。
7. **Metric：**前瞻候选验证比例、FDR/阴性表现、方向一致性、物理 CN/阵列长度误差、不可判定比例；以独立材料为统计单位，报告每家族证据。
8. **Ground truth：**双侧单拷贝锚定的超长 ONT/光学图谱或 fiber-FISH 可对特定位点的长度提供独立层级；ddPCR/qPCR 可支持家族总量，普通 FISH/CENH3 主要支持定位/身份而非精确缺失 bp。错位、污染、提取差异或供体不匹配未排除时记 `unresolved`；新 T2T 仍是组装代理，不独立证明物理 copy。
9. **Development or held-out：**Ey15-2、Macadamia、YSD56 与已公开查看过的候选为开发/上下文；真正 prospective 材料需新建来源和盲化清单，分子与物种级留出。
10. **Pass/fail criterion：**先导建议至少 20 个已知 family×donor 描述单元、预注册正负层及阳性/匹配阴性各至少 3 个可物理核验位点，但 2 物种×2 供体最多只有 4 个生物簇，**不足以估计总体 FDR、零误认率或强方法学优势**。正式准确率/恢复率的独立供体和物种数须据开发期簇间差异及目标精度做功效计算后冻结；不足时只给完整个案、分母、失败和描述性区间。前瞻生物发现需全部候选预测先锁定，至少一个预定位置/方向假说获独立物理证据才可称该位点获支持，不能推广到所有家族。物理阳性须锚定长度区间高于组装、阴性与组装相容；否则该位点失败或未解决。ONT 去新需另过精确模拟+真实真值，未过时仅作交叉平台核查。
11. **Estimated compute/data requirement：**按 0.5-Gb 基因组×4 供体×HiFi/ONT/PCR-free WGS 各 30×，原始碱基约 180 Gb，另需约 0.5–1 TB 或按真实压缩率/中间文件重算的工作空间；测序及 FISH/ddPCR 金额须当地报价。可考虑 [Medicago A17/R108 公开项目](https://ngdc.cncb.ac.cn/gsa/browse/CRA027262) 作流程先导，但该项目已被仓库列为待入组，单株/提取链和标准 PCR-free WGS 资格尚不闭合，不能直接叫严格外部留出。T7 完整性门通过前不下载。
12. **Publication relevance：**新方法是否解决真实问题的必要示范；不宜以十文库运行数量代替。

### G5：生产级效率、恢复与可复现性

1. **Current claim：**流式扫描、Rust 并行和阶段级指纹复用可处理真实子集/部分 Gb 输入。
2. **Current evidence：**Mo17 1.129 Gb 与 Morex 1.170 Gb 完成诊断；完整 Ey15/Macadamia 定量分别约 9,495/11,293 s；T7 掉盘、K30076 损坏、历史比较超时留下明确失败状态。
3. **Weakness：**无 `discover`/`quantify` 阶段内部 checkpoint；`tandemx run` 只在已完成步骤按输出和 SHA-256 指纹复用，中断的单步重跑。T7 重挂不证明文件完整；9/14 K30076 第二块哈希变化、gzip 失败、随后外盘 EIO/卸载，块 3 未开始。现有 profiler 可采样存活进程树 RSS/进程数、scratch 和累计 CPU，但缺 I/O 字节/速率与逐子进程完整高水位；`wait4` 直属子进程 RSS 不等于整树峰值。尚无独立 10 Gb 与完整 wheat/barley/rye 的全流程重复资源曲线。
4. **Required experiment：**先对源文件作独立 SHA-256、gzip/FASTQ 语义和来源检查，再用新目录有限写入—读回证明外盘稳定。固定机器、OS、CPU/内存、盘位、线程、构建、命令、输入哈希/条数/bp；停止背景重任务，记录并交错随机化工具顺序和冷/暖缓存轨。用固定完整读的嵌套输入在 ~100 Mb、~1 Gb、~10 Gb 各至少三次新目录隔离重复，完整谷物若不能三次则单次只作为可行性 smoke。分别报告发现、聚类、定量、定位/比较、报告和端到端。对阶段约 30%/70% 注入 SIGKILL、EIO、盘满、输入/收据损坏，核对恢复前后规范化稳定产品哈希、重算量和临时盘清理。
5. **Required algorithm change：**性能原型中添加内部原子分块检查点、输入/参数版本指纹、临时文件上限与可重入输出；仅在科学输出逐字节一致及故障注入通过后并入生产。
6. **Comparator：**同输入同线程/整机预算下的 TideHunter、TRF、SRF 完整流水线及适用 assembly 工具；只比较共同任务，流程外模块分表。
7. **Metric：**wall/user/system CPU、存活进程树合计 RSS 高水位及采样精度、PID/线程数、读写字节与 IOPS、压缩/解压 CPU、峰值/最终 temp、bp/s 与 reads/s、1/2/4/8 线程扩展、技术失败/超时/信号、恢复重做比例、规范化稳定产品哈希；三次技术重复不是三个生物样本。不可把各 PID 独立 RSS 峰值求和或最大单一 PID 作为进程树峰值。
8. **Ground truth：**资源仪器与文件哈希；科学准确度另按 G1–G4 真值评分。
9. **Development or held-out：**资源调优使用开发输入；最终冻结机器与命令后用未用于优化的真实大输入复测。不得把不同主机或并发任务拼成 scaling 曲线。
10. **Pass/fail criterion：**至少三次隔离重复（预先声明极大全量运行的一次可行性例外）；四级输入运行状态和全过程资源完整。研究门建议 RAM 峰值 ≤可用物理内存 70%、scratch 高水位 ≤预留配额 70%、运行前可用空间 ≥预测 scratch 峰值 1.5 倍；这些是容量规划待硬件实测后冻结，绝非已通过阈值。增线程使中位 wall 变差 >10% 或稳定输出改变判失败。内部检查点需只重算未提交块，原子收据含输入指纹、计数及哈希；SIGKILL/盘满/EIO/损坏收据必须非零退出、保留失败证据、不接受半成品，恢复结果哈希等于清洁运行。未通过不宣称生产级。
11. **Estimated compute/data requirement：**现有单次 Morex 1.17 Gb/1,909.5 s 若假定同类读分布与完全线性，10 Gb discover 约 4.5 小时，300/480 Gb 的 10/16-Gb 基因组 30× 全库约 5.7/9.1 天；这只是容量估算，未包括全读 quantify、组装扫描、I/O 或非线性开销，不能作为实测速度。10× 发现预算也不是生物饱和阈值。三工具/三重复可能需多 TB 工作空间与多日机时；按真实压缩字节和 scratch 乘数预检。
12. **Publication relevance：**是可复现软件发布的必要工程门，不能代替科学新意。

## 3. 只保留两个核心方法候选，以及一个共享真值框架

**候选 M1，优先级 1：可辨识且校准的家族丰度推断。** 总假说是读段级竞争归属/拒判或共享特征的可识别分解能够在同源背景下改善家族丰度及区间覆盖。先让两条小型研究原型在同一开发协议中与 A0、A1 对打；只保留胜出且可解释的一条。不能让不可辨识家族靠正则化或先验变出窄区间，也不能用 Macadamia 新旧代理的 158,602 bp 作训练标签。M1 若未超过普通映射或拒判过高则停止。

**候选 M2，优先级 2，条件性：读段支持的局部单体/HOR 结构。** 只在自身逐拷贝结构真值和读段/组装来源能支撑关联、且开发试验显示共同 assembly HOR 任务至少具竞争力后推进；首要独特终点是对错误组装的 read-supported 结构矛盾。该端点须与相同 reads+assembly 输入的映射审计基线竞争，assembly-only 结构工具只在共同结构任务比较。M2 的 exact-HOR 试验不以 M1 丰度模型通过为前提。若无多个独立材料的结构证据，停在候选关系表。

**共享框架 B1：精确受控 array contraction/deletion。** 属于评估与数据资源贡献，不算第三个新算法。必须优先建立，以避免在旧/新参考代理上继续评价自己。只有在生成器、编辑账本、评分器和独立再现均通过后，它才能作为未来 `audit-benchmark` 命令的候选；第一阶段不新增公共命令。

现阶段明确不推进：A3 第三轮、把已有 median 加置信标签改名成概率模型、无消融的 NNLS/EM 堆叠、仅为论文加入 AI、`unitFinder` 重试、为优胜结果临时改阈值、先做大量 GUI/图形改造。

## 4. 下一阶段实验序列、拆分与冻结规则

1. **P0：来源和同终点协议。** 保存目前生产源码与失败实验；建立 2026 竞品版本表，CENdetectHOR pipeline 和 library 分别固定 commit/0.0.6，TideCluster 现行版与历史本地 1.21.2 分开。先核实新 T2T/HiFi/ONT 的同供体、抽提、单倍型、相位及校验和；T7 仍须读写稳定性与每文件哈希门，未通过不下载/跑大计算。
2. **P1：真值与评分先行。** 独立实现 B1 的编辑账本和盲化评分协议（本轮仅设计）。开发材料先 10 Mb exact、再约 100–150 Mb 半合成；拆分单位是整物种/供体/家族同源簇。预注册零表示、等长非重复删除、扩张、donor-swap 等负控。签署主终点/最小有意义差异和 power 之后封存至少两个植物物种的 held-out，绝不重复揭盲。
3. **P2：M1 有界原型对打。** 固定目录条件先比较共享特征分解、逐读段竞争/拒判、A0 和普通映射；再端到端未知目录。每项有独立 branch/worktree、假说、冻结输入/参数、before/after 表、消融、聚焦回归、短决策说明。开发集失败即停止分支，失败和 timeout 进入 `docs/negative_results/`。
4. **P3：M2 可行性门。** 不先做完整 HOR 产品。先在 exact 结构与 CENdetectHOR human/Arabidopsis positive control 证明逐拷贝/顺序评分可复现；共同 assembly HOR 终点与现行工具平等比较。再测 read-supported 组装结构差，使用同输入的常规读段映射审计基线和去 reads 消融；若无增益或独立分子不足，停止 M2。
5. **P4：一次性留出与真实竞争者。** 将完成的 M1（及通过门槛的 M2）冻结后，对新植物材料一次性运行全部已宣布工具、参数、资源上限与分母；输出包含 `ok`、合法空目录、技术失败、超时、存储失败、拒判、N/A。human/Col-CEN 结构比较使用与 CENdetectHOR 文献可核数据；它们不替代未使用植物留出。
6. **P5：自然组装塌缩与前瞻生物实验。** 对不同覆盖/组装器的 assemblies 做独立真值检查，与人工编辑结果分层；新材料先锁**全部**家族/坐标/方向、匹配规则和阴性清单，再用独立长分子、PCR-free WGS、ddPCR/qPCR 或 FISH 盲化核验。所有预注册预测及失败进入分母；独立供体才是生物重复，技术平台或池样本不能冒充。失败不删。
7. **P6：隔离性能、可恢复发布和论文重构。** 通过科学门后补齐 ~100 Mb→1 Gb→10 Gb→完整谷物资源/故障实验；清洁安装、容器、完整输出定义与正式 release 分别审核。仅此后按“方法→精确真值→当前竞品→校准审计→植物发现→扩展性”重写论文，不维护 v5 原 Result 顺序的惯性。

**总体停止门：**若 B1 编辑账本不能独立复现，停止**受控组装缺失量**主张，但 M1 的独立合成丰度真值和 M2 的 exact copy/order 真值仍可各自评价；若 M1 不胜过强映射基线或无法校准，停止“准确丰度/缺失量” Methods 主张；若 M2 在共同结构终点不具竞争力且同输入读段审计无增益，停止 HOR 主张；若 prospective 假说无独立支持，稿件只报告失败/有限审计；若 full cereal 与恢复实验不通过，不宣称生产级。任何一项失败都不因既往投入而改写分母或隐藏结果。

## 5. 后续发布前需要重新确认的独立事项

- 将 `docs/comparator_matrix.md` 和 README 中过时的 competitor/版本与 MVP 口径更新，但**本轮先保留原样作为历史档案**。
- `docs/current_status.md` 的旧时间段包含已被 9 月 14 日检查点覆盖的“正在运行”；新实验不得引用被覆盖状态。
- 可复现 benchmark 必须保存完整 input/source/executable hash、配置、机器/线程、全部输出及空/失败/拒判行；模型成功才推进生产代码与文稿。
- 该审计没有实际安装或运行 2026 CENdetectHOR、TideCluster 新版或 HiCAT-human；文献/README 能力与本地同机性能必须分开。

## 6. 独立严苛审稿后的修订记录

最后一轮独立审稿专门检查本方案，而非替现有结果背书。它指出：原草案把 2–4 个供体和 20 个 family×donor 混为足够统计样本；容许“可辨识”家族或“最佳”基线在留出结果后选择；在读段+组装结构矛盾上直接和仅有组装输入的工具比较；把注入编辑量偷换成物理拷贝真值；以 B1 失败错误地阻断 M1/M2 独立真值实验；把 SRF 原生重复单位目录视作现成家族聚类。这些均已在 G1–G4、专项竞品协议及总体停止门中改为预先锁定的分母/基线、按独立材料做功效设计、同输入比较和分层真值。**这些修订是方案修正，尚无新实验结果，因此没有一条因此变成“已通过”。**

## 7. 审计覆盖与尚未验证的范围

本轮核对了仓库规则与[当前状态总账](current_status.md)、[9 月 16 日专家现状盘点](tandemx_expert_status_20260916_zh.md)、生产发现/定量/层级/管线实现、相关单元与集成测试、[性能说明](performance.md)、已冻结的模拟/SRF/旧新组装/正交读段证据摘要及[稿件 v5](../paper/0910/manuscript_v5.md)。文献侧核对官方 CENdetectHOR、SRF、TideCluster、HiCAT-human 等项目与相应方法原文，并把官方能力与尚未本地复现的能力区分。仓库中还有大量历史运行收据、原始比较数据和图稿；本轮**没有逐文件复验全部仓库内容或 T7 上的大文件**。2026 竞品没有安装在同一机器重跑，新增真值材料没有下载，资源预测没有取代测量；这些是下一阶段明确的未完成工作。
