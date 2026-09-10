# 豆科下一阶段数据入组决策表（2026-09-10）

本表以已核实的逐样本映射为准，详见 [豆科 T2T/近 T2T 数据清单](legume_t2t_dataset_inventory_20260910.md)。入组指允许进入 TandemX 的小规模 read-sampling 和 assembly-aware 试跑，不等于立即下载完整 FASTQ/BAM。排序优先考虑：

1. assembly、raw reads 与材料的可追溯匹配；
2. 倍性和 haplotype 复杂度；
3. 是否有论文支持的 T2T/端粒/着丝粒和已知 tandem-repeat 背景；
4. raw 数据量和预估计算成本。

原始 bases 与压缩文件 bytes 分开记录；`估算` 不应写成仓库已报告值。`T2T`、`T2T-claimed/partial` 和 `near-T2T` 按原清单定义使用。

## 结论：按当前“寻找未报道 TR”的目标，YSD56 之后先进入花生二倍体

YSD56 作为野大豆 T2T 流程基线保留。若问题是尽快发现具有跨材料价值的未报道 TR，下一主线应按用户指定的花生资源推进 **V14167（AA）→ K30076（BB）→ 一个四倍体花生**。两份二倍体祖先基因组先建立 A/B 亚基因组的已知与候选 repeat 基线，再进入四倍体，才能区分已有祖先家族、亚基因组特异变体和四倍化后出现的候选结构。首个四倍体科学优先选 S245；若受下载与计算预算限制，先选 HN51。

ZH13 仍是最合适的大豆阳性对照，而不是当前 novel-discovery 主线的下一材料：它是二倍体栽培大豆，已有 CentGm91/92/273/413/444 的明确背景，适合检查 YSD56 结果是否能复现已知家族和筛查规则，但对“未报道家族”的独立信息增益低于花生。A17 是首个 Medicago T2T 材料；苜蓿 ZM4 为四倍体且仍有 8 个 gaps，应在 A17 和至少一个四倍体花生之后进入。

当前推荐顺序为：**完成 YSD56 10.96x 候选复现 → V14167 小样本 → K30076 小样本 → S245（或成本优先的 HN51）→ A17 → ZM4**。ZH13、HJD 可作为阳性/低成本流程控制穿插，但不替代花生主线。

## 分阶段入组表

| 阶段/顺序 | 材料（物种、倍性） | donor/材料匹配证据 | T2T 与 TR 背景 | 已核实 raw 规模 | 预计计算成本 | 决策 |
|---|---|---|---|---|---|---|
| 基线 0 | YSD56（*Glycine soja*，2x） | assembly `GCA_040083835.1` 与 `PRJNA1095640/SRP502474` 同为 YSD56，材料级匹配明确 | 论文报告 0 gap、20 染色体/40 端粒；TRF 标签 `trf91/92/182/183/184/273/276`；Hi-C 量有 23.84/123.84 Gb 冲突 | Illumina 91.05 + HiFi 44.19 + ONT 240.85 = **376.09 Gb** core bases；计 Hi-C 则约 **399.93 或 499.93 Gb**，Hi-C 待复核 | 高；已作为基线，不重复下载 | 保留为参考和回归样本 |
| 阳性对照 | ZH13（*Glycine max*，2x） | `PRJCA015269/CRA010060`、GWH `GWHBWDJ00000000.1`；BioSample 标为 Zhonghuang 13/monoisolate；library 级链仍需入组前复核 | 论文报告 39/40 端粒、全部着丝粒，`T2T-claimed/archive-complete`；CentGm91/92/273/413/444 背景 | 三类核心 DNA raw 约 238.7 Gb bases（HiFi 约 67×、ONT 约 117×、BGI 约 53×，覆盖度换算） | 中 | **大豆阳性对照首选**；同属大豆、二倍体、重复背景最清楚，但不取代花生 novel-discovery 主线 |
| 低成本流程控制 | HJD（*Vigna unguiculata*，2x） | `PRJCA044301` 与 HJD assembly `GWHHKIK00000000.1` 在同一研究；具体 DNA 提取链入组前复核 | 论文报告 11 条 gap-free 染色体/22 端粒；当前清单未闭合命名 satellite 家族 | HiFi 37.33 + ONT 61.01 + Hi-C 32.12 = 130.46 Gb bases | 低至中 | **首选跨物种低成本控制**；若先验证流程吞吐量，可排在 ZH13 前 |
| Medicago 第一材料 | A17（*Medicago truncatula*，2x） | `PRJCA042220/CRA027262`；assembly `GWHGEXC00000000.1` 与 A17 BioSample/raw 在同一项目；单株 DNA 链未额外公开 | 论文报告 complete/gap-free T2T；当前清单未闭合独立 satellite accession | HiFi 估算 39.86 + ONT 67.72 Gb；项目总量 246.99 GB（含 R108 等数据） | 中 | **建议早期入组**；小基因组、T2T、raw 关系清楚，适合跨豆科 repeat discovery |
| 花生 1（下一发现材料） | V14167（*Arachis duranensis*，AA，2x） | `GCA_054824555.1`；BioSample isolate V14167；`SRR33996330/96240/96449` 按材料闭合；论文为 20 株 accession-level pooled material | 论文称六个花生基因组 T2T，V14167 有 20 端粒；重复丰富并有着丝粒/端粒分析，但当前清单未固定命名 TR catalog | HiFi 211.520 + ONT 108.019 + Illumina 131.114 = **450.653 Gb** | 中至高 | **花生首个推荐**；AA 二倍体、装配/三类 raw 映射完整，可先建立 A 亚基因组基线 |
| 花生 2 | K30076（*Arachis ipaensis*，BB，2x） | `GCA_059453495.1`；BioSample isolate K30076；`SRR33996329/6606/96322` 按材料闭合；同样是 pooled material | 论文 T2T、20 端粒；重复丰富，具体命名 TR catalog 尚未在当前清单闭合 | HiFi 175.311 + ONT 105.865 + Illumina 74.345 = **355.521 Gb** | 中至高 | **与 V14167 成对入组**；BB 二倍体可在四倍体前建立 parental repeat 对照，raw 量比 V14167 低 |
| 6 | YP4（*Phaseolus vulgaris*，2x） | `GCA_051622355.1/PRJNA1072282`；论文称山西单株取样，材料级 donor 较强；逐 run modality 尚未完全拆出 | 论文称 T2T，但只有 9/11 条 pseudomolecule、20/22 端粒；应标 `T2T-claimed/partial-T2T-evidence`；无固定命名 TR catalog | HiFi 31.75 + ONT 177.04 + Hi-C 144.79 = **353.58 Gb** core bases | 中至高 | 作为二倍体部分 T2T 对照；排在 HJD/A17 后，避免把部分 T2T 误当 gap-free 真值 |
| 7 | S245（*Arachis hypogaea*，AABB，4x） | `GCA_054824585.1`；BioSample isolate LaiyangSilihong；`SRR33996311/96595/96310` 按材料闭合；20 株 pooled accession-level material | 论文 T2T、40 端粒；异源四倍体重复丰富，可用于 A/B 亚基因组与 abundance deficit；未固定命名 TR catalog | HiFi 316.441 + ONT 198.190 + Illumina 154.720 = **669.351 Gb** | **很高** | 科学优先的首个四倍体花生；需在二倍体 parent pilot 通过后进入 |
| 8 | HN51（*A. hypogaea*，AABB，4x） | `GCA_054824525.1`；Huayu23；`SRR33996516/96327/96182` 按材料闭合；pooled material | 论文 T2T、40 端粒；重复丰富；未固定命名 TR catalog | HiFi 309.594 + ONT 146.376 + Illumina 51.203 = **507.173 Gb** | 高 | **预算受限时四倍体第一选择**；raw 低于 S245，材料名称和 run 映射完整 |
| 9 | HN873（*A. hypogaea*，AABB，4x） | `GCA_054824565.1`；Chedouzi；`SRR33996436/96231/96193` 按材料闭合；pooled material | 论文 T2T、40 端粒；地方品种背景；未固定命名 TR catalog | HiFi 282.991 + ONT 178.740 + Illumina 55.668 = **517.399 Gb** | 高 | 在 HN51/S245 之后扩展；先复用已验证参数 |
| 10 | S83（*A. hypogaea*，AABB，4x） | `GCA_054824515.1`；YunnanRainbow；`SRR33996477/96210/96584`；assembly/long-read 匹配，但 Illumina 当前 reported bases=0 | 论文 T2T、40 端粒；短读段缺失使 abundance/QV 交叉验证受限 | HiFi 316.600 + ONT 237.065 = **553.665 Gb**；Illumina 0/`metadata_blocked` | 高至很高 | 后置；先补齐 short-read payload 或明确不使用 Illumina polishing 结论 |
| 11 | ZM4（*Medicago sativa* cv. Zhongmu-4，4x、四 haplotype） | 新装配与旧 `PRJCA004062/CRA005190`、ONT `PRJCA030790`、新项目 `PRJCA041059` 跨批次；cultivar-level 匹配，个体/donor identity 未完全证明 | 新装配 **near-T2T**，仍有 8 gaps；旧 GWH `GWHBECI00000000` 仅 chromosome-level；重复丰富但当前清单未固定独立 TR catalog | 旧项目报告 267.60 GB；新装配按 34× HiFi + 35× ONT 估算约 216 Gb bases，实际新项目 bytes 待核实 | **高至很高**；四 haplotype 和 gap 处理增加内存/索引成本 | 在一个二倍体和一个四倍体花生通过后入组；适合测试多倍体和 abundance deficit，不适合作为最早 smoke test |
| 12 | R108（*Medicago littoralis*，2x） | `PRJCA042220/CRA027262`；assembly `GWHGEXD00000000.1`、BioSample 与 raw 同一项目；单株 DNA 链未额外公开 | 论文 complete/gap-free T2T；无固定命名 TR catalog | HiFi 估算 47.03 + ONT 40.51 Gb；部分 BAM/Hi-C bytes 已报告，项目总量 246.99 GB（与 A17 合计） | 中 | A17 通过后成对加入，主要用于材料/物种比较，不先于 A17 |
| 13 | FC6（*V. unguiculata*，2x） | assembly `GWHHKIL00000000.1`；raw 来自前一研究，当前论文未闭合原始项目 accession；donor chain `metadata_blocked` | 论文报告 11 条 gap-free 染色体/22 端粒；无固定命名 TR catalog | HiFi 28.59 + ONT 67.27 + Hi-C 24.10 = 119.96 Gb bases，但 accession 尚未闭合 | 低至中（技术上）；证据风险高 | HJD 通过且 raw accession 补齐后再入组，不能用作第一 cowpea 样本 |

## 六个花生材料的内部决策

六个材料全部是同一研究的论文报告 T2T assembly，且 assembly isolate 与 HiFi/ONT/Illumina run 已按材料闭合到 accession 层面；共同限制是论文称每个 accession 种植 20 株并收集幼叶，因此这是 accession-level pooled donor，而非已证明的单株 donor。四倍体 A/B subgenome 也不能拆写成两个独立 donor。

在二倍体 parent 阶段，先 **V14167 → K30076**：V14167 的 AA 祖先身份和较小 1.178-Gb assembly 适合建立 A 基线；K30076 的 BB assembly 1.485 Gb，但三类 raw 合计 355.5 Gb，低于 V14167，可作为第二个 parent。两者通过后再进入四倍体。

四倍体中，**S245 是科学优先、HN51 是成本优先**。S245 的 HiFi/ONT/Illumina 都充分且为栽培四倍体，适合第一个 AABB 试验，但约 669 Gb core raw 成本最高；HN51 约 507 Gb、run 映射同样闭合，若内存或存储受限应先用 HN51。HN873 随后用于地方品种扩展；S83 最后，因为当前 Illumina run reported bases=0，不能承担完整的 read-based abundance 交叉验证。

## 入组闸门和计算策略

| 闸门 | 进入条件 | 失败处理 |
|---|---|---|
| G0 metadata-only | BioSample/assembly isolate、raw modality、release 状态、reported bases/bytes、md5/sha256 均记录；不下载大文件 | 标为 `metadata_blocked`，不进入 full-run |
| G1 小规模 read pilot | 每个材料先抽取少量 HiFi/ONT/Illumina reads；验证 FASTQ/BAM 格式、长度分布、adapter/低质量比例和 donor-consistency | 只保留可追溯的 pilot 结果，不能宣称全数据结论 |
| G2 assembly-aware pilot | 检查 FASTA 总长、N 数、染色体数、端粒/着丝粒报告和 repeat 输入版本；严格保留 `T2T`/`partial`/`near-T2T` 标签 | 装配异常只作为诊断，不替换公开 accession 或论文状态 |
| G3 full enrollment | G1/G2 通过，且存储、索引和峰值内存预算明确；先 long reads，再按需要加入 Illumina；Hi-C/Pore-C 默认不进 discovery/quantify 核心输入 | 回退到 pilot 或只运行 read-first 模块 |

计算预算应按 **raw volume × ploidy/haplotype complexity × 是否需要 assembly indexing** 估计。四倍体花生和四 haplotype ZM4 即使单次文件量相近，也应按更高内存和更多索引/比较任务预算。所有表中成本是相对等级，不是已经实测的运行时间。

## 目前不应提前承诺的结论

- 论文称 T2T 不等于 TandemX 可以把所有卫星阵列视为完整真值；尤其 YP4 只有部分染色体 T2T 证据，ZM4 仍有 8 gaps。
- accession-level donor match 不等于单株 DNA 提取和 library chain 已证明；花生六材料、ZM4、A17/R108 均要保留这一边界。
- 新装配不能自动作为绝对 copy-number truth。TandemX 输出应使用 read–assembly abundance deficit 或 estimated under-representation 等描述。
- 花生四倍体和 ZM4 的 A/B 或四 haplotype 不能在没有相位/同源定位证据时当作独立重复样本。
