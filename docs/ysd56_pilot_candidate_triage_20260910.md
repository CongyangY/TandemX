# YSD56 pilot candidate triage（只读审计）

**日期：** 2026-09-10  
**输入：** `/Volumes/T7/Codex/TandemX/results/legume_ysd56_novel_tr_pilot_v1_20260910/run/discover/`  
**样本：** *Glycine soja* YSD56；assembly `GCA_040083835.1`；HiFi `SRR28726931`；本次使用嵌套的 `sample_003.fastq.gz`。  
**边界：** 只读检查既有 `families.tsv`、`monomers.fa`、`family_similarity.tsv` 和 `family_hierarchy.tsv`。本结果是 pilot candidate triage，任何行均不得称为 novel；不改变算法、不重跑 discovery。

## 输入规模和审计结果

本次 discovery 已处理 65,106 条 reads、1,099,009,791 bp，产生 1,474 个 operational families。`family_similarity.tsv` 有 2,774 个 emitted related rows，其中 2,700 个为非 distinct 关系；`family_hierarchy.tsv` 有 2,700 条边：807 条为 `putative_period_multiple`，1,893 条为 `unresolved_related_or_partial`。运行审计记录了 1,085,601 个理论 family pair，其中 1,080,878 个在 k-mer gate 被剪枝；这说明表中的“没有关系”只能解释为当前 gate 下未发出关系，不能解释为生物学独立。

| 检查 | 观察结果 | 解释 |
|---|---:|---|
| families | 1,474 | 代表性 consensus；`cluster_identity=0.95`，不是已命名 family |
| similarity emitted | 2,774 | 仅保存非 distinct 关系 |
| hierarchy edges | 2,700 | 2,396 forward、304 reverse |
| putative period multiple | 807 | 549 条 `multiple_error=0`，258 条为近整数倍；均带 `heuristic_period_multiple_not_validated_hor` |
| unresolved related/partial | 1,893 | 相关、部分重复或其他同源解释未解决 |
| reverse orientation similarity | 308 | 代表最佳局部比对方向；不等于生物学反向排列 |
| likely redundant pairs | 74 | 来自现有 `redundant_candidate=true` 规则；未执行 collapse |
| exact cyclic+reverse-complement groups | 0 | 对 `monomers.fa` 做全旋转和反向互补 canonical key 后，没有两个 family 共享完全相同的 key |
| low-complexity families | 24 | 应保留记录，但暂停生物学解释 |

输入文件的运行 receipt 还给出：`families.tsv` SHA-256 为 `482d6daa6dcde1bf1869b7ddd85224e1ff34e76b150ed3f1a447e1523c4daf01`，`monomers.fa` 为 `c4af546f6f7648c35ea547510dcaed16fff34bc3755f3b29c85a9a34dd6da7ce`，`family_similarity.tsv` 为 `432c93c8b8160ca5f8e2dfc7e9776941da28223257bf22b7c84834ba5e759d9d`，`family_hierarchy.tsv` 为 `aef0e1e354bd42eafa48ed1d8d7286c427e8144515b9b2373dbee4c969c115d6`。

## 周期倍数观察

周期倍数关系高度集中在少数长度层级。最常见的 shorter lengths 为 92 bp、91 bp 和 48 bp；其中 92-bp shorter family 的 outgoing edges 最多，单个 family 的计数为：TXF000214 43、TXF000013 42、TXF001226 41、TXF000126 33、TXF000823 31、TXF000478 28。91-bp shorter hubs 包括 TXF000075（23）、TXF001189（18）、TXF000112（14）和 TXF000608（12）；48-bp hubs 包括 TXF000006（12）、TXF000052（11）和 TXF001020（10）。这些是“周期关系审计入口”，不是 91/92-bp 已知 soybean repeat 的序列确认，也不是 HOR。

较长端的关系同样聚集在 144、184、192、276/368/552、736/828/920 bp 等长度。示例包括 TXF001167（552 bp，24 个 incoming edges）、TXF000623（144 bp，16）、TXF000344（736 bp，16）、TXF000108（552 bp，15）和 TXF000328（920 bp，15）。应在后续序列级核对中同时保留 shorter unit 和 longer consensus；不能只保留较长者而把 shorter family 当作噪声。

`multiple_error=0` 只说明报告长度是整数倍；即使同时有 `local_identity=1` 和 `local_overlap_fraction_shorter=1`，也可能是重复共识、部分单元或读长/周期表示造成的关系。只有多单元次序、边界和独立长读长或组装证据一致时，才可将其升级为 `HOR_candidate`；本 pilot 没有完成该升级。

## 反向互补、旋转和冗余

`monomers.fa` 的 canonical 检查为每条 consensus 计算所有 cyclic rotations，并对序列及 reverse complement 取最小表示。没有发现 exact cyclic/反向互补重复，说明 1,474 条 consensus 没有完全相同的规范化序列；它不能排除近似重复。

现有 local audit 报告 308 条 `orientation=reverse`，其中 304 条进入 hierarchy，4 条属于 `likely_redundant`。这里的 `orientation` 是最佳局部比对方向，已包含循环相位搜索；所以应记为 `reverse_orientation_related`，不能写成“反向重复家族”。

74 条 `likely_redundant` pair 组成 32 个连通分量（24 个二节点、4 个三节点，以及大小为 4、5、9、10 的各 1 个分量）。最大两个分量分别包括：

- 32–56 bp 的低复杂度簇：TXF000007、TXF000024、TXF000069、TXF000080、TXF000159、TXF000313、TXF000929、TXF001120、TXF001145；
- 32–64 bp 的另一低复杂度簇：TXF000092、TXF000100、TXF000174、TXF000235、TXF000276、TXF000358、TXF000587、TXF000757、TXF000876、TXF000898。

这些 cluster 应作为 redundancy hold，并保留代表选择所需的 support/span/sequence 证据。现有 pilot 没有 collapse，故不能把 component 数量当作最终 family 数量。

## 透明 triage 表

紧凑 TSV 已写入：

`/Volumes/T7/Codex/TandemX/results/legume_ysd56_novel_tr_pilot_v1_20260910/ysd56_pilot_candidate_triage.tsv`

该表共 83 行，按以下可复核规则选取：

1. support reads 和 support span 排名前 15 的 abundance anchors；
2. outgoing `putative_period_multiple` 最多的前 20 个 shorter hubs；
3. incoming `putative_period_multiple` 最多的前 15 个 longer hubs；
4. `likely_redundant` 图中最大的 3 个连通分量的全部成员；
5. reverse-orientation local relation 按 identity、overlap 和 shared-k-mer fraction 排名前 10 对的全部端点。

表中 `period_multiple_as_shorter_edges` 与 `period_multiple_as_longer_edges` 是 hierarchy 中的关系计数，`exact_multiple_incidence` 是该 family 作为 shorter 或 longer 端参与 `multiple_error=0` 边的次数，`reverse_orientation_related_pairs` 和 `likely_redundant_pairs` 来自 emitted similarity rows。`putative_targets_top`、`putative_parents_top` 只列最接近整数倍/局部 identity 较高的最多 8 个端点，用于人工复核，不代表完整邻接表。

`triage_class` 的含义如下：

| 类别 | 进入条件 | 后续动作 |
|---|---|---|
| `abundance_anchor` | support 排名前 15，且当前表没有 emitted relation | 作为高覆盖入口；先做已知 repeat/序列库排除 |
| `period_multiple_hub_shorter` | shorter 端 outgoing edges ≥8 | 检查 91/92/48 等父周期、子周期和替代解释 |
| `period_multiple_hub_longer` | longer 端 incoming edges ≥8 | 检查是否是多条 shorter consensus 的共同上位表示 |
| `period_multiple_review` | 参与至少一条 putative period edge | 保留所有 parent/child，不做 HOR 命名 |
| `reverse_orientation_review` | 至少一条 reverse-orientation relation | 做反向互补和 cyclic phase 的序列复核 |
| `redundancy_hold` | 进入 likely-redundant component | 暂停去重，先保留代表和被替代序列 |
| `low_complexity_hold` | `low_complexity_flag=true` | 暂停生物学解释，单独记录简单重复风险 |

因此，`triage_class` 是工作顺序，不是 novelty、着丝粒功能或 HOR 结论。TSV 中的 family ID、MD5、长度和 support 值均应回指原始输入；不能把 TSV 作为新的发现结果源。

## 后续检查边界

下一步若要把 pilot 候选推进到 Result 6，必须在不改变本次 frozen 输出的前提下：

1. 用完整的 YSD56 10.957x nested sample 做描述性复现；nested reads 不是独立重复；
2. 对候选 consensus 做 exact、反向互补、cyclic rotation、短序列敏感比对及低复杂度敏感性分析；
3. 加入 soybean/豆科已知库（包括 `trf91`、`trf92` 等既有命名）后再更新 `known_or_homologous` 判断；本次 triage 没有做这一步；
4. 将候选定位到 `GCA_040083835.1`，检查 assembly gap、centromere boundary、rDNA、telomere、TE 和 unplaced contig；
5. 对 abundance、assembly span 和平台 concordance 分开报告；不能把 read--assembly deficit 当作 binary collapse，也不能把新组装自动当作 copy-number truth；
6. 只有在序列、连续阵列、位置和独立证据均通过后，才可使用 `previously_unreported_candidate` 或更强的功能性措辞。

在完成上述步骤前，本表行只能保留为 `preliminary_candidate_not_novel`、
`unresolved_related_or_hor` 或相应的 `known_family_variant`/`metadata_blocked`
状态；后续有直接序列证据的行可被排除为已知类别。例如，后续固定 rDNA
screen 已将 `TXF000367` 记为 `rDNA_related_known_repeat`，而 `TXF000708`
仍为 `preliminary_candidate_not_novel`。本 triage 本身没有任何 `novel`
结论，也不能把有限库的无匹配解释为 novelty proof。
