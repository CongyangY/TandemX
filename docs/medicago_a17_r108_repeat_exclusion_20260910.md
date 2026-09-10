# Medicago A17/R108 已报道重复排除清单（2026-09-10）

## 目的与证据边界

这是一个面向 Medicago A17/R108 T2T 候选分析的紧凑、只读排除清单。它只登记论文或权威项目页面明确报道过的重复名称、长度、材料适用性和功能/位置证据。它不把未取回的序列补成共识，也不把“本次没有序列命中”写成 novel。

机器可读文件位于 `benchmarks/inputs/medicago_a17_r108_known_repeats_20260910/`：

- `source_units.tsv`：8 个已报道重复单位及证据/不确定性字段。
- `medicago_a17_r108_known_repeats_v1.fa`：有意为零序列记录的 FASTA 注释文件；没有任何推断序列可供比对。
- `receipt.json`：来源、范围、零大型数据下载和 unresolved 清单。

本审计没有下载 T2T 组装、reads 或其他大型文件，也没有改变 TandemX 算法。由于当前可追溯论文/项目页面没有给出可安全复用的 GenBank/RefSeq 或补充 FASTA accession，序列文件保持空记录是有意的负结果。

## A17/R108 的材料适用性

2025 T2T 研究将 A17 写作 *Medicago truncatula* Jemalong A17，将 R108 写作 *Medicago littoralis* R108；NGDC 项目 PRJCA042220 同时列出两种物种。2004 论文使用的是 *Medicago* accession R108-1 的旧称。因此，MtR1/MtR2/MtR3 的 R108-1 FISH 证据可以作为历史材料证据，但在把它映射到 2025 R108 T2T 序列时必须保留 `R108_taxon_identity_requires_reconciliation`。不能把 R108 当作 A17 的第二份序列，也不能把 A17 的 sequence applicability 自动转移给 R108。

## 已报道单位

| repeat | 报道长度 | A17 | R108 | 可复用状态 | 主要证据与限制 |
| --- | ---: | --- | --- | --- | --- |
| MtR1 | 166 bp | pericentromeric | 2004 paper 报道 R108-1 缺失 | 仅名称/长度排除 | PMID 15480726；论文还列出 AQ841077 marker/BAC context，但它不是已验证的 MtR1 monomer 序列 |
| MtR2 | 183 bp | pericentromeric | 2004 paper 报道 R108-1 缺失 | 仅名称/长度排除 | PMID 15480726；未取得直接序列 accession |
| MtR3 | 166 bp | functional centromere | 2004 paper 报道 R108-1 各二价体 centromeric | 仅名称/长度排除 | PMID 15480726；FISH 与阵列跨度不能替代序列 accession |
| CentM168 | 168 bp | 主要 centromeric satellite | R108 centromere 几乎全为该 repeat | 仅名称/长度排除 | PMID 40714838；CENH3 enrichment 是功能线索，精确序列仍 unresolved |
| CentM183 | 183 bp | species-specific satellite，外围 | 未在公开摘要中报告 | 仅名称/长度排除 | PMID 40714838；不能把未报道当作缺失 |
| CentM51 | 未报告 | unresolved | unresolved | 仅名称排除 | 2025 研究官方摘要列为 chromosome-specific satellite，材料/染色体归属未取回 |
| CentM515 | 未报告 | unresolved | unresolved | 仅名称排除 | 同上 |
| CentM287 | 未报告 | unresolved | unresolved | 仅名称排除 | 同上 |

“仅名称/长度排除”表示候选可在报告表中标成已报道标签，不能进行序列相似性排除。即使候选周期为 166/168/183 bp，也不能仅凭周期判为同一家族或 HOR。

## 给 Result 6 的使用规则

1. 先保留完整候选表，再添加 `reported_repeat_label`。只有当直接序列、GenBank/RefSeq accession 或作者补充 FASTA 可复核时，才进入 exact、reverse-complement、rotation 和宽松相似性排除。
2. 对本清单的空 FASTA，工具输入应记录 `sequence_not_retrieved`，而不是“无已知序列”。候选如果仅与 MtR3 或 CentM168 的长度相同，仍是 `candidate_novel_not_confirmed` 或 `known_family_unresolved`，不能称 novel。
3. 若未来得到序列，保留 accession/version、原始方向、monomer 边界、identity、alignment coverage、e-value/score 和检索日期。旋转或反向互补匹配要报告变换，不要覆盖原始序列。
4. 单元长度的整数倍只能触发 `heuristic_period_multiple_not_validated_hor`。只有稳定多单元次序、阵列边界和独立长读长/T2T支持，才可升级为 `novel_array_structure_candidate`；新序列家族和功能性 centromere 仍需分开判定。
5. CentM168/CentM183 的 CENH3、FISH、repeat-density 和 assembly location 字段要分开保存。没有直接 CENH3/FISH 证据的候选使用 `assembly_predicted_centromere` 或 `centromere_inferred_not_CENH3`。

## 必须保留的 unresolved 标签

本批次至少保留：`sequence_not_retrieved`、`reported_length_only`、`taxon_scope_unresolved`、`chromosome_assignment_unresolved`、`R108_taxon_identity_requires_reconciliation`、`R108_is_M_littoralis`、`CENH3_unvalidated`、`no_HOR_evidence`、`candidate_novel_not_confirmed`。这些标签描述证据边界，不是算法失败，也不是新颖性结论。

## 来源

- Kulikova et al. 2004, “Satellite repeats in the functional centromere and pericentromeric heterochromatin of *Medicago truncatula*”, *Chromosoma*, DOI [10.1007/s00412-004-0315-3](https://doi.org/10.1007/s00412-004-0315-3), [PMID 15480726](https://pubmed.ncbi.nlm.nih.gov/15480726/)。
- Shen et al. 2025, “Two complete telomere-to-telomere Medicago genomes reveal the landscape and evolution of centromeres”, *Molecular Plant*, DOI [10.1016/j.molp.2025.07.016](https://doi.org/10.1016/j.molp.2025.07.016), [PMID 40714838](https://pubmed.ncbi.nlm.nih.gov/40714838/)。
- NGDC [PRJCA042220](https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA042220)：A17/R108 项目及物种范围；本任务不下载其中组装或 reads。
- 中国科学院官方研究摘要（2025-08-07）：列出 CentM51、CentM515、CentM287，[官方页面](https://english.genetics.cas.cn/news_/researchnews/202508/t20250807_1049276.html)。

本文件及其机器可读附件没有把任何单位宣称为 novel；缺失 sequence locator 的项将在获得直接 accession 或作者补充序列后重新审计。
