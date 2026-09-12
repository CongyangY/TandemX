# YSD56 TXF000367/TXF000708 序列复杂度与潜在伪阳性审计

**日期：** 2026-09-10  
**范围：** 只读检查 YSD56 pilot 的 `families.tsv`、`monomers.fa`、`candidate_monomers.fa`、`candidate_reads.tsv` 和 `monomer_membership.tsv`。输入目录为 `/Volumes/T7/Codex/TandemX/results/legume_ysd56_novel_tr_pilot_v1_20260910/run/discover/`。  
**边界：** 不上传序列、不改算法、不重跑 discovery；以下结果只描述序列复杂度和局部退化风险，不能称为 novel、着丝粒重复或 HOR。

## 结论

TXF000367 是 335 bp 的高一致性 family consensus，碱基组成接近均衡，Shannon entropy 为 1.9746 bit/base（理论上限 2），没有显著低复杂度判据；但其 5′ 端有一个 16-bp poly-A run。TXF000708 是 785 bp consensus，A/T 各 260 bp、GC fraction 0.3376，entropy 为 1.9103 bit/base；它的短 motif 重复和压缩性比 TXF000367 更强，但仍没有接近纯单一 motif 的全长周期。

两条序列之间没有现有 audit 输出的 related pair，也没有 exact cyclic/反向互补相同序列。它们共享全部 3-mers、114/126 与 114/120 个 canonical 4-mers，但只共享 14 个 8-mers，且不共享 11-mers。这是短 motif 退化的预期现象，不能作为同一 family、同源或生物学联系的证据。与各自其余 1,473 个本地 monomer 的 11-mer 比较中，TXF000367 的最高共享仅为 2/320 个目标 unique 11-mers，TXF000708 的最高共享为 4/734；没有看到支持全长关系的短 motif 证据。

## family 与 read-level 支持

| 指标 | TXF000367 | TXF000708 | 审计解释 |
|---|---:|---:|---|
| monomer length | 335 bp | 785 bp | `families.tsv` 的 dominant period |
| consensus MD5 | `8f2bc6414dc11839d9ea814edc427707` | `76c737f3773ed281915f26212b49ad46` | 同时用于回指 `monomers.fa` |
| support reads | 16 | 8 | pilot operational support，不是 abundance truth |
| support span | 217,921 bp | 94,397 bp | family 汇总的 read-derived span |
| mean identity | 0.9945 | 0.9933 | cluster 内 consensus 一致性 |
| confidence | high | high | 仍带 `uncalibrated_confidence` |
| membership rows / unique reads | 16 / 16 | 9 / 8 | TXF000708 的 read `SRR28726931.2021985` 有两个 candidate 区间 |
| membership edit distance | 14 个 0、1 个 2、1 个 6 | 7 个 0、1 个 1、1 个 17 | 支持 consensus 稳定，但 TXF000708 有一个较远成员 |
| candidate periods | 335 bp ×14；334、330 各 1 | 785 bp ×8；768 bp ×1 | candidate period 与 family consensus 存在轻微偏移 |
| candidate repeat span | 5,180–22,399 bp | 4,321–16,344 bp | 约 15.509–67.876 和 5.626–20.820 estimated units |
| candidate sequence variants | 3 个 sequence；14/16 为同一 MD5 | 3 个 sequence；7/9 为同一 MD5 | 不是 16 或 9 个独立序列证据 |

`candidate_reads.tsv` 中的高 score（TXF000367 为 0.9865–0.9987，TXF000708 为 0.9771–0.9998）说明当前 read-to-consensus 匹配较好；它没有解决完整 array、材料特异性或功能位置问题。TXF000708 的重复 read ID 和 768-bp/785-bp period 差异应在后续 full-depth 复核中保留。

## 复杂度指标

指标由 `monomers.fa` 中的 family consensus 计算；压缩比为 zlib level 9 压缩字节数/原始序列长度，数值越低表示越容易压缩。低复杂度标签不以单一指标决定。

| 指标 | TXF000367 | TXF000708 | 解读 |
|---|---:|---:|---|
| A/C/G/T counts | 108/86/65/76 | 260/162/103/260 | TXF000708 明显 A/T 偏高 |
| GC fraction | 0.4507 | 0.3376 | 与 `families.tsv` 一致 |
| Shannon entropy | 1.9746 bit/base | 1.9103 bit/base | TXF000367 接近四碱基均衡；TXF000708 略低 |
| normalized entropy | 0.9873 | 0.9552 | entropy/2 |
| maximum homopolymer | A×16 | A×8 | TXF000367 的 5′ poly-A 更长 |
| zlib-9 compressed bytes / ratio | 139 / 0.4149 | 260 / 0.3312 | TXF000708 的整体重复性/可压缩性更高 |
| low-complexity flag | false | false | 这是 discovery 输出标签，不替代本次审计 |

TXF000367 的最长 run 是 5′ 端 `A×16`，其余较长 run 主要为 A×5、G×4、T×4 或 C/G×3。TXF000708 的最长 run 是 `A×8`，但 A/T 短 run 更密集，且序列中多次出现 `AACCGTCG`、`TGAAGAAGGAA`、`CAACAACCGTC` 等局部 motif。

## 低阶 motif periodicity

对每个 period (p) 计算 circular shift match fraction，即 `s[i] == s[(i+p) mod length]` 的比例；只报告 (p=2)–80，避免把单核苷酸组成误当作 period 1。以每条序列的单核苷酸组成计算的随机组成 baseline 为 TXF000367 0.2590、TXF000708 0.2792。

| 序列 | 最高非平凡 period（match fraction） | 其次候选 | 结论 |
|---|---|---|---|
| TXF000367 | p=10: 0.3433 | p=12: 0.3343；p=3/4: 0.3224 | 高于组成 baseline，但远低于稳定 tandem array 应有的近完全匹配；更像局部 motif/共识结构信号 |
| TXF000708 | p=43: 0.3669 | p=13: 0.3592；p=25: 0.3567；p=3: 0.3452 | 存在较强局部重复背景，但没有单一低阶 period 主导全长序列 |

补充的 exact canonical k-mer 统计如下。`repeat_window_fraction` 是落在出现至少两次的 canonical k-mer 中的窗口比例；短 k 值受到 A/T 和有限字母表的强烈影响。

| k | TXF000367 unique/窗口；repeat-window fraction | TXF000708 unique/窗口；repeat-window fraction |
|---:|---:|---:|
| 3 | 32/333；1.0000 | 32/783；1.0000 |
| 4 | 126/332；0.8825 | 120/782；0.9693 |
| 5 | 239/331；0.4743 | 300/781；0.8335 |
| 6 | 298/330；0.1636 | 495/780；0.5756 |
| 8 | 320/328；0.0274 | 678/778；0.2237 |
| 11 | 320/325；0.0185 | 734/775；0.0955 |

因此 TXF000708 的短 motif 退化风险高于 TXF000367；TXF000367 的主要局部伪阳性风险集中在 poly-A 起点和较短的 3–6-mer，而不是全长低复杂度。两者都不能只凭 period 峰升级为 HOR。

## 与彼此及本地候选的短 motif 关系

比较使用每个序列的 exact canonical k-mer set（序列和 reverse complement 取 canonical；不做外部数据库检索）。结果为：

| k | TXF000367 unique | TXF000708 unique | 两者 shared | TXF000367 fraction | TXF000708 fraction | 解释 |
|---:|---:|---:|---:|---:|---:|---|
| 3 | 32 | 32 | 32 | 1.0000 | 1.0000 | 仅说明 3-mer 字母表已饱和 |
| 4 | 126 | 120 | 114 | 0.9048 | 0.9500 | 短 motif 退化明显 |
| 5 | 239 | 300 | 147 | 0.6151 | 0.4900 | 仍受短字母表和 A/T 组成影响 |
| 6 | 298 | 495 | 92 | 0.3087 | 0.1859 | 不能视为 family-level relation |
| 8 | 320 | 678 | 14 | 0.0438 | 0.0206 | 低共享 |
| 11 | 320 | 734 | 0 | 0.0000 | 0.0000 | 无共享 canonical 11-mer |

与其余本地 families 的最高共享（按目标序列 shared/target-unique 计）为：

- TXF000367：TXF001291 共享 2/320（0.00625），TXF001300 3/320（0.00938），TXF000640 2/320（0.00625）；
- TXF000708：TXF001290 4/734（0.00545），TXF000421 5/734（0.00681），TXF000464 4/734（0.00545）。

在 (k=6) 时，一些 30–50 bp 的短 families 可以覆盖其自身全部 6-mer，但只占 TXF000367 的约 3.7–4.7% 或 TXF000708 的约 0.4–0.6%；这是“短 motif 被长序列包含”的退化关系，不是全长同源。现有 `family_similarity.tsv` 和 `family_hierarchy.tsv` 中没有 TXF000367 或 TXF000708 的 emitted related row；两条序列之间的最佳 ungapped local match 只有 28 bp、identity 0.5357，也未达到当前 family audit 的最小有效 overlap 语境。

对 monomer 序列做全 cyclic rotation 和 reverse-complement canonical key 后，TXF000367 与 TXF000708 均没有 exact duplicate group。上述本地短 motif 检查也没有替代 GenBank/RefSeq/ENA/NGDC 或作者库的已知 repeat 排除；当前结论只能限定在这个 YSD56 pilot catalogue。

## 伪阳性分层与处理建议

| 审计对象 | 当前观察 | 建议状态 |
|---|---|---|
| TXF000367 全长序列 | entropy 高、无 low-complexity flag；有 A×16 起点和弱 p=10 信号 | `candidate_sequence_review`; 保留 `homopolymer_edge_risk`、`no_HOR_evidence` |
| TXF000708 全长序列 | entropy 略低、zlib ratio 更低、4–6-mer 重复更多；无 low-complexity flag | `candidate_sequence_review`; 保留 `short_motif_degeneracy_risk`、`no_HOR_evidence` |
| TXF000708 read support | 8 unique reads、9 membership rows；一个 read 有两个区间；一个成员 edit distance 17 | `copy_number_unresolved`; `duplicate_read_interval` |
| TXF000367/708 彼此关系 | 无 emitted related pair；11-mer 不共享；仅 3–5-mer 高共享 | `unresolved_related_or_partial` 不应赋予；记录 `short_motif_only_no_family_relation` |
| 与本地其他 candidates | 11-mer 最高仅占目标 unique set 的 0.625% 或 0.681% | `local_catalogue_no_strong_short_motif_match`; 不等于外部数据库无匹配 |
| biological interpretation | 后续全组装定位给出 chr13/chr19 长阵列；rDNA 排除已将 TXF000367 识别为 5S-related，TXF000708 仍缺 CENH3/FISH、TE 全库和全数据库检索 | TXF000367=`rDNA_related_known_repeat`; TXF000708=`preliminary_candidate_not_novel` |

后续整合定位和 rDNA 排除已把队列缩小到 TXF000708；TXF000367 因 119/119-bp 5S rDNA 完全匹配而退出。TXF000708 已在 10.957x nested 深度精确复现。用户提供的 `core_nt` BLASTN 结果进一步显示其前 172/785 bp（21.91%）与多条公开 *Glycine max* 序列有显著局部同源，但没有近全长命中；因此状态更新为 `no_near_full_length_core_nt_match_local_Glycine_max_homology_present`。它仍需更完整的 soybean/legume 与 TE 库检索、独立 ONT/Illumina 支持及必要的实验验证。在这些证据完成前，保持 `preliminary_candidate_not_novel`，不使用 novel、centromeric 或 HOR 命名。
