# 豆科 tandem/satellite repeat novelty 审计

**审计日期：** 2026-09-10  
**用途：** 为 TandemX 后续候选发现提供可追溯的豆科背景、排除规则和发现优先级；不预先建立 manuscript 的第六个 Result。  
**范围：** 以论文、GenBank/RefSeq、ENA、NGDC 以及可追溯的重复序列数据库或作者重复库为证据。此次审计只读，不下载大规模测序数据，不改变 TandemX 算法。

## 结论先行

花生（*Arachis hypogaea*）是当前最适合做跨材料、跨亚基因组的 discovery 审计的对象：2026 年论文给出了二倍体祖先和 4 个四倍体材料的 T2T 组装，同时已有 B-c、H-b 以及 centromere monomer 的历史证据。优先问题应是 **同一已知家族在不同材料/At-Bt 亚基因组中的阵列组成、边界和结构是否产生未命名的变体**，而不是把旧家族重新命名为新家族。

苜蓿（*M. sativa*）有清楚的 E180（约 185–189 bp）和 MtR3 证据，但本次核对到的 *M. sativa* 组装仍是 chromosome/draft 级，不能把搜索结果直接解释为完整着丝粒新家族。*Medicago truncatula* A17 与 *M. littoralis* R108 已有 T2T 和 CentM168/CentM183 证据，可作为 *Medicago* 的已知阳性对照；它们不是 *M. sativa*。

对“发现潜力”的判断必须分为三种：新序列家族、新的阵列结构（例如 HOR 或单元次序）以及新的功能性着丝粒证据。数据库或论文中未检出只能记为 `sequence_not_retrieved`、`metadata_blocked` 或 `candidate_novel_not_confirmed`，不能记为“novel”。

## 物种和材料优先级

| 优先级 | 物种/材料 | 组装和可追溯锚点 | 已知重复证据 | Result 6 可检验问题 | 主要限制 |
|---|---|---|---|---|---|
| 最高 | *A. hypogaea*；*A. duranensis* V14167、*A. ipaensis* K30076、四倍体 S245、HN873、HN51、S83 | 2026 T2T 论文；NCBI BioProject PRJNA1259200、NGDC PRJCA026588；原始数据 PRJNA1259747；资源页 omicsplant.cn/peanut 与 figshare 数据集 | B-c 约 115 bp；H-b-Ah 317 bp，GenBank KF957858；论文报告 peanut 的 116-bp repeat，但其角色/命名在正文中不清楚 | 比较 At/Bt 及材料间着丝粒阵列的 monomer、变体、单元次序、阵列长度和边界；检查旧 B-c/H-b 是否覆盖候选 | 2026 文献的着丝粒单元命名和 116-bp 句子需要回到补充数据核对；功能证据和 HOR 不能仅由周期峰推定 |
| 高 | *Lablab purpureus* cv. Linghu | 2025 T2T 论文；GenBank GCA_051403735.1、SRA SRP583761 | 11 条染色体均报道 centromere，CentroMiner 预测卫星和逆转座子富集；摘要未给出稳定的家族名/monomer | 在完整着丝粒上下文中建立作者库之外的 monomer 与位置证据，并寻找可复现的材料特异变体 | “centromere”主要为计算预测，需 CENH3/FISH 或其他正交证据 |
| 高 | *Astragalus membranaceus* | 2025/2026 T2T 论文（GigaScience，DOI 10.1093/gigascience/giaf117） | 8 条染色体、16 个端粒；centromere 平均约 2.949 Mb，tandem repeat 约占 centromere 54.49%；未核对到命名 monomer/HOR | 研究充分组装但重复命名不足的材料，优先查明是否存在跨染色体共享家族或染色体特异家族 | 本次审计尚未从论文页确认 GenBank/RefSeq accession；应先补齐 accession/version |
| 高 | *Sesbania cannabina* | 2023 T2T 论文，PMID 37897613 | 报道约 2.087-Gb 四倍体组装，含 centromeric regions；未核对到家族名 | 在异源四倍体中区分 homeolog-specific 与共享 TR | accession、亚基因组和作者重复库需要回查，当前为 `metadata_blocked` |
| 中高 | *Vigna unguiculata* HJD、FC6 | 2025/2026 T2T 论文；两材料、11 条染色体 | CEN455（455 bp）为已知 FISH satellite；完整 CEN455 阵列覆盖 7/11 centromeres，其余 4 条为部分阵列或按 TR 密度/片段推断 | 比较 CEN455 阵列断裂、边界和材料间结构；检验剩余 4 条 centromere 是否由另一家族组成 | 只有部分位置有 CEN455 支持；其余应标 `centromere_inferred_not_CENH3` |
| 中高 | *M. sativa*（四倍体栽培材料及 PI464715） | 2020 参考基因组 PRJNA540215；PI464715 为染色体级但非 T2T；NGDC GWHFIBZ00000000.1/PRJCA033031 为 draft chromosome-level | E180 185–189 bp；MtR3 在多个 *Medicago* 物种/材料有保守和着丝粒 FISH 证据 | 重点做 E180 家族内变体、阵列扩张/收缩、亚基因组分化和非着丝粒位置，而不是宣称新 centromeric family | E180 PCR/FISH 跨种结果表明存在广泛但位置多变的家族；FISH 分辨率约 4 Mb，不能证明小阵列位置 |
| 中 | *Phaseolus vulgaris* YP4 | 2025 组装，BioProject PRJNA1072282；11 条染色体，9 条 T2T pseudomolecule | 预测 11 个 centromere，报告 11 个 tandem units；7 个 unit 按 85% CD-HIT 聚类；平均 tandem repeat 约 46.99% | 结构新颖性和单位组合值得筛查，尤其是未被 FISH/ChIP 支持的阵列 | 作者明确要求 FISH/ChIP 确认 authentic centromere，不能直接称功能性 satellite |
| 中 | *Vigna radiata* IPM02-03 | 2025 chromosome-scale、非 T2T 组装 | CEN177、CEN174；CEN174 被描述为 CEN177 阵列的起始/连接单元；chr9 缺乏预期信号 | 检查 CEN174/177 的 junction、chr9 和 unplaced scaffold；可产生阵列结构候选 | TRASH/组装证据为主，功能性验证不足；`assembly_gap_or_partial` |
| 低（新家族）/高（流程对照） | *Glycine soja* YSD56、*G. max* | 2025 T2T/高完整度资源；*G. soja* 论文 DOI 10.1038/s41597-025-05741-y | trf91、trf92、trf182、trf273、trf276、trf183、trf184；soybean 已有 CentGm91/92/273/413/444 | 用于检出已知家族、反向互补、period multiple 和位置排除流程 | 新家族优先级低；CentGm/91、92 等既有命名使假阳性风险高 |

## 已报道的重复家族和功能边界

下表只记录本次能由论文或序列记录直接追溯的内容。缺失的 accession、单元序列或功能实验不填猜测值。

| 物种/家族 | monomer 或重复单位 | 位置/功能证据 | 可追溯来源 | 审计解释 |
|---|---:|---|---|---|
| *A. hypogaea* B-c | 约 115 bp | 2012 Cot-1 文库；主要在 B 染色体，类似着丝粒 satellite，但非 CENH3 功能证明 | Zhang et al. 2012，PMID 22797674，DOI 10.1159/000339455 | `repeat_role_ambiguous`；本次未确认稳定 GenBank accession |
| *A. hypogaea* H-b-Ah | 317 bp | FISH 22 个信号（18 个 A-genome、4 个 B-genome），多为 A-genome pericentromeric；GenBank KF957858 | 2016 New Phytologist，DOI 10.1111/nph.13999；[KF957858](https://www.ncbi.nlm.nih.gov/nuccore/KF957858) | 反复命名/相似序列不等于 CENH3 功能；B-c 和 H-b 必须分开登记 |
| *A. hypogaea* 2026 T2T repeat | 116 bp（论文正文此处角色不清） | Lastz 用于核对 telomere 位置的 repeat 句子；不能据此确定为 CentO | [2026 peanut T2T paper](https://www.nature.com/articles/s41588-026-02577-z) | `repeat_role_ambiguous`、`sequence_not_retrieved`；等待补充表/作者库核对 |
| *M. sativa* E180 | 185–189 bp；共识约 189 bp；历史报道约 1.8×10^5 copies/genome | 37/40 *Medicago* taxa PCR 阳性；FISH 可见 subtelomeric、intercalary、proximal/pericentromeric 多种位置 | [1993 alfalfa paper](https://pubmed.ncbi.nlm.nih.gov/8349123/)；[2012 Medicago study](https://pmc.ncbi.nlm.nih.gov/articles/PMC3286279/)；历史 accession L07381、L20583、L20584、X64425–X64428、X68165–X68166、Z22861 | E180 是已知家族；跨种扩增不代表保守的 centromere 功能；FISH 分辨率约 4 Mb |
| *Medicago* MtR3 | 论文确认保守重复，但未在本审计中把单一长度当作统一 monomer | 在 *M. truncatula*, M. minima, M. edgeworthii, M. ruthenica, M. caerulea, M. sativa 和四倍体 *M. falcata* 有 centromeric FISH；其他测试物种无信号 | 2017 paper，DOI 10.1007/s13258-017-0556-1 | `taxon_scope_unresolved`；序列存在与 FISH/功能位置须分开 |
| *M. truncatula* MtR1/MtR2/MtR3 | MtR3 的 centromeric array 约 450 kb 至 >1 Mb；MtR1/MtR2 为 pericentromeric | A17 中 MtR3 为功能 centromere repeat；R108-1 也见于所有 centromeres，MtR1/MtR2 缺失 | [2004 centromere paper](https://pubmed.ncbi.nlm.nih.gov/15480726/) | *Medicago* 已知阳性对照，不能转移命名到 *M. sativa* |
| *M. truncatula* / *M. littoralis* CentM168/CentM183 | CentM168 168 bp；CentM183 183 bp | 2025 T2T 研究：A17 主要为 CentM168 + species-specific CentM183；R108 几乎全为 CentM168；CentM168 有 CENH3 enrichment，CentM183 在外围 | [2025 Medicago T2T paper](https://pubmed.ncbi.nlm.nih.gov/40714838/)；NGDC PRJCA042220/CRA027262 | 这是 *Medicago* 的 T2T 功能参照，不是 *M. sativa* 证据 |
| *V. unguiculata* CEN455 | 455 bp | FISH/组装关联；两 T2T 材料中 7/11 染色体有完整阵列，其余部分或推断 | cowpea T2T paper；已知 satellite paper DOI 10.1038/s42003-020-01507-x | `centromere_inferred_not_CENH3` 适用于仅由 TR 密度/片段推断的染色体 |
| *V. radiata* CEN177/CEN174 | 177 bp、174 bp | TRASH 识别；CEN174 被解释为 CEN177 array 的起始/连接单元；chr9 signal 缺失 | [2025 mungbean assembly](https://pmc.ncbi.nlm.nih.gov/articles/PMC11842788/) | 已知家族；junction 可作为 structure candidate，不能称新家族 |
| *G. soja* / *Glycine* CentGm | CentGm91/92/273/413/444；SB92 accession U11026 | CENH3-associated satellite 证据和 T2T 比较资源 | [pan-centromere study](https://academic.oup.com/plcell/article/35/4/1527/6973384)；[wild soybean T2T record](https://www.nature.com/articles/s41597-025-05741-y) | 流程已知阳性，任何相同/近似候选均应降为 `not_novel` 或 `known_family_variant` |
| *Cicer arietinum* CaSat1/CaSat2/CaRep1/2 | CaSat1 约 162–168 bp；CaSat2 100 bp；CaRep 为 Ty3-gypsy-like pericentric repeat | CaSat1 邻近 18S rDNA；CaSat2 为各染色体 pericentric heterochromatin | [chickpea repeat paper](https://pubmed.ncbi.nlm.nih.gov/10344208/)；GenBank AJ005998–AJ006007 | 已知豆科对照，不是 T2T 新发现 |

Fabeae 的广泛研究在 14 个物种中找到 64 个 centromeric satellite families，monomer 从 33 到 2,979 bp；多数为物种特异，也有跨物种 superfamily，且同源序列不一定都在 CENH3 中富集。[该研究](https://academic.oup.com/mbe/article/37/8/2341/5817320)说明“物种特异”“长度独特”与“功能性着丝粒”是三个独立命题。*Vicia faba* 的 23 个新 satellite 中仅 7 个由 CENH3 证据支持为 centromeric；*Pisum sativum* 的 FabTR-10 还显示 459-bp 与 1,975-bp 结构变体。它们适合作为严格排除与功能验证的参考，不能替代目标豆科材料的证据。

## `novel` 分层判据

输出中建议同时保留 `novelty_tier` 和 `novelty_status`，不要用一个布尔值代替证据链。

| 层级 | 含义 | 允许的措辞 | 必需证据 |
|---|---|---|---|
| N0 | 与已命名家族有 exact/近似匹配，或只是把旧 accession 重新聚类 | `not_novel`、`known_family` | 序列比对、反向互补检查、来源 accession 和版本 |
| N1 | 未找到 exact match，但与已知 family 有明显同源关系；可能是 subfamily/变体 | `known_family_variant`、`candidate_novel` | 多拷贝共识、相似性敏感性分析、周期和材料范围；不能直接命名新 family |
| N2 | monomer 已知，但发现新的阵列结构、单元次序、junction 或材料特异扩张/收缩 | `novel_array_structure_candidate` | 完整/近完整阵列、长读长或 T2T 组装支持、父家族与子结构同时报告；HOR 需独立证据 |
| N3 | 在预先登记的数据库、论文重复库和同属/近缘材料中未形成 exact/近似匹配的候选 sequence family | `candidate_novel_sequence_family` | 可复现共识、多拷贝周期、阵列跨度/拷贝证据、来源和 assembly version 完整；仍需正交验证 |
| N4 | N3 候选同时有功能性 centromere 证据，且在独立材料或实验中可复现 | `validated_novel_centromeric_family` | CENH3 ChIP/CUT&Tag、FISH/定位或等价正交证据，序列和坐标可复核 |

`N3` 不能从“数据库没有命中”单独得到；数据库覆盖范围、检索参数、序列是否实际取回、组装是否含完整阵列都必须进入结果。仅有周期峰、RepeatMasker 类别或 centromere repeat density 时，最多为 `candidate_novel`/`assembly_predicted_centromere`。

## 序列、周期和位置排除流程

1. **登记来源。** 对每个候选保存物种、材料、亚基因组、assembly accession/version、染色体、坐标、gap/N 状态、来源论文、作者重复库和检索日期。没有 accession 或版本时使用 `accession_unresolved`/`metadata_blocked`。
2. **建立已知库。** 合并 GenBank/RefSeq、ENA、NGDC、PlantSat（旧数据库，需重新核实可用性）、作者补充库，以及有明确 accession 的 Repbase/REXdb/RepeatExplorer 结果。保留原始名称、别名、monomer 序列和文献，而不是只保留一个聚类 ID。
3. **序列相似性排除。** 依次做 exact、反向互补、短序列敏感比对和较宽松相似性比对；报告 identity、coverage、alignment orientation 和 e-value/score。相似性阈值应做敏感性分析，不把单一 80% 或 90% 阈值写成生物学定律。对低复杂度、AT-rich、短 monomer 和断裂组装分别标记。
4. **确认 tandem 证据。** 记录 monomer 共识如何生成、period 分布、支持拷贝数、连续阵列跨度、方向、GC/entropy 和跨 reads/contigs 的一致性。简单 TR、satellite、LTR 内部重复和 organelle-derived repeat 分开分类。
5. **检验 period multiple/HOR。** 检查候选 period 是否为已知 period 的整数倍，并同时报告父周期、子周期、单元次序和替代解释。周期比值本身只能是 `heuristic_period_multiple_not_validated_hor`；只有重复的多单元模式、稳定边界以及长读长/T2T 支持时才可写 `HOR_candidate`，有独立实验或多材料复现后才可写 validated HOR。
6. **定位排除。** 将候选坐标与 centromere boundary、端粒、rDNA、基因密度、TE、pericentromere 和 unplaced scaffold 分开比较。`repeat_dense`、CentroMiner 预测或 assembly gap 不能独立证明着丝粒功能；没有 CENH3/FISH 时写 `assembly_predicted_centromere` 或 `centromere_inferred_not_CENH3`。
7. **跨材料和亚基因组核对。** 在花生中必须分别检查 At/Bt 及二倍体祖先；在四倍体苜蓿、*Sesbania* 等材料中区分 homeolog-shared、homeolog-specific 和样本/单倍型特异序列。保留“同家族不同位置”和“相似序列不同功能”的可能性。
8. **正交验证。** 对 N2/N3 候选优先使用独立 assembly、Illumina/ONT/HiFi 覆盖、FISH、CENH3 ChIP/CUT&Tag 或独立材料。读段支持 abundance 不能替代位置功能；新 assembly 也不应自动当作 copy-number truth。

## 必须保留的 uncertainty 标签

建议在 Result 6 的候选表中至少保留以下字段或等价标签：

`candidate_novel_not_confirmed`、`known_family_variant`、`sequence_not_retrieved`、`accession_unresolved`、`metadata_blocked`、`assembly_gap_or_partial`、`sample_or_haplotype_mismatch`、`copy_number_unresolved`、`repeat_role_ambiguous`、`assembly_predicted_centromere`、`centromere_inferred_not_CENH3`、`CENH3_unvalidated`、`FISH_resolution_limited`、`no_HOR_evidence`、`heuristic_period_multiple_not_validated_hor`、`taxon_scope_unresolved`、`organellar_origin_candidate`、`unplaced_contig_location`。

其中 `sequence_not_retrieved` 表示本次没有取得可比序列，`not_novel` 才表示有足够的 exact/近似匹配证据；二者不可互换。`centromere_inferred_not_CENH3` 表示位置来自重复密度、CentroMiner、组装结构或文献推断，不表示已证明 CENH3 结合。`organellar_origin_candidate` 特别适用于 *Medicago polymorpha* 中由线粒体 rps10-like 片段扩增形成的 2,157、1,064、987、971 和 587 bp repeat motifs；它是有趣的核内 tandem repeat 起源候选，但不能直接称 centromeric satellite。

## 推荐的发现顺序

1. **先做花生六材料/两个祖先的统一库。** 这是证据最完整且跨亚基因组问题最清楚的对象。第一轮目标是复现 B-c/H-b/已知 centromeric repeat，并量化其变体与阵列结构；第二轮才筛选 N1–N3 候选。该分析在完成验证前只形成独立候选报告，不进入 manuscript 正文。
2. **随后做 Lablab、Astragalus 和 Sesbania。** 三者的共同价值是 T2T 或接近 T2T、已有 centromere 区域、但论文中缺少可直接复用的家族命名。它们是“组装完整而重复注释不足”的高潜力材料。先补齐 GenBank/RefSeq/NGDC accession 和作者库，再开始 novelty 判定。
3. **用 cowpea 做结构验证。** CEN455 已知且部分染色体仍是推断位置，适合测试“已知 monomer + 新 array architecture”和“未被 CEN455 覆盖的 centromere”两条分支。
4. **苜蓿用于家族变体和位置多样性。** E180、MtR3 已知，且 *M. sativa* 非 T2T；可产出可靠的 variant/position 结果，但 N3/N4 新 centromeric family 的优先级低于花生和上述 T2T 豆科。
5. **Phaseolus 和 Vigna radiata 用于中等优先级候选。** 两者已有 11 个或 2 个主要 TR 单元，但仍缺少统一功能验证；重点放在结构候选和不确定性保留。
6. **Glycine、chickpea 和 Fabeae 作为已知阳性/负对照。** 先确保流程能把 CentGm、CaSat、Fabeae 旧家族识别为 N0 或 N1，再解释其他豆科候选。

本审计没有把任何物种或序列宣布为“已发现的新家族”。后续只有在 accession/version、可复核序列、连续 tandem 证据、位置证据和独立验证均达到相应层级后，才可在论文中使用 `candidate novel` 或 `validated novel centromeric family`。

## 主要来源

- *A. hypogaea* 重复/FISH：Zhang et al. 2012，PMID [22797674](https://pubmed.ncbi.nlm.nih.gov/22797674/)；2016 New Phytologist，DOI [10.1111/nph.13999](https://nph.onlinelibrary.wiley.com/doi/10.1111/nph.13999)。
- 花生 T2T：Nature Genetics 2026，DOI [10.1038/s41588-026-02577-z](https://www.nature.com/articles/s41588-026-02577-z)；NCBI PRJNA1259200；NGDC PRJCA026588。
- E180：1993，PMID [8349123](https://pubmed.ncbi.nlm.nih.gov/8349123/)；2012，PMC [3286279](https://pmc.ncbi.nlm.nih.gov/articles/PMC3286279/)。
- *Medicago* centromere：2004，PMID [15480726](https://pubmed.ncbi.nlm.nih.gov/15480726/)；2025 T2T，PMID [40714838](https://pubmed.ncbi.nlm.nih.gov/40714838/)；NGDC PRJCA042220/CRA027262。
- PlantSat：Macas et al. 2002，DOI [10.1093/bioinformatics/18.1.28](https://academic.oup.com/bioinformatics/article/18/1/28/243554)。旧站点可用性和当前内容需在真正分析前复核。
- 豆科比较：Fabeae，DOI [10.1093/molbev/msaa090](https://academic.oup.com/mbe/article/37/8/2341/5817320)；*Vicia faba*，DOI [10.1038/s41598-018-24196-3](https://www.nature.com/articles/s41598-018-24196-3)；*Pisum sativum*，PLOS Genetics [10.1371/journal.pgen.1010633](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1010633)。
- 其他 T2T/高完整度资源：*L. purpureus* PMID [41233345](https://pubmed.ncbi.nlm.nih.gov/41233345/)；*P. vulgaris* PMID [40366866](https://pubmed.ncbi.nlm.nih.gov/40366866/)；*G. soja* DOI [10.1038/s41597-025-05741-y](https://www.nature.com/articles/s41597-025-05741-y)；*V. radiata* PMC [11842788](https://pmc.ncbi.nlm.nih.gov/articles/PMC11842788/)。
