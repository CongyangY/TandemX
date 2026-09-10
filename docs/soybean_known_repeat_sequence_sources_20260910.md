# 大豆已报道卫星重复的公开序列来源核查（2026-09-10）

本表只记录可由公开原始论文或主数据库核验的序列来源。`resolved_direct` 表示存在明确的 GenBank accession，可通过 NCBI FASTA 端点取得；`family_source` 表示有同一家族的克隆序列，但该记录不能安全地等同于某一个长度标签；`unresolved` 表示论文报告了重复名称或 TRF 标签，但在本次核查范围内没有发现对应的独立序列 accession。没有下载基因组或原始 reads；仅读取了小型 GenBank 记录和论文/数据库元数据。

## 1. 可直接取得的 GenBank 序列

NCBI 的统一 FASTA 端点为：

```text
https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=<ACCESSION.VERSION>&db=nuccore&report=fasta&retmode=text
```

| 重复/历史名称 | accession（记录长度） | 原始记录中的物种/材料与注释 | 状态与可用性 |
|---|---|---|---|
| SB92 | [U11026.1](https://www.ncbi.nlm.nih.gov/nuccore/U11026.1)（183 bp） | *Glycine max* BSR-101；clone `pNAU-401`；记录注释为 SB92 tandem satellite，183 bp 为 XbaI dimer clone | **resolved_direct**；[FASTA](https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=U11026.1&db=nuccore&report=fasta&retmode=text)。该记录来自 [Vahedian et al. 1995（PMID 8541510）](https://pubmed.ncbi.nlm.nih.gov/8541510/)，是最清楚的 SB92/U11026 对应关系。 |
| CentGm-1（家族序列；PNAS 2023 的长度标签映射需谨慎） | [OP605952.1](https://www.ncbi.nlm.nih.gov/nuccore/OP605952.1)（1,070 bp） | *G. max* clone `CentGm-1`；`repeat_region` 注释为 `CentGm-1` tandem/centromeric repeat；[Tek et al. 2024（PMID 39353738）](https://pubmed.ncbi.nlm.nih.gov/39353738/)，Sanger clone | **resolved_direct / family_source**；[FASTA](https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=OP605952.1&db=nuccore&report=fasta&retmode=text)。这是包含多个变异单体的克隆片段，不应截取整条记录后直接当作一个单体。 |
| CentGm-1 原始 ChIP-clone 集合 | [AB536705.1](https://www.ncbi.nlm.nih.gov/nuccore/AB536705.1)（176 bp）、[AB536707.1](https://www.ncbi.nlm.nih.gov/nuccore/AB536707.1)（178 bp）、[AB536708.1](https://www.ncbi.nlm.nih.gov/nuccore/AB536708.1)（73 bp）、[AB536710.1](https://www.ncbi.nlm.nih.gov/nuccore/AB536710.1)（299 bp）、[AB536712.1](https://www.ncbi.nlm.nih.gov/nuccore/AB536712.1)（168 bp） | *G. max* clones G48/G56/G61/G71/G87；NCBI `repeat_region` 注释为 `CentGm-1` | **resolved_direct / family_source**；每条 accession 均可按上述 NCBI FASTA 端点取得。原始论文为 [Tek et al. 2010（PMID 20204495）](https://pubmed.ncbi.nlm.nih.gov/20204495/)。 |
| CentGm-4（约 413 bp 家族） | [OP605953.1](https://www.ncbi.nlm.nih.gov/nuccore/OP605953.1)（904 bp） | *G. max* clone `CentGm-4`；`repeat_region` 注释为 `CenGm-4` tandem/centromeric repeat；[Tek et al. 2024（PMID 39353738）](https://pubmed.ncbi.nlm.nih.gov/39353738/) | **resolved_direct / family_source**；[FASTA](https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=OP605953.1&db=nuccore&report=fasta&retmode=text)。记录含多个变异单体/片段，整条 904 bp 不是单个 413-bp 单体。 |
| CentGm-4 原始 ChIP-clone 集合 | [AB536703.1](https://www.ncbi.nlm.nih.gov/nuccore/AB536703.1)（432 bp）、[AB536704.1](https://www.ncbi.nlm.nih.gov/nuccore/AB536704.1)（361 bp）、[AB536706.1](https://www.ncbi.nlm.nih.gov/nuccore/AB536706.1)（147 bp）、[AB536711.1](https://www.ncbi.nlm.nih.gov/nuccore/AB536711.1)（356 bp）、[AB536714.1](https://www.ncbi.nlm.nih.gov/nuccore/AB536714.1)（178 bp） | *G. max* clones G29/G30/G51/G75/G93；NCBI `repeat_region` 注释为 `CentGm-4`（部分记录为截短片段） | **resolved_direct / family_source**；每条 accession 均可按 NCBI FASTA 端点取得。原始论文为 [Tek et al. 2010（PMID 20204495）](https://pubmed.ncbi.nlm.nih.gov/20204495/)。AB536713.1 为 G92a retrotransposon + G92b/CentGm-4 混合克隆，不能作为干净的 CentGm-4 单体来源。 |

## 2. 名称和 accession 的对应边界

* Gill et al. 2009 的原始研究把 CentGm-1 定义为约 92 bp 单体、CentGm-2 定义为约 91 bp 单体；Findley et al. 2010 的图注也明确将 92-bp probe 标为 CentGm-1、91-bp probe 标为 CentGm-2。论文和数据库入口：[PMC2773056](https://pmc.ncbi.nlm.nih.gov/articles/PMC2773056/)、[PMC2907198](https://pmc.ncbi.nlm.nih.gov/articles/PMC2907198/)。
* Liu et al. 2023 的 PNAS 论文则按长度命名 `CentGm91`, `CentGm92`, `CentGm273`, `CentGm413`, `CentGm444`，并在正文写为 `CentGm91 = CentGm-1`, `CentGm92 = CentGm-2`, `CentGm413 = CentGm-4`；[PMC10589659](https://pmc.ncbi.nlm.nih.gov/articles/PMC10589659/)。这与上述历史的 91/92 命名顺序不一致。因此，TandemX 输入准备时应以实际序列比对确认 `CentGm91`/`CentGm92`，不能仅按名字把 OP605952 或 U11026 强行归入其中一个长度标签。
* NCBI 的 `esearch` 对精确词 `CentGm-2` 未返回独立 nucleotide record；Gill et al. 2009 的 CentGm-2 序列来自 WGS clone 搜索/提取（其 WGS 数据范围包括 `CL866971–CL866979`, `CL866987–CL868441` 等），论文没有给出一个可单独、无歧义地代表 CentGm-2 单体的 GenBank accession。故 `CentGm92`（或按历史命名的 CentGm-2）目前标为 **unresolved**，只能将上述 WGS 数据作为进一步检索边界。

## 3. CentGm273 和 CentGm444

PNAS 2023 明确称 `CentGm273` 与 `CentGm444` 为该研究中新发现的两类重复；正文及关联数据只给出 ZH13/27 个材料的比对、FISH/注释结果，研究数据提交在 NGDC GSA [CRA009445](https://ngdc.cncb.ac.cn/gsa/browse/CRA009445)，并没有列出独立的 GenBank nucleotide accession。两者目前均为 **unresolved**：可从该研究的 ZH13 组装/附表按坐标提取候选单体，但这属于 assembly-derived sequence，需要另行保存提取坐标和版本，不能把 CRA009445（ChIP-seq/input 数据）当作 repeat FASTA accession。

## 4. YSD56 的 TRF 标签

Lian et al. 2025 的 YSD56 论文报告 TRF 扫描得到 `trf91`, `trf92`, `trf182`, `trf273`, `trf276`, `trf183`, `trf184` 七个高周期重复，并指出 `trf91`/`trf92` 富集区对应候选着丝粒；[PMC12350797](https://pmc.ncbi.nlm.nih.gov/articles/PMC12350797/)，方法段明确说明这些是 TRF 识别的标签，而非已提交的独立 repeat record。

| YSD56 标签 | 独立 repeat accession | 可追溯的公开序列来源 | 状态 |
|---|---|---|---|
| trf91, trf92, trf182, trf183, trf184, trf273, trf276 | 未发现 | YSD56 GenBank assembly [GCA_040083835.1](https://www.ncbi.nlm.nih.gov/assembly/GCA_040083835.1/)，assembly/annotation 另见 [Figshare 10.6084/m9.figshare.26004370.v2](https://doi.org/10.6084/m9.figshare.26004370.v2)；原始 reads 项目为 [PRJNA1095640](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1095640)，run 集合 SRP502474 | **unresolved / assembly-derived**。论文的补充文件为 [41597_2025_5741_MOESM1_ESM.docx](https://pmc.ncbi.nlm.nih.gov/articles/instance/12350797/bin/41597_2025_5741_MOESM1_ESM.docx)，但主文没有为七个 TRF 标签提供独立 FASTA accession。使用时应从固定版本 GCA_040083835.1 按论文/附表坐标提取并记录坐标；不能把 `trfXXX` 当作 NCBI accession。 |

## 5. 对 TandemX 数据准备的直接结论

1. 可立即作为公开小型序列输入的是 `U11026.1`（SB92，183-bp dimer）、`OP605952.1` 及 AB536705/707/708/710/712（CentGm-1 家族），以及 `OP605953.1` 及 AB536703/704/706/711/714（CentGm-4 家族）。
2. `CentGm91`/`CentGm92` 的长度标签与历史 `CentGm-1`/`CentGm-2` 命名存在顺序冲突；在未做序列比对前，只能写 family-level source，不能声称一一对应。
3. `CentGm273`、`CentGm444` 和 YSD56 的七个 `trf` 标签目前没有被核验为独立 GenBank repeat accession。它们的可复现来源是固定版本组装（或研究附表中的 assembly 坐标），不应伪造 FASTA accession。
