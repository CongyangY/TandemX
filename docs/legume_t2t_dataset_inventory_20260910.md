# 豆科 T2T / 近 T2T 数据清单（2026-09-10）

本清单用于为 TandemX 的读段优先、装配感知分析选择公开豆科数据。核查范围是论文及其原始数据仓库页面；截至 2026-09-10 未下载 FASTQ、FASTA 或 Figshare 大文件。`T2T` 仅在原始论文明确使用 telomere-to-telomere、gap-free/无 gap，且给出端粒或着丝粒证据时使用；`chromosome-level` 单独记录，不能升级为 T2T。

## 推荐顺序

| 优先级 | 数据组 | TandemX 价值 | 主要限制 |
|---|---|---|---|
| P0 | 北大农研院花生：V14167、K30076、S245、HN873、HN51、S83 | 同一研究覆盖二倍体祖先和四倍体栽培材料；论文报告 0 gap、端粒和着丝粒分析，重复序列比例高；可做 A/B 亚基因组与多材料比较 | raw run 与 assembly 已按材料名闭合，但四倍体缺少独立 donor/提取链的公开证明；S83 的短读段记录为 released 但 reported bases=0 |
| P0 | 苜蓿 ZM4（Zhongmu-4）近 T2T | 四倍体、四单倍型、重复丰富；新装配报告仅剩 8 个未解析 gap，适合检验复杂多倍体中的 repeat discovery 和 abundance deficit | 严格说是 near-T2T，不是 gap-free；新旧数据项目和不同测序批次的个体级 donor identity 未完全公开证明 |
| P1 | 豇豆 HJD（可并入 FC6） | 两个材料均为 11 条 gap-free 染色体、22 个端粒；HJD 的 HiFi/ONT/Hi-C 量和 assembly accession 可追溯 | FC6 原始数据来自前一研究，当前论文未在数据声明中给出其原始项目 accession；论文正文的 HJD 组装长度与 Table 1 有约 25 Mb 差异 |
| P1 | 普通菜豆 YP4 | 公开 GenBank 完整基因组和 NCBI BioProject；11 条无 gap 染色体、20/22 个端粒、9 条 T2T pseudomolecule；数据量适中，适合先跑端到端 | 论文标题/结论称 T2T，但作者同时报告只有 9/11 条染色体达到 T2T；应按“部分染色体 T2T”使用 |
| P1 | ZH13（中黄13）大豆 | 完整 GWH assembly，HiFi、ONT、BGI 短读段均有覆盖信息；39/40 端粒和所有着丝粒被报道；单材料、项目级元数据较清楚 | “gap-free/T2T”来自论文与 GWH 的完整基因组记录，仍应在本地下载后检查 FASTA gap、端粒边界及 haplotype 表示 |
| P2 | YSD56 野大豆 | 对抗大豆胞囊线虫 X12 的野生种质；GenBank + SRA 双重公开；20 染色体、40 端粒、0 gap，且有 Illumina/HiFi/ONT 三类验证读段 | 不是栽培大豆；Hi-C 体积在论文不同位置写成 23.84 Gb 与 123.84 Gb，下载前必须以 SRA 文件元数据为准 |

## 逐数据核实

### 1. 北大农研院花生六材料（P0）

**来源和 accession。** 原始研究为 Bian et al., *Nature Genetics* 2026，论文 DOI 为 [10.1038/s41588-026-02577-z](https://doi.org/10.1038/s41588-026-02577-z)，全文及补充材料见 [PMC13175896](https://pmc.ncbi.nlm.nih.gov/articles/PMC13175896/)。装配数据声明为 NCBI [PRJNA1259200](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1259200) 和 NGDC [PRJCA026588](https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA026588)；论文将“genome assemblies and annotations, transcriptome and PacBio Iso-Seq reads”以及相关全基因组数据统一指向 NCBI 原始数据项目 [PRJNA1259747](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1259747)。装配和注释的作者发布镜像为 [Figshare 10.6084/m9.figshare.28883636](https://doi.org/10.6084/m9.figshare.28883636)，Figshare 页面显示总下载约 3.53 GB（装配/注释资源，不等于 raw reads）。

| 材料 | 物种/倍性 | 论文装配长度 | 公开 assembly accession | raw HiFi / ONT / Illumina accession | donor-matched 判断 |
|---|---|---:|---|---|---|
| V14167 | *Arachis duranensis*，AA 二倍体祖先 | 1.178 Gb | [GCA_054824555.1](https://www.ncbi.nlm.nih.gov/assembly/GCA_054824555.1/)；逐样本 raw 见下表 | PRJNA1259747；SRR33996330/96240/96449 | 材料匹配：是。论文称来自广西农科院的 V14167；单株/叶片池化细节需按补充表确认 |
| K30076 | *A. ipaensis*，BB 二倍体祖先 | 1.485 Gb | [GCA_059453495.1](https://www.ncbi.nlm.nih.gov/assembly/GCA_059453495.1/)；逐样本 raw 见下表 | PRJNA1259747；SRR33996329/6606/96322 | 材料匹配：是；个体级 donor identity 未单独给出 |
| S245（Laiyang Silihong） | *A. hypogaea*，异源四倍体 | 2.572 Gb | [GCA_054824585.1](https://www.ncbi.nlm.nih.gov/assembly/GCA_054824585.1/)；逐样本 raw 见下表 | PRJNA1259747；SRR33996311/6595/96310 | 材料匹配：是；论文称 S245 为莱阳四粒型材料 |
| HN873（Chedouzi） | *A. hypogaea*，异源四倍体、地方品种 | 2.622 Gb | [GCA_054824565.1](https://www.ncbi.nlm.nih.gov/assembly/GCA_054824565.1/)；逐样本 raw 见下表 | PRJNA1259747；SRR33996436/6231/6193 | 材料匹配：是；论文称材料来自四川南充 |
| HN51（Huayu23） | *A. hypogaea*，异源四倍体、栽培品种 | 2.625 Gb | [GCA_054824525.1](https://www.ncbi.nlm.nih.gov/assembly/GCA_054824525.1/)；逐样本 raw 见下表 | PRJNA1259747；SRR33996516/6327/6182 | 材料匹配：是；论文称为山东花生研究所选育的 Huayu23 |
| S83（Yunnan Rainbow Peanut） | *A. hypogaea*，异源四倍体 | 2.620 Gb | [GCA_054824515.1](https://www.ncbi.nlm.nih.gov/assembly/GCA_054824515.1/)；逐样本 raw 见下表 | PRJNA1259747；SRR33996477/6210/6584 | 材料匹配：是；论文以 S83 作为云南彩色种皮材料 |

**装配和 raw-read 证据。** 论文 Table 1 报告六个装配均为 0 gap；V14167/K30076 各有 20 个端粒，四个四倍体材料各有 40 个端粒。论文称生成六个 T2T peanut genomes，并同时使用 PacBio HiFi、ONT ultra-long、Hi-C 和 Illumina polishing。总 HiFi 为约 1.42 Tb、平均约 116×；ONT 为平均约 70×。六个装配长度加总约 13.10 Gb，因此按论文平均覆盖度估算，ONT 原始碱基量约 0.92 Tb；这是覆盖度换算值，不是仓库压缩文件大小。Illumina 用于基因组大小估计和 polishing，但本次可访问文本未给出六材料 Illumina raw 总量，记为 `unresolved`，不要把 521 个群体重测序材料的量混入装配 raw 数据。

**donor-matched 和风险。** 论文说明每个 accession 种植 20 株并收集幼叶，说明至少在 accession/材料层面是 donor-matched 的 pooled material；这不是“一个单倍体供体”的证据。四倍体材料的 A/B subgenome 不能当作两个独立 donor。建议先选 S245（栽培四倍体）和 V14167（A 亚基因组祖先）做 pilot，再扩展其余四个材料。

**第二轮逐样本核实（NCBI/ENA，2026-09-10）。** NCBI Assembly 对 `PRJNA1259200` 返回 6 个 assembly 记录；NCBI BioSample 的 `isolate` 和 `data_category` 将它们映射到六个材料。NCBI SRA `runinfo` 与 ENA `read_run` 的当前（2025 年提交）记录中，raw 样本名使用 `<accession>_hifi`, `<accession>_ont`, `<accession>_wgs`, `<accession>_hic`，因此可以按材料名一一对应。下表的 bases 来自 NCBI SRA runinfo，bytes 是 ENA `fastq_bytes` 各文件之和（压缩文件大小；`--` 表示仓库记录为 0/空，不代表可以安全下载）。每行保留最新的 `SAMN490…/SRR339…` 记录；同一材料还有旧的 `SAMN412…/SRR289…` 重复提交，bases 相同。

| 材料（assembly isolate） | 组装物种；assembly BioSample / accession / FASTA 文件 | raw BioSample（HiFi / ONT / short-read） | HiFi run：bases / bytes | ONT run：bases / bytes | Illumina-like short-read run：bases / bytes | 同一材料判断 |
|---|---|---|---:|---:|---:|---|
| V14167 | *A. duranensis*；SAMN48356820 / [GCA_054824555.1](https://www.ncbi.nlm.nih.gov/assembly/GCA_054824555.1/) / `GCA_054824555.1_ASM5482455v1_genomic.fna.gz` | SAMN49012632 / SAMN49012638 / SAMN49012650 | SRR33996330：211.520 Gb / 161.392 GB | SRR33996240：108.019 Gb / 99.089 GB | SRR33996449（DNBSEQ）：131.114 Gb / 97.367 GB | accession 名、BioSample isolate 和论文材料一致；材料级是，单株 donor 链未公开 |
| K30076 | *A. ipaensis*；SAMN48356821 / [GCA_059453495.1](https://www.ncbi.nlm.nih.gov/assembly/GCA_059453495.1/) / `GCA_059453495.1_ASM5945349v1_genomic.fna.gz` | SAMN49012633 / SAMN49012639 / SAMN49012651 | SRR33996329：175.311 Gb / 131.425 GB | SRR33996606：105.865 Gb / 99.787 GB | SRR33996322（DNBSEQ）：74.345 Gb / 59.310 GB | accession 名、BioSample isolate 和论文材料一致；材料级是，单株 donor 链未公开 |
| S245（Laiyang Silihong） | *A. hypogaea*；SAMN48356822 / [GCA_054824585.1](https://www.ncbi.nlm.nih.gov/assembly/GCA_054824585.1/) / `GCA_054824585.1_ASM5482458v1_genomic.fna.gz` | SAMN49012634 / SAMN49012640 / SAMN49012652 | SRR33996311：316.441 Gb / 247.160 GB | SRR33996595：198.190 Gb / 179.338 GB | SRR33996310（DNBSEQ）：154.720 Gb / 116.200 GB | assembly isolate=`LaiyangSilihong`；raw sample 名=S245；材料级是，单株 donor 链未公开 |
| HN873（Chedouzi） | *A. hypogaea*；SAMN48356823 / [GCA_054824565.1](https://www.ncbi.nlm.nih.gov/assembly/GCA_054824565.1/) / `GCA_054824565.1_ASM5482456v1_genomic.fna.gz` | SAMN49012635 / SAMN49012641 / SAMN49012653 | SRR33996436：282.991 Gb / 138.051 GB | SRR33996231：178.740 Gb / 161.521 GB | SRR33996193（DNBSEQ）：55.668 Gb / 42.992 GB | assembly isolate=`Chedouzi`；raw sample 名=HN873；材料级是，单株 donor 链未公开 |
| HN51（Huayu23） | *A. hypogaea*；SAMN48356824 / [GCA_054824525.1](https://www.ncbi.nlm.nih.gov/assembly/GCA_054824525.1/) / `GCA_054824525.1_ASM5482452v1_genomic.fna.gz` | SAMN49012636 / SAMN49012642 / SAMN49012654 | SRR33996516：309.594 Gb / 150.772 GB | SRR33996327：146.376 Gb / 131.956 GB | SRR33996182（DNBSEQ）：51.203 Gb / 38.098 GB | assembly isolate=`Huayu23`；raw sample 名=HN51；材料级是，单株 donor 链未公开 |
| S83（Yunnan Rainbow） | *A. hypogaea*；SAMN48356825 / [GCA_054824515.1](https://www.ncbi.nlm.nih.gov/assembly/GCA_054824515.1/) / `GCA_054824515.1_ASM5482451v1_genomic.fna.gz` | SAMN49012637 / SAMN49012643 / SAMN49012655 | SRR33996477：316.600 Gb / 235.576 GB | SRR33996210：237.065 Gb / 219.860 GB | SRR33996584（DNBSEQ）：0 / --（ENA 和 SRA runinfo 均无 payload bases） | assembly isolate=`YunnanRainbow`；raw sample 名=S83；材料级是，short-read 当前记录应标 `metadata_blocked` |

六个 assembly accession 的 `PRJNA1259200[BioProject]` 检索、assembly 摘要和 BioSample XML 可分别在 [NCBI Assembly E-utilities](https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=assembly&term=PRJNA1259200%5BBioProject%5D&retmode=json)、[Assembly summaries](https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=assembly&id=34785121,31198251,31198241,31198221,31198211,31198201&retmode=json) 和 [BioSample XML](https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=biosample&id=48356820,48356821,48356822,48356823,48356824,48356825&rettype=xml&retmode=text) 复核；raw 的项目入口是 [PRJNA1259747](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1259747)，逐 run 字段来自 [ENA Portal API](https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&query=study_accession%3D%22PRJNA1259747%22&fields=run_accession%2Csample_accession%2Cinstrument_platform%2Cbase_count%2Cfastq_bytes&format=tsv&limit=10000) 和 [NCBI SRA runinfo](https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc=PRJNA1259747)。

逐样本映射解决了“六材料 raw run 与 assembly 是否一一对应”的问题：可以按上述 6 个材料直接组织 TandemX 输入。它没有解决两个更高层的不确定性：论文称每个 accession 种植 20 株并收集幼叶，因此 raw 与 assembly 更可能是 accession-level pooled material；NCBI BioSample 的 cultivar 字段对二倍体为 `Wild`、四倍体为 `Cultivar`，真正品系名称在 `isolate`/sample name 中，不能把泛化字段误写为 cultivar 证明。Illumina-like run 实际是 DNBSEQ WGS；Hi-C 另有 `<accession>_hic` run，S245 和 S83 当前记录 reported bases=0，不应混入 short-read core 输入。

### 2. 苜蓿 Zhongmu-4（ZM4）近 T2T（P0）

**来源和 accession。** 新研究为 Sod et al., *Plant Communications* 2026，DOI [10.1016/j.xplc.2025.101691](https://doi.org/10.1016/j.xplc.2025.101691)，全文见 [PMC13084066](https://pmc.ncbi.nlm.nih.gov/articles/PMC13084066/)。论文将新 raw sequencing data 指向 NGDC [PRJCA041059](https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA041059)，将 assembly/annotation 指向 [Figshare 10.6084/m9.figshare.30768716](https://doi.org/10.6084/m9.figshare.30768716)。需要拆开三个 raw 来源：后续方法论文明确说明 ZM4 的 HiFi 与 Hi-C 来自 [PRJCA004062 / CRA005190](https://ngdc.cncb.ac.cn/gsa/browse/CRA005190)，ONT ultra-long 来自 [PRJCA030790](https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA030790)，Pore-C 来自 PRJCA041059。旧版 chromosome-level assembly 为 GWH [GWHBECI00000000](https://ngdc.cncb.ac.cn/gwh/Assembly/21836/show)。

| 材料 | 装配状态 | assembly accession | raw accession（按 modality） | 论文质量/大小 | donor-matched 判断 |
|---|---|---|---|---|---|
| *Medicago sativa* cv. Zhongmu-4 | **near-T2T**；4 个 haplotype，最终约 3,136 Mb，仍有 8 个 unresolved gaps；不能写成 gap-free T2T | 新装配：Figshare 10.6084/m9.figshare.30768716；旧版：GWHBECI00000000 | 新研究：PRJCA041059；旧/基础 HiFi、Illumina、Hi-C：PRJCA004062 / CRA005190；ONT UL：PRJCA030790（来源论文的 accession 交叉说明） | 新研究初始组装使用约 34× HiFi（>15 kb）和 35× ONT UL（>100 kb）；旧 GWH 组装 2,563,638,311 bp、32 条 pseudo-chromosomes、明确标为 chromosome-level | 材料匹配： cultivar-level 是；个体级 donor/extraction identity 未完全证明。旧研究从一株 5 年生 ZM4 植株取样，新研究沿用 ZM4 名称但应保留跨批次不确定性 |

**逐样本 accession 表（Medicago A17/R108，PRJCA042220）。** NGDC 的 [BioProject 页面](https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA042220) 和 [GSA CRA027262](https://ngdc.cncb.ac.cn/gsa/browse/CRA027262) 将两个材料的 assembly、BioSample、HiFi、ONT 和 Illumina/Hi-C 实验明确列出。GWH 的 [A17 assembly 页面](https://ngdc.cncb.ac.cn/gwh/Assembly/98366/show) 给出 `GWHGEXC00000000.1`、494,467,734 bp、HiFi 80.64×、ONT 137.02×；GWH 高级搜索结果给出 R108 的 [GWHGEXD00000000.1](https://ngdc.cncb.ac.cn/gwh/Assembly/98367/show)、415,266,007 bp、HiFi 113.23×、ONT 97.54×。GSA 项目总量为 246.99 GB；run-level file size 目前只在逐 run 页面明确暴露了部分记录，未把 coverage 换算值冒充 reported bases。

| 材料 / cultivar | assembly BioSample；assembly accession / file | BioSample：HiFi / ONT / Illumina Input / Hi-C | HiFi run；ONT run；Illumina run；Hi-C run | reported bases / bytes | 是否同一材料 |
|---|---|---|---|---|---|
| A17 (*Medicago truncatula*) | SAMC5586315；[GWHGEXC00000000.1](https://ngdc.cncb.ac.cn/gwh/Assembly/98366/show) / `GWHGEXC00000000.1.genome.fasta.gz` | SAMC5545447 / SAMC5545448+SAMC5545449 / SAMC5545454 / SAMC5545450 | CRR1954560.bam；CRR1954561.fq.gz + CRR1954562.fq.gz；CRR1954567_r1.fq.gz + CRR1954567_r2.fq.gz；CRR1954563_r1.fq.gz + CRR1954563_r2.fq.gz | GWH coverage（估算 bases）：HiFi 39.86 Gb、ONT 67.72 Gb；GSA 项目总 246.99 GB；逐 run bytes `metadata_blocked` | 是；GSA BioSample cultivar=A17，assembly 与四类实验均在同一 PRJCA042220/CRA027262；单株 DNA 链未额外提供 |
| R108 (*Medicago littoralis*，页面将其写作 R108/R10R) | SAMC5586316；[GWHGEXD00000000.1](https://ngdc.cncb.ac.cn/gwh/Assembly/98367/show) / `GWHGEXD00000000.1.genome.fasta.gz` | SAMC5545455 / SAMC5545456 / SAMC5545461 / SAMC5545457 | CRR1954568.bam；CRR1954569.fq.gz；CRR1954574_r1.fq.gz + CRR1954574_r2.fq.gz；CRR1954570_r1.fq.gz + CRR1954570_r2.fq.gz | GWH coverage（估算 bases）：HiFi 47.03 Gb、ONT 40.51 Gb；CRR1954568.bam **23,929.88 MB**；Hi-C CRR1954570 两文件 **22,239.22 + 24,528.92 MB**；项目总 246.99 GB；其余 run bytes `metadata_blocked` | 是；GSA BioSample cultivar=R108，assembly 与四类实验均在同一 PRJCA042220/CRA027262；单株 DNA 链未额外提供 |

这里的 Illumina run 是 `Input`（CENH3/染色质实验的 input control），不是独立的 whole-genome polishing library；如果 TandemX 只需要 DNA short reads，应把该列留作可选验证而不是核心输入。A17 的两个 ONT run 属同一 ONT ultra-long 样本，R108 为一个 ONT run；HiFi 文件在 NGDC 以 BAM 归档，不能直接把 BAM 压缩字节数当作碱基量。[CRA027262 页面](https://ngdc.cncb.ac.cn/gsa/browse/CRA027262) 是当前最完整的 accession 关系来源，A17/R108 raw run 因此已从 `metadata_blocked` 解析到具体 CRR accession。

该项目对应 Yi et al., *Molecular Plant* 2025，DOI [10.1016/j.molp.2025.07.016](https://doi.org/10.1016/j.molp.2025.07.016)，标题和 NGDC 项目描述均明确使用 **complete, gap-free T2T**。因此 A17/R108 可列为论文报告的 T2T；GWH 也将两个 assembly 标为 `Complete genome`。这与 ZM4 的 near-T2T 记录分开，不能因同属 *Medicago* 而混用术语。

**ZM4 的逐样本核实。** 旧 raw 项目 [PRJCA004062 / CRA005190](https://ngdc.cncb.ac.cn/gsa/browse/CRA005190) 的样本链已完全公开：`SAMC466290`（Zhongmu4_genome1，PacBio long read；CRR330257–CRR330260）、`SAMC466291`（Zhongmu4_genome2，Illumina genome short read；CRR330261）、`SAMC466292`（Hi-C；CRR330262–CRR330266）和 `SAMC466293`（RNA；CRR330267–CRR330269）。仓库报告该项目 **22 files / 267.60 GB**；其中 genome short-read 两文件为 **35,697.66 + 38,763.44 MB**，Hi-C 文件为 9,723.70/10,368、11,709.74/12,523.47、9,893.30/10,731.62、10,684.60/11,583.11、11,426.92/12,195.18 MB。PacBio 四个 run 的文件名已列出，但逐文件大小尚未在当前可访问页面闭合，记为 `metadata_blocked`。这些记录和 BioSample 页面都标注 cultivar Zhongmu-4、leaf，说明 **material-level matched**；旧 ZM4 与新 PRJCA041059/PRJCA030790 数据是否来自同一株、同一 DNA 提取，仍不能仅凭 cultivar 名称证明。

**预估下载量。** 旧 [CRA005190 页面](https://ngdc.cncb.ac.cn/gsa/browse/CRA005190) 报告 22 个文件、267.60 GB，包含 genome short reads、Hi-C 和长读段。新装配按 3.136 Gb × 34×/35× 换算，约为 106.6 Gb HiFi 和 109.8 Gb ONT UL 原始碱基；这两个值是估算，不代表 PRJCA041059 压缩文件总量。若只做 TandemX discovery/quantify，先请求或列出 accession 元数据，避免直接取完整 Hi-C/Pore-C。

**推荐用途。** ZM4 是最适合测试“多倍体 + haplotype-resolved + 高重复”的主数据，但应同时保留 `near_T2T` 标签和 8 个 gap 的坐标。旧 GWHBECI00000000 只能作为 chromosome-level 对照，不能作为 T2T truth。

### 3. 豇豆 HJD / FC6（P1）

原始研究为 Wei et al., “Complete telomere-to-telomere genomes of cowpea…”，*Horticulture Research* 2026，DOI [10.1093/hr/uhaf359](https://doi.org/10.1093/hr/uhaf359)，全文见 [PMC13102514](https://pmc.ncbi.nlm.nih.gov/articles/PMC13102514/)。其数据声明给出 HJD 的 NGDC 项目 [PRJCA044301](https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA044301) 和辅助资源 [Figshare 10.6084/m9.figshare.29878346](https://doi.org/10.6084/m9.figshare.29878346)。NGDC Genome Warehouse 检索可见 HJD 的 assembly accession [GWHHKIK00000000.1](https://ngdc.cncb.ac.cn/gwh/search/advanced/result?query_box=GWHHKIK00000000.1)，FC6 为 [GWHHKIL00000000.1](https://ngdc.cncb.ac.cn/gwh/search/advanced/result?query_box=GWHHKIL00000000.1)。

| 材料 | 物种/类型 | assembly accession / 装配证据 | raw reads accession 与量 | donor-matched 判断 |
|---|---|---|---|---|
| HJD | *Vigna unguiculata*，粮用自交系 | GWHHKIK00000000.1；论文 Table 1 列 548,737,058 bp、11 contigs、22 telomeres、0 gap；正文另写去除细胞器后 523.32 Mb，存在内部长度矛盾 | PRJCA044301；HiFi 37.33 Gb（约 73×）、ONT UL 61.01 Gb（约 119×）、Hi-C 32.12 Gb | 是，当前论文的 HJD assembly 与 HJD raw 属同一研究材料；仍应以 BioSample/experiment 元数据确认具体 DNA 提取 |
| FC6 | *V. unguiculata*，蔬菜型材料 | GWHHKIL00000000.1；论文 Table 1 列 542,521,633 bp、11 contigs、22 telomeres、0 gap | 论文说明取自前一研究：HiFi 28.59 Gb（约 57×）、ONT 67.27 Gb（约 133×）、Hi-C 24.10 Gb；当前论文数据声明未给 FC6 的 raw project accession，记为 `metadata_blocked` | 材料匹配到 FC6；raw 数据来自前一研究，跨论文的 donor identity 需单独核查 |

HJD 和 FC6 均可称为论文报告的 **complete T2T/gap-free assembly**；这与普通 chromosome-level 有明确区别。HJD 的 core raw 约 130.46 Gb 长读段，连同 Hi-C 约 162.58 Gb；FC6 core raw 约 119.96 Gb，但其 accession 尚未闭合。HJD 适合先跑，FC6 适合作为第二材料和 subspecies comparison。

### 4. 普通菜豆 YP4（P1）

原始研究为 Wang et al., *GigaScience* 2025，DOI [10.1093/gigascience/giaf001](https://doi.org/10.1093/gigascience/giaf001)，全文见 [PMC12077395](https://pmc.ncbi.nlm.nih.gov/articles/PMC12077395/)。GenBank/GWH 可核实到 assembly [GCA_051622355.1](https://ngdc.cncb.ac.cn/gwh/ncbi_assembly/3182086/show)，BioSample 为 SAMN41276676，原始项目为 NCBI [PRJNA1072282](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1072282)。

YP4 是 *Phaseolus vulgaris* 品种 Pinjinyun No. 4（红芸豆/菜豆类型）；论文称从山西太原一个单株取叶片做基因组 DNA。装配长度 560,297,700 bp，11 条染色体，0 gap，BUSCO 99.5%，QV 54.86；但端粒仅检出 20/22，作者给出 9 条 T2T pseudomolecule。因此本清单将其标记为 `T2T-claimed / partial-T2T-evidence`，不能把 11 条染色体全部宣称为 T2T。

| raw modality | accession | 论文报告量 | donor-matched |
|---|---|---:|---|
| PacBio HiFi | PRJNA1072282（逐 run 未在本次页面核实） | 31.75 Gb，约 56.67× | 是，基因组 DNA 来自该 YP4 单株 |
| ONT ultra-long | PRJNA1072282（逐 run 未核实） | 177.04 Gb，约 315.97× | 是，论文称用于同一 YP4 assembly；具体 BioSample 关联待元数据确认 |
| Hi-C | PRJNA1072282（逐 run 未核实） | 144.79 Gb，约 258.42× | 材料级匹配，具体库/叶片关联待核实 |
| RNA-seq（可选） | PRJNA1072282 | 118.68 Gb clean reads | 组织来自 YP4 材料，非 TandemX core raw |

核心基因组 raw 约 353.58 Gb（未压缩碱基量）；含 RNA-seq 约 472.26 Gb。该数据组适合验证较小二倍体基因组上的高覆盖 ONT/HiFi 读取，但下载前应从 NCBI SRA run table 拆出三种 DNA modality。

### 5. ZH13（中黄13）大豆（P1）

原始研究为 Wang et al., “The T2T genome assembly of soybean cultivar ZH13 and its epigenetic landscapes”，论文全文入口 [ScienceDirect](https://www.sciencedirect.com/science/article/pii/S1674205223003179)。NGDC BioProject [PRJCA015269](https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA015269) 标记为 *Glycine max*、monoisolate，并列有 HiFi、ONT、WGS 和 Hi-C 样本；GWH assembly 页面为 [GWHBWDJ00000000.1](https://ngdc.cncb.ac.cn/gwh/Assembly/66216/show)，BioSample SAMC1127443，GSA assembly raw project [CRA010060](https://ngdc.cncb.ac.cn/gsa/browse/CRA010060)，GSA raw run 记录中可见 CRR705247。

论文报告 ZH13-T2T 装配 1,007,237,669 bp、contig N50 48.76 Mb、39/40 个端粒、所有着丝粒，BUSCO 99.6%、QV 59.46；GWH 记录将其作为 Complete genome。严格使用时仍需本地检查 FASTA gap 和端粒两端，标签建议为 `T2T-claimed / archive-complete`。

| raw modality | accession | 覆盖/下载量估计 | donor-matched |
|---|---|---:|---|
| PacBio HiFi | PRJCA015269 / CRA010060；GWH 页面列 CRR705247，run-level modality 需确认 | 约 67×，按 1.007 Gb 估计约 67.5 Gb bases | 是，BioSample 标注 Zhonghuang 13，项目为 monoisolate；个体 DNA 提取链仍需 metadata 验证 |
| ONT | PRJCA015269 / CRA010060 | 约 117×，约 117.8 Gb bases | 同上 |
| BGI DNBSEQ Illumina | PRJCA015269 / CRA010060 | 约 53×，约 53.4 Gb bases | 同上 |
| Hi-C / 表观组 / RNA | PRJCA015269 下相应 BioSamples；不要混入 core DNA raw | 逐 run 量未在本次页面闭合 | 组织/实验匹配需按用途拆分 |

三类 core DNA raw 合计约 238.7 Gb bases，是较合适的 soybean pilot。ZH13 与后续其他大豆 T2T assembly 比较时，不能把不同 cultivar 当作 donor-matched truth。

### 6. YSD56 野大豆（P2）

原始研究为 Lian et al., *Scientific Data* 2025，DOI [10.1038/s41597-025-05741-y](https://doi.org/10.1038/s41597-025-05741-y)，全文见 [PMC12350797](https://pmc.ncbi.nlm.nih.gov/articles/PMC12350797/)。GenBank assembly 为 [GCA_040083835.1](https://www.ncbi.nlm.nih.gov/assembly/GCA_040083835.1/)，SRA/BioProject 为 [PRJNA1095640](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1095640)，研究数据集索引为 SRP502474。

YSD56 是对大豆胞囊线虫 X12 抗性稳定的 *Glycine soja* 野生种质。最终装配 1,008,523,555 bp、20 条染色体、40 个端粒、0 gap，BUSCO 99.7%、QV 52.16；论文明确使用 T2T reference genome。论文报告 Illumina 91.05 Gb、HiFi 44.19 Gb、ONT ultra-long 240.85 Gb，并在不同位置把 Hi-C 写成 23.84 或 123.84 Gb。2026-09-10 复核的 ENA run table 则给出 Illumina 92.173 Gb、HiFi 44.193 Gb、ONT 71.019 Gb、Hi-C 123.845 Gb。ONT 差异尤其大，可能涉及过滤或统计层级，但公开记录不足以闭合；实际分析必须以具体 run 文件的 checksum 和 observed bases 为准。

| raw modality | accession | ENA run-level bases（论文报告） | donor-matched |
|---|---|---:|---|
| PacBio HiFi | `SRR28726931` | 44.193 Gb（44.19 Gb） | 是，run 和 assembly 均关联 `SAMN40909152` / YSD56 leaf |
| ONT ultra-long | `SRR28726929` | 71.019 Gb（240.85 Gb；`source_conflict`） | 是，run 和 assembly 均关联 `SAMN40909152`；统计量不闭合 |
| Illumina WGS | `SRR28726932` | 92.173 Gb（91.05 Gb） | 是，run 和 assembly 均关联 `SAMN40909152` |
| Hi-C | `SRR28726930` | 123.845 Gb（23.84 或 123.84 Gb） | 材料级匹配；不进入 discovery/abundance 输入 |

ENA 当前三类 WGS genomic reads（HiFi、ONT、Illumina；不计 Hi-C/Iso-Seq）合计 207.385 Gb。该合计是 run metadata，不是本地已验证总量；目前只有 HiFi 完整文件通过本地 checksum 和 FASTQ 审计。YSD56 适合做野生大豆重复/着丝粒案例和抗性相关区域背景，但不要把它当作 ZH13 的 donor-matched replicate。

## 明确排除或降级的候选

- **一般 chromosome-level 苜蓿 ZM4 旧版 GWHBECI00000000**：可用于装配对照，但页面明确写 `Draft genome in chromosome level`，不能写成 T2T；新近 T2T 论文的 Figshare 资源应优先。
- **经典 chickpea ICC4958、普通菜豆 G19833/BAT93 等旧参考**：适合作为 fragmented/chromosome-level baseline，但本次目标是 T2T/近 T2T raw-read workflow，不列入主推荐组。YP4 论文将 BAT93 accession 列为 GCA_001517995.1，并将 G19833 作为旧比较装配。
- **HAAS216 野大豆（Scientific Data 2026）**：论文摘要报告 1,018.17 Mb、20 染色体、20 centromeres/40 telomeres、gap-free T2T，但本次核查尚未从其公开正文/NCBI 页面闭合 BioProject、assembly 和 modality-specific raw accession，因此保留为后续候选，不把它放入可立即运行的主清单。

## 下载前最小核查清单

1. 对每个 BioProject 导出 BioSample–Experiment–Run 表，确认 HiFi、ONT/ultra-long、Illumina 和 Hi-C 的 modality，不按文件名猜测。
2. 核对 raw reads 的 release 状态、文件大小、md5/sha256 和压缩格式；先只取 metadata，不下载大文件。
3. 下载前为每个材料保存 assembly accession、raw accession、BioSample accession、论文 DOI 和仓库页面快照。
4. 首轮只选择 S245、V14167、ZM4、HJD、YP4、ZH13、YSD56 中的 1–2 个材料做小规模 read sampling；完整 Hi-C/Pore-C 不属于 TandemX discovery/quantify 的第一步输入。
5. 对新装配先运行 assembly gap、端粒候选、染色体数和 FASTA 总长检查；任何论文称谓都不能替代本地 QC。

## 证据等级和术语约定

- `T2T / gap-free（论文报告）`：论文和仓库记录都支持无 gap/端粒证据；例如 HJD、FC6、YSD56，以及花生六材料论文的声明。
- `T2T-claimed / partial-T2T-evidence`：论文使用 T2T，但端粒证据只覆盖部分染色体；例如 YP4 的 9/11 条 pseudomolecule。
- `near-T2T`：仍有明确 unresolved gaps；例如 ZM4 新装配的 8 个 gap。
- `chromosome-level`：染色体锚定或 pseudo-chromosome 完整，但没有足够端粒/无 gap 证据；例如旧版 ZM4。
- `donor-matched` 只表示 assembly 与 raw reads 在 accession/BioSample 或论文材料层面对应；若没有明确的单株、DNA 提取和 library chain，写 `material-matched`，不写“同一 donor 已证明”。
