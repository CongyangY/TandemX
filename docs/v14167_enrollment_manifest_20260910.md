# V14167 首个花生 enrollment manifest 说明

配置文件：[benchmarks/configs/peanut_v14167_enrollment_manifest.json](../benchmarks/configs/peanut_v14167_enrollment_manifest.json)。该 manifest 是 **metadata-only freeze**：截至 2026-09-10 未下载 FASTQ/BAM/FASTA，也未运行 TandemX 分析。

V14167（*Arachis duranensis*，AA 二倍体）被选作花生首个入组材料。公开装配为 `GCA_054824555.1`，论文报告 1.178 Gb、0 gap 和 20 个端粒。逐样本 raw 映射为：PacBio HiFi `SRR33996330`（211.520 Gb bases）、ONT `SRR33996240`（108.019 Gb）和 DNBSEQ WGS `SRR33996449`（131.114 Gb）。这些 run 分别关联 V14167 的 raw BioSamples `SAMN49012632`、`SAMN49012638` 和 `SAMN49012650`；assembly BioSample 为 `SAMN48356820`。研究项目入口为 [PRJNA1259747](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1259747)，装配入口为 [GCA_054824555.1](https://www.ncbi.nlm.nih.gov/assembly/GCA_054824555.1/)，论文 DOI 为 [10.1038/s41588-026-02577-z](https://doi.org/10.1038/s41588-026-02577-z)。

donor-match 采用 `A2_material_level_pooled`：assembly isolate、材料名和三类 raw run 在 accession/BioSample 层面闭合，但论文说明每个 accession 取 20 株幼叶，公开信息没有证明单株 DNA 提取和每个 library 的单一 donor 链。因此 manifest 允许 V14167 作为材料级 matched 输入，不把它写成单株 donor 已证明。

首轮 pilot 目标是相对 assembly 的 HiFi 5×、ONT 2×、DNBSEQ 2×，约 10.61 Gb bases；每种 modality 使用独立但固定的 seed（914167、914168、914169），做两个重复。该深度和 seed 是建议值，未被消耗，也没有生成抽样文件。Hi-C 不进入首轮 discovery/quantify 核心输入。完整入组必须等两个 pilot 重复都通过 accession、BioSample、平台、压缩格式、非零 payload、observed base count、checksum、读长摘要和资源预算检查。

硬停止条件覆盖 accession/material 不匹配、受限或空 payload、平台/格式错误、observed bases 与 reported bases 偏差超过 10% 且无来源解释、checksum/文件大小改变、paired-end 不完整、assembly accession/checksum 不符，以及 pilot 资源超预算 20% 以上。assembly 的本地 N 数、总长或序列数异常属于 review stop，需要保留原始报告并降级状态；不能静默修改 manifest。所有停止规则都只针对数据入组和 QC，不把组装当作绝对 copy-number truth。

该配置不改变代码和 manuscript；它只冻结 V14167 的公开 accession、证据等级、建议 pilot 参数、输入校验字段和停止规则，供后续实际下载前审阅。
