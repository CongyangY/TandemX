# 本轮决策：A3 开发候选未显示足够整体优势；保留精确计数优化

基线固定为 `81827c3`。准确度和效率分别在独立 Git 分支上运行，未合并到
正式开发分支。正式 manuscript v3、生产 quantify 及 k-mer 实现与 baseline
逐字节一致。Bioconda、Zenodo、tagged release 继续暂停。

## 准确度：两轮结束，未进入 held-out

比较的是两个开发源种子、十类条件、2/5/10x 覆盖，共 60 个数据集的
120 个 family-condition 行；不是 120 个独立生物学样本。
全部家族包括歧义和零预测均保留。指标的分母和人工 source truth 见
`phase_gate_protocol_20260910.md`。A3 synthetic FPR 的分母、truth 和阈值
也与历史 classifier FPR 不可直接比较；完整调和见
`fpr_scope_reconciliation_20260910.md`。

| 方法 | 开发 MARE | indel recall | clean MARE | background excess bp | synthetic FPR |
| --- | ---: | ---: | ---: | ---: | ---: |
| A0 frozen diagnostic k-mer | 0.20697 | 不提供区间 | 0.03015 | 37,757 | 0.0667 |
| A1 competitive alignment occupancy | 0.05683 | 0.90895 | 0.02486 | 0 | 0.1167 |
| A2 strict periodic context | 0.24549 | 0.10788 | 0.04569 | 0 | 0.1167 |
| A3 第一轮 local re-phasing | 0.08224 | 0.89816 | 0.04569 | 0 | 0.1167 |
| A3 第二轮 short-gap inclusion | 0.07469 | 0.91656 | 0.04569 | 0 | 0.1167 |

本开发夹具中的正向观察是：局部 re-phasing 将 indel recall 从 strict A2 的
0.10788 提升至第二轮 A3 的 0.91656；A3 在 `background_homology` 条件的
background excess 为 0，而 A0 为 37,757 bp。第二轮仅纳入已接受相位链
内部的短间隙，所有阈值保持不变。

未获支持的假设：A3 的丰度准确度优于普通 competitive mapping。最终 MARE
仍比 A1 高约 31.4%，平均计时约为 A0 的 2.77 倍，clean 条件误差也高于
A0。因此这两轮开发结果没有显示足以支持继续扩展或进入独立验证的整体
优势。历史 `binary_FPR_at_most_002` gate 因其 FPR 与历史 classifier FPR
的定义不可通约而撤回；原始 receipt、代码和数值均保留不改。

消融还显示：sequence-only 与单独 phase 步骤 MARE 同为 0.1984；加入 phase
coverage 后为 0.0699，加入严格 recurrence/order 后反而为 0.2309；re-phasing
恢复到 0.0822。specificity weight 的有/无在此开发集没有改变结果。
因此不能把所有组件都描述为有效贡献，更不能用术语组合代替方法增益。

**本轮不再扩展 A3 accuracy branch：不增加 HMM、EM、detector、repair 或
第三轮修改。** 未生成、未查看 final held-out 数据；没有独立验证成功率可
报告。这是一个尚未显示足够整体优势的开发候选，不是否定所有 phase-aware
abundance 方法。SRF、mm2-ivh、MCS、TideHunter 等先例见
`novelty_matrix_phase_gate_20260910.md`；本轮不命名新算法、不声称首创。

## 效率：真实输入仍有收益，GB 级尚未验证

每个真实输入上旧/新各三次串行、交错、fresh-process 测量。计时覆盖完整
quantify 调用（含 read parsing），不含解释器启动、导入及额外源文件哈希。

| 输入 | 旧/新 wall 中位数(s) | 倍数 | 旧/新 CPU 中位数(s) | 旧/新 peak RSS 中位数(MB) |
| --- | --- | ---: | --- | --- |
| 玉米 Mo17，11,680,888 bases | 15.045 / 4.170 | 3.61 | 14.982 / 4.150 | 63.65 / 60.70 |
| Ey15，113,336,152 bases | 199.995 / 93.914 | 2.13 | 198.494 / 92.576 | 204.91 / 194.15 |

每个输入的六次 `copy_number.tsv` SHA 完全一致。两组最终输出树大小分别
为 103,755 和 588,166 bytes，新旧相同；这是保留输出大小，不是临时磁盘
峰值。RSS 的小幅差异不足以支持显著节省内存，尤其没有达到 40% 减少标准。

Mo17 是 FASTA、Ey15 是 FASTQ.gz，且 catalogue 不同；两者不是控制其余变量
的单一 scaling 曲线。计时期间同机曾有较小的 accuracy 计算，故记录的是
该主机状态下的结果，而非完全隔离主机测量。Ey15 配置中的 genome_size
1.76 Gb 是本次 performance-only 参数：显式 haploid_depth=1 使其不参与
归一化，不能当作 Ey15 的基因组大小或生物学丰度分析。

效率阶段约 21 分钟完成，用户给定约 30 分钟上限。1 Gbase 旧/新各一次
预计至少还需约 27.5 分钟，无法装入剩余预算，故没有启动。完整 wheat
workload、GB 级重复性能、分配热点 profile 和 fused multi-k 均未执行。
之前的错误元数据中止尝试单独保留，不算有效测量，也不算算法失败。

**保留 baseline 中已有的精确 rolling-code 优化，定位为实现优化。**
100-Mb 级有 >=2x 和输出一致的证据；不能据此填写 GB 级或最终论文性能数字。
没有新的生产算法变更需要合并。

## 归档、验证与下一步边界

大文件、每一条件和原始命令保存在 T7：

- `results/phase_coherence_gate_development_r1_20260910`
- `results/phase_coherence_gate_development_r2_20260910`
- `results/rolling_target_count_scaling_gate_v1_20260910`

仓库小型归档在 `benchmarks/evidence/phase_and_scaling_gate_20260910/`，包括
协议、全部方法摘要、效率原始重复值、无效尝试、源路径/哈希和 baseline
逐字节核验。A0 最终评分直接重读 native `estimated_bp`；与第一轮重构
数值的最大差异为 0 bp，未改变任何判定。第二轮实际执行源码快照另存 T7。

完整回归：769 passed、2 skipped；之后新增的短间隙边界检查所在测试文件
13 passed。效率独立 worktree 38 个聚焦测试通过；Rust parity 因该 worktree
不含 editable native extension 未执行，不能称其通过。主 checkout 的完整
回归包含原有 native 测试环境。

所有子任务已收尾。下一阶段只保留正式 SRF/ordinary-mapping 统一终点比较
的准备和既有证据投稿整理；需要最终性能数字时再做明确排期的 GB 级重复。
本轮不启动这些后续工作、不修改 Abstract 或核心 claims，也不进行发布。
