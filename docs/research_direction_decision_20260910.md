# 三份外部讨论后的研究决策（2026-09-10）

## 结论

这些建议有助于确定问题，但不能作为现成算法方案或期刊录用预测。
当前最值得检验的是：家族代表中的局部共享序列能否产生高诊断
k-mer 深度，而完整串联上下文并不支持相同的家族丰度。
TXF000695 提供动机和已消费的开发案例，不是这一机制已被证明的证据。
本阶段先做机制可辨识性与公平基线设计，再决定是否改变生产估计器。

主任务承担科学决策和验收；Terra 完成竞争方法、输入和接口核查；
Luna 完成有明确边界的稿件编辑。并行任务不能自行选择阈值、消费新验证集
或把一种合理解释写成已验证机制。

## 对建议的取舍

| 建议 | 决策 | 理由 |
| --- | --- | --- |
| 减少防御性写作 | 立即落实 | 合并重复限定、移出内部施工记录；关键证据限制仍与结果相邻 |
| 增加 SRF 比较 | 补齐正式共同终点 | 已有四个成功端到端 pilot 和 32 次开发调用，不能写成从未比较；仍缺正式 HOR-aware 丰度验证 |
| 增加 Merqury/KAT、TandemTools、AniAnn’s | 按任务定位 | 分别对应全局 k-mer 评估、已组装阵列评估、assembly 注释；不强行套用家族发现准确率 |
| 用 read occupancy 替换 k-mer | 先作为透明基线 | SRF 已有 mapped-bp abundance，单独采用不是创新，也存在比对和模板偏差 |
| k-mer 与 occupancy 双证据一致即 high confidence | 不采用自动升级规则 | 共享 reads、catalogue 和归一化，误差可能相关；一致不等于独立验证 |
| 上下文纯度或稀有 seed 权重 | 进入受控机制实验 | catalogue 内稀有不等于全基因组特异；用同一 k-mer 招募和验证会形成循环证据 |
| 多代表、softmax、EM | 暂不并入主算法 | 多代表可能改变家族定义；分数 softmax 不是校准后验；EM 不能消除不可辨识性 |
| read bootstrap 置信区间 | 限定为条件采样误差 | 不能覆盖提取、平台、归一化、家族定义等系统偏差 |
| 继续修复组装以提高论文档次 | 不继续原恢复路线 | 当前有零个验证成功的 repeat-representation gain；需要新的独立假设才能重开 |

## 第一项方法实验的定义

目标量先限定为给定 operational catalogue 的 sampled read repeat-base
fraction。对 read i、family f 的接受区间先求并集，然后在不同 family
重叠处保留歧义。一个 read 可以包含不重叠的不同家族，每个 base 只能计一次。
总分母是全部输入 read bases，包括没有命中的 reads。命中 read 比例不是该量。
只有在明确采样与 genome-size 假设下，fraction × genome size 才转换为 bp 丰度。

开发对照必须分开检验：

1. 冻结的诊断 k-mer 估计器。
2. 同一冻结 catalogue 的普通竞争比对 interval-union occupancy。
3. 在 2 的基础上加入串联上下文过滤，保留未分配和歧义 bases。

保持 catalogue 不变可以隔离估计误差。另做完整 de novo 流程评估才可反映
发现误差；提供 truth catalogue 的分析不能替代端到端比较。
SRF 使用自己的 native motifs/HORs，是端到端竞争方法，不能被强迫使用
TandemX 的 catalogue 后仍称为 SRF 原生工作流。

受控机制因素包括：完整 tandem arrays、仅共享短片段的非串联背景、
近缘家族、真实混合家族 read、单位变异和测序错误。因素分别改变，避免同时
更换 detector、代表序列、归一化和 assignment 后无法解释收益来自哪里。
背景负对照必须参与 false-assigned-bp 评分。上下文证据需要与招募证据区分，
任何学习到的参数只使用开发分割；真实已观察家族不能成为新 holdout。

主要终点是 family-level repeat-bp error、错误分配的背景 bp、真实 repeat-bp
recall、歧义/不可估计比例及全部 eligible families 的覆盖率。不能通过多拒绝
困难家族得到较好的条件平均误差而称为全面提升。资源测量包括完整步骤的
时间和峰值内存；不预先承诺更快。各项结果按独立生成的源基因组汇总，
不能把同一基因组的多个条件当成独立生物学重复。

## 开始新验证前的门槛

- 先完成普通 occupancy 和 SRF native abundance 的公平评分接口。
  HOR/monomer 对应允许循环、反向互补及合法周期倍数，但跨家族 mosaic
  不能任意归给最有利的单一家族；歧义和未匹配输出必须保留。
- 在已消费开发数据上确认可检出条件和成本，再锁定参数。现有 SRF runner
  固定 k=151，不能把不存在的 k=101 参数命令写进预注册或直接运行。
- 开发通过后才定义新的独立 family/genome 分割、seed 注册、数值验收门槛、
  资源上限和一次性执行清单；本文件是研究决策，不是已完成的预注册。
- 5801–5803、6401–6403 已消费；7401–7403 拒用，3201–3203 仍保留。
  历史状态中的“untouched”不得覆盖后续消费记录。
- 若上下文方法不能优于普通 occupancy，或收益仅来自排除更多家族，就停止
  该创新主张。若只有共同偏差、参数不稳或资源成本不合理，也停止升级。
  阴性结果保留；冻结生产算法及旧验证结论保持原样。

## 本轮实际交付与未完成事项

已完成三项审计、SRF 执行接口核查以及 manuscript_v3 编辑稿。
v3 是同一证据的叙事调整，完整恢复负结果另存补充文件；它不代表算法性能
提升或论文科学等级已经上升。论文主图、完整补充包和引用仍须完成核对。
初始审计时未运行新 benchmark。用户随后明确要求同步开展准确度、速度和内存
研究：已完成结果一致的 Python 计数优化和独立的上下文机制开发实验，见
`performance_opportunity_audit_20260910.md` 与
`context_prototype_development_20260910.md`。未消费新 holdout、未安装新依赖，
生产估计公式不变；新上下文原型未接入生产。

后续顺序：评分接口与开发负对照 → 普通 occupancy 基线 → 单因素上下文
消融 → 独立验证 → 根据真实收益决定算法和稿件主张。重要设计改变、
基线完成、实测改善或失败、验证门槛结果均及时汇报。
Genome Research 是否合适取决于可推广的方法收益或更充分的独立生物学发现；
不能从软件功能数量、术语包装或三份外部评价推导出来。

Bioconda、Zenodo、正式 tagged release 继续暂停，unitFinder 永久停止。

依据：`competitor_gap_audit_20260910.md`、
`occupancy_input_audit_20260910.md`、
`srf_benchmark_execution_contract_20260910.md`、
`manuscript_revision_audit_20260910.md`；
[SRF 官方方法和 abundance 工作流](https://github.com/lh3/srf)。
