# Novelty and claim audit

The relevant standard includes existing Genome Research methods, not only simple
repeat scanners. SRF was published in Genome Research in 2023 and already reports
satellite repeat units/HORs from accurate reads or assemblies, including cases of
assembly under-representation ([primary paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC10760446/)).
A read-first interface or a plot of reads versus assemblies therefore does not,
on its own, establish TandemX's methodological novelty.

## Current evidence and unproven hypotheses

| Candidate contribution | Actual current evidence | Required next evidence |
| --- | --- | --- |
| Indel-aware, multiple-array read discovery | Independent Python/Rust kernels; corrected diagnostic failures; 16 simulated scenarios | Wider unseen sequence/error distributions, matched output granularity, comparative speed/memory scaling |
| Explicit operational monomer resolution | Candidate sequences, fixed representatives, assignment bounds and alternatives; a related-unit sensitivity sweep | Fragmentation/merging errors over diverse variant structure; sensitivity independent of simulator labels |
| Joint abundance and assembly-under-representation inference | Toy downstream commands and an SRF workflow comparator exist | Controlled coverage/copy number/collapse truth, bias and interval coverage, matched real plant reads/assemblies |
| Probe candidate specificity | Toy scoring exists | Independent sequence-based off-target evaluation, uncertainty and documented experimental evidence when available |
| Optional learned confidence | Not implemented; no AI novelty claim | Grouped held-out evaluation, transparent non-ML calibration baseline, ablations, calibration and domain shift |

This table expresses research hypotheses, not completed results. Any AI component
must improve a predefined accuracy/calibration/cost endpoint and remain optional
when its calibration domain does not cover the input. Adding a pretrained model
or a conversational interface is not an evaluated contribution.

The mature-paper gate remains open: real biological ground truth, independent
sample replication, fully paired comparator experiments, statistical uncertainty,
complete multipanel main/supplement figures and an install/reuse audit are still
needed. Perfect recovery on diagnostic seeds does not satisfy this gate.
