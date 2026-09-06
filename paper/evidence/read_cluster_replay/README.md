# Read-cluster uncertainty development replay

Committed source 6f29db5bc3a77db834239d38e91925f80e6b6018 replayed the same 27
development abundance inputs (81 family conditions), preserving source and input
hashes. All executions completed; no planted error rate or read-coordinate truth
entered the estimator. The original baseline remains separate and unchanged.

Only 33/81 intervals met the support guard; 24/33 included truth. The other 48
intervals are missing, not [0,0]. At 20x, coverage among seven available intervals
was 7/7 without substitutions, 6/7 at 0.1% and 4/7 at 1% substitutions. At 5x it
was 3/4, 3/4 and 1/4 respectively. No 1x condition had enough effective support.
These small, correlated development counts are not a 95% calibration claim.

At 20x/1% substitutions mean signed CN error remains -17.97%, with
estimator-minus-sampling-oracle -19.74%. Read-cluster uncertainty alone does not
correct error attenuation or biological unit divergence. The reference model
therefore remains separate from the default CLI. Results include every interval
state and all point-estimate errors, without selecting only available intervals.
