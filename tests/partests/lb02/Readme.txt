This test checks basic parallel load balancing. Three sub-domains.
The workload is artificially perturbed on partition 0 in the first 
solution step to enforce load imbalance and force load balancing 
engine to redistribute the load. The amount of work being balanced
is, however, dependent on actual computing environment.

This test requires parallel oofem build with petcs and parmetis modules.
