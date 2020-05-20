
Here we used clusters with the following parameters:
* Filter cells > 80% MT
* BBKNN batch removal with default parameters
* Clustering with 0.5 resolution AND 1.0

The aim of this parameter sweep was to identify sparsity l1 and l2.

Best clustering with 1.0 resolution:
```bash
Best: 0.993051 using {'activation': 'softmax', 'loss': 'categorical_crossentropy', 'optimizer': 'sgd', 'sparsity_l1__activity': 0.0001, 'sparsity_l1__bias': 1e-10, 'sparsity_l1__kernel': 0.0001, 'sparsity_l2__activity': 0.0, 'sparsity_l2__bias': 1e-06, 'sparsity_l2__kernel': 0.0}
```

Notes:
* Here we tested L1 and L2 sparsity for activation, bais, and kernel. The best setttings where when spasity set the same across all variables.
* Not shown: we did a test of CP10K vs ln(CP10k+1) and the latter worked slightly better (max accuracy 0.98 vs 0.99).
