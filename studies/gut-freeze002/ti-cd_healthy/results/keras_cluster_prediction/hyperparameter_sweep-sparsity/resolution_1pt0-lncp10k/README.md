Best clustering:
```bash
Best: 0.993051 using {'activation': 'softmax', 'loss': 'categorical_crossentropy', 'optimizer': 'sgd', 'sparsity_l1__activity': 0.0001, 'sparsity_l1__bias': 1e-10, 'sparsity_l1__kernel': 0.0001, 'sparsity_l2__activity': 0.0, 'sparsity_l2__bias': 1e-06, 'sparsity_l2__kernel': 0.0}
```

However, keeping a simple model where l1 sparsity all the same and l2 spasity is 0 is basically just as good (mean test score = 0.9928568). We will use these settings since it is simpler:
```R
dat = read.csv("keras_model-grid_search-grid_result.tsv", sep = "\t")
subset(dat, param_sparsity_l1__activity == 1e-04 & param_sparsity_l1__bias == 1e-04 & param_sparsity_l1__kernel == 1e-04 & param__sparsity_l2__kernel == 0.0 & param__sparsity_l2__activity == 0.0 & param__sparsity_l2__bias == 0.0)
```

Notes:
* Here we tested L1 and L2 sparsity for activation, bais, and kernel. The best setttings where not too far off when spasity set the same across all variables.
* We used ln(CP10k+1)
