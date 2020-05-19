Best clustering:
```bash
Best: 0.989152 using {'activation': 'softmax', 'loss': 'categorical_crossentropy', 'optimizer': 'sgd', 'sparsity_l1__activity': 0.1, 'sparsity_l1__bias': 0.0, 'sparsity_l1__kernel': 0.0001, 'sparsity_l2__activity': 0.0, 'sparsity_l2__bias': 0.0, 'sparsity_l2__kernel': 0.0}
```

Notes:
* Here we tested L1 and L2 sparsity for activation, bais, and kernel.
* We used CP10k+1 ... using ln(CP10+1) gives a slight boost to 0.99
