
Here we used clusters with the following parameters:
* Filter cells > 80% MT
* BBKNN batch removal with default parameters
* Clustering with 0.5 resolution

The aim of this parameter sweep was to identify sparsity l1 and l2.

```bash
Best: 0.925153 using {'activation': 'softmax', 'loss': 'categorical_crossentropy', 'optimizer': 'sgd', 'sparsity_l1': 0.0001, 'sparsity_l2': 0.0}
```
