
Here we used clusters with the following parameters:
* Filter cells > 80% MT
* BBKNN batch removal with default parameters
* Clustering with 0.5 resolution

The aim of this parameter sweep was to identify set activation, loss, optimizer.

```bash
Best: 0.928782 using {'activation': 'softmax', 'loss': 'categorical_crossentropy', 'optimizer': 'sgd', 'sparsity_l1': 0.001}
```
