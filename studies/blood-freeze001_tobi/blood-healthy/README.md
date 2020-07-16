# Single-cell RNAseq workflow for healthy blood

Data is 0h, 24h and 48h. This workflow was modelled on the existing workflow found [https://github.com/andersonlab/sc_nextflow/tree/master/studies/gut-freeze002/ti-cd_healthy](here).

## Steps taken
I copied the `example_run.sh` script but modified the variables to use `PROJECT_HOME` instead of `HOME`. `PROJECT_HOME` can be defined as the user runs the script e.g. 
```
PROJECT_HOME=/nfs/users/nfs_o/oa3/2020-06-11_rotation-blood-work/ bash example_run.sh
```
Then for each of the following options I made these changes:
```
--file_paths_10x [Created samples file with 3 samples, full names and paths follow by just number suffixes as short names]
--file_metadata [Unmodified apart from change to PROJECT_HOME path]
--file_sample_qc [Left as is but am skeptical about the values for pct_mito and n_genes]
--genes_exclude_hvg [Left as is]
--genes_score [Left as is, should probably customise for blood data?]
--output_dir [Left as is]
-params-file [Removed most marker gene datasets]
```

