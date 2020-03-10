
# Notes

We use conda for version control for command line programs and python pacakges.

Conda, however, does not work well with R packages. Instead, we use conda to install a minimal R environment and then use the `renv` package for R package versioning.

```R
# initialize renv
renv::init()

# update renv after intalling any packages
renv::settings$snapshot.type("simple")
renv::snapshot()

# add renv.lock, .Rprofile, and renv/activate.R to git
```
