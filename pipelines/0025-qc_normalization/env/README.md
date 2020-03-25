
# Notes

## General

We use conda for version control for command line programs and python pacakges.

## R

Conda, however, does not work well with R packages. One could use conda to install a base version of R. Alternatively, one can pull from from a rocker image (e.g., `rocker/r-ver:3.6.2` instead of `buildpack-deps:curl` or a minimal install like `alpine` ).

For R package management we use `renv`.

```R
# initialize renv
renv::init()

# update renv after intalling any packages
renv::settings$snapshot.type("simple")
renv::snapshot()

# add renv.lock, .Rprofile, and renv/activate.R to git
```
