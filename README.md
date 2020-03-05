
# Description

Bits of code from Leland's personal projects.

This repo is pulled from a default template for workflows.

# Workflow template setup

## lib

* The `lib` directory contains general libraries that may be referenced by multiple workflows, for instance cromwell configs and python configs.


## modules

* Each module is a full analysis. Think of it like the heading of a methods section in a paper. For instance if this were genetic summary statistics workflow, a module might be "fine-mapping" that does both conditional and credible set analysis. Another module may be "colocalization".

* Modules may have numbers prior to their name (e.g., `example_module_1` to `0025-example_module_1`). These numbers do not mean anything, but merely used to keep modules in their general order of execution. These are optional.

* A module consists of :
1. A workflow, for instance Cromwell (see `0025-example_module_1`) or Snakemake (see `0033-example_module_2`).
2. A `scripts` directory with *all* scripts referenced by that workflow (unless a general lib script is called).
3. A `docs` directory that contains a documentation of the default parameters written in a style that is publishable as methods in a paper (including citations). Within the `docs` directory there may be a `reference` with any additional reference materials.
4. An `example_runtime_setup` directory contains files that give an example of actual config files and any other files used to run the pipeline.


## studies

* A studies directory should either exist within the workflow repo or be a separate repo that has the same name as the workflow repo, but with `studies` appended to it (e.g. `template-workflow` becomes `template-workflow-studies`).
* If there is a standard set of plots that will always look the same way, a module should generate such plots. Otherwise, all code to analyze the results of a module run should be in the `studies` directory. For instance if this were genetic summary statistics workflow, `studies` may contain a `t2d` directory and a `weight` directory.
* Within a study is either an Jupyter notebook (either python or R kernel) or an R markdown file. Nearly all plots / analysis of the results of running the various modules should be done in the notebook / markdown file.
* A study may also contain a scripts directory with scripts to aggregate data for a one off analysis (if the analysis is going to be repeated, consider making a new module or adding it to an existing module) or for special plots that cannot be done in the notebook / markdown file.


# New workflow reminders

- [ ] Documentation
- [ ] Environment version control
- [ ] Module version control
- [ ] Git branches
- [ ] Code review


# Documentation

Be sure to document your code!


# Environment version control

Analysis environment is controlled using conda. Each module should have an `environment.yml` file with all of the packages used. If a required package or library is missing from conda (and therefore not in the `environment.yml`), it should be noted in the `README.md` of the module.

```bash
conda env export --no-builds | grep -v prefix | grep -v name > environment.yml
```


# Module version control

Each module within this workflow uses [bumpversion](https://pypi.org/project/bumpversion) for automatic [semantic versioning](https://semver.org).

```bash
# bump the appropriate increment
bumpversion patch --verbose --dry-run
bumpversion minor --verbose --dry-run
bumpversion major --verbose --dry-run

# commit with tags
git push --tags
```


# GitHub forks

Forking the repository allows developers to work independently while retaining well-maintained code on the master fork. For instructions on how to fork, follow the [Fork a repo](https://help.github.com/en/articles/fork-a-repo) instructions.

After forking the repo, clone the repo to your local desktop:

```bash
# to use SSH
git clone git@github.com:<username>/template-workflow.git

# to use Https
git clone https://github.com/<username>/template-workflow.git
```

This creates a replica of the remote repository on your local desktop. *Note*: When you create your local repository, it will also make a local clone of the remote repository (typically as ```origin```). So, your local master branch would simply be ```master```. But, your remote master branch will be ```origin/master```. You can also add multiple remote repositories. For instance, let us say our main repository is under the remote repository ```my_repo```. We will want to add it as a remote repository, so we can fetch the most up-to-date code. You could add it by:

```bash
# Add the my_repo remote repo to your local desktop -- this will allow you to pull and push to branches on the my_repo repository
git remote add my_repo git@github.com:my_repo/template-workflow.git
```


# Git branches

Branching is how git actually tracks code development. For more information, see the [Git Branch Tutorial](https://www.atlassian.com/git/tutorials/using-branches) on Atlassian. If you want to add a new feature, module, or fix a bug, a common work flow would look like this:

```bash
# Update your local copy of the master branch to make sure you are getting the most up-to-date code
git fetch my_repo master

# Create a new branch based on your feature from the recently pulled master branch
git checkout -b sc-atac-seq my_repo/master
```

As you develop, you want to commit your work to your branch, so you don't lose it all if something happens!

```bash
# Confirm we're on the right branch
git branch

# Add all your work to be tracked (Note: there are many ways to add specific files, etc. See https://git-scm.com/docs/git-add for more information). The following command adds everything in your currently directory.
git add .

# Commit your work to the branch with a message describing what's in the commit
git commit -m "Created the scATAC-seq pipeline!"

# When we are done with our feature, or just want to push it to GitHub for safe keeping:
git push origin HEAD # HEAD pushes everything up to the most recent commit
```


# Code review

Create a [GitHub Pull Request](https://help.github.com/en/articles/creating-a-pull-request). A PR allows other developers a chance to go through and comment on lines of code they believe can be improved. In addition, it will tell you if the code you are trying to merge into the ```my_repo``` branch actually conflicts with code that already exists in the branch, so you don't overwrite someone else's work.

Once another developer approves the PR, you have the go-ahead to merge your code! Congrats, you finished your feature!

*Note*: There are some cases where you may just want to push directly to the my_repo fork, thereby avoiding code reviews. For instance, if you're working on a one-off project that you want people to be able to see, but no one else is necessarily working on, you can always push directly to the branches on my_repo fork. Or, you could also still go through the steps of a PR, but simply merge your own code without CR. 
