# The _bmers_ package

The _bmers_ package (Bayesian Mixed Effects Regressions fit with Stan) fits Bayesian mixed effects regressions with weakly informative priors and maximal random effects structures in _Stan_ with the No U-Turn Sampler (through the _RStan_ interface).  The maximal random effects structure is determined automatically, and only the fixed effects structure and the random effects grouping factors need to be specified.  Continuous variables are automatically scaled, and factor contrasts are automatically set; neither of these needs to be done beforehand.  NA values in both the fixed effects and random effects are supported through the use of sum contrasts with NAs set to zero in the model matrix.

## Installation

Prior to installing *bmers*, the *rstan* package must already be installed and functioning (see instructions here: <https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started>), and the *devtools* package must also already be installed; to install *devtools* call `install.packages(devtools)`.  Then, to install *bmers*, call:

&nbsp;

`devtools::install_github("CDEager/bmers")`

&nbsp;

A vignette is provided which shows how to use the major functions **bmer**, **build_bmer_model** and **fit_bmer_build**, and also demonstrates some of the basic querying functions for **bmerBuild** and **bmerFit** objects.

&nbsp;

**NOTE:** If you fit a model using **bmer** or **fit_bmer_build** which requires a new **stanmodel** to be compiled (if warmup/sampling does not begin for a minute or so after you enter the **bmer** / **fit_bmer_build** call, then *rstan* is creating a new **stanmodel**), and then try to fit another model which uses this same **stanmodel** in the same R session, the session will crash.  I will add a patch that fixes this.  For the time being, after you fit a new type of model (one which requires the compilation of a new **stanmodel**), save your workspace, close the current R session, and reopen a new session.  From there you should be able to fit that type of model as many times as you want within the same R session without issue, and without needing to recompile that type of model.

&nbsp;
