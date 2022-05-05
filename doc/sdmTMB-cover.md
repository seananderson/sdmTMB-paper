---
name: Sean C. Anderson
title: Research Scientist
opening: Dear Editorial Board
closing: Sincerely
signedname: Sean Anderson on behalf of coauthors
email: sean.anderson@dfo-mpo.gc.ca
fontsize: 10pt
topmargin: -20pt
textheightextra: 65pt
margin: 1.50in
output:
  pdf_document:
    template: template.tex
    keep_tex: true
---

We are pleased to submit our manuscript "sdmTMB: an R package for fast, flexible, and user-friendly generalized linear mixed effects models with spatial and spatiotemporal random fields" for consideration as Research Article in *Methods in Ecology and Evolution*.

Ecological data are frequently collected over space and often also at discreet points in time.
Such data often suffer from spatial and/or spatiotemporal correlation.
An increasingly common way to model these data is through generalized linear mixed effects models (GLMMs) with Gaussian random fields.
In addition to accounting for spatial correlation, random fields can represent the combined effects of latent variables and estimate meaningful quantities describing the spatial correlation.

In this manuscript, we introduce the R package sdmTMB for fitting GLMMs with Gaussian random fields for spatial and or spatiotemporal latent variables.
The package uses the SPDE (stochastic partial differential equation) approximation to Gaussian random fields for computational efficiency, as popularized by the INLA R package.
However, unlike INLA, sdmTMB uses the TMB R package to calculate the marginal log likelihood, integrating over random effects with the Laplace approximation for rapid model fitting.
It is also possible to pass a fitted model object to the tmbstan package to conduct full Bayesian inference with Stan.

We believe sdmTMB fills an important niche in statistical ecology: a fast and flexible R package for fitting spatial/spatiotemporal random field GLMMs for users familiar with the frequently used lme4, glmmTMB, or mgcv packages.
In our experience, existing tools, although powerful, prove challenging for many applied ecologists to use.
We provide comparisons with the most closely related packages VAST and INLA/inlabru in the appendices.
In addition to user friendliness, we show that sdmTMB is faster than non-TMB-based alternatives and offers a combination of features not available in other packages such as penalized smoothers, break-point effects, and anisotropy (where correlation decays faster in one direction than another).

Our manuscript has not been previously published, nor is it under consideration elsewhere.
We have posted a preprint version of our manuscript on bioRxiv at <https://doi.org/10.1101/2022.03.24.485545>.
We look forward to constructive feedback from your reviewers and editors.
