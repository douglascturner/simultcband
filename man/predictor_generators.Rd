\name{predictor_generators}
\alias{predictor_generators}
\title{Predictor generators (pgens)}
\description{Documentation for creating a custom predictor generator (pgen)}
\details{
simultcband ships with pgens that will fit the needs
of many users. If you need to create a custom pgen in order to use
an underlying predictor not implemented by simultcband, this help file
gives the information for how to do that, as well as the intuition of
the design.

The goal of the pgen is to create a function that in turn
returns the predictions and optionally standard errors. The reason for this
pgen layer design is so that some prediction parameters used
in the bootstraps (e.g., smoothing parameters, number and location of knots)
can be functions of the original data and thus fixed across the bootstraps.  As
an example, for a Generalized Additive Model, suppose that the generalized
cross-validation (GCV) smoothing parameter (sp) is 0.3. Then, there are two
predictor functions that would produce the same predictions of the original
data consistent with this prediction:
  (1) the predictor that uses GCV.
  (2) the predictor with sp fixed to 0.3.
This framework allows for, e.g., the smoothing parameter to be chosen using
GCV, and then fixed throughout the bootstraps.

A predictor generator returns a *function* that does the following:

- Takes as arguments "data", "newdata", "se.fit", and "use_boot_weights".

- Returns a list with two elements: "fit" and "se.fit".

If a predictor does not know how to produce standard errors, it can set the
returned se.fit to NA; this way, if a user attempts to use a cband method that
requires standard errors, a user-friendly error will be given. Alternatively,
the predictor can give an error directly if the se.fit arg is TRUE.

When "use_boot_weights" is TRUE, the data contain a variable named
"WEIGHT_FOR_BOOT". The name `WEIGHT_FOR_BOOT` can thus be used for what is
usually the "weight" argument of the underlying predictor. Note that there is
usually non-standard evaluation of this argument (e.g., in the case of
mgcv::gam() and lm()). See the documentation (e.g., ?lm) of the predictor for
more information.
}
