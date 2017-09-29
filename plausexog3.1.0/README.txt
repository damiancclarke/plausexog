===============================================================================
README.txt for plausexog code version 3.1.0

yyyy-mm-dd: 2014-08-16 [Original]
yyyy-mm-dd: 2017-09-29 [Update]
contact: damian.clarke@usach.cl
===============================================================================

This folder contains the ado plausexog.ado and the help file plausexog.hlp whi-
ch run Conley et al's (2012) Plausibly Exogenous inference routines in Stata.  
The heart of this code comes from Hansen's original code, however this has been
augmented to include additional error capture, functionality with missing data
and factor variables/wildcards, automated graphing, and additional documentati-
on.

From version 2.0.0 onwards, and additional feature has been added to use any di-
stribution for the LTZ approch instead of assuming a Gaussian prior (as describ-
ed in Conley et al. 2012, page 265). This requires the user written ado mvtnorm 
for simulating multivariate normal draws.

