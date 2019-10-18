# The code below can be used to reproduce the results of
# An information-theoretic approach for selecting arms in clinical trials
# by P. Mozgunov and T. Jaki (2019)

# wde-main-code.R
Implements the WE-design (using either Select-the-Best or Randomization Rules) that uses the leading term of the asymptotics

# Phase-2-Run.R
Provides the results of the WE Design for the trial with binary endpoints
and corresponds to Section 4.3

# wde-numerical-code.R
Implements the WE-design with the minimum efficacy bound and the numerical intergration for the information gain
computation.

# Phase-2-Min-Efficacy-Run.R
Provides the results of the WE Design with the minimum efficacy bounds for the trial with binary endpoints and corresponds
to Section 4.4

# wde-co-primary-code.R
Implements the WE-design (using either Select-the-Best or Randomization Rules) for the co-primary efficacy endpoints (using either Select-the-Best or Randomization Rules) 
for the trial with two binary endpoints and uses the leading terms of the asymptotics

# Phase-2-Co-Primary-Run.R
Provides the results for WE design for co-primary endpoints and corresponds to Section 5.

