__References:__
- Results contained in [Universal location of Yang-Lee edge singularity in classic O⁡(𝑁) universality classes](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.107.116013).
- More details are contained in my [thesis](https://repository.lib.ncsu.edu/items/ab3dd7a0-9955-4f1b-8be9-adb25e4cfcfd).
- See also [Universal Location of the Yang-Lee Edge Singularity in O⁡(𝑁) Theories](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.125.191602) for first computation which uses a different method.

__Contents:__

1. ONvars_YZ
  - Outputs list of Taylor expansion variables for given truncations. Imported by following files.
2. ON_YLE_BetaFunctions
  - Starting from the $O(\partial^2)$ SDE expressions, computes the simplified 3D flows of Taylor expansion coefficients in the YLE expansion for both $\rho$ and $\phi$ parametrizations.
  - Computed locally. Output to file and results are used for later computations.
  - Arbitrary d (dimension) beta functions computations/expressions not used.
3. Scaling_sol_solver
  - Using flows from ON_YLE_BetaFunctions, computes 3D scaling solutions for desired range of $N$ and regulator scale parameter $\alpha$.
  - Determines PMS (principle of minimal sensitivity), $\alpha_{\text{PMS}}$, locations of exponents.
  - $\Delta$ PMS location is used determination of critical amplitudes from integrated flows.
  - Computed locally. Outputs to file and results are used to determine initial conditions for integrating flows via thermal perturbations at fixed renormalized mass (YLE flow).
4. HPC_YZ_delta
  - Integrates $\rho$ flows to compute exponent $\delta$ and ampltidue $B_c$ using corrections to scaling.
  - Implemeneted on NCSU's high performance cluster (HPC).
5. HPC_YZ_mcdata_once, HPC_YZ_gamma
  - Integrates $\rho$ flows to compute the renormalized mass for various $t>0$ which tunes to vanishing (sufficiently small) extneral field $H \approx 0$ for a following computation of $\gamma$ and $C_2^+$.
  - Implemeneted on NCSU's high performance cluster (HPC).
6. HPC_YZ_YLE_R40 (not cleaned up)
  - Integrates $\rho$ flows for $t>0$ and renormalized mass $m_R^2 \approx 0$ into the complex plane using a modified set of partially (R) renormalized parameters down to $t=-40$.
  - Utilizes a switch to $\phi$ flows (transformed to purely real expressions) which approach YLE fixed point. Extracts relevant YLE location/amplitude parameters.
  - Implemeneted on NCSU's high performance cluster (HPC).
