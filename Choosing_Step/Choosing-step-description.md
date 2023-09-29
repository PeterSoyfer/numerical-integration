The task was to solve Cauchy problem for the ODE system of harmonic oscillator from 'Euler_failure' folder with Dormand-Prince, now choosing the optimal integration step $h$
instead of setting it voluntarily to $0.1, 0.01$ or $0.001$. Given fixed values of acceptable error $tol = 10^{-7}, 10^{-9}$ and $10^{-11}$ and time periods
$T = 100 \pi, \dots , 10^5 \pi$, parameters of 'numerical-vs-analytical difference' $\left|\widetilde{x}(T) - x(T)\right|$ and $\left|\widetilde{z}(T) - z(T)\right|$,
number of integration steps $N_{\text{steps}}$ and the actual value of integration step $h$. Global error and Runge numbers $R_x, R_z$ were also calculated.

