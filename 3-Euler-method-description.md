'Euler-method-1.cpp' is the first look at the brilliant numerical integration method known as Euler method (or rather first-order Runge-Kutta method), which, however, is ought to be treated accordingly to its peculiarities.

   a) The task was to solve numerically the following Cauchy problem: $\dot{y} = 2t, \quad y_0 = y(t_0) = 0, \quad t \in \left[ 0 , 10 \right].$
   Three even partitions were used, with steps $h = 0.1, 0.01, 0.001$. Given the analytical solution $y(t) = t^2$, one can easily calculate the error for each step and see if the method works at all and where the error does maximise. Another *gnuplot* graph shows it clearly.

   b) The second part of the task was to solve another Cauchy problem: $y' = f'(x),\quad y(x_0) = f(x_0),\quad x \in \left[ x_0 , X \right] $, where $f$ is taken from the previous program, and $\left[ x_0 , X \right]$ is chosen so that the function is differentiable (and therefore integrable) on the entire segment. I picked $x_0 = 0.5$ and $X = 0.9$ since the function is discontionuous at $x = 1$ . This fact certainly has its effect on the quality of the numerical approximation, so we can see the error $R_1$ with the naked eye.
