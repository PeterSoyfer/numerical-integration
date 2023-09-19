# Numerical Integration and Optimisation Theory

The first cluster of programs represents basics of numerical analysis:
1. The first piece of code demonstrates instability occuring whilst calculating the integral $I_n = \int\limits_{0}^{1} \frac{x^n}{x + 6} dx = \frac{1}{n} - 6 I_{n-1}$ by a direct recurrent method, given that $I_0 = \log\left(\frac{7}{6}\right)$. Namely, calculating $I_{31}$ in such way gives us a cumulative error of order $\approx (-5)\cdot 10^{7}$, which in comparison to the real exact value of $I_{31} \approx 0+$ clearly persuades a courteous reader that this method makes no sense in our case.
