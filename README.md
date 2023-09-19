# Numerical Integration and Optimisation Theory

The first cluster of programs represents basics of numerical analysis:
1. The first piece of code demonstrates instability occuring whilst calculating the integral $I_n = \int\limits_{0}^{1} \frac{x^n}{x + 6} dx$ by a direct recurrent method, given that $I_0 = \log\left(\frac{7}{6}\right)$. Namely, calculating $I_{31}$ in such way gives us a cumulative machine error of order $\approx (-6)^{31}$, which in comparison to the real exact value of $I_{31} \approx 0+$ clearly persuades a courteous reader that this method makes no sense in our case.
