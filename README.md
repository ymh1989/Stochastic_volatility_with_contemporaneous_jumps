## Stochastic volatility with contemporaneous jumps ##

### Introduction ###
The Black–Scholes model [1] concerns with the problems of option pricing and has achieved great success, especially in stock option. However, many empirical studies indicate the drawbacks of the Black–Scholes model. First of all, the most apparent violated assumption is that the volatility of underlying assets is constant. In order to complement the property of changing volatility, several variations of stochastic volatility(SV) models, for example the Heston model [2], are regarded as more appropriate than classical Black–Scholes model. Next, classical Black–Scholes model can not reflect the property of discontinuous jumps in the asset price. Hence, the jump model generally provide more efficient and accurate valuation of derivatives. [3-4] Nevertheless, both SV models and jump models are not totally describe empirical cases. For example, SV models have trouble describing the crashes in 1987 or 2008. Jump models, on the other hand, generally can explain the crashes using appropriate parameters which make a negative jump. However, jump models cannot explain the fluctuation of implied volatility after the crashes. For these reasons, models which include the properties of not only SV, but also jump models are proposed.

In this repo, I have implemented numerical solution of stochastic volatility with contemporaneous jumps(SVCJ) model[5-6] by using finite difference method. Especially, I would like to show how boundary condition affects the solution. First, I applied linear boundary condition. Next, hybrid boundary condition is used. 

###About SVCJ###
- `PDF` file which is described about SVCJ is attached in `code` folder.

### Implementation ###
- `MATLAB` codes and figures are uploaded.
- Operator spliting method(OSM) is used. 


###Reference###

\[1\] Black, Fischer, and Myron Scholes. "The pricing of options and corporate liabilities." The journal of political economy (1973): 637-654.

\[2\] Heston, Steven L. "A closed-form solution for options with stochastic volatility with applications to bond and currency options." Review of financial studies 6.2 (1993): 327-343.

\[3\] Kremer, Joseph W., and Rodney L. Roenfeldt. "Warrant pricing: jump-diffusion vs. Black-Scholes." Journal of Financial and Quantitative Analysis 28.02 (1993): 255-272.

\[4\] Tankov, Peter. Financial modelling with jump processes. Vol. 2. CRC press, 2003.

\[5\] Eraker, Bjørn. "Do stock prices and volatility jump? Reconciling evidence from spot and option prices." The Journal of Finance 59.3 (2004): 1367-1403.

\[6\] Zhang, Ying-Ying, et al. "Quadratic finite element and preconditioning for options pricing in the SVCJ model." Journal of Computational Finance, Forthcoming (2011).