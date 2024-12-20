# Computational Mathematics & Computational Finance :
This is a repository where I put all my stuff from the course of Computational Finance, AY 2024/2025, Polytecnic University of Milan (Politecnico di Milano). In this course we learn computational methods to simulate Lévy processes and Stochastic Volatility processes, to numerically approximate solutions of PDEs, we also focus on Monte Carlo simulations, simulations based on Characteristic Function, Fourier Transform, etc.

<p align="center">
    <img src="./Summary.png" alt="drawing" width="800"/>
</p>


## Content
| Folder | What's inside ? |
|:-------|:------------|
|Lecture 01 - Unrealistic B&S| Two ways to show why B&S is unrealistic & why we need better models|
|Lecture 04 - Merton & Kou Models| Simulation of jump times & Merton and Kou models for $(S_t)_t$|
|Lecture 07 - Carr Madan (CM)| Pricing, and calibration of $\sigma$ on market data, using CM formula |
|Lecture 08 - Monte Carlo| Some MC simulations with variance reduction techniques|
|Lecture 09 - Monte Carlo, bis| Other MC simulations for path dependent option pricing, with variance reduction techniques|
|Lecture 10 - FFT| Convolution methods for pricing path dependent options (extension of CM)|
|Lecture 10 - MC American| MC simulations for American Option pricing, based on L&S article|
|Lecture 11 - 12 - PDE BS| Finite difference on B&S PDE for pricing (Euler Explicit / Implicit & Theta Method)|
|Lecture 13 - 14 - PDE AM & PIDE| B&S PDE for American option pricing & PIDE for pricing under Lévy|
|Lecture 15 - 2dPDE| |
|Lecture 16 - Heston Model| Heston model implementations for stochastic volatility (Euler & Andersen's article schemes) |
|Tests| Just some random simulations |
|Portfolio Management| Everything related to the $\textit{"Portfolio Management"}$ module of the CF course|
|Financial Engineering| Financial Engineering labs|

## Portfolio Management, Ginevra Angelini

### Returns :
- Total return : $H_{t,\tau} = \frac{P_t}{P_{t-\tau}}$ : $invariant$ in the equity market ;
- Linear return : $L_{t,\tau} = \frac{P_t}{P_{t-\tau}} - 1$ ;
- Compounded return : $C_{t,\tau} = ln(\frac{P_t}{P_{t-\tau}})$.

### Estimators :
Imagine our invariants $X_t$ are distributed following $N(\mu, \Sigma)$.
- The location estimator is the sample mean : $\hat \mu [ i_T ] = \frac{1}{T} \sum_{t=1}^{T} x_t$ ;
- The dispersion estimator is the sample covariance matrix : $\hat \Sigma [i_T] = \frac{1}{T} \sum_{t=1}^{T} (x_t - \hat \mu)(x_t - \hat \mu)'$.

### Evaluating allocations :
- The allocation (the nb of units bought for each securities) is represented by the N-dimensional vector $\alpha$ ;
- At the time of investment decision, the value of the portfolio is : $w_T(\alpha) = \alpha'p_T$ ;
- At the investment horizon $\tau$, the portfolio is a one-dimensional random variable : $W_{T+\tau}(\alpha) = \alpha'P_{T+\tau}$ ;
- The investor has one or more $objectives$ $\Psi$, namely quantities that the investor perceives as beneficial and therefore desires in the largest possible amounts. For example, it can be the absolute wealth : $\Psi_\alpha = W_{T+\tau}$.

### Index of satisfaction :
- We can summarize all the features of a given allocation $\alpha$ into one single number $S$ that indicates the respective degree of satisfaction : $\alpha \rightarrow S(\alpha)$. See the properties, slide $35/57$ in lecture $1$.
- See the notion of expected utility.

### Building strategies :
- A strategy is a set of investment choices based on a determined information set $I_t$, function of the information set : $S(t) = f(I(t))$ ;
- Let's assume $N$ is the number of possible investment assets, a single choice could be defined as a signal $s_t^i$ : $S(t) = (s^1_t, s^2_t, ..., s^N_t)$ ;
- Equity curve : $X_\tau = \Pi_{t=1}^\tau (1+\sum_{i=1}^Ns^i_{t-1}r^i_t)$ ;
- Annual return : $AnnRet = (\frac{X_T}{X_0})^{T/250} - 1$ ;
- Annual volatility : $AnnVol = std(\frac{X_t}{X_{t-1}}-1)$ ;
- Maximum drawback : $DD = min\{\frac{X_t}{X_{max}} - 1 : t = 1, ..., T \}$ ;
- Sharpe ratio : $Sharpe = \frac{AnnRet - RiskFree}{AnnVol}$ : how much extra return you receive for the extra volatility you endure for holding a riskier asset ;
- Calmar ratio : $Calmar = \frac{AnnRet - RiskFree}{MaxDD}$ : like Sharpe, but instead of using volatility to assess risk, it uses the maximum drawback.
