# Computational_Finance
A repo where I can practice for the course Computational Finance

## Content
| Folder | What's inside ? |
|:-------|:------------|
|Unrealistic B&S| Two ways to show why B&S is unrealistic & why we need better models |
|Tests| Just some random simulations |
|Portfolio Management| Everything related to this module of the Computational Finance course |

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