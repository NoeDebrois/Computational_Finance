function [c,ceq] = nonlinConstrVar(x, Ret, pval, tgtVaR)
    retPf = x'*Ret';
    VaR = quantile(retPf, pval);
    c = VaR-tgtVaR;
    ceq = [];
end