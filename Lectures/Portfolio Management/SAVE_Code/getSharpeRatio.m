function s = getSharpeRatio(x, ExpRet, V, rf)
    s = (x'*ExpRet'-rf)/sqrt(x'*V*x);
end