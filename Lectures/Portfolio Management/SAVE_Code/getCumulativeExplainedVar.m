function CumExplainedVar = getCumulativeExplainedVar_old(latent, n)
    ExplainedVar = latent(1:n)/sum(latent);
    CumExplainedVar = sum(ExplainedVar);
end