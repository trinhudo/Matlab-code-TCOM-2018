function out = Xi(kappa,alpha,K,lrsi)
out = kappa * (K+1)*exp(-K)/lrsi...
    *exp(K*(K+1)/2/lrsi/(alpha+(K+1)/lrsi))...
    /sqrt((alpha+(K+1)/lrsi)*K*(K+1)/lrsi)...
    *whittakerM(-1/2, 0, K*(K+1)/lrsi/(alpha+(K+1)/lrsi));
end