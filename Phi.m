function out = Phi(nu,kappa,alpha,K,lrsi)
out = 0;
for ll=0:50
    out = out + kappa*(K+1)*exp(-K)/lrsi/(factorial(ll)^2)...
        *((K*(K+1)/lrsi)^ll)...
        *((alpha+(K+1)/lrsi)^(-ll-1))...
        *gamma(ll+1)*gammainc((alpha+(K+1)/lrsi)*nu,ll+1);
end
end

