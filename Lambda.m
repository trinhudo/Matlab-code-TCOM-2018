function out = Lambda(rho,kappa,alpha,K,a3,lSF,lNF,lrsi)
out = 1 - Xi(kappa,alpha,K,lrsi)*Omega(rho,a3,lSF,lNF);
end