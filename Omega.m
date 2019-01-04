function out = Omega(rho,a3,lSF,lNF)
out = 1/lSF*exp(-rho/a3/lNF)/(rho/lNF + 1/lSF);
end
