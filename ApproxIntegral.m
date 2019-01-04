function out = ApproxIntegral(mu,chi,xi)
%
A1 = exp(-mu*chi)/chi;
A2 = xi*igamma(0,mu*chi);
A3 = 0;
for uu=2:49
    B1 = ((-1)^uu)*(xi^uu)/(factorial(uu));
    B21 = exp(-(mu*chi));
    B22 = 0;
    for vv = 1:(uu-1)
        temp = (factorial(vv-1))*((-chi)^(uu-vv-1))/...
            ((factorial(uu-1))*(mu^vv));
        B22 = B22 + temp;
    end
    B23 = ((-chi)^(uu-1))/(factorial(uu-1))*(ei(-mu*chi));
    B2 = B21*B22-B23;
    A3 = A3 + B1*B2;
end
out = A1-A2+A3;