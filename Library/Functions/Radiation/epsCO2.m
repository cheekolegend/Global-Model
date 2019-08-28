function eps_c = epsCO2(T,pCO2,Lm)
%EPSCO2 calculates the total emissivity of CO2 at T > 400 K
    %T = gas temperature (K)
    %pCO2 = partial pressure of CO2 (Pa)
    %Lm = mean beam length (m)

    t = T/1000;
    pCO2 = pCO2/100000; %convert to bar
    Lm = Lm*100; %convert to cm
    lamda = log10(pCO2*Lm); %path length, bar*cm
    
    %Coefficients for CO2 gas
    c = [-3.9781     2.7353     -1.9882      0.31054     0.015719    
         1.9326     -3.5932      3.7247     -1.4535      0.20132
         -0.35366    0.61766    -0.84207     0.39859    -0.063356
         -0.080181   0.31466    -0.19973     0.046532   -0.0033086];
    
    %Essentially each number in a row is a coefficient for a polynomial of
    %form a = c0 + c1*t^1 + c2*t^2 + ...
    for i = 1:4
        sumC = 0;
        for j = 1:5
            sumC = sumC + c(i,j)*t^(j-1); %j-1 since Leckner starts index at 0
        end
        a(i) = sumC;
    end
    
    sumA = 0;
    for i = 1:3
        sumA = sumA + a(i)*lamda^(i-1); %i-1 since Leckner starts index at 0
    end
    
    eps_c = exp(sumA);
end