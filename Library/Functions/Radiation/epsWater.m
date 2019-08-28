function eps_w = epsWater(T,pH2O,Lm)
%EPSWATER calculates the total emissivity of water vapour at T > 400 K
    %T = gas temperature (K)
    %pH2O = partial pressure of water (Pa)
    %Lm = mean beam length (m)

    t = T/1000;
    pH2O = pH2O/100000; %convert to bar
    Lm = Lm*100; %convert to cm
    lamda = log10(pH2O*Lm); %path length, bar*cm
    
    %Coefficients for water vapour
    c = [-2.2118    -1.1987     0.035596
          0.85667    0.93048    -0.14391
         -0.10838   -0.17156    0.045915];
    
    %Essentially each number in a row is a coefficient for a polynomial of
    %form a = c0 + c1*t^1 + c2*t^2
    for i = 1:3
        sumC = 0;
        for j = 1:3
            sumC = sumC + c(i,j)*t^(j-1); %j-1 since Leckner starts index at 0
        end
        a(i) = sumC;
    end
    
    sumA = 0;
    for i = 1:3
        sumA = sumA + a(i)*lamda^(i-1); %i-1 since Leckner starts index at 0
    end
    
    eps_w = exp(sumA);
end