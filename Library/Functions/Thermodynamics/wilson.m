function gam = wilson(mixture,x,T)
    R = 1.98721; %cal/mol/K, from Aspen HYSYS
    n   = length(mixture.species);
    gam  = zeros(n,1);
    
    i = 1; j = 1; k = 1;
    sum1 = 0; sum2 = 0; sum3 = 0;
    x = x(:);
    
    molarVolume = getWilsonVol(mixture);
    coefficients = getWilsonCoeffs(mixture);
    coeff = coefficients;

    for i = 1:n
        for j = 1:n
            coeff(i, j) = (molarVolume(j)/molarVolume(i))*...
                exp(-1.*coefficients(i, j)/R/T);
        end
    end
    
    for i = 1:n
        sum1 = 0; sum2 = 0;
        for j = 1:n
            sum3 = 0;
            for k = 1:n
                sum3 = sum3 + x(k, 1)*coeff(j, k);
            end
            sum1 = sum1 + x(j, 1)*coeff(i, j);
            sum2 = sum2 + x(j, 1)*coeff(j, i)/sum3;
        end
        sum1 = log(sum1);
        gam(i) = exp(1 - sum1 - sum2);
    end
end
