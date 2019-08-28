function T_bubble = bubbleT(mixture,x,P)
%bubbleT Calculates bubble point temperature for a mixture
    %P: Ambient pressure, Pa
    %T_guess: Initial guess for bubble point temperature, K
    %x: fuel mole fractions
    %species: species array
    
    %Algorithm: 
    %1. Guess temperature
    %2. Calculate Psat of each component at T_guess
    %3. Calculate P; Pcalc = x(1)Psat(1)gamma1 + x(2)Psat(2)gamma2
    %4. If Pcalc != P, guess again
    
    nSpecies = length(mixture.species);
    Tguess = 0;    
    Psat = zeros(nSpecies,1);

    for i = 1:nSpecies
        Tguess = Tguess + mixture.species(i).bubbleTemp;
    end
    Tguess = Tguess/nSpecies + 273.15;
    %Tguess = [min(mixture.bubbleTemps)+273,max(mixture.bubbleTemps)+273];
    
%     T1 = 298;
%     T = 323;
%     Pcal = 0;
%     while abs(Pcal-P)>=1
%         if (Pcal>P)
%             T=T1+((T-T1)/2);
%         elseif (Pcal<P)
%             T=T+((T-T1)/2);
%         end
%         P1=mixture.species(1).vapourPressure(T);
%         P2=mixture.species(2).vapourPressure(T);
%         Pcal=P1*x(1) + P2*x(2);
%     end
%     T_bubble = T;
    

    objective = @(T) P - Pcalc(T);
    T_bubble = fzero(objective,Tguess); 
    
    function res = Pcalc(T)
       gamma = wilson(mixture,x,T);
       y = zeros(nSpecies,1);
       for i = 1:nSpecies
           Psat(i) = mixture.species(i).vapourPressure(T);
           y(i) = (x(i)*Psat(i)*gamma(i))/P;
       end
       res = sum(y)*P;
       
       %Check if imagniary number...
       if isreal(res)==0
           res = P;
       end
    end
end
    


