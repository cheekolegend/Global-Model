function sol = mdotAir(mixture,mdot_fuel,yO2)
%mdotAir calculates air entrainment into the flame (kg
%   v = stoichiometric coefficient of [fuel, O2, N2, CO2, H2O]
%   yO2 = limiting mole fraction of oxygen
%   X = reaction conversion

%Ensure mdot_fuel is a column vector
if isrow(mdot_fuel) == (true)
    mdot_fuel = mdot_fuel';
end

    function y = f(n_air,i)
        n_HC = 1/mixture.species(i).molarMass; %mols fuel per kg fuel
        n_N2 = 0.79*n_air;
        n_O2 = 0.21*n_air;
        n_CO2 = mixture.species(i).stoichCoeff(3)/mixture.species(i).molarMass;
        n_H2O = mixture.species(i).stoichCoeff(4)/mixture.species(i).molarMass;
        
        X = mixture.species(i).combustionEfficiency;
        y = n_air - (yO2*(n_HC*(1-X) + (n_O2+n_N2+0.5*n_H2O)*X) + (n_CO2+0.5*n_H2O-n_O2)*X)/(0.21 - yO2);
        
        if mixture.species(i).molarMass == Water().molarMass
            y = n_air;
        end
    end

    n = zeros(length(mixture.species),1);
    for i = 1:length(mixture.species)
        solution = @(nAir) f(nAir,i);
        n(i) = fzero(solution,10);
    end
    sol = sum(n.*mdot_fuel)*Air().molarMass;
end

