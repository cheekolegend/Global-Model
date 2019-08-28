function sol = mdotAir2(mixture,mdot_fuel,yO2)
%mdotAir calculates air entrainment into the flame (kg
%   v = stoichiometric coefficient of [fuel, O2, N2, CO2, H2O]
%   yO2 = limiting mole fraction of oxygen
%   x = mass fraction of fuel mixture
%   X = reaction conversion

%Ensure mdot_fuel is a column vector
if isrow(mdot_fuel) == (true)
    mdot_fuel = mdot_fuel';
end

%Get average molar mass
x = zeros(length(mdot_fuel),1);
for i = 1:length(x)
    x(i) = mdot_fuel(i)/sum(mdot_fuel);
end
x_mol = massToMoles(mixture,x)';
Mavg = sum(x_mol.*mixture.molarMasses);
X = averageProperty(mixture.combustionEfficiencies,mdot_fuel);

stoichCoeffs = mixture.stoichCoeffs;
    function y = f(n_air)
        n_HC = 1/Mavg; %mols fuel per kg fuel
        n_N2 = 0.79*n_air;
        n_O2 = 0.21*n_air;
        n_CO2 = sum(stoichCoeffs(:,3)'./mixture.molarMasses.*x'); %Works.
        n_H2O = sum(stoichCoeffs(:,4)'./mixture.molarMasses.*x');
             
        y = n_air - (yO2*(n_HC*(1-X) + (n_O2+n_N2+0.5*n_H2O)*X) + (n_CO2+0.5*n_H2O-n_O2)*X)/(0.21 - yO2);  
    end

    %solution = @(nAir) f(nAir);
    n = fzero(@f,10);
        
    sol = n*sum(mdot_fuel)*Air().molarMass;
end

