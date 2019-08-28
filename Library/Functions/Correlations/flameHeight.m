function L = flameHeight(fuel_mixture,mdot_fuel,D)
%flameHeight calculates the average flame height using Heskestad's
%correlation.

if iscolumn(mdot_fuel) == (true)
    mdot_fuel = mdot_fuel';
end

P = 101325; %Pa
g = 9.81; %m s-2
T_a = Air().referenceTemperature;
rho_a = Air().gasDensity(P,T_a);
cp_a = Air().gasHeatCapacity(T_a);

Hc = fuel_mixture.heatsCombustion;
Xa = fuel_mixture.combustionEfficiencies;
Q = sum(mdot_fuel.*Hc.*Xa);
Qstar = Q/(rho_a*cp_a*T_a*power(g*D,0.5)*D);
L = D*(-1.02 + 3.7*power(Qstar,2/5));

%Check for negative flame height in case of pure inert fuel
if L <= 0
    L = 1e-6;
end
end