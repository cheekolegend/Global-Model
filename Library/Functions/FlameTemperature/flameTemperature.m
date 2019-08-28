function res = flameTemperature(fuel_mixture,mdot_fuel,reactants,mdotReact,products,mdotProd,T0,Ts,Tf)
%flameTemperature calculates the flame temperature for a multicomponent
%fuel. The energy balance is:
%Energy in = Energy out
%Heat of combustion = energy to heat to Tf + energy to evaporate
%fuel_mixture: evaporating fuel mixture object
%mdot_fuel: evaporation rate of individual fuel components, kg/s
%reactants: combustion reactants mixture object
%mdotReact: mass flow rate of individual reactants, kg/s
%products: combustion products mixture object
%mdotProd: mass flow rate of inidividual combustion products, kg/s
%T0: initial fuel temperature, K
%Ts: pool surface temperature, K
%Tf: flame temperature, K

%Convert all variables to row vectors
if iscolumn(mdot_fuel) == (true)
    mdot_fuel = mdot_fuel';
end
Hc = fuel_mixture.heatsCombustion;
Hv = fuel_mixture.latentHeatsVap(Ts);
Xa = fuel_mixture.combustionEfficiencies;
Xr = fuel_mixture.radiationFracs;
    res = sum(mdotReact.*reactants.gasHeatCapacityIntegralsTo298(Ts)) + ...
            sum(mdotProd.*products.gasHeatCapacityIntegrals(Tf)) - ...
            sum(mdot_fuel.*(Hc.*(Xa-Xr))) + ...
            sum(mdot_fuel.*(Hv+fuel_mixture.liquidHeatCapacityIntegrals(T0,Ts)));
end

