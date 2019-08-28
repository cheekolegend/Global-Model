function test_flameTemperature()

%Miscellaneous startup tasks
clear;
addpath(genpath('C:\Users\aaron\OneDrive - Dalhousie University\G-Drive\MASc\Research\Fire\Modelling\Global Model 2.0\Global Model 2.0')); %adds all subfolders to path

%Declare species
ethanol = Ethanol();
water = Water();
species = [ethanol,water];
massFracs = [1,0];

mixture = Mixture(species,massFracs);

%Creating reactants list for use in flame temperature summation
O2N2 = [Oxygen(),Nitrogen()]; O2N2massFracs = [0,0]; %mass fracs don't matter
reactantList = horzcat(species,O2N2);
reactantMassFracs = horzcat(massFracs,O2N2massFracs);
reactants = Mixture(reactantList,reactantMassFracs);

%Creating products list for use in flame temperature summation
O2N2CO2H2O = [Oxygen(),Nitrogen(),CarbonDioxide(),Water()]; O2N2CO2H2OmassFracs = [0,0,0,0]; %mass fracs don't matter
productList = horzcat(species,O2N2CO2H2O);
productMassFracs = horzcat(massFracs,O2N2CO2H2OmassFracs);
products = Mixture(productList,productMassFracs);

T0 = 298;
yO2 = 0.1;
Ts = 350;
mdot_fuel = [1,0];
mdot_air = mdotAir(mixture,mdot_fuel,yO2);
mdotReact = [mdot_fuel,Air().xO2*mdot_air,Air().xN2*mdot_air];
mdotProd = mdotProducts(mixture,mdot_fuel,mdot_air);

%Tf = 1600;
% Q_etOH = 0.08*Ethanol().gasHeatCapacityIntegral(Tf)*Ethanol().molarMass;
% Q_O2 = 13.1*Oxygen().gasHeatCapacityIntegral(Tf)*Oxygen().molarMass;
% Q_N2 = 274.5*Nitrogen().gasHeatCapacityIntegral(Tf)*Nitrogen().molarMass;
% Q_CO2 = 39.9*CarbonDioxide().gasHeatCapacityIntegral(Tf)*CarbonDioxide().molarMass;
% Q_H2O = 59.9*Water().gasHeatCapacityIntegral(Tf)*Water().molarMass;
% Q = [Q_etOH;Q_O2;Q_N2;Q_CO2;Q_H2O];

Tflame = fzero(@func,1500)

%Flame temperature being overestimated somehow... 4th sum is okay. Gas heat
%capacity integrals are being overestimated, so should Tf be lower?
function res = func(Tf)
    res = sum(mdotReact.*reactants.gasHeatCapacityIntegralsTo298(Ts)) + ...
            sum(mdotProd.*products.gasHeatCapacityIntegrals(Tf)) - ...
            sum(mdot_fuel.*(mixture.heatsCombustion.*(mixture.combustionEfficiencies-mixture.radiationFracs))) + ...
            sum(mdot_fuel.*(mixture.latentHeatsVap(Ts)+mixture.liquidHeatCapacityIntegrals(T0,Ts)));
end
end
