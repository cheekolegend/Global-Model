%mols O2/mol etOH = 3
%mols O2/kg etOH = 3/MW_etOH = 65.2
%kg O2/kg etOH = 2.1
%If 1 kg etOH burns completely in 100 kg of air, then O2 out should equal:
%(100 kg air)*(0.233 kg O2/kg air) - (2.1 kg O2/kg etOH)*(1 kg etOH) = 21.2 kg O2 out 
%We get slightly more O2 out since only 92% of the etOH combusts

%Miscellaneous startup tasks
clear;
addpath(genpath(pwd)); %adds all subfolders to path

%Declare species
ethanol = Ethanol();
heptane = Heptane();
species = [ethanol,heptane];
massFracs = [0.4,0.6];

mixture = Mixture(species,massFracs);

mdot = [1;1];
mdot_air = 100;

mdot_out = mdotProducts(mixture,mdot,mdot_air)

%Test partial pressure
%Creating products list for use in flame temperature summation
species = [Ethanol(),Heptane()];
O2N2CO2H2O = [Oxygen(),Nitrogen(),CarbonDioxide(),Water()]; O2N2CO2H2OmassFracs = [0,0,0,0]; %mass fracs don't matter
productList = horzcat(species,O2N2CO2H2O);
productMassFracs = horzcat(massFracs,O2N2CO2H2OmassFracs);
products = Mixture(productList,productMassFracs);

Tf = 1500;
Ts = 360;
D = 0.3;
Qrad = radiation(mixture,mdot,products,mdot_out,Tf,Ts,D)

