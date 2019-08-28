%Miscellaneous startup tasks
clear;
addpath(genpath(pwd)); %adds all subfolders to path

%Declare species
ethanol = Ethanol();
water = Water();
species = [ethanol,water]; nSpecies = length(species);
massFracs = [1,0];
mixture = Mixture(species,massFracs);

yO2_min = 0.1;
mdot_fuel = [1;1]; %Column vector
mdot_air1 = mdotAir(mixture,mdot_fuel,yO2_min)
mdot_air2 = mdotAir2(mixture,mdot_fuel,yO2_min)

%mdot_air increases as water content goes up because of the assumption that
%X = average(X_a)

%Air entrainment is the same using mdotAir and mdotAir2 for Ethanol-Heptane
%mixtures.