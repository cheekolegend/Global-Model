%Use the excel spreadsheet to confirm the prediction of mdotAir
%Miscellaneous startup tasks
clear;
addpath(genpath(pwd)); %adds all subfolders to path

%Declare species
ethanol = Heptane();
water = Water();
species = [ethanol,water]; nSpecies = length(species);
massFracs = [0.4,0.6];
mixture = Mixture(species,massFracs);

yO2 = 0.12;
mdot_fuel = [0.1002;0]; %Column vector
mdot_air = mdotAir(mixture,mdot_fuel,yO2)