%test averageProperty()

%Miscellaneous startup tasks
clear;
addpath(genpath(pwd)); %adds all subfolders to path

%Declare species
ethanol = Ethanol();
water = Water();
species = [ethanol,water]; nSpecies = length(species);
massFracs = [1,1]; %Doesn't matter

mixture = Mixture(species,massFracs);

mdot_fuel = [1,1];
fv = averageProperty(mixture.sootVolFracs,mdot_fuel)
