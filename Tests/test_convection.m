%Miscellaneous startup tasks
clear;
addpath(genpath(pwd)); %adds all subfolders to path

%Declare species
ethanol = Ethanol();
heptane = Heptane();
species = [ethanol,heptane]; nSpecies = length(species);
massFracs = [0.99,0.01];
mixture = Mixture(species,massFracs);

D = 0.3;
A = pi*D^2/4;
mdotPUA = 0.015;
Tf = 1500;
Ts = 360;

y = [1,0];
mdot_fuel = [mdotPUA*A,0.0001];

Q = convection(mixture,mdot_fuel,y,Tf,Ts,D)