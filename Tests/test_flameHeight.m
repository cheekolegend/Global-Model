%Test Flame Height
%Miscellaneous startup tasks
clear;
addpath(genpath('C:\Users\aaron\OneDrive - Dalhousie University\G-Drive\MASc\Research\Fire\Modelling\Global Model 2.0\Global Model 2.0')); %adds all subfolders to path

%Declare species
ethanol = Ethanol();
heptane = Heptane();
species = [ethanol,heptane];
massFracs = [1,1]; %Doesn't matter here

mixture = Mixture(species,massFracs);

mdot_fuel = [1,0]; %kg
D = 0.3;
L = flameHeight(mixture,mdot_fuel,D)
