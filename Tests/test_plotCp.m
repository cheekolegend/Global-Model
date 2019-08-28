%Miscellaneous startup tasks
clear;
addpath(genpath('C:\Users\aaron\OneDrive - Dalhousie University\G-Drive\MASc\Research\Fire\Modelling\Global Model 2.0\Global Model 2.0')); %adds all subfolders to path

%Declare species
ethanol = Ethanol();
water = Water();
species = [ethanol,water];
massFracs = [1,0];

mixture = Mixture(species,massFracs);

T = linspace(298,1500,100);
cpEthanol = Ethanol().gasHeatCapacityIntegral(T)*Ethanol().molarMass;
cpOxygen = Oxygen().gasHeatCapacityIntegral(T)*Oxygen().molarMass;
cpNitrogen = Nitrogen().gasHeatCapacityIntegral(T)*Nitrogen().molarMass;
cpCarbonDioxide = CarbonDioxide().gasHeatCapacityIntegral(T)*CarbonDioxide().molarMass;
cpWater = Water().gasHeatCapacityIntegral(T)*Water().molarMass;

%plot(T,cpCarbonDioxide)
plot(T,cpEthanol,T,cpOxygen,T,cpNitrogen,T,cpCarbonDioxide,T,cpWater)
legend('EtOH','O2','N2','CO2','H2O');