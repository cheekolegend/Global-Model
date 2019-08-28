%Inputs
species1 = Ethanol();
species2 = Water();
V1 = [57.2;60.35;58.84;56.94;52.63;51.5;46.87;42.19;35.36;25.55;16.81;14.42]; %mL
V2 = 100*ones(length(V1),1)-V1;

%Calculations
T = 20+273; %K
m1 = V1*species1.liquidDensity(T);
m2 = V2*species2.liquidDensity(T);
x1 = m1./(m1 + m2);
x2 = ones(length(x1),1) - x1;