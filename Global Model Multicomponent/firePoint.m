%What is the concentration of alcohol where the fire will extinguish?
clear;
addpath(genpath(pwd));

P = 101325; %Pa
yO2 = 0.233; %mass frac
T_infty = 298.15; %K

species1 = Ethanol();
species2 = Water();
mixture = Mixture([species1,species2],[1,0]);

%Defining liquid fraction range
n = 100;
x1 = linspace(0,1,n);
x2 = 1 - x1;
x = [x1;x2];
x = x';

%Storage variables
x_mol = x;
y = x;
B = zeros(n,1);

%Mixture propertis
Xr_mix = mixture.radiationFracs;
r_mix = mixture.airFuelRatios; r_mix(1) = 2.08; %oxygen to fuel mass ratio
Hc_mix = mixture.heatsCombustion;

for i = 1:n
    x_mol(i,:) = massToMoles(mixture,x(i,:));
    y(i,:) = solveY(mixture,x_mol(i,:),P);
    Ts = bubbleT(mixture,x_mol(i,:),P);
    cp_mix = mixture.gasHeatCapacityIntegrals(Ts); %from 298 to Ts K
    Hv_mix = mixture.latentHeatsVap(Ts);
    B(i) = (yO2*(sum(y(i,:).*Hc_mix.*(1-Xr_mix))/sum(y(i,:).*r_mix)) - sum(y(i,:).*cp_mix))/sum(y(i,:).*Hv_mix);
end

a = [0,1];
b = [0,1];
%plot(x(:,1),y(:,1),a,b)
plot(x(:,1),B)