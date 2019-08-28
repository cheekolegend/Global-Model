%Inputs
species1 = Ethanol();
species2 = Water();
V1 = 700; %mL
V2 = 300; %mL

%Calculations
T = 20+273; %K
V1 = V1/1000/1000; %m3
V2 = V2/1000/1000;
m1 = V1*species1.liquidDensity(T);
m2 = V2*species2.liquidDensity(T);
x1 = m1/(m1 + m2);
x2 = 1 - x1;

%Outputs
x = [x1,x2]
m = m1 + m2

