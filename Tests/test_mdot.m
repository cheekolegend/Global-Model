%test mdot
ethanol = Ethanol();
water = Water();
species = [ethanol]; nSpecies = length(species);
massFracs = [1];

mixture = Mixture(species,massFracs);

%Creating reactants list for use in flame temperature summation
O2N2 = [Oxygen(),Nitrogen()]; O2N2massFracs = [0,0]; %mass fracs don't matter
reactantList = horzcat(species,O2N2);
reactantMassFracs = horzcat(massFracs,O2N2massFracs);
reactants = Mixture(reactantList,reactantMassFracs);

%Creating products list for use in flame temperature summation
O2N2CO2H2O = [Oxygen(),Nitrogen(),CarbonDioxide(),Water()]; O2N2CO2H2OmassFracs = [0,0,0,0]; %mass fracs don't matter
productList = horzcat(species,O2N2CO2H2O);
productMassFracs = horzcat(massFracs,O2N2CO2H2OmassFracs);
products = Mixture(productList,productMassFracs);

%Ambient conditions
P = 101325; %Pa

%Tray dimensions
D = 0.3; %m
H = 0.05; %m
A = pi*D^2/4; %m2
V_cont = A*H; %container volume

%Entrainment assumptions
yO2 = 0.1; %limiting oxygen mol fraction in fire
X = 0.9; %reaction conversion

%Testing variables
Tf = 1500;
Ts = 78+273;
fv = mixture.species(1).sootVolFrac;
mdot_fuel = [0.001,0];
func = @(mdot) mdot - ((convection(mixture,mdot,Tf,Ts,D) + radiation(mixture,mdot,Tf,Ts,D,fv) - reradiation(Ts,D))/(sum(mixture.species(1).latentHeatVap(Ts)*mixture.x(1))));
%res = fzero(func,mdot_fuel)

Qconv = @(mdot) convection(mixture,mdot,Tf,Ts,D);
Qrad = @(mdot) radiation(mixture,mdot_fuel,Tf,Ts,D,fv);

x = linspace(0.0001,0.1,100);
plot(x,func(x))

% figure(1);
% plot(x,Qconv(x));
% figure(2);
% plot(x,Qrad(x));
