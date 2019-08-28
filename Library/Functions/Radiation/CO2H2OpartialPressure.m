function [pCO2,pH2O] = CO2H2OpartialPressure(products,mdot_prod)
%CO2H2OpartialPressure outputs the partial pressure of a gas stream in Pa
%products - combustion products mixture object
%mdot_prod - mass flow rate of individual combustion products, kg/s

P = 101325;
ndot_out = mdot_prod./products.molarMasses;
p = zeros(1,length(products.species));

%Ideal gas assumption: mole frac = volume frac ~ partial pressure
for i = 1:length(products.species)
    p(i) = ndot_out(i)/sum(ndot_out)*P;
end
pCO2 = p(end-1);
pH2O = p(end);
end