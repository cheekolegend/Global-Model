function m = mdotProducts(mixture,mdot_fuel,mdot_air)
%mdotProducts returns an array of mass flow rates for [fuel1,fuel2,...,O2,N2,CO2,H2O]
%mixture is the fuel mixture that is evaporating
%mdot_fuels is the mass flowrate of the fuel mixture evaporating
%mdot_air is the total mass flow rate of air entrained
%All entries must be column vectors!
air = Air();
O2 = Oxygen();
CO2 = CarbonDioxide();
H2O = Water();

coefficients = mixture.stoichCoeffs;
%Must be row vectors
vfuel = coefficients(1:end,1)';
vO2 = coefficients(1:end,2)';
vCO2 = coefficients(1:end,3)';
vH2O = coefficients(1:end,4)';
Xa = mixture.combustionEfficiencies;
molarMasses = mixture.molarMasses;

%Convert all variables to row vectors
if iscolumn(mdot_fuel) == (true)
    mdot_fuel = mdot_fuel';
end

m_fuel = zeros(1,length(mixture.species));
for i = 1:length(mixture.species)
    m_fuel(i) = mdot_fuel(i) - mdot_fuel(i)*Xa(i);
end
m_O2 = mdot_air*air.xO2-O2.molarMass*sum(vO2./(vfuel.*molarMasses).*mdot_fuel.*Xa);
m_N2 = mdot_air*air.xN2;
m_CO2 = CO2.molarMass*sum(vCO2./(vfuel.*molarMasses).*mdot_fuel.*Xa);
m_H2O = H2O.molarMass*sum(vH2O./(vfuel.*molarMasses).*mdot_fuel.*Xa);
m = [m_fuel,m_O2,m_N2,m_CO2,m_H2O];
end