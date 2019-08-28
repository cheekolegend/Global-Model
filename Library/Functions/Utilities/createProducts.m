function products = createProducts(fuel_mixture)
%CREATEPRODUCTS generates a products list of object type mixture().
%It contains the combustion reactants: [fuel1,fuel2,...,O2,N2,CO2,H2O]
%The list is used in the flame temperature calculation.
species = fuel_mixture.species;
massFracs = fuel_mixture.x;

%Creating products list for use in flame temperature summation
O2N2CO2H2O = [Oxygen(),Nitrogen(),CarbonDioxide(),Water()]; O2N2CO2H2OmassFracs = [0,0,0,0]; %mass fracs don't matter
productList = horzcat(species,O2N2CO2H2O);
productMassFracs = horzcat(massFracs,O2N2CO2H2OmassFracs);
products = Mixture(productList,productMassFracs);
end