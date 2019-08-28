function reactants = createReactants(fuel_mixture)
%CREATEREACTANTS generates a reactant of object type mixture().
%It contains the combustion reactants: [fuel1,fuel2,...,O2,N2]
%The list is used in the flame temperature calculation.
species = fuel_mixture.species;
massFracs = fuel_mixture.x;

O2N2 = [Oxygen(),Nitrogen()]; O2N2massFracs = [0,0]; %mass fracs don't matter
reactantList = horzcat(species,O2N2);
reactantMassFracs = horzcat(massFracs,O2N2massFracs);
reactants = Mixture(reactantList,reactantMassFracs);
end