function x_mass = molesToMass(mixture,x_mol)
molarMasses = mixture.molarMasses;
%All variables must be column vectors
if iscolumn(molarMasses) ~= (true)
    molarMasses = molarMasses';
end
if iscolumn(x_mol) ~= (true)
    x_mol = x_mol';
end
    x_mass = zeros(length(mixture.species),1);
    for i = 1:length(mixture.species)
        x_mass(i) = x_mol(i).*mixture.species(i).molarMass/(sum(x_mol.*molarMasses));
    end
end