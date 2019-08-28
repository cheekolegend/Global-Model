function x_mol = massToMoles(mixture,x)
molarMasses = mixture.molarMasses;
%All variables must be column vectors
if iscolumn(molarMasses) ~= (true)
    molarMasses = molarMasses';
end
if iscolumn(x) ~= (true)
    x = x';
end
    x_mol = zeros(length(mixture.species),1);
    for i = 1:length(mixture.species)
        x_mol(i) = x(i)./mixture.species(i).molarMass/(sum(x./molarMasses));
    end
end