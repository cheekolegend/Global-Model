function res = speciesInMixture(species,mixture)
%SPECIESINMIXTURE returs true if the designated species in the mixture
%object
%species: desired species
%mixture: mixture object
    if any(ismember(species.molarMass,mixture.molarMasses))
        res = (true);
    else
        res = (false);
    end
end