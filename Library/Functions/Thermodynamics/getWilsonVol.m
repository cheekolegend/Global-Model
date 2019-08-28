function V = getWilsonVol(mixture)

nSpecies = length(mixture.species);

if nSpecies == 2
    if speciesInMixture(Ethanol(),mixture) == (true) && speciesInMixture(Water(),mixture) == true
        EthanolWater_WilsonVol
    elseif speciesInMixture(Ethanol(),mixture) == (true) && speciesInMixture(Heptane(),mixture) == true
        EthanolHeptane_WilsonVol
    elseif speciesInMixture(Ethanol(),mixture) == (true) && speciesInMixture(Hexane(),mixture) == true
        EthanolHexane_WilsonVol
    elseif speciesInMixture(Ethanol(),mixture) == (true) && speciesInMixture(isoPropanol(),mixture) == true
        EthanolisoPropanol_WilsonVol
    elseif speciesInMixture(Propanol(),mixture) == (true) && speciesInMixture(Water(),mixture) == true
        PropanolWater_WilsonVol
    end
end

if nSpecies == 3
    if speciesInMixture(Ethanol(),mixture) == (true) && speciesInMixture(isoPropanol(),mixture) == (true) && speciesInMixture(Water(),mixture) == true
        EthanolisoPropanolWater_WilsonVol
    end
end


end