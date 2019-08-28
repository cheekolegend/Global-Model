function mwk = getWilsonCoeffs(mixture)
%getWilsonCoeffs obtains the binary interaction coefficients for a given
%mixture. The interaction coefficient matrix must be coded by hand.

nSpecies = length(mixture.species);

%Binary mixtures
if nSpecies == 2
    if speciesInMixture(Ethanol(),mixture) == (true) && speciesInMixture(Water(),mixture) == true
        EthanolWater
    elseif speciesInMixture(Ethanol(),mixture) == (true) && speciesInMixture(Heptane(),mixture) == true
        EthanolHeptane
    elseif speciesInMixture(Ethanol(),mixture) == (true) && speciesInMixture(Hexane(),mixture) == true
        EthanolHexane
    elseif speciesInMixture(Ethanol(),mixture) == (true) && speciesInMixture(isoPropanol(),mixture) == true
        EthanolisoPropanol
    elseif speciesInMixture(Propanol(),mixture) == (true) && speciesInMixture(Water(),mixture) == true
        PropanolWater
    end
end

%Ternary mixtures
if nSpecies == 3
    if speciesInMixture(Ethanol(),mixture) == (true) && speciesInMixture(isoPropanol(),mixture) == (true) && speciesInMixture(Water(),mixture) == true
        EthanolisoPropanolWater
    end
end  

end