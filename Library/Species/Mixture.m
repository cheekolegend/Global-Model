classdef Mixture
    properties
    species = [];
    x = [];
    x_mol = [];
    criticalTemp = [];
    criticalPres = [];
    molarMass = [];
    standardEnthalpy = [];
    referenceTemperature = [];
    referencePressure = [];
    heatCombustion = [];
    combustionEfficiency = [];
    radiationFrac = [];
    sootVolFrac = [];
    %stoichCoeff = []; %Stoichiometric coefficients for combustion (fuel,O2,CO2,H2O)
    bubbleTemp = [];
    airFuelRatio = [];
    end
    methods
        function obj = Mixture(list,x)
            obj.species = list;
            obj.x = x;
            
            %obj.x_mol = ...
            totalMoles = 0;
            for i = 1:length(obj.species)
                totalMoles = totalMoles + obj.x(i)/obj.species(i).molarMass;
            end
            for i = 1:length(obj.species)
                obj.x_mol(i) = obj.x(i)/obj.species(i).molarMass/totalMoles;
                obj.criticalTemp(i) = obj.species(i).criticalTemp;
                obj.criticalPres(i) = obj.species(i).criticalPres;
                obj.molarMass(i) = obj.species(i).molarMass;
                obj.standardEnthalpy(i) = obj.species(i).standardEnthalpy;
                obj.referenceTemperature(i) = obj.species(i).referenceTemperature;
                obj.referencePressure(i) = obj.species(i).referencePressure;
                obj.bubbleTemp(i) = obj.species(i).bubbleTemp;
            end
        end
        
        function y = bubbleTemps(obj)
            y = length(obj.species);
            for i = 1:length(obj.species)
                y(i) = obj.species(i).bubbleTemp;
            end
        end
        
        function y = liquidMolarVolumes(obj,T)
            y = zeros(1,length(obj.species));
            for i = 1:length(obj.species)
                y(i) = obj.species(i).liquidMolarVolume(T);
            end
        end
        
        function y = liquidDensities(obj,T)
            y = zeros(1,length(obj.species));
            for i = 1:length(obj.species)
                y(i) = obj.species(i).liquidDensity(T);
            end
        end
        
        function y = heatsCombustion(obj)
            y = length(obj.species);
            for i = 1:length(obj.species)
                y(i) = obj.species(i).heatCombustion;
            end
        end
        
        function y = latentHeatsVap(obj,T)
            y = length(obj.species);
            for i = 1:length(obj.species)
                y(i) = obj.species(i).latentHeatVap(T);
            end
        end
        
        function y = airFuelRatios(obj)
            y = length(obj.species);
            for i = 1:length(obj.species)
                y(i) = obj.species(i).airFuelRatio;
            end
        end
        
        function y = radiationFracs(obj)
            y = length(obj.species);
            for i = 1:length(obj.species)
                y(i) = obj.species(i).radiationFrac;
            end
        end
        
        function y = combustionEfficiencies(obj)
            y = length(obj.species);
            for i = 1:length(obj.species)
                y(i) = obj.species(i).combustionEfficiency;
            end
        end 
        
        function y = sootVolFracs(obj)
            y = length(obj.species);
            for i = 1:length(obj.species)
                y(i) = obj.species(i).sootVolFrac;
            end
        end 
        
        function y = mdotAirs(obj,yO2)
            y = length(obj.species);
            for i = 1:length(obj.species)
                y(i) = 0.029*ndotAir(obj.species(i),yO2);
            end
        end  
        
        function y = stoichCoeffs(obj)
            y = length(obj.species);
            for i = 1:length(obj.species)
                for j = 1:length(obj.species(i).stoichCoeff)
                    y(i,j) = obj.species(i).stoichCoeff(j);
                end
            end
        end 
        
        function y = molarMasses(obj)
            y = length(obj.species);
            for i = 1:length(obj.species)
                y(i) = obj.species(i).molarMass;
            end
        end 
        
        function y = gasHeatCapacityIntegrals(obj,T)
            y = length(obj.species);
            for i = 1:length(obj.species)
                y(i) = obj.species(i).gasHeatCapacityIntegral(T);
            end
        end 
        
        function y = gasHeatCapacityIntegralsTo298(obj,T)
            y = length(obj.species);
            for i = 1:length(obj.species)
                y(i) = obj.species(i).gasHeatCapacityIntegralTo298(T);
            end
        end 
        
        function y = liquidHeatCapacityIntegrals(obj,T0,T)
            y = length(obj.species);
            for i = 1:length(obj.species)
                y(i) = obj.species(i).liquidHeatCapacityIntegral(T0,T);
            end
        end         
   end
end