classdef SpeciesAir < matlab.mixin.Heterogeneous
    %SPECIES Summary of this class goes here
    %   Detailed explanation goes here
    properties
        criticalTemp
        criticalPres
        molarMass
        standardEnthalpy
        referenceTemperature
        referencePressure
        heatCombustion
        combustionEfficiency
        radiationFrac
        sootVolFrac
        stoichCoeff %Stoichiometric coefficients for combustion (fuel,O2,CO2,H2O)
        bubbleTemp
    end
    methods
       function obj = SpeciesAir()
           obj = obj@matlab.mixin.Heterogeneous();
       end
    end
    methods (Abstract)
        %Returns the VP in Pa given T in K for T()
        vapourPressure(T);
        %Returns the liquid density in kg/m3 given T in K
        liquidDensity(T);
        %Returns the ideal gas density in kg/m3 given T in K and P in Pa
        gasDensity(P,T);
        %Returns the liquid heat capacity in J/kg.K given T in K for T()
        liquidHeatCapacity(T);
        %Returns the integral of the liquid heat capacity in J/kg.K given T in K for T()
        liquidHeatCapacityIntegral(T0,T);
        %Returns the ideal gas heat capacity in J/kg.K given T in K for T()
        gasHeatCapacity(T);
        %Returns the integral of the ideal gas heat capacity in J/kg given T in K
        gasHeatCapacityIntegral(T);
        %Returns the integral of the ideal gas heat capacity in J/kg given T in K
        gasHeatCapacityIntegralTo298(T);
        %Returns the latent heat of vapourization in J/kg given T in K for T()
        latentHeatVap(T);
        %Returns the liquid phase enthalpy in J/kmol given T in K
        liquidEnthalpy(T);
        %Returns the gas phase enthalpy in J/kmol given T in K
        gasEnthalpy(T);
        %Returns the gas phase dynamic viscosity in Pa*s given T in K
        gasViscosity(T);
        %Returns the gas phase thermal conductivity in W/m/K given T in K
        gasThermalConductivity(T);
    end
    
end

