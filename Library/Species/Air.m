classdef Air < SpeciesAir
    %Air
    %   All properties are in SI units unless otherwise specified: kg, m3,
    %   K, J, Pa
    
    properties
        xO2
        xN2
        
    end
    methods
        function obj = Air()
            obj.criticalTemp = 0;
            obj.criticalPres = 0;
            obj.molarMass = 28.97/1000; %kg/mol
            obj.standardEnthalpy = 0; %J/kg
            obj.referenceTemperature = 298.15;
            obj.referencePressure = 101325;
            obj.heatCombustion = 0;
            obj.combustionEfficiency = 0;
            obj.radiationFrac = 0; 
            obj.sootVolFrac = 0; 
            obj.stoichCoeff = [0,0,0,0];
            obj.xO2 = 0.233;
            obj.xN2 = 1-obj.xO2;
        end
        %Returns the VP in Pa given T in K for T(277-396)
        function y = vapourPressure(obj, T)
            y = 0; %Antoine takes T in Celcius
        end
        %Returns the liquid density in kg/m3 given T in K (182-540)
        function y = liquidDensity(obj,T)
            y = 0;
        end        
        %Returns the ideal gas density in kg/m3 given T in K and P in Pa
        function y = gasDensity(obj, P, T)
            y = P*obj.molarMass/8.314/T;
        end
        %Returns the liquid heat capacity in J/kg.K given T in K for T()
        function y = liquidHeatCapacity(obj, T)
            y = 0;
        end
        function y = liquidHeatCapacityIntegral(obj, T0, T)
            y = 0;
        end
        %Returns the ideal gas heat capacity in J/kg.K given T in K for T(150-2000)
        function y = gasHeatCapacity(obj, T)
            %Smith, Van Ness (2005)
            R = 8.314; %J/mol/K
            A = 3.355;
            B = 0.575e-3;
            C = 0;
            D = 0.016e5;
            y = R*(A+B*T+C*power(T,2)+D*power(T,-2))/obj.molarMass;
        end
        %Returns the integral of the ideal gas heat capacity in J/kg given T in K
        function y = gasHeatCapacityIntegral(obj, T)
            R = 8.314; %J/mol/K
            A = 3.280;
            B = 0.593e3;
            C = 0;
            D = 0.04e-5;
            y = 0;
        end
        function y = gasHeatCapacityIntegralTo298(obj, T)
            R = 8.314; %J/mol/K
            A = 3.280;
            B = 0.593e3;
            C = 0;
            D = 0.04e-5;
            y = 0;
        end
        %Returns the latent heat of vapourization in J/kg given T in K for T()
        function y = latentHeatVap(obj, T)
            y = 0;
        end
        %Returns the liquid phase enthalpy in J/kmol given T in K
        function y = liquidEnthalpy(obj, T)
            y = gasEnthalpy(obj,T)-latentHeatVap(obj,T);
        end
        %Returns the gas phase enthalpy in J/kmol given T in K
        function y = gasEnthalpy(obj, T)
            y = obj.standardEnthalpy+gasHeatCapacityIntegral(obj,T);
        end
        %Returns the gas phase dynamic viscosity in Pa*s given T in K (100-800)
        function y = gasViscosity(obj,T)
            y = 1.458e-6*T^1.5/(T + 110.4);
        end
        %Returns the gas phase thermal conductivity in W/m/K given T in K (100-1000)
        function y = gasThermalConductivity(obj,T)
            y = 0.002624*T^1.5/(T + 245.4*10^(-12/T));
        end
    end
end

