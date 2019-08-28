classdef Oxygen < Species
    %1-Propanol, n-propyl alcohol
    %   All properties are in SI units unless otherwise specified: kg, m3,
    %   K, J, Pa
    
    properties

    end
    methods
        function obj = Oxygen()
            obj.criticalTemp = 154.6;
            obj.criticalPres = 50.43*100000;
            obj.molarMass = 31.999/1000; %kg/mol
            obj.standardEnthalpy = 0; %J/kg
            obj.referenceTemperature = 298.15;
            obj.referencePressure = 101325;
            obj.heatCombustion = 0;
            obj.combustionEfficiency = 0;
            obj.radiationFrac = 0; 
            obj.sootVolFrac = 0; 
            obj.stoichCoeff = [0,0,0,0];
            obj.bubbleTemp = 0;
        end
        %Returns the VP in Pa given T in K for T()
        function y = vapourPressure(obj, T)
            y = 0; %Antoine takes T in Celcius
        end
        %Returns the liquid density in kg/m3 given T in K ()
        function y = liquidDensity(obj,T)
            y = 0;
        end
        %Returns the molar volume in m3/mol given T in K
        function y = liquidMolarVolume(obj,T)
            y = 0;
        end
        %Returns the ideal gas density in kg/m3 given T in K and P in Pa
        function y = gasDensity(obj, P, T)
            y = P*obj.molarMass/8.314/T;
        end
        %Returns the liquid heat capacity in J/kg.K given T in K for T(173.15-483)
        function y = liquidHeatCapacity(obj, T)
            y = 0;
        end
        function y = liquidHeatCapacityIntegral(obj, T0,T)
            y = 0;
        end
        %Returns the ideal gas heat capacity in J/kg.K given T in K for T(150-1500)
        function y = gasHeatCapacity(obj, T)
            %Smith, Van Ness (2005)
            R = 8.314; %J/mol/K
            A = 3.639;
            B = 0.506e-3;
            C = 0;
            D = -0.227e5;
            y = R*(A+B*T+C*power(T,2)+D*power(T,-2))/obj.molarMass;
        end
        %Returns the integral of the ideal gas heat capacity in J/kg given T in K
        function y = gasHeatCapacityIntegral(obj, T)
            R = 8.314; %J/mol/K
            A = 3.639;
            B = 0.506e-3;
            C = 0;
            D = -0.227e5;
            y = R*(A*(T-obj.referenceTemperature)+B/2.*(power(T,2)-power(obj.referenceTemperature,2))+C/3.*(power(T,3)-power(obj.referenceTemperature,3))+(-1)*D*(power(T,-1)-power(obj.referenceTemperature,-1)))/obj.molarMass;
        end
        function y = gasHeatCapacityIntegralTo298(obj, T)
            R = 8.314; %J/mol/K
            A = 3.639;
            B = 0.506e-3;
            C = 0;
            D = -0.227e5;
            y = R*(A*(obj.referenceTemperature-T)+B/2.*(power(obj.referenceTemperature,2)-power(T,2))+C/3.*(power(obj.referenceTemperature,3)-power(T,3))+(-1)*D*(power(obj.referenceTemperature,-1)-power(T,-1)))/obj.molarMass;
        end
        %Returns the latent heat of vapourization in J/kg given T in K for T(182.56-540.26)
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
    end
end

