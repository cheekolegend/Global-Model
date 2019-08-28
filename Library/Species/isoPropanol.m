classdef isoPropanol < Species
    %isoPropanol
    %   All properties are in SI units unless otherwise specified: kg, m3,
    %   K, J, Pa
    %   April 26, 2019 11:17 AM: All properties updated.
    %   CombustionEfficiency and sootVolFrac still need updating
    properties

    end
    methods
        function obj = isoPropanol()
            obj.criticalTemp = 508.3;
            obj.criticalPres = 47.62*100000;
            obj.molarMass = 60.096/1000; %kg/mol
            obj.standardEnthalpy = 0; %J/kg
            obj.referenceTemperature = 298.15;
            obj.referencePressure = 101325;
            obj.heatCombustion = 33.1*1000000; %SFPE pg. 560
            obj.combustionEfficiency = calcCombustionEfficiency(0.118); %Smoke point from https://pdfs.semanticscholar.org/68b2/3dc6348fe21ff26bdf89b5b4d713c8558337.pdf
            obj.radiationFrac = calcRadiationFraction(0.118); 
            obj.sootVolFrac = 0.07e-6; %Assumed the same as ethanol
            obj.stoichCoeff = [1,9/2,3,4];
            obj.bubbleTemp = 82;
        end
        %Returns the VP in Pa given T in K for T(185.28-508.3)
        function y = vapourPressure(obj, T)
            %Perry and Green, 1997 https://www.accessengineeringlibrary.com/browse/perrys-chemical-engineers-handbook-eighth-edition/p200139d89972_48001#p200139d89982_55
            A = 96.094;
            B = -8575.4;
            C = -10.292;
            D = 1.6665e-17;
            E = 6;
            y = exp(A + B/T + C*log(T) + D*T^E);
        end
        %Returns the liquid density in kg/m3 given T in K (146-536)
        function y = liquidDensity(obj,T)
            %Yaws' Thermophysical Properties of Chemicals and Hydrocarbons (Electronic Edition)
            A = 0.27300;
            B = 0.27093;
            C = 508.31;
            n = 0.24300;
            y = A*power(B,-(power(1-T/C,n)))*1000;
        end
        %Returns the molar volume in m3/mol given T in K
        function y = liquidMolarVolume(obj,T)
            y = obj.molarMass/obj.liquidDensity(T);
        end
        %Returns the ideal gas density in kg/m3 given T in K and P in Pa
        function y = gasDensity(obj, P, T)
            y = P*obj.molarMass/8.314/T;
        end
        %Returns the liquid heat capacity in J/kg.K given T in K for T(173-513)
        function y = liquidHeatCapacity(obj, T)
            % Yaws' Handbook of Thermodynamic Properties for Hydrocarbons and Chemicals
            A = 873.1467177;
            B = -10.05925883;
            C = 0.045744611;
            D = -8.37662E-05;
            E = 5.5582E-08;
            y = (A+B*T+C*power(T,2)+D*power(T,3)+E*power(T,4))/obj.molarMass;
        end
        function y = liquidHeatCapacityIntegral(obj,T0,T)
            A = 873.1467177;
            B = -10.05925883;
            C = 0.045744611;
            D = -8.37662E-05;
            E = 5.5582E-08;
            y = (A*(T-T0)+B/2.*(power(T,2)-power(T0,2))+C/3.*(power(T,3)-power(T0,3))+D/4.*(power(T,4)-power(T0,4))+E/5.*(power(T,5)-power(T0,5)))/obj.molarMass;
        end
        %Returns the ideal gas heat capacity in J/kg.K given T in K for T(150-1500)
        function y = gasHeatCapacity(obj, T)
            % Yaws' Handbook of Thermodynamic Properties for Hydrocarbons and Chemicals
            A = 48.7560898;
            B = -0.096083222;
            C = 0.00142278;
            D = -2.85799E-06;
            E = 2.68889E-09;
            F = -1.24099E-12;
            G = 2.25613E-16;
            y = (A+B*T+C*power(T,2)+D*power(T,3)+E*power(T,4)+F*power(T,5)+G*power(T,6))/obj.molarMass;
        end
        %Returns the integral of the ideal gas heat capacity in J/kg given T in K
        function y = gasHeatCapacityIntegral(obj, T)
            A = 48.7560898;
            B = -0.096083222;
            C = 0.00142278;
            D = -2.85799E-06;
            E = 2.68889E-09;
            F = -1.24099E-12;
            G = 2.25613E-16;
            y = (A*(T-obj.referenceTemperature)+B/2.*(power(T,2)-power(obj.referenceTemperature,2))+C/3.*(power(T,3)-power(obj.referenceTemperature,3))+D/4.*(power(T,4)-power(obj.referenceTemperature,4))+E/5.*(power(T,5)-power(obj.referenceTemperature,5))+F/6.*(power(T,6)-power(obj.referenceTemperature,6)+G/7.*(power(T,7)-power(obj.referenceTemperature,7))))/obj.molarMass;
        end
        function y = gasHeatCapacityIntegralTo298(obj, T)
            A = 48.7560898;
            B = -0.096083222;
            C = 0.00142278;
            D = -2.85799E-06;
            E = 2.68889E-09;
            F = -1.24099E-12;
            G = 2.25613E-16;
            y = (A*(obj.referenceTemperature-T)+B/2.*(power(obj.referenceTemperature,2)-power(T,2))+C/3.*(power(obj.referenceTemperature,3)-power(T,3))+D/4.*(power(obj.referenceTemperature,4)-power(T,4))+E/5.*(power(obj.referenceTemperature,5)-power(T,5))+F/6.*(power(obj.referenceTemperature,6)-power(T,6)+G/7.*(power(obj.referenceTemperature,7)-power(T,7))))/obj.molarMass;
        end         
        %Returns the latent heat of vapourization in J/kg given T in K for T(182.56-540.26)
        function y = latentHeatVap(obj, T)
            %Yaws https://www.sciencedirect.com/science/article/pii/B9780815515968500122
            A = 58.9824;
            n = 0.3260;
            y = A*power((1-T/obj.criticalTemp),n)*1000/obj.molarMass;
        end
        %Returns the liquid phase enthalpy in J/kmol given T in K
        function y = liquidEnthalpy(obj, T)
            y = gasEnthalpy(obj,T)-latentHeatVap(obj,T);
        end
        %Returns the gas phase enthalpy in J/kmol given T in K
        function y = gasEnthalpy(obj, T)
            y = obj.standardEnthalpy+gasHeatCapacityIntegral(obj,T);
        end
        function y = airFuelRatio(obj)
            y = obj.stoichCoeff(2)*32e-3/(obj.stoichCoeff(1)*obj.molarMass)/0.231;
        end
    end
end

