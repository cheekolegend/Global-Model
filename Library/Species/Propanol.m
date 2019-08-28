classdef Propanol < Species
    %1-Propanol, n-propyl alcohol
    %   All properties are in SI units unless otherwise specified: kg, m3,
    %   K, J, Pa
    
    properties

    end
    methods
        function obj = Propanol()
            obj.criticalTemp = 536.8;
            obj.criticalPres = 51.75*100000;
            obj.molarMass = 60.096/1000; %kg/mol
            obj.standardEnthalpy = 0; %J/kg
            obj.referenceTemperature = 298.15;
            obj.referencePressure = 101325;
            obj.heatCombustion = 31.3*1000000; %SFPE pg. 560
            obj.combustionEfficiency = calcCombustionEfficiency(0.118); %Smoke point from https://pdfs.semanticscholar.org/68b2/3dc6348fe21ff26bdf89b5b4d713c8558337.pdf
            obj.radiationFrac = calcRadiationFraction(0.118); 
            obj.sootVolFrac = 0.07e-6; %Assumed the same as ethanol
            obj.stoichCoeff = [1,9/2,3,4];
            obj.bubbleTemp = 97;
        end
        %Returns the VP in Pa given T in K for T(185.28-508.3)
        function y = vapourPressure(obj, T)
            %Perry and Green, 1997
            A = 76.964;
            B = -7623.8;
            C = -7.4924;
            D = 5.94E-18;
            E = 6;
            y = exp(A + B/T + C*log(T) + D*T^E);
        end
        %Returns the liquid density in kg/m3 given T in K (146-536)
        function y = liquidDensity(obj,T)
            %Yaws Hydrocarbons and Chemicals
            A = 0.27500;
            B = 0.26976;
            C = 536.71;
            n = 0.24940;
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
            %Yaws Hydrocarbons and Chemicals
            A = 407.9200837;
            B = -4.626366627;
            C = 0.025163289;
            D = -5.77018E-05;
            E = 5.19894E-08;
            y = (A+B*T+C*power(T,2)+D*power(T,3)+E*power(T,4))/obj.molarMass;
        end
        function y = liquidHeatCapacityIntegral(obj,T0,T)
            A = 407.9200837;
            B = -4.626366627;
            C = 0.025163289;
            D = -5.77018E-05;
            E = 5.19894E-08;
            y = (A*(T-T0)+B/2.*(power(T,2)-power(T0,2))+C/3.*(power(T,3)-power(T0,3))+D/4.*(power(T,4)-power(T0,4))+E/5.*(power(T,5)-power(T0,5)))/obj.molarMass;
        end
        %Returns the ideal gas heat capacity in J/kg.K given T in K for T(150-1500)
        function y = gasHeatCapacity(obj, T)
            %Yaws Hydrocarbons and Chemicals 2009
            A = 63.79880716;
            B = -0.240449787;
            C = 0.001829791;
            D = -3.42665E-06;
            E = 3.1388E-09;
            F = -1.43091E-12;
            G = 2.58595E-16;
            y = (A+B*T+C*power(T,2)+D*power(T,3)+E*power(T,4)+F*power(T,5)+G*power(T,6))/obj.molarMass;
        end
        %Returns the integral of the ideal gas heat capacity in J/kg given T in K
        function y = gasHeatCapacityIntegral(obj, T)
            A = 63.79880716;
            B = -0.240449787;
            C = 0.001829791;
            D = -3.42665E-06;
            E = 3.1388E-09;
            F = -1.43091E-12;
            G = 2.58595E-16;
            y = (A*(T-obj.referenceTemperature)+B/2.*(power(T,2)-power(obj.referenceTemperature,2))+C/3.*(power(T,3)-power(obj.referenceTemperature,3))+D/4.*(power(T,4)-power(obj.referenceTemperature,4))+E/5.*(power(T,5)-power(obj.referenceTemperature,5))+F/6.*(power(T,6)-power(obj.referenceTemperature,6)+G/7.*(power(T,7)-power(obj.referenceTemperature,7))))/obj.molarMass;
        end
        function y = gasHeatCapacityIntegralTo298(obj, T)
            A = 63.79880716;
            B = -0.240449787;
            C = 0.001829791;
            D = -3.42665E-06;
            E = 3.1388E-09;
            F = -1.43091E-12;
            G = 2.58595E-16;
            y = (A*(obj.referenceTemperature-T)+B/2.*(power(obj.referenceTemperature,2)-power(T,2))+C/3.*(power(obj.referenceTemperature,3)-power(T,3))+D/4.*(power(obj.referenceTemperature,4)-power(T,4))+E/5.*(power(obj.referenceTemperature,5)-power(T,5))+F/6.*(power(obj.referenceTemperature,6)-power(T,6)+G/7.*(power(obj.referenceTemperature,7)-power(T,7))))/obj.molarMass;
        end         
        %Returns the latent heat of vapourization in J/kg given T in K for T(182.56-540.26)
        function y = latentHeatVap(obj, T)
            %Yaws https://www.sciencedirect.com/science/article/pii/B9780815515968500122
            A = 70.1792;
            n = 0.451;
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

