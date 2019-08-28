classdef Heptane < Species
    %Heptane
    %   All properties are in SI units unless otherwise specified: kg, m3,
    %   K, J, Pa
    
    properties

    end
    methods
        function obj = Heptane()
            obj.criticalTemp = 540.26;
            obj.criticalPres = 27.40*100000;
            obj.molarMass = 100.204/1000; %kg/mol
            obj.standardEnthalpy = -187780/obj.molarMass; %J/kg
            obj.referenceTemperature = 298.15;
            obj.referencePressure = 101325;
            obj.heatCombustion = 44.4e6;
            obj.combustionEfficiency = 0.92; %Hamins (1999)
            obj.radiationFrac = 0.30; %Hamins (1999)
            obj.sootVolFrac = 0.50e-6; %Hamins (1999)
            obj.stoichCoeff = [1,11,7,8];
            obj.bubbleTemp = 98.42;
        end
        %Returns the VP in Pa given T in K for T(182.57-540.2)
        function y = vapourPressure(obj, T)
            %Perry and Green, 1997
            A = 87.829;
            B = -6996.4;
            C = -9.8802;
            D = 7.2099E-06;
            E = 2;
            y = exp(A + B/T + C*log(T) + D*T^E);
        end
        %Returns the liquid density in kg/m3 given T in K (182-540)
        function y = liquidDensity(obj,T)
            %Yaws Hydrocarbons and Chemicals
            A = 0.23200;
            B = 0.25649;
            C = 540.26;
            n = 0.27910;
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
        %Returns the liquid heat capacity in J/kg.K given T in K for T(198.15-520)
        function y = liquidHeatCapacity(obj, T)
            %Yaws Hydrocarbons and Chemicals
            A = 1251.285246;
            B = -13.34709676;
            C = 0.061092865;
            D = -0.000119248;
            E = 8.64746e-08;
            y = (A+B*T+C*power(T,2)+D*power(T,3)+E*power(T,4))/obj.molarMass;
        end
        function y = liquidHeatCapacityIntegral(obj, T0, T)
            A = 1251.285246;
            B = -13.34709676;
            C = 0.061092865;
            D = -0.000119248;
            E = 8.64746e-08;
            y = (A*(T-T0)+B/2.*(power(T,2)-power(T0,2))+C/3.*(power(T,3)-power(T0,3))+D/4.*(power(T,4)-power(T0,4))+E/5.*(power(T,5)-power(T0,5)))/obj.molarMass;
        end
        %Returns the ideal gas heat capacity in J/kg.K given T in K for T(150-1500)
        function y = gasHeatCapacity(obj, T)
            %Yaws Hydrocarbons and Chemicals
            A = 99.51463088;
            B = -0.27476053;
            C = 0.002878833;
            D = -5.30451E-06;
            E = 4.65818E-09;
            F = -2.02447E-12;
            G = 3.48821E-16;
            y = (A+B*T+C*power(T,2)+D*power(T,3)+E*power(T,4)+F*power(T,5)+G*power(T,6))/obj.molarMass;
        end
        %Returns the integral of the ideal gas heat capacity in J/kg given T in K
        function y = gasHeatCapacityIntegral(obj, T)
            A = 99.51463088;
            B = -0.27476053;
            C = 0.002878833;
            D = -5.30451E-06;
            E = 4.65818E-09;
            F = -2.02447E-12;
            G = 3.48821E-16;
            y = (A*(T-obj.referenceTemperature)+B/2.*(power(T,2)-power(obj.referenceTemperature,2))+C/3.*(power(T,3)-power(obj.referenceTemperature,3))+D/4.*(power(T,4)-power(obj.referenceTemperature,4))+E/5.*(power(T,5)-power(obj.referenceTemperature,5))+F/6.*(power(T,6)-power(obj.referenceTemperature,6)+G/7.*(power(T,7)-power(obj.referenceTemperature,7))))/obj.molarMass;
        end
        function y = gasHeatCapacityIntegralTo298(obj, T)
            A = 99.51463088;
            B = -0.27476053;
            C = 0.002878833;
            D = -5.30451E-06;
            E = 4.65818E-09;
            F = -2.02447E-12;
            G = 3.48821E-16;
            y = (A*(obj.referenceTemperature-T)+B/2.*(power(obj.referenceTemperature,2)-power(T,2))+C/3.*(power(obj.referenceTemperature,3)-power(T,3))+D/4.*(power(obj.referenceTemperature,4)-power(T,4))+E/5.*(power(obj.referenceTemperature,5)-power(T,5))+F/6.*(power(obj.referenceTemperature,6)-power(T,6)+G/7.*(power(obj.referenceTemperature,7)-power(T,7))))/obj.molarMass;
        end
        %Returns the latent heat of vapourization in J/kg given T in K for T(182.56-540.26)
        function y = latentHeatVap(obj, T)
            %Yaws https://www.sciencedirect.com/science/article/pii/B9780815515968500122
            A = 49.7300;
            n = 0.3860;
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
            y = obj.stoichCoeff(2)*Oxygen().molarMass/(obj.stoichCoeff(1)*obj.molarMass)/0.231;
        end
    end
end

