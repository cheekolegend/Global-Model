classdef Ethanol < Species
    %Heptane
    %   All properties are in SI units unless otherwise specified: kg, m3,
    %   K, J, Pa
    
    properties

    end
    methods
        function obj = Ethanol()
            obj.criticalTemp = 513.92;
            obj.criticalPres =61.48e5;
            obj.molarMass = 46.07/1000; %kg/mol
            obj.standardEnthalpy = 0; %J/kg
            obj.referenceTemperature = 298.15;
            obj.referencePressure = 101325;
            obj.heatCombustion = 26.8e6;
            obj.combustionEfficiency = 0.92; %Hamins (1999)
            obj.radiationFrac = 0.23; %Hamins (1999)
            obj.sootVolFrac = 0.07e-6; %Hamins (1999)
            obj.stoichCoeff = [1,3,2,3]; %fuel,O2,CO2,H2O
            obj.bubbleTemp = 78.37;
        end
        %Returns the VP in Pa given T in K for T(159.05-513.92)
        function y = vapourPressure(obj, T)
            %Perry and Green, 1997
            A = 74.475;
            B = -7164.3;
            C = -7.327;
            D = 3.13E-06;
            E = 2;
            y = exp(A + B/T + C*log(T) + D*T^E);
        end
        %Returns the liquid density in kg/m3 given T in K (182-540)
        function y = liquidDensity(obj,T)
            %Yaws Hydrocarbons and Chemicals
            A = 0.27600;
            B = 0.27668;
            C = 516.25;
            n = 0.23670;
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
        %Returns the liquid heat capacity in J/kg.K given T in K for T(173.15-483)
        function y = liquidHeatCapacity(obj, T)
            %Yaws Hydrocarbons and Chemicals
            A = 238.3080385;
            B = -2.380642326;
            C = 0.013317056;
            D = -3.19962E-05;
            E = 3.15051E-08;
            y = (A+B*T+C*power(T,2)+D*power(T,3)+E*power(T,4))/obj.molarMass;
        end
        function y = liquidHeatCapacityIntegral(obj,T0,T)
            A = 238.3080385;
            B = -2.380642326;
            C = 0.013317056;
            D = -3.19962E-05;
            E = 3.15051E-08;
            y = (A*(T-T0)+B/2.*(power(T,2)-power(T0,2))+C/3.*(power(T,3)-power(T0,3))+D/4.*(power(T,4)-power(T0,4))+E/5.*(power(T,5)-power(T0,5)))/obj.molarMass;
        end
        %Returns the ideal gas heat capacity in J/kg.K given T in K for T(150-1500)
        function y = gasHeatCapacity(obj, T)
            %Yaws Hydrocarbons and Chemicals 2009
            A = 53.77704843;
            B = -0.207301605;
            C = 0.001422271;
            D = -2.62982E-06;
            E = 2.39649E-09;
            F = -1.08912E-12;
            G = 1.96393E-16;
            y = (A+B*T+C*power(T,2)+D*power(T,3)+E*power(T,4)+F*power(T,5)+G*power(T,6))/obj.molarMass;
        end
        %Returns the integral of the ideal gas heat capacity in J/kg given T in K
        function y = gasHeatCapacityIntegral(obj, T)
            A = 2.7091E+01;
            B = 1.1055E-01;
            C = 1.0957E-04;
            D = -1.505E-07;
            E = 4.6601E-11;
            F = 0;
            G = 0;
            y = (A*(T-obj.referenceTemperature)+B/2.*(power(T,2)-power(obj.referenceTemperature,2))+C/3.*(power(T,3)-power(obj.referenceTemperature,3))+D/4.*(power(T,4)-power(obj.referenceTemperature,4))+E/5.*(power(T,5)-power(obj.referenceTemperature,5))+F/6.*(power(T,6)-power(obj.referenceTemperature,6)+G/7.*(power(T,7)-power(obj.referenceTemperature,7))))/obj.molarMass;
        end
        function y = gasHeatCapacityIntegralTo298(obj, T)
            A = 2.7091E+01;
            B = 1.1055E-01;
            C = 1.0957E-04;
            D = -1.505E-07;
            E = 4.6601E-11;
            F = 0;
            G = 0;
            y = (A*(obj.referenceTemperature-T)+B/2.*(power(obj.referenceTemperature,2)-power(T,2))+C/3.*(power(obj.referenceTemperature,3)-power(T,3))+D/4.*(power(obj.referenceTemperature,4)-power(T,4))+E/5.*(power(obj.referenceTemperature,5)-power(T,5))+F/6.*(power(obj.referenceTemperature,6)-power(T,6)+G/7.*(power(obj.referenceTemperature,7)-power(T,7))))/obj.molarMass;
        end                
        %Returns the latent heat of vapourization in J/kg given T in K for T(182.56-540.26)
        function y = latentHeatVap(obj, T)
            %Yaws https://www.sciencedirect.com/science/article/pii/B9780815515968500122
            A = 60.8036;
            n = 0.3800;
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
        %Returns mass air to fuel ratio http://ethanolpro.tripod.com/id213.html
        function y = airFuelRatio(obj)
            y = obj.stoichCoeff(2)*32e-3/(obj.stoichCoeff(1)*obj.molarMass)/0.231;
        end
    end
end

