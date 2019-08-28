classdef Water < Species
    %Heptane
    %   All properties are in SI units unless otherwise specified: kg, m3,
    %   K, J, Pa
    
    properties

    end
    methods
        function obj = Water()
            obj.criticalTemp = 647.1;
            obj.criticalPres = 220.55*100000;
            obj.molarMass = 18.015/1000; %kg/mol
            obj.standardEnthalpy = 0; %J/kg
            obj.referenceTemperature = 298.15;
            obj.referencePressure = 101325;
            obj.heatCombustion = 0;
            obj.combustionEfficiency = 0.000001; %Hamins (1999)
            obj.radiationFrac = 0; %Hamins (1999)
            obj.sootVolFrac = 0; %Hamins (1999)
            obj.stoichCoeff = [1,0,0,0];
            obj.bubbleTemp = 100;
        end
        %Returns the VP in Pa given T in K for T(273.16-647.13)
        function y = vapourPressure(obj, T)
            %Perry and Green, 1997
            A = 73.649;
            B = -7258.2;
            C = -7.3037;
            D = 4.17E-06;
            E = 2;
            y = exp(A + B/T + C*log(T) + D*T^E);
        end
        %Returns the liquid density in kg/m3 given T in K (273-648)
        function y = liquidDensity(obj,T)
            %http://ddbonline.ddbst.de/DIPPR105DensityCalculation/DIPPR105CalculationCGI.exe?component=Water
            A = 0.14395;
            B = 0.0112;
            C = 649.727;
            D = 0.05107;
            y = A/(B^(1+(1-T/C)^D));
        end
        %Returns the molar volume in m3/mol given T in K
        function y = liquidMolarVolume(obj,T)
            y = obj.molarMass/obj.liquidDensity(T);
        end
        %Returns the ideal gas density in kg/m3 given T in K and P in Pa
        function y = gasDensity(obj, P, T)
            y = P*obj.molarMass/8.314/T;
        end
        %Returns the liquid heat capacity in J/kg.K given T in K for T(273-615)
        function y = liquidHeatCapacity(obj, T)
            %https://www.accessengineeringlibrary.com/browse/chemical-properties-handbook/c9780070734012ch03
            A = 92.053;
            B = -3.9953E-02;
            C = -2.1103E-04;
            D = 5.3469E-07;
            E = 0;
            y = (A+B*T+C*power(T,2)+D*power(T,3)+E*power(T,4))/obj.molarMass;
        end
        function y = liquidHeatCapacityIntegral(obj, T0, T)
            A = 92.053;
            B = -3.9953E-02;
            C = -2.1103E-04;
            D = 5.3469E-07;
            E = 0;
            y = (A*(T-T0)+B/2.*(power(T,2)-power(T0,2))+C/3.*(power(T,3)-power(T0,3))+D/4.*(power(T,4)-power(T0,4))+E/5.*(power(T,5)-power(T0,5)))/obj.molarMass;
        end
        %Returns the ideal gas heat capacity in J/kg.K given T in K for T(100-1500)
        function y = gasHeatCapacity(obj, T)
            %https://www.accessengineeringlibrary.com/browse/chemical-properties-handbook/c9780070734012ch02
            A = 33.933;
            B = -6.4186E-03;
            C = 2.9906E-05;
            D = -1.7825E-08;
            E = 3.6934E-12;
            F = 0;
            G = 0;
            y = (A+B*T+C*power(T,2)+D*power(T,3)+E*power(T,4)+F*power(T,5)+G*power(T,6))/obj.molarMass;
        end
        %Returns the integral of the ideal gas heat capacity in J/kg given T in K
        function y = gasHeatCapacityIntegral(obj, T)
            A = 33.933;
            B = -6.4186E-03;
            C = 2.9906E-05;
            D = -1.7825E-08;
            E = 3.6934E-12;
            F = 0;
            G = 0;
            y = (A*(T-obj.referenceTemperature)+B/2.*(power(T,2)-power(obj.referenceTemperature,2))+C/3.*(power(T,3)-power(obj.referenceTemperature,3))+D/4.*(power(T,4)-power(obj.referenceTemperature,4))+E/5.*(power(T,5)-power(obj.referenceTemperature,5))+F/6.*(power(T,6)-power(obj.referenceTemperature,6)+G/7.*(power(T,7)-power(obj.referenceTemperature,7))))/obj.molarMass;
        end
        function y = gasHeatCapacityIntegralTo298(obj, T)
            A = 33.933;
            B = -6.4186E-03;
            C = 2.9906E-05;
            D = -1.7825E-08;
            E = 3.6934E-12;
            F = 0;
            G = 0;
            y = (A*(obj.referenceTemperature-T)+B/2.*(power(obj.referenceTemperature,2)-power(T,2))+C/3.*(power(obj.referenceTemperature,3)-power(T,3))+D/4.*(power(obj.referenceTemperature,4)-power(T,4))+E/5.*(power(obj.referenceTemperature,5)-power(T,5))+F/6.*(power(obj.referenceTemperature,6)-power(T,6)+G/7.*(power(obj.referenceTemperature,7)-power(T,7))))/obj.molarMass;
        end   
        %Returns the latent heat of vapourization in J/kg given T in K for T(273-373)
        function y = latentHeatVap(obj, T)
            %https://www.accessengineeringlibrary.com/browse/chemical-properties-handbook/c9780070734012ch05
            A = 52.053;
            n = 0.321;
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
            y = 0;
        end
    end
end

