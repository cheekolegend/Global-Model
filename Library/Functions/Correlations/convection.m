function Qconv = convection(mixture,mdot_fuel,y,Tf,Ts,D)
%CONVECTION calculates convective heat transfer to the fuel surface using a
%stagnant film model.

if iscolumn(mdot_fuel) == (true)
    mdot_fuel = mdot_fuel';
end

if iscolumn(y) == (true)
    y = y';
end
    P = 101325;
    n = length(mixture.species);
    Hc = mixture.heatsCombustion;           Hc = sum(Hc.*y);
    Xa = mixture.combustionEfficiencies;    Xa = sum(Xa.*y);
    Xr = mixture.radiationFracs;            Xr = sum(Xr.*y);
    r = mixture.airFuelRatios;              r = sum(r.*y);
    
    A = pi*D^2/4;
    
    %Parameters for Rayleigh number
    g = 9.81; %gravitational acceleration
    T = (Tf + Ts)/2; %Temperature to evaluate air properties at
    nu = Air().gasViscosity(T)/Air().gasDensity(P,T); %kinematic viscosity, m2/s
    alpha = Air().gasThermalConductivity(T)/(Air().gasDensity(P,T)*Air().gasHeatCapacity(T)); %thermal diffusivity, m2/s
    Pr = nu/alpha; %Prandtl number
    Gr = 1/T*g*D^3*(Ts - Air().referenceTemperature)/(nu)^2; %Grashof number
    Ra = Pr*Gr; %Rayleigh number
    
    if Ra < 1e7
        Nu = 0.54*Ra^(1/4);
        h = Nu*Air().gasThermalConductivity(T)/D;
    else
        Nu = 0.15*Ra^(1/3);
        h = Nu*Air().gasThermalConductivity(T)/D;
    end
    
    Cpa = Air().gasHeatCapacity(T);
    Y = sum(mdot_fuel)*(Cpa/h)/A;
    
%     %Exclude Water
%     if speciesInMixture(Water(),mixture) == (true)
%         Hc = mixture.species(1).heatCombustion;
%         Xa = mixture.species(1).combustionEfficiency;   
%         Xr = mixture.species(1).radiationFrac;          
%         r = mixture.species(1).airFuelRatio; 
%         Y = mdot_fuel(1)*(Cpa/h)/A;
%     end
    
    Q = A*h/Cpa*(Hc*(Xa-Xr)*r/Xa - Cpa*(Ts-Air().referenceTemperature))*Y/(exp(Y)-1);
    %Q = zeros(length(mixture.species),1);
%     for i = 1:length(mixture.species)
%         %Ignore inert fuels
%         if Hc(i) == 0 || mdot_fuel(i) < 10^-8
%             Q(i) = 0;
%         else
%             if Y(i)/exp(Y(i)) == 0 
%                 warning('mdot(i)/(exp(mdot(i))-1) == Inf, reduce mdot(i)');
%             end
%             Q(i) = A*h/Cpa*(Hc(i)*(Xa(i)-Xr(i))*r(i)/Xa(i) - Cpa*(Ts-Air().referenceTemperature))*Y(i)/(exp(Y(i))-1);
%         end
%     end
    Qconv = Q;
end