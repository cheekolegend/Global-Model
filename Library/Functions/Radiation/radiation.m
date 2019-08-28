function Qrad = radiation(mixture,mdot_fuel,products,mdot_prod,Tf,Ts,D)
%RADIATION calculates radiation heat transfer from a flame of
%height L to the surface of a liquid pool of diameter D.
%mixture: evaporating fuel mixture ojbect
%mdot_fuel: fuel evaporation rate of each species, kg/s
%products: combustion products mixture object
%mdot_prod: mass flow rate of all combustion products, kg/s
%Tf: flame temperature, K
%Ts: pool surface temperature, K
%D: pool diameter, m

    A = pi*D^2/4;
    sig = 5.67e-8; %W m-2 K-4
    p = 0.08; %surface reflectivity of liquid fuels
    L = flameHeight(mixture,mdot_fuel,D);
    fv = averageProperty(mixture.sootVolFracs,mdot_fuel);
    F = viewFactor(D,L);
    
    %Ensure view factor is <= 1
    if F > 1
        F = 1;
    end
    
    %Emissivity of soot calculation
    n = 3.49; s = 2.17; %infrared optical averaged constants
    c1 = 36*pi*fv*n^2*s/((n^2 - (n*s)^2 + 2)^2 + 4*n^4*s^2); 
    c2 = 0.014388; %Plank's second constant, m K
    kappa = 3.6*c1*Tf/c2; %effective soot emission parameter
    Af = pi*D*(D/2 + L); %flame area
    V = pi*D^2/4*L; %flame volume
    Lm = 3.6*V/Af; %mean beam length
    eps_soot = 1 - exp(-kappa*Lm);
    
    %Emmissivities for H2O and CO2
    [pCO2,pH2O] = CO2H2OpartialPressure(products,mdot_prod);
    eps_w = epsWater(Tf,pH2O,Lm);
    eps_c = epsCO2(Tf,pCO2,Lm);
    
    Qrad = F*sig*(Tf^4 - Ts^4)*A*(1 - p)*(1 - (1 - eps_soot)*(1 - eps_w)*(1 - eps_c)); %W
end