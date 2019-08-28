function [mixture,t,y] = main()
%Miscellaneous startup tasks
clear;
addpath(genpath(pwd)); %adds all subfolders to path

%User Inputs
species1 = Ethanol();
species2 = isoPropanol();
species3 = Water();
species = [species1,species2,species3]; nSpecies = length(species);
massFracs = [0.5,0.25,0.25];
mass = 1; %kg
%Time span
t_end = 1000; %s
dt = 1; %s
tSpan = 0:dt:t_end;

%Creating mixture objects
mixture = Mixture(species,massFracs);
reactants = createReactants(mixture);
products = createProducts(mixture);

%Ambient conditions
P = 101325; %Pa

%Tray dimensions
D = 0.33; %m
H = 0.05; %m
A = pi*D^2/4; %m2

%Entrainment assumptions
yO2 = 0.1; %limiting oxygen mol fraction in fire

%Initial conditions
T0 = 293.15;
mdot_0 = 0.001; %kg/s
m_0 = mass; %kg
Tf_0 = 1500; %K
x_0 = zeros(nSpecies,1);
for i = 1:nSpecies
    x_0(i) = mixture.x(i);
end
y_0 = solveY(mixture,massToMoles(mixture,x_0),P);
mi_0 = m_0*x_0;

%Indexing for ease of use in writing equations
m_1 = 4; %m1 is the fourth variable. The list order is found in var0
m_end = m_1 + nSpecies - 1;
y_1 = m_end + 1;
y_end = y_1 + nSpecies - 1;
x_1 = y_end + 1;
x_end = x_1 + nSpecies - 1;

%Solving algorithm
var0 = [mdot_0;m_0;Tf_0;mi_0;y_0;x_0];
M = zeros(length(var0),length(var0));
index0 = 4; %Index to start listing equations dependent on number of species
    %Define mass matrix
    M(2,2) = 1; %dm/dt = -mdot has a derivative
    j = 0;
    for a = m_1:m_end
        M(index0+j,a) = 1; %dm_i/dt = -mdot*y_i has a derivative
        j = j+1;
    end
opt = odeset('Mass', M, 'InitialStep', 1e-2, 'AbsTol',1e-8,'Events', @myEvent);
[t,y] = ode15s(@func,tSpan,var0,opt);

    function res = func(~,var)
        res = zeros(length(var),1);
        mdot = var(1);
        m = var(2);
        Tf = var(3);
        
        %Declaring miscellaneous variables
        x_mol = massToMoles(mixture,var(x_1:x_end));
        y_mol = massToMoles(mixture,var(y_1:y_end));
        Ts = bubbleT(mixture,x_mol,P);
        gamma = wilson(mixture,x_mol,Ts);
        mdot_fuel = mdot*var(y_1:y_end);
        mdot_air = mdotAir(mixture,mdot_fuel,yO2);
        mdotReact = [mdot_fuel',Air().xO2*mdot_air,Air().xN2*mdot_air]; %Row vector
        mdotProd = mdotProducts(mixture,mdot_fuel,mdot_air);
        
        %Equations
        res(1) = mdot*sum((mixture.latentHeatsVap(Ts)' + mixture.liquidHeatCapacityIntegrals(T0,Ts)').*var(x_1:x_end)) - convection(mixture,mdot_fuel,var(y_1:y_end),Tf,Ts,D) - radiation(mixture,mdot_fuel,products,mdotProd,Tf,Ts,D) + reradiation(Ts,D);
        res(2) = -mdot; %dm/dt = -mdot
        res(3) = flameTemperature(mixture,mdot_fuel,reactants,mdotReact,products,mdotProd,T0,Ts,Tf);

        res_index = index0;
        
        %mdot_i = -mdot*y_i
        j = 0;
        for i = res_index:res_index+nSpecies-1
            res(i) = -mdot*var(y_1+j);
            j = j+1;
        end
        res_index = i+1;

        %0 = y_i*P - x_i*Psat_i
        j = 0;
        for i = res_index:res_index+nSpecies-2
            res(i) = y_mol(j+1)*P - x_mol(j+1)*gamma(j+1)*mixture.species(j+1).vapourPressure(Ts);
            j = j+1;
        end
        res_index = i+1;
        
        %0 = x_i*m - m_i (x1 = m1/m)
        j = 0;
        for i = res_index:res_index+nSpecies-2
            res(i) = var(x_1+j)*m - var(m_1+j);
            j = j+1;
        end        
        
        res(end-1) = 1 - sum(var(y_1:y_end));
        res(end) = 1 - sum(var(x_1:x_end));
        
        %Check to prevent negative mass fractions
        for i = 0:nSpecies-1
            if var(x_1 + i) <= 0
                var(x_1 + i) = 0;
                var(m_1 + i) = 0;
            end
        end       
    end
    function [value,isterminal,direction] = myEvent(~,y)
        value = y(2); %the value we want to be zero
        isterminal = 1; %Halt integration
        direction = []; %The zero can be approached from either direction
    end
end


