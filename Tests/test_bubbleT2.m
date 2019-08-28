species1 = Ethanol();
species2 = Heptane();
species = [species1,species2]; nSpecies = length(species);
massFracs = [0.8204,1-0.8204];
mass = 0.7771; %kg

x = 0:0.01:1;
T = zeros(length(x),1);
for i = 1:length(x)
    T(i) = bubbleT(mixture,[x(i),1-x(i)],101325);
end
plot(x,T)