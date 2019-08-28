clear; clc;
addpath(genpath(pwd)); %adds all subfolders to path

species1 = Ethanol();
species2 = Heptane();
species = [species1,species2]; nSpecies = length(species);
mixture = Mixture(species,[1,1]);

P = 101325;

xtemp = (0:0.01:1);
x = zeros(length(xtemp),2); x(:,1) = xtemp; x(:,2) = 1 - xtemp;
y = zeros(length(xtemp),2);
Ts = zeros(length(xtemp),1);
gamma = zeros(length(xtemp),2);

for i = 1:length(x)
    Ts(i) = bubbleT(mixture,x(i,:),P);
    gamma(i,:) = wilson(mixture,x(i,:),Ts(i));
    y(i,1) = x(i,1)*gamma(i,1)*mixture.species(1).vapourPressure(Ts(i))/P;
    y(i,2) = x(i,2)*gamma(i,2)*mixture.species(2).vapourPressure(Ts(i))/P;
end

figure(1);
plot(x(:,1),Ts,'-o',y(:,1),Ts,'-*');
figure(2);
plot(x(:,1),x(:,1),'-',x(:,1),y(:,1),'-*');
figure(3);
plot(x(:,1),gamma(:,1),'-',x(:,1),gamma(:,2),'-*');