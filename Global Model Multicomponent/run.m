clear;
addpath(genpath(pwd)); %adds all subfolders to path

%Run
[mixture,t,res] = main();

%Indexing for output
m_1 = 4;
nSpecies = length(mixture.species); %length(res(1,m_1:end))/3;

m_end = m_1 + nSpecies - 1;
y_1 = m_end + 1;
y_end = y_1 + nSpecies - 1;
x_1 = y_end + 1;
x_end = x_1 + nSpecies - 1;

%Data
mdot_tot = res(:,1)*1000; %g/s
m_tot = res(:,2)*1000; %g/s
Tf = res(:,3);

m = res(:,m_1:m_end)*1000; %g
y = res(:,y_1:y_end);
x = res(:,x_1:x_end);
mdot = mdot_tot.*y;

Q = zeros(length(t),length(mixture.species));
Hc = mixture.heatsCombustion;
for i = 1:nSpecies
    Q(:,i) = mdot(:,i)*1000.*y(:,i)*Hc(i);
end
Q_tot = sum(Q,2); %sum along row

T = table(t,mdot_tot,Q_tot,mdot,m_tot,m,x,y,Tf);

figure(1)
plot(t,mdot_tot)
xlabel('Time (s)'); ylabel('Burning Rate (g/s)');

figure(2)
plot(t,Q_tot)
xlabel('Time (s)'); ylabel('HRR (W)');