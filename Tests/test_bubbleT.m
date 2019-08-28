mixture = Mixture([Ethanol(),Water()],[1,1]);
x_mol = [0.99,0.01];
P = 101325;

Ts = bubbleT(mixture,x_mol,P)