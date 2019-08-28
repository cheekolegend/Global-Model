function y = solveY(mixture,x_mol,P)
    Ts = bubbleT(mixture,x_mol,P);
    gamma = wilson(mixture,x_mol,Ts);
    
    y_mol = zeros(length(mixture.species),1);
    for i = 1:length(mixture.species)
        y_mol(i) = x_mol(i)*mixture.species(i).vapourPressure(Ts)*gamma(i)/P;
    end
    y = molesToMass(mixture,y_mol);
end