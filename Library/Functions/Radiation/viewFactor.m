function F = viewFactor(D,L)
%VIEWFACTOR calculates the view factors for a cylinder. If the top, middle,
%and bottom are labelled 1,2,3, then this function outputs F = F13 + F23. 
%D: fire diameter (m)
%L: flame height (m)
A1 = pi*D^2/4;
A2 = pi*D*L;
A3 = A1;
R = D/2;
Ri = R/L;
Rj = R/L;
S = 1 + (1+Rj^2)/Ri^2;

%View factors
F13 = 0.5*(S - (S^2-4*(Rj/Ri)^2)^0.5);
F12 = 1 - F13; %summation rule (F11 = 0)
F21 = A1/A2*F12; %reciprocity rule
F23 = F21; %symmetry rule
F = F13 + F23;
end