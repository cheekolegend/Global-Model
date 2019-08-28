function Xa = calcCombustionEfficiency(Lsp)
%CALCRADIATIONEFFICIENCY estimates the radiation fraction using the
%correlation of Tewarson (1986) (see pg. 1212 in SFPE Handbook 5th ed)

%Lsp: laminar smoke point of pure fuel (m)

Xa = 1.15*Lsp^0.1;
end