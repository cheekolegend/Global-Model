function Xr = calcRadiationFraction(Lsp)
%CALCRADIATIONFRACTION estimates the radiation fraction using the
%correlation of Tewarson (1986) (see pg. 1212 in SFPE Handbook 5th ed)

%Lsp: laminar smoke point of pure fuel (m)

Xr = 0.41 - 0.85*Lsp;
end