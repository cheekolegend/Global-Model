function Qrerad = reradiation(Ts,D)
%reradiation(ambient,Ts,D) calculates reradiation from a flame against a
%liquid pool surface
    sig = 5.67e-8;
    A = pi*D^2/4;
    Qrerad = sig*A*(Ts^4 - Air().referenceTemperature^4);
end