function res = averageProperty(mixture_property,mdot_fuel)
    %Takes a mixture property and fuel mass flow rate as input
    %Outputs mass averaged property
    
    y = zeros(1,length(mixture_property));
    for i = 1:length(mixture_property)
        y(i) = mdot_fuel(i)/sum(mdot_fuel);
    end
    
    res = sum(y.*mixture_property);
end