function sum = sumEq(y,equation)
%sumEq computes a sum for y_i*equation(a_i,b_i,...)
%  

sum = 0;
for i = 1:length(y)
    sum = sum + y(i)*equation(i);
end

end

