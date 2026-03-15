function y = FranckCondonElementRatio(n,lambda)
% Gn/G0 at eq

n = round(n);
y = lambda.^(2*n)./(factorial(n));%2*

end