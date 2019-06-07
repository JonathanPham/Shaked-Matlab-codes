tol=10^-8;
func=@(x) exp(-x)-x^2;
der=@(x)  -exp(-x)-2*x;
x=1/2;
while abs(func(x))>tol
    x=x-func(x)/der(x);
end