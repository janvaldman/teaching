syms x f

f=sqrt(x^2+exp(x))+cos(x^2+exp(x));
Df=diff(f);
DDf=diff(f,2);

hf = matlabFunction(f);
hDf = matlabFunction(Df);
hDDf = matlabFunction(DDf);


a=-3; b=2;
xx=linspace(a,b,1000);
plot(xx,hf(xx))
hold on
plot(xx,hDf(xx));
plot(xx,hDDf(xx));

hold off
enlarge_axis

legend('f(x)','Df(x)','DDf(x)')








