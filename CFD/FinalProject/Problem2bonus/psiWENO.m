function psi = psiWENO(a,b,c,d)
eps=1e-6;
IS0=13*(a-b)^2+3*(a-3*b)^2;
IS1=13*(b-c)^2+3*(b+c)^2;
IS2=13*(c-d)^2+3*(3*c-d)^2;

a0=(eps+IS0)^(-2);
a1=6*(eps+IS1)^(-2);
a2=3*(eps+IS2)^(-2);

w0=a0/(a0+a1+a2);
w2=a2/(a0+a1+a2);

psi = (a-2*b+c)*w0/3+(w2-0.5)*(b-2*c+d)/6;
end