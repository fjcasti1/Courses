function T=calculateT(S,t)
tol=1e-13;
% Use trapezoidal integration
    T=0;
    a=find(abs(t-10)<tol);
    b=find(abs(t-20)<tol)-1;
    for i=a:b
        T=T+0.5*(t(i+1)-t(i))*(S(i+1)+S(i));
    end
    T=T/10;
end