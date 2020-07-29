function [t,w] = modEuler(dydt,tspan,y0,N)

h = diff(tspan)/N;
t = tspan(1) + h*(0:N)';
w = zeros(length(y0),N+1);
w(:,1) = y0(:).';

for i = 1:N
    k1 = h*dydt(t(i)    ,w(:,i)      );
    k2 = h*dydt(t(i)+h, w(:,i)+k1);
    w(:,i+1) = w(:,i) + (k1 + k2)/2;
end
w = w.';
end