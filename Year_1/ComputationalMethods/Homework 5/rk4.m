function [t,w] = rk4(dydt,tspan,y0,N)

h = diff(tspan)/N;
t = tspan(1) + h*(0:N)';
w = zeros(length(y0),N+1);
w(:,1) = y0(:).';

for i = 1:N
    k1 = h*dydt(t(i)    ,w(:,i)      );
    k2 = h*dydt(t(i)+h/2, w(:,i)+k1/2);
    k3 = h*dydt(t(i)+h/2, w(:,i)+k2/2);
    k4 = h*dydt(t(i)+h  , w(:,i)+k3  );
    w(:,i+1) = w(:,i) + (k1 + 2*k2 + 2*k3 + k4)/6;
end
w = w.';
end