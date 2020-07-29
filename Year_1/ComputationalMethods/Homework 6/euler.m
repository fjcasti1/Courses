function [t,w] = euler(dydt,tspan,y0,N)

h = diff(tspan)/N;
t = tspan(1) + h*(0:N)';
w = zeros(length(y0),N+1);
w(:,1) = y0(:).';

for i = 1:N
    k1 = h*dydt(t(i)    ,w(:,i)      );
	w(:,i+1) = w(:,i) + k1;
end
w = w.';
end