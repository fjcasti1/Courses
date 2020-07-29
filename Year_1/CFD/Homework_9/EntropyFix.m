function phi=EntropyFix(y)
eps=0.1;
phi=abs(y).*(abs(y)>=eps)+((y.^2+eps^2)./(2*eps)).*(abs(y)<eps);
end