function [A,B,C]=inletCoefficients(uavg,a,b,time,startTime)
% Calculates the coefficients of the parabola of the inlet 
% velocity profile
if time>=startTime
    avector=[0 b^2 0.5*(b^2-a^2)-(b^3-a^3)/(3*a)];
    bvector=[a^2-b^2 b b-a-(b^3-a^3)/(3*a^2)];
    cvector=[a-b 1 0];
    dvector=[0 0 uavg*(b-a)];
    coeffs=GaussTriSol(avector,bvector,cvector,dvector);
    A=coeffs(1);
    B=coeffs(2);
    C=coeffs(3);
else
    A=0;
    B=0;
    C=0;
end

end