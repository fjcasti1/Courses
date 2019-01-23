function phileft = phi_left(t)
if 0<=t && t<=0.25
    phileft=1;
elseif 0.25<t && t<=0.5
    phileft=0;
elseif 0.5<t && t<=1
    phileft=(1-cos(4*pi*t))/2;
elseif t>1
    phileft=0;
else
    error('t negative')
end
end