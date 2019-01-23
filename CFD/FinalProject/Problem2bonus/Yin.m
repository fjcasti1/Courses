function Yin(time)
global Yin1; global Yin2;
global Yin3; global period;
global nperiod;

if time>=nperiod*period
    Yin1=abs(Yin1-1);
    Yin2=abs(Yin2-1);
    Yin3=abs(Yin3-1);
    nperiod=nperiod+1;
end
end