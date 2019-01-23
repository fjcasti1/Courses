clear all; close all; format long; clc
path='DataSavedSteady/GCIdata_Steady/';
DATA=dir('DataSavedSteady/GCIdata_Steady/*.mat');

% For normal simulation M up to 512
% % sortDATA(1)=DATA(2);
% % sortDATA(2)=DATA(4);
% % sortDATA(3)=DATA(6);
% % sortDATA(4)=DATA(1);
% % sortDATA(5)=DATA(3);
% % sortDATA(6)=DATA(5);

% For steady state simulation M up to 256
sortDATA(1)=DATA(2);
sortDATA(2)=DATA(4);
sortDATA(3)=DATA(5);
sortDATA(4)=DATA(1);
sortDATA(5)=DATA(3);

Tu5=nan(length(DATA),7);
Tv5=nan(length(DATA),7);
TY5=nan(length(DATA),7);
TR5=nan(length(DATA),7);

Tuend=nan(length(DATA),7);
Tvend=nan(length(DATA),7);
TYend=nan(length(DATA),7);
TRend=nan(length(DATA),7);

put5Vector=zeros(length(DATA),1);
pvt5Vector=zeros(length(DATA),1);
pYt5Vector=zeros(length(DATA),1);
Rt5Vector=zeros(length(DATA),1);
% put10Vector=zeros(length(DATA),1);
% pvt10Vector=zeros(length(DATA),1);
% pYt10Vector=zeros(length(DATA),1);
% Rt10Vector=zeros(length(DATA),1);
putendVector=zeros(length(DATA),1);
pvtendVector=zeros(length(DATA),1);
pYtendVector=zeros(length(DATA),1);
RtendVector=zeros(length(DATA),1);

% This because I didn't save M and N (IDIOT)
Mv=[16 32 64 128 256];
Nv=Mv/2;
for i=1:length(Mv)
M=Mv(i);
N=Nv(i);

Tu5(i,1)=M;  Tu5(i,2)=N;
Tv5(i,1)=M;  Tv5(i,2)=N;
TY5(i,1)=M;  TY5(i,2)=N;
TR5(i,1)=M;  TR5(i,2)=N;

Tuend(i,1)=M;  Tuend(i,2)=N;
Tvend(i,1)=M;  Tvend(i,2)=N;
TYend(i,1)=M;  TYend(i,2)=N;
TRend(i,1)=M;  TRend(i,2)=N;
end

for i=1:length(DATA)
    name=sortDATA(i).name; 
    txt=[path, name]; load(txt)
    
    put5Vector(i)=put5;
    pvt5Vector(i)=pvt5;
    pYt5Vector(i)=pYt5;
    Rt5Vector(i)=Rt5;
    
    putendVector(i)=putend;
    pvtendVector(i)=pvtend;
    pYtendVector(i)=pYtend;
    RtendVector(i)=Rtend;
    
    if i>=3
        [Tu5,Tv5,TY5,TR5]=updateTable(put5Vector,pvt5Vector,pYt5Vector,Rt5Vector,i,Tu5,Tv5,TY5,TR5);
        [Tuend,Tvend,TYend,TRend]=updateTable(putendVector,pvtendVector,pYtendVector,RtendVector,i,Tuend,Tvend,TYend,TRend);
    end
end
Tu5=array2table(Tu5,'VariableNames',{'M','N','p','u_h0','GCI12','GCI23','beta'});
Tv5=array2table(Tv5,'VariableNames',{'M','N','p','v_h0','GCI12','GCI23','beta'});
TY5=array2table(TY5,'VariableNames',{'M','N','p','Y_h0','GCI12','GCI23','beta'});
TR5=array2table(TR5,'VariableNames',{'M','N','p','R_h0','GCI12','GCI23','beta'});

Tuend=array2table(Tuend,'VariableNames',{'M','N','p','u_h0','GCI12','GCI23','beta'});
Tvend=array2table(Tvend,'VariableNames',{'M','N','p','v_h0','GCI12','GCI23','beta'});
TYend=array2table(TYend,'VariableNames',{'M','N','p','Y_h0','GCI12','GCI23','beta'});
TRend=array2table(TRend,'VariableNames',{'M','N','p','R_h0','GCI12','GCI23','beta'});
%% Save GCI results
txt='GCIresults_Steady';
save(txt,'Tu5','Tv5','TY5','TR5','Tuend','Tvend','TYend','TRend')