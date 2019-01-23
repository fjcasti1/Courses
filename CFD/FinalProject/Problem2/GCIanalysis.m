clear all; close all; format long; clc
CFL=0.8;
% path=['RESULTS_CFL0',num2str(10*CFL),'/'];
path=['RESULTS_Design01/'];
DATA=dir([path,'GCIdata*.mat']);

for i=1:length(DATA)
        data=load([path DATA(i).name]);
        auxM(i)=data.M;
        auxN(i)=data.N;
        auxT(i)=data.T;
end
M=sort(auxM);
N=sort(auxN);
for i=1:length(DATA)
    T(i)=auxT(auxM==M(i));
end

tableMatrix=zeros(length(DATA),3);
tableMatrix(:,1)=M;
tableMatrix(:,2)=N;
tableMatrix(:,3)=T;


tableMatrix=updateTable(tableMatrix);
table=array2table(tableMatrix,'VariableNames',{'M','N','T','p','T_h0','GCI12','GCI23','beta'});

%% Save GCI results
% In mat format
txt=[path 'GCIresults_CFL0' num2str(10*CFL)];
save(txt,'table')
% In txt format
writetable(table,txt)