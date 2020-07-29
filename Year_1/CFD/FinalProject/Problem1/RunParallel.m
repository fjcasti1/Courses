% Run different simulations in parallel
% M=[16,32,64,128,256,512];
% N=[8,16,32,64,128,256];
% CFL=0.8;
% movie=0;
% 
% parfor (i=1:length(M),4)
%     Problem1(M(i),N(i),CFL,movie)
% end

M=256;
N=128;
CFL=0.8;
movie=1;

parfor (i=1:length(M),4)
    Problem1(M(i),N(i),CFL,movie)
end