% Run different simulations in parallel
N=[16,32,64,128,256,512,1024];
M=N/2;
CFL=0.8;
movie=[0 0 0 1 0 0 0];

design=1;

parfor (i=1:length(M),4)
    Problem2(M(i),N(i),CFL,movie(i),design)
end

