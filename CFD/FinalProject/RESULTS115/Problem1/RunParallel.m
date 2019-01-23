% Run different simulations in parallel
M=[16,32,64,128,256,512,1024];
N=M/2;
CFL=0.8;
movie=0;

parfor (i=1:length(M),4)
    Problem1(M(i),N(i),CFL,movie)
end
