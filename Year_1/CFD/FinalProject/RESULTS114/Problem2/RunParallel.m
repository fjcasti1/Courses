% Run different simulations in parallel
N=[16,32,64,16,32,64]
M=N/2;
CFL=[0.8 0.8 0.8 0.9 0.9 0.9];
movie=0;

design=1;

parfor (i=1:length(M),4)
    Problem2(M(i),N(i),CFL(i),movie,design)
end

disp(' ')
disp(' ')
disp('FIRST PARFOR LOOP COMPLETED')
disp(' ')
disp(' ')

N=[256,256,512,512,1024,1024]
M=N/2;
CFL=[0.8 0.9 0.8 0.9 0.8 0.9];
M=N/2;
movie=0;

parfor (i=1:length(M),4)
    Problem2(M(i),N(i),CFL(i),movie,design)
end


disp(' ')
disp(' ')
disp('SECOND PARFOR LOOP COMPLETED')
disp(' ')
disp(' ')

