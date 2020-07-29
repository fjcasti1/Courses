% Run different simulations in parallel
N=[16 32 64 16 32 64];
M=N/2;
CFL=[0.8 0.8 0.8 0.9 0.9 0.9];
movie=0;

design=1;

parfor (i=1:length(M),4)
    Problem2bonus(M(i),N(i),CFL(i),movie,design)
end

disp(' ')
disp(' ')
disp('FIRST PARFOR LOOP COMPLETED ')
disp(' ')
disp(' ')

N=[128 256 128 256];
M=N/2;
CFL=[0.8 0.8 0.9 0.9];
movie=0;

design=1;

parfor (i=1:length(M),4)
    Problem2bonus(M(i),N(i),CFL(i),movie,design)
end

disp(' ')
disp(' ')
disp('SECOND PARFOR LOOP COMPLETED ')
disp(' ')
disp(' ')

N=[512 512];
M=N/2;
CFL=[0.8 0.9];
movie=0;

design=1;

parfor (i=1:length(M),4)
    Problem2bonus(M(i),N(i),CFL(i),movie,design)
end

disp(' ')
disp(' ')
disp('THIRD PARFOR LOOP COMPLETED ')
disp(' ')
disp(' ')

