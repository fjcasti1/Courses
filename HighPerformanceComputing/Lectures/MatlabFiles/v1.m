A = randn(1000);
B = randn(1000);
C = zeros(1000);
tic;
for k=1:1000
    for j=1:1000
        C(j,k) = A(j,k) + B(j,k);
    end
end
toc