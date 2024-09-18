function dev = dev(A)

I = eye(size(A,1));
dev = A - (1/3) * trace(A) * I;