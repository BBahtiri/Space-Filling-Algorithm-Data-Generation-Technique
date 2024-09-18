function [a0,vf_ratio] = orientation_vector(A)

[a0,lambda] = eig(A);
vf_ratio    = diag(lambda);
