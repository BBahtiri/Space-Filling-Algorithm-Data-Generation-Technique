function [output] = toMat_full(A)
output = zeros(3,3);
output(1,1) = A(1);
output(2,1) = A(2);
output(3,1) = A(3);
output(1,2) = A(4);
output(2,2) = A(5);
output(3,2) = A(6);
output(1,3) = A(7);
output(2,3) = A(8);
output(3,3) = A(9);
end

