function Y = tensorproduct_AX(A, X)
%tensorproduct - Description
%
% Syntax: Y = tensorproduct_AX(A, X)
% A is a [3x3] matrix, X is a [:, :, 3] tensor, Y = AX, which means that taking
% X for a [3x1] vector like this [x_1; x_2; x_3], and Y will be a [:, :, 3]
% tensor with the same shape like X.
    f = @(a) (a(1) * X(:, :, 1) + a(2) * X(:, :, 2) + a(3) * X(:, :, 3));
    lambda = @(idx) f(A(idx, :));
    Y = cat(3, lambda(1), lambda(2), lambda(3));
end