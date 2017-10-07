%% Numerical Solution
ts = [0, 0.01, 0.02, 0.06, 0.1];
dx = 0.01; dy = 0.01; dt = 0.0001;
[U, X, Y] = solver_ex1(dx, dy, dt);
plotSlice_ex1(ts, dt, U, X, Y);
%% Code: Solver of Function
% 
%   function [U, X_, Y_] = solver_ex1( dx, dy, dt )
%       A1 = [0, 4/3, -2/3; 1, -2/3, -5/3; -1, -4/3, -1/3];
%       A2 = [5/3, 1, 5/3; -1/3, -1, -7/3; 1/3, -1, 1/3];
%       X = round(1 / dx);
%       Y = round(1 / dy);
%       T = round(1 / dt);
%       x = 1; y = 1; t = 0.1;
%       Mx = x * X; My = y * Y; N = t * T;
%       v1_ = @(X, Y) (sin(pi * X) .* sin(pi / 2 * Y));
%       v2_ = @(X, Y) (sin(2 * pi * X));
%       v3_ = @(X, Y) (cos(2 * pi * X) .* cos(2 * pi * Y));
%       b0_ = @(X, Y) (cat(3, v1_(X, Y), v2_(X, Y), v3_(X, Y)));
%       A1_s = A1 * A1; A2_s = A2 * A2; 
%       rx = dt / dx;
%       ry = dt / dy;
%       x_ = 0 : dx : x;
%       y_ = 0 : dy : y;
%       t_ = 0 : dt : t;
%       [X_, Y_] = meshgrid(x_, y_);
%       U = zeros([size(X_), 3, N + 1]);
%       U(:, :, :, 1) = b0_(X_, Y_);
%       % U_t = zeros(size(U(2 : end - 1, 2 : end - 1,:, 1)));
%       x__ = x_(2 : end - 1); y__ = y_(2 : end - 1);
%       dy_0 = ry / 2 * speye(length(y__));
%       dy_1 = ry / 2 * ones(1, length(y__) - 1);
%       ryDy1 = sparse(diag(dy_1, -1)) + sparse(diag(-dy_1, 1));
%       ryDy2 = -2 * dy_0 + sparse(diag(dy_1, 1)) + sparse(diag(dy_1, -1));
%       ryDy1(1, 1) = -ryDy1(1, 2); ryDy1(end, end) = -ryDy1(end, end - 1);
%       ryDy2(1, 1) = -ryDy2(1, 2); ryDy2(end, end) = -ryDy2(end, end - 1);
%       dx_0 = rx / 2 * speye(length(x__));
%       dx_1 = rx / 2 * ones(1, length(x__) - 1);
%       rxDx1 = sparse(diag(dx_1, -1)) + sparse(diag(-dx_1, 1));
%       rxDx2 = -2 * dx_0 + sparse(diag(dx_1, 1)) + sparse(diag(dx_1, -1));
%       rxDx1(1, 1) = -rxDx1(1, 2); rxDx1(end, end) = -rxDx1(end, end - 1);
%       rxDx2(1, 1) = -rxDx2(1, 2); rxDx2(end, end) = -rxDx2(end, end - 1);
%       % half, U_h := U_0 + A1 $ [U_0 * Lx1] + A1^2 $ [U_0 * Lx2]
%       % next, U_t := U_h + A2 $ [Ly1 * U_h] + A2^2 $ [Ly2 * U_h]
%       Lx1 = rxDx1; Lx2 = rx * rxDx2;
%       Ly1 = ryDy1; Ly2 = ry * ryDy2;
%       coef_x = @(u, Lx) cat(3, u(:, :, 1) * Lx, u(:, :, 2) * Lx, u(:, :, 3) * Lx);
%       coef_y = @(u, Ly) cat(3, Ly * u(:, :, 1), Ly * u(:, :, 2), Ly * u(:, :, 3));
%       coef_x1 = @(u) coef_x(u, Lx1); coef_x2 = @(u) coef_x(u, Lx2);
%       coef_y1 = @(u) coef_y(u, Ly1); coef_y2 = @(u) coef_y(u, Ly2);
%       for kk = 1 : N
%           U__ = U(:, :, :, kk);
%           U_ = U__(2 : end - 1, 2 : end - 1, :);
%           tensor_A1_Lx_1 = -tensorproduct_AX(A1, coef_x1(U_));
%           tensor_A1_Lx_2 = tensorproduct_AX(A1_s, coef_x2(U_));
%           U_t = U_ + tensor_A1_Lx_1 + tensor_A1_Lx_2;
%           tensor_A2_Ly_1 = -tensorproduct_AX(A2, coef_y1(U_t));
%           tensor_A2_Ly_2 = tensorproduct_AX(A2_s, coef_y2(U_t));
%           U_t = U_t + tensor_A2_Ly_1 + tensor_A2_Ly_2;
%           U(2 : end - 1, 2 : end - 1, :, kk + 1) = U_t;
%           U(:, end, 1, kk + 1) = sin(t_(kk + 1));
%           U(1, :, 1, kk + 1) = sin(2 * t_(kk + 1));
%           U(end, :, 1, kk + 1) = sin(pi * x_);
%           U(:, 2, 1, kk + 1) = U(:, 1, 1, kk + 1);
%           U(:, 1, 2, kk + 1) = sin(4 * t_(kk + 1));
%           U(:, end - 1, 2, kk + 1) = U(:, end, 2, kk + 1);
%           U(2, :, 2, kk + 1) = U(1, :, 2, kk + 1);
%           U(end - 1, :, 2, kk + 1) = U(end, :, 2, kk + 1);
%           U(:, 1, 3, kk + 1) = cos(2 * t_(kk + 1)) .* cos(2 * pi * y_);
%           U(end, :, 3, kk + 1) = cos(2 * t_(kk + 1)) .* cos(2 * pi * x_);
%           U(:, end - 1, 3, kk + 1) = U(:, end, 3, kk + 1);
%           U(2, :, 3, kk + 1) = U(1, :, 3, kk + 1);
%       end
%   end
% 
%% Code: Value of Slices
% 
%   function plotSlice_ex1( t_slices, dt, u, X_, Y_ )
%       t_slices_id = round(t_slices / dt) + 1;
%       for ii = 1 : length(t_slices)
%           figure; hold on;
%           v = u(:, :, :, t_slices_id(ii));
%           subplot(1, 3, 1); mesh(X_, Y_, v(:, :, 1));
%           view(3); xlabel('x'); ylabel('y'); zlabel('v_1(numerical)');
%           subplot(1, 3, 2); mesh(X_, Y_, v(:, :, 2));
%           view(3); xlabel('x'); ylabel('y'); zlabel('v_2(numerical)');
%           subplot(1, 3, 3); mesh(X_, Y_, v(:, :, 3));
%           view(3); xlabel('x'); ylabel('y'); zlabel('v_3(numerical)');
%           suptitle(sprintf('t = %.2f', t_slices(ii)));
%       end
%   end
% 
%% Code: tensor product Y = AX
% 
%   function Y = tensorproduct_AX(A, X)
%   %tensorproduct - Description
%   %
%   % Syntax: Y = tensorproduct_AX(A, X)
%   % A is a [3x3] matrix, X is a [:, :, 3] tensor, Y = AX, which means that taking
%   % X for a [3x1] vector like this [x_1; x_2; x_3], and Y will be a [:, :, 3]
%   % tensor with the same shape like X.
%       f = @(a) (a(1) * X(:, :, 1) + a(2) * X(:, :, 2) + a(3) * X(:, :, 3));
%       lambda = @(idx) f(A(idx, :));
%       Y = cat(3, lambda(1), lambda(2), lambda(3));
%   end
% 