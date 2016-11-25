%% HW 4.4.8
%% Initial-boundary-value Problem
% * $v_t = \nu (v_{xx} + v_{yy}), (x, y) \in (0, 1) \times (0, 1), t > 0$
% * $v(0, y, t) = v(1, y, t) = 0, y \in [0, 1], t \geq 0$
% * $v(x, 0, t) = v(x, 1, t) = 0, x \in [0, 1], t \geq 0$
% * $v(x, y, 0) = \sin \pi x \sin 2 \pi y, (x, y) \in [0, 1] \times [0, 1]$
%%%
%% Numerical Solution(Wrong Mesh Grid)
%%%
ts = [0.06, 0.1, 0.9];
dx = 0.05; dy = 0.05; dt = 0.01;
[u, X, Y] = solver_4_4_8(dx, dy, dt);
plotSlice_4_4_3(ts, dt, u, X, Y);
%% Numerical Solution(Correct Mesh Grid)
%%%
ts = [0.06, 0.1, 0.9];
dx = 0.005; dy = 0.005; dt = 0.001;
[u, X, Y] = solver_4_4_8(dx, dy, dt);
plotSlice_4_4_3(ts, dt, u, X, Y);

%% Code: Solver of Functionn
%%%
%
%   function [u, X_, Y_] = solver_4_4_8( dx, dy, dt )
%       f_ = @(X, Y)(sin(pi * X) .* sin(2 * pi * Y));
%       nu = 1;
%   
%       X = round(1 / dx);
%       Y = round(1 / dy);
%       T = round(1 / dt);
%       x = 1; y = 1; t = 1;
%       Mx = x * X; My = y * Y; N = t * T;
%   
%       rx = nu * dt / dx / dx;
%       ry = nu * dt / dy / dy;
%   
%       x_ = 0 : dx : x;
%       y_ = 0 : dy : y;
%       [X_, Y_] = meshgrid(x_, y_);
%       u = zeros([size(X_), N + 1]);
%       u(:, :, 1) = f_(X_, Y_);
%   
%       u_t = zeros(size(u(2 : end - 1, 2 : end - 1, 1)));
%       u_h = zeros(size(u_t));
%   
%       x__ = x_(2 : end - 1); y__ = y_(2 : end - 1);
%       [X__, Y__] = meshgrid(x__, y__);
%   
%       A_0 = ry * speye(length(y__));
%       A_1 = ry * ones(1, length(y__) - 1);
%       A = -2 * A_0 + sparse(diag(A_1, 1)) + sparse(diag(A_1, -1));
%   
%       B_0 = rx * speye(length(x__));
%       B_1 = rx * ones(1, length(x__) - 1);
%       B = -2 * B_0 + sparse(diag(B_1, 1)) + sparse(diag(B_1, -1));
%   
%       coef_Ix = speye(length(x__));
%       coef_Iy = speye(length(y__));
%   
%       coef_mrdx2 = coef_Ix - B;
%       coef_prdy2 = coef_Iy + A;
%       coef_mrdy2 = coef_Iy - A;
%       coef_inv_mrdx2 = inv(coef_mrdx2);
%   
%       for kk = 1 : N
%           u__ = u(:, :, kk);
%           u_ = u__(2 : end - 1, 2 : end - 1);
%           u_star = coef_prdy2 * u_ * coef_inv_mrdx2;
%           coef_Right_Hand = u_star - A * u_;
%           u_t = mldivide(coef_mrdy2, coef_Right_Hand);
%           u(2 : end - 1, 2 : end - 1, kk + 1) = u_t;
%       end
%   end
%   
%%%
%% Code: Value of Slices
%%%
%   
%   function plotSlice_4_4_3( t_slices, dt, u, X_, Y_ )
%       t_slices_id = round(t_slices / dt) + 1;
%       for ii = 1 : length(t_slices)
%           figure; hold on;
%           title(sprintf('t = %.2f', t_slices(ii)));
%           mesh(X_, Y_, u(:, :, t_slices_id(ii)));
%           view(3); xlabel('x'); ylabel('y'); zlabel('u(numerical)');
%       end
%   end
%
%%%
