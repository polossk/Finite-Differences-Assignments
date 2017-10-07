%% Numerical Solution: Peaceman-Rachford scheme
ts = [0.06, 0.1, 0.9];
dx = 0.05; dy = 0.05; dt = 0.01;
[u1, X1, Y1] = solver_PR(dx, dy, dt);
plotSlice_ex2(ts, dt, u1, X1, Y1);
%% Numerical Solution: FTCS scheme(Mx = My = 20, dt = 0.01)
% *Warning: This config won't converge.*
[u2, X2, Y2] = solver_FTCS(dx, dy, dt);
plotSlice_ex2(ts, dt, u2, X2, Y2);
%% Numerical Solution: FTCS scheme(Mx = My = 20, dt = 0.0005)
% _Notice: This config will converge._
[u2, X2, Y2] = solver_FTCS(dx, dy, 0.0005);
plotSlice_ex2(ts, dt, u2, X2, Y2);
%% Code: Solver of Peaceman-Rachford scheme
%%%
%
%   function [u, X_, Y_] = solver_PR( dx, dy, dt )
%       f_ = @(X, Y)(sin(pi * X) .* sin(2 * pi * Y));
%       nu = 1;
%       X = round(1 / dx);
%       Y = round(1 / dy);
%       T = round(1 / dt);
%       x = 1; y = 1; t = 1;
%       Mx = x * X; My = y * Y; N = t * T;
%       hrx = nu * dt / dx / dx / 2;
%       hry = nu * dt / dy / dy / 2;
%       x_ = 0 : dx : x;
%       y_ = 0 : dy : y;
%       [X_, Y_] = meshgrid(x_, y_);
%       u = zeros([size(X_), N + 1]);
%       u(:, :, 1) = f_(X_, Y_);
%       u_t = zeros(size(u(2 : end - 1, 2 : end - 1, 1)));
%       u_h = zeros(size(u_t));
%       x__ = x_(2 : end - 1); y__ = y_(2 : end - 1);
%       [X__, Y__] = meshgrid(x__, y__);
%       A_0 = hry * speye(length(y__));
%       A_1 = hry * ones(1, length(y__) - 1);
%       A = -2 * A_0 + sparse(diag(A_1, 1)) + sparse(diag(A_1, -1));
%       B_0 = hrx * speye(length(x__));
%       B_1 = hrx * ones(1, length(x__) - 1);
%       B = -2 * B_0 + sparse(diag(B_1, 1)) + sparse(diag(B_1, -1));
%       coef_Ix = speye(length(x__));
%       coef_Iy = speye(length(y__));
%       coef_prdx2 = coef_Ix + B;
%       coef_mrdx2 = coef_Ix - B;
%       coef_prdy2 = coef_Iy + A;
%       coef_mrdy2 = coef_Iy - A;
%       coef_inv_mrdx2 = inv(coef_mrdx2);
%       coef_inv_mrdx2_times_prdx2 = coef_inv_mrdx2 * coef_prdx2;
%       for kk = 1 : N
%           u__ = u(:, :, kk);
%           u_ = u__(2 : end - 1, 2 : end - 1);
%           coef_Right_Hand = coef_prdy2 * u_ * coef_inv_mrdx2_times_prdx2;
%           u_t = mldivide(coef_mrdy2, coef_Right_Hand);
%           u(2 : end - 1, 2 : end - 1, kk + 1) = u_t;
%       end
%   end
%
%%%
%% Code: Solver of FTCS scheme
%%%
%
%   function [u, X_, Y_] = solver_FTCS( dx, dy, dt )
%       f_ = @(X, Y)(sin(pi * X) .* sin(2 * pi * Y));
%       nu = 1;
%       X = round(1 / dx);
%       Y = round(1 / dy);
%       T = round(1 / dt);
%       x = 1; y = 1; t = 1;
%       Mx = x * X; My = y * Y; N = t * T;
%       rx = nu * dt / dx / dx;
%       ry = nu * dt / dy / dy;
%       x_ = 0 : dx : x;
%       y_ = 0 : dy : y;
%       [X_, Y_] = meshgrid(x_, y_);
%       f = f_(X_, Y_);
%       u = zeros([size(X_), N + 1]);
%       u(:, :, 1) = f;
%       for kk = 1 : N
%           for ii = 2 : My
%               for jj = 2 : Mx
%                   u(ii, jj, kk + 1) = (1 - 2 * (rx + ry)) * u(ii, jj, kk) ...
%                       + ry * (u(ii + 1, jj, kk) + u(ii - 1, jj, kk)) ...
%                       + rx * (u(ii, jj + 1, kk) + u(ii, jj - 1, kk));
%               end
%           end
%       end
%   end
%
%%%
%% Code: Value of Slices
%%%
%   
%     function plotSlice_ex2( t_slices, dt, u, X_, Y_ )
%         t_slices_id = round(t_slices / dt) + 1;
%         for t_slice_id = t_slices_id
%             figure; hold on;
%             title(sprintf('t = %.2f', t_slice));
%             mesh(X_, Y_, u(:, :, t_slices_id(ii)));
%             view(3); xlabel('x'); ylabel('y'); zlabel('u(numerical)');
%         end
%     end
%
%%%
