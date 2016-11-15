%% HW 4.2.1
%% Initial-boundary-value Problem
% * $v_t = \nu (v_{xx} + v_{yy}), (x, y) \in (0, 1) \times (0, 1), t > 0$
% * $v(0, y, t) = v(1, y, t) = 0, y \in [0, 1], t \geq 0$
% * $v(x, 0, t) = v(x, 1, t) = 0, x \in [0, 1], t \geq 0$
% * $v(x, y, 0) = \sin \pi x \sin 2 \pi y, (x, y) \in [0, 1] \times [0, 1]$
%%%
%% Numerical Solution
%%%
% ts = [0.06, 0.1, 0.9];
ts = [0.06, 0.1, 0.9];
% 1
% dx = 0.05; dy = 0.05; dt = 0.0005;
dx = 0.05; dy = 0.05; dt = 0.0005;
[u, X, Y] = solver_4_2_1(dx, dy, dt);
plotSlice_4_2_1(ts, dt, u, X, Y);
% 2
% dx = 0.05; dy = 0.05; dt = 0.0001;
dx = 0.05; dy = 0.05; dt = 0.0001;
[u, X, Y] = solver_4_2_1(dx, dy, dt);
plotSlice_4_2_1(ts, dt, u, X, Y);
%% Code: Solver of Functionn
%%%
%
%   function [u, X_, Y_] = solver_4_2_1( dx, dy, dt )
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
%       f = f_(X_, Y_);
%       u = zeros([size(X_), N + 1]);
%       u(:, :, 1) = f;
%   
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
%     function plotSlice_4_2_1( t_slices, dt, u, X_, Y_ )
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
