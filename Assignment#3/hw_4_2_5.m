%% HW 4.2.5
%% Initial-boundary-value Problem
% * $v_t = \nu (v_{xx} + v_{yy}) + F(x, y, t), (x, y) \in R, t > 0$
% * $v(x, y, t) = g(x, y, t), (x, y) \in \partial R, t > 0$
% * $v(x, y, 0) = f(x, y), (x, y) \in \bar{R}$
% where
% * $R = (0, 1) \times (0, 1)$
% * $F(x, y) = 0, f = 0, \nu = 1$
% * $g(0, y, t) = \sin t \sin \pi y, g(x, 0, t) = \sin t \sin \pi x$
% * $g(1, y, t) = g(x, 1, t) = 0$
%%%
%% Numerical Solution
%%%
% ts = [1.5, 4.5, 6.25];
ts = [1.5, 4.5, 6.25];
% dx = 0.05; dy = 0.05; dt = 0.0005;
dx = 0.05; dy = 0.05; dt = 0.0005;
[u, X, Y] = solver_4_2_5(dx, dy, dt, ts);
%% Code: Solver of Functionn
%%%
%
%   function [u, X_, Y_] = solver_4_2_5( dx, dy, dt, slices )
%       f_ = @(X, Y)( zeros(size(X)) );
%       F_ = @(X, Y, t)( zeros(size(X)) );
%       g_0t = @(W, t)( sin(t) .* sin(pi * W) );
%       g_1t = @(W, t)( zeros(size(W)) );
%       nu = 1;
%       slices_id = round(slices / dt);
%   
%       X = round(1 / dx);
%       Y = round(1 / dy);
%       x = 1; y = 1;
%       Mx = x * X; My = y * Y;
%       
%       rx = nu * dt / dx / dx;
%       ry = nu * dt / dy / dy;
%   
%       x_ = 0 : dx : x;
%       y_ = 0 : dy : y;
%       [X_, Y_] = meshgrid(x_, y_);
%       
%       u = zeros(size(X_));
%       u_t = zeros(size(X_));
%   
%       A_0 = ry * ones(1, length(y_));
%       A_1 = ry * ones(1, length(y_) - 1);
%       A = -2 * diag(A_0) + diag(A_1, 1) + diag(A_1, -1);
%   
%       B_0 = rx * ones(1, length(x_));
%       B_1 = rx * ones(1, length(x_) - 1);
%       B = -2 * diag(B_0) + diag(B_1, 1) + diag(B_1, -1);
%   
%       for id = 1 : max(slices_id)
%           t = id * dt;
%           u(:, 1) = reshape(g_0t(y_, t), size(u(:, 1)));
%           u(1, :) = reshape(g_0t(x_, t), size(u(1, :)));
%           u(end, :) = zeros(size(u(end, :)));
%           u(:, end) = zeros(size(u(:, end)));
%           u_t = u + A * u + u * B;
%           if (any(id == slices_id))
%               plotSlice_4_2_3(u_t, X_, Y_, dt, id);
%           end
%           u = u_t;
%       end
%   end
%
%%%
%% Code: Value of Slices
%%%
%   
%   function plotSlice_4_2_3(u, X_, Y_, dt, id)
%       figure; hold on;
%       title(sprintf('t = %.2f', id * dt));
%       mesh(X_, Y_, u);
%       view(3); xlabel('x'); ylabel('y'); zlabel('u(numerical)');
%   end
%
%%%
