%% HW 4.2.3
%% Initial-boundary-value Problem
% * $v_t = \nu (v_{xx} + v_{yy}) + F(x, y, t), (x, y) \in R, t > 0$
% * $v(x, y, t) = g(x, y, t), (x, y) \in \partial R, t > 0$
% * $v(x, y, 0) = f(x, y), (x, y) \in \bar{R}$
% where
% * $R = (0, 1) \times (0, 1)$
% * $F(x, y) = \sin t \sin 2 \pi x \sin 4 \pi y$
% * $g = 0, f = 0, \nu = 0.5$
%%%
%% Numerical Solution
%%%
% ts = [1.5, 4.5, 6.25];
ts = [1.5, 4.5, 6.25];
% dx = 0.05; dy = 0.05; dt = 0.0005;
dx = 0.05; dy = 0.05; dt = 0.0005;
[u, X, Y] = solver_4_2_3(dx, dy, dt, ts);
%% Code: Solver of Functionn
%%%
%
%   function [u, X_, Y_] = solver_4_2_3( dx, dy, dt, slices )
%       f_ = @(X, Y)( zeros(size(X)) );
%       F_ = @(X, Y, t)( sin(t) .* sin(2 * pi * X) .* sin(4 * pi * Y) );
%       nu = 0.5;
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
%       u_fin = f_(X_, Y_);
%   
%       x__ = x_(2 : end - 1); y__ = y_(2 : end - 1);
%       [X__, Y__] = meshgrid(x__, y__);
%       u = zeros(size(X__));
%       u_t = zeros(size(X__));
%   
%       A_0 = ones(1, length(y__));
%       A_1 = ones(1, length(y__) - 1);
%       A = -2 * diag(A_0) + diag(A_1, 1) + diag(A_1, -1);
%   
%       B_0 = ones(1, length(x__));
%       B_1 = ones(1, length(x__) - 1);
%       B = -2 * diag(B_0) + diag(B_1, 1) + diag(B_1, -1);
%   
%       for id = 1 : max(slices_id)
%           u_t = u + A * u + u * B + F_(X__, Y__, id * dt);
%           if (any(id == slices_id))
%               u_fin(2 : end - 1, 2 : end - 1) = u_t;
%               plotSlice_4_2_3(u_fin, X_, Y_, dt, id);
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
