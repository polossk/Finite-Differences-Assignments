%% HW 4.2.6
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
% dx = 0.01; dy = 0.01; dt = 0.00002;
dx = 0.01; dy = 0.01; dt = 0.00002;
[u, X, Y] = solver_4_2_5(dx, dy, dt, ts);
