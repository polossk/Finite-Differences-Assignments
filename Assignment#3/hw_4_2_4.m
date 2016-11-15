%% HW 4.2.4
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
% dx = 0.01; dy = 0.01; dt = 0.00002;
dx = 0.01; dy = 0.01; dt = 0.00002;
[u, X, Y] = solver_4_2_3(dx, dy, dt, ts);
