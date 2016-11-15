%% HW 4.2.8
%% Neumann-boundary-value Problem
% * $v_t = \nu (v_{xx} + v_{yy}), (x, y) \in (0, 1) \times (0, 1), t > 0$
% * $v(x, 1, t) = v(x, 0, t) = 0, x \in [0, 1], t > 0$
% * $v(0, y, t) = 0, y \in [0, 1], t > 0$
% * $v_x(1, y, t) = - \pi e ^{-5 \nu \pi^2 t}, y \in [0, 1]. t \geq 0$
% * $v(x, y, 0) = \sin \pi x \sin 2 \pi y, (x, y) \in [0, 1] \times [0, 1]$
% where
% * \nu = 1$
%%%
%% Numerical Solution
%%%
% ts = [0.06, 0.1, 0.9];
ts = [0.06, 0.1, 0.9];
% 1
% dx = 0.05; dy = 0.05; dt = 0.0005;
dx = 0.05; dy = 0.05; dt = 0.0005;
[u, X, Y] = solver_4_2_7(dx, dy, dt, ts);
% 2
% dx = 0.025; dy = 0.025; dt = 0.0001;
dx = 0.025; dy = 0.025; dt = 0.0001;
[u, X, Y] = solver_4_2_7(dx, dy, dt, ts);