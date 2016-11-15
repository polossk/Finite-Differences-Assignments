%% HW 4.2.2
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
% dx = 0.1; dy = 0.05; dt = 0.0005;
dx = 0.1; dy = 0.05; dt = 0.0005;
[u, X, Y] = solver_4_2_1(dx, dy, dt);
plotSlice_4_2_1(ts, dt, u, X, Y);
% 2
% dx = 0.1; dy = 0.05; dt = 0.0001;
dx = 0.1; dy = 0.05; dt = 0.0001;
[u, X, Y] = solver_4_2_1(dx, dy, dt);
plotSlice_4_2_1(ts, dt, u, X, Y);