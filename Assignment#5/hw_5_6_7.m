%% HW 5.6.7
%% Initial-boundary-value Problem
% * $v_t - v_x = 0, x \in (0, 1), t > 0$
% * $v(x, 0) = 1, x \in [0.4, 0.6]$
% * $v(x, 0) = 0, x \in [0, 0.4) \cup (0.6, 1]$
% * $v(0, t) = v(1, t), t \geq 0$
% * One Analytic Solution: $v = v_0(|\{ x + t \}|), v_0 = v(x, 0)$
%%%
ts = [0.0, 0.5, 1.0, 2.0];
dx = 0.01; dt = 0.008;
%% Numerical Solution #1: FTFS Difference Scheme
%%%
[u, v, X, T] = solver_5_6_7_A(dx, dt, 2);
plotSlice_5_6_1(ts, dt, u, v, X, T);
%% Numerical Solution #2: Lax-Wendroff Difference Scheme
%%%
[u, v, X, T] = solver_5_6_7_B(dx, dt, 2);
plotSlice_5_6_1(ts, dt, u, v, X, T);
%% Numerical Solution #3: Lax-Friedrichs Difference Scheme
%%%
[u, v, X, T] = solver_5_6_7_C(dx, dt, 2);
plotSlice_5_6_1(ts, dt, u, v, X, T);
%%%
ts = [0.0, 0.5, 1.0, 2.0];
dx = 0.001; dt = 0.0008;
%% Numerical Solution #1: FTFS Difference Scheme
%%%
[u, v, X, T] = solver_5_6_7_A(dx, dt, 2);
plotSlice_5_6_1(ts, dt, u, v, X, T);
%% Numerical Solution #2: Lax-Wendroff Difference Scheme
%%%
[u, v, X, T] = solver_5_6_7_B(dx, dt, 2);
plotSlice_5_6_1(ts, dt, u, v, X, T);
%% Numerical Solution #3: Lax-Friedrichs Difference Scheme
%%%
[u, v, X, T] = solver_5_6_7_C(dx, dt, 2);
plotSlice_5_6_1(ts, dt, u, v, X, T);

%% Code: FTFS Difference Scheme
%%%
%
%   function [u, v, X_, T_] = solver_5_6_7_A( dx, dt, t )
%       f_ = @(x) ((x >= 0.4) & (x <= 0.6));
%       g_ = @(x) abs(rem(x, 1));
%       v_ = @(X, T) f_(g_(X + T)) * 1;
%       nu = -1;
%   
%       X = round(1 / dx);
%       T = round(1 / dt);
%       x = 1; Mx = x * X; N = t * T;
%   
%       rx = nu * dt / dx;
%   
%       x_ = 0 : dx : x;
%       t_ = 0 : dt : t;
%       [X_, T_] = meshgrid(x_, t_);
%       u = zeros(size(X_));
%       v = v_(X_, T_);
%       u(1, :) = f_(x_);
%       u_t = zeros(size(x_'));
%   
%       A_0 = (1 + rx) * speye(length(x_));
%       A_1 = -rx * ones(1, length(x_) - 1);
%       A = A_0 + sparse(diag(A_1, 1)); A(end, 1) = -rx;
%   
%       for kk = 1 : N
%           u_t = (A * u(kk, :)')';
%           u_t(end) = u_t(1);
%           u(kk + 1, :) = u_t;
%       end
%   end
%
%%%
%% Code: Lax-Wendroff Difference Scheme
%%%
%
%   function [u, v, X_, T_] = solver_5_6_7_B( dx, dt, t )
%       f_ = @(x) ((x >= 0.4) & (x <= 0.6));
%       g_ = @(x) abs(rem(x, 1));
%       v_ = @(X, T) f_(g_(X + T)) * 1;
%       nu = -1;
%   
%       X = round(1 / dx);
%       T = round(1 / dt);
%       x = 1; Mx = x * X; N = t * T;
%   
%       rx = nu * dt / dx;
%   
%       x_ = 0 : dx : 1;
%       t_ = 0 : dt : t;
%       [X_, T_] = meshgrid(x_, t_);
%       u = zeros(size(X_));
%       v = v_(X_, T_);
%       u(1, :) = f_(x_);
%       u_t = zeros(size(x_'));
%   
%       D_0 = (1 - rx * rx) * speye(length(x_));
%       D_1 = ones(1, length(x_) - 1);
%       A = rx * (rx + 1) * sparse(diag(D_1, -1)) + rx * (rx - 1) * sparse(diag(D_1, 1));
%       A = 0.5 * A + D_0;
%       A(1, end) = rx * (rx + 1) / 2; A(end, 1) = rx * (rx - 1) / 2;
%   
%       for kk = 1 : N
%           u_t = (A * u(kk, :)')';
%           u_t(end) = u_t(1);
%           u(kk + 1, :) = u_t;
%       end
%   end
%
%%%
%% Code: Lax-Friedrichs Difference Scheme
%%%
%
%   function [u, v, X_, T_] = solver_5_6_7_C( dx, dt, t )
%       f_ = @(x) ((x >= 0.4) & (x <= 0.6));
%       g_ = @(x) abs(rem(x, 1));
%       v_ = @(X, T) f_(g_(X + T)) * 1;
%       nu = -1;
%   
%       X = round(1 / dx);
%       T = round(1 / dt);
%       x = 1; Mx = x * X; N = t * T;
%   
%       rx = nu * dt / dx;
%   
%       x_ = 0 : dx : x;
%       t_ = 0 : dt : t;
%       [X_, T_] = meshgrid(x_, t_);
%       u = zeros(size(X_));
%       v = v_(X_, T_);
%       u(1, :) = f_(x_);
%       u_t = zeros(size(x_'));
%   
%       D_d = (1 + rx) / 2 * ones(1, length(x_) - 1);
%       D_u = (1 - rx) / 2 * ones(1, length(x_) - 1);
%       A = sparse(diag(D_d, -1)) + sparse(diag(D_u, 1));
%       A(1, end) = (1 + rx) / 2; A(end, 1) = (1 - rx) / 2;
%   
%       for kk = 1 : N
%           u_t = (A * u(kk, :)')';
%           u_t(end) = u_t(1);
%           u(kk + 1, :) = u_t;
%       end
%   end
%
%%%