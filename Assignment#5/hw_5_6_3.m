%% HW 5.6.3
%% Initial-boundary-value Problem
% * $v_t - v_x = 0, x \in (0, 1), t > 0$
% * $v(x, 0) = \sin^{40} \pi x, x \in [0, 1]$
% * $v(0, t) = v(1, t), t \geq 0$
% * One Analytic Solution: $v = \sin^{40} (\pi (x + t))$
%%%
%% Numerical Solution
%%%
ts = [0.0, 0.12, 0.2, 0.8];
% #1: $M = 20, \Delta t = 0.04$
dx = 0.05; dt = 0.04;
[u, v, X, T] = solver_5_6_3(dx, dt, 1);
plotSlice_5_6_1(ts, dt, u, v, X, T);
% #2: $M = 100, \Delta t = 0.008$
dx = 0.01; dt = 0.008;
[u, v, X, T] = solver_5_6_3(dx, dt, 1);
plotSlice_5_6_1(ts, dt, u, v, X, T);
% #3: $t_slice = [0, 5, 10, 20]$
ts = [0, 5, 10, 20];
[u, v, X, T] = solver_5_6_3(dx, dt, 20);
plotSlice_5_6_1(ts, dt, u, v, X, T);

%% Code: Solver of Functionn
%%%
%
%   function [u, v, X_, T_] = solver_5_6_3( dx, dt, t )
%       f_ = @(x)(sin(pi * x) .^ 80);
%       v_ = @(X, T)(sin(pi * (X + T)) .^ 40);
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