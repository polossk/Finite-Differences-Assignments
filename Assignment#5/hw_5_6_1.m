%% HW 5.6.1
%% Initial-boundary-value Problem
% * $v_t + v_x = 0, x \in (0, 1), t > 0$
% * $v(x, 0) = \sin^{80} \pi x, x \in [0, 1]$
% * $v(0, t) = v(1, t), t \geq 0$
% * One Analytic Solution: $v = \sin^{80} (\pi (x - t))$
%%%
%% Numerical Solution
%%%
ts = [0.0, 0.12, 0.2, 0.4];
% #1: $M = 20, \Delta t = 0.01$
dx = 0.05; dt = 0.01;
[u, v, X, T] = solver_5_6_1(dx, dt);
plotSlice_5_6_1(ts, dt, u, v, X, T);
% #2: $M = 100, \Delta t = 0.002$
dx = 0.01; dt = 0.002;
[u, v, X, T] = solver_5_6_1(dx, dt);
plotSlice_5_6_1(ts, dt, u, v, X, T);

%% Code: Solver of Functionn
%%%
%
%   function [u, v, X_, T_] = solver_5_6_1( dx, dt )
%       f_ = @(x)(sin(pi * x) .^ 80);
%       v_ = @(X, T)(sin(pi * (X - T)) .^ 80);
%       nu = 1;
%   
%       X = round(1 / dx);
%       T = round(1 / dt);
%       x = 1; t = 1;
%       Mx = x * X; N = t * T;
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
%       A_0 = (1 - rx) * speye(length(x_));
%       A_1 = rx * ones(1, length(x_) - 1);
%       A = A_0 + sparse(diag(A_1, -1)); A(1, end) = rx;
%   
%       for kk = 1 : N
%           u_t = (A * u(kk, :)')';
%           u_t(1) = u_t(end);
%           u(kk + 1, :) = u_t;
%       end
%   end
%
%%%
%% Code: Value of Slices
%%%
%   
%   function plotSlice_5_6_1( t_slices, dt, u, v, X_, T_ )
%       figure; hold on;
%       title('Analytic Solution');
%       surf(X_, T_, v); shading flat; view(3);
%       xlabel('x'); ylabel('t'); zlabel('v(analytic)');
%       figure; hold on;
%       title('Numerical Solution');
%       surf(X_, T_, u); shading flat; view(3);
%       xlabel('x'); ylabel('t'); zlabel('u(numerical)');
%       t_slices_id = round(t_slices / dt) + 1;
%       for ii = 1 : length(t_slices)
%           figure; hold on;
%           title(sprintf('t = %.2f', t_slices(ii)));
%           plot(X_(1, :), u(t_slices_id(ii), :), 'b-');
%           plot(X_(1, :), v(t_slices_id(ii), :), 'r--');
%           xlabel('x'); ylabel('y');
%           legend('numerical', 'analytic');
%       end
%   end
%
%%%
