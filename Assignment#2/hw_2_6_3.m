%% HW 2.6.3
%% Eqution with Lower Order Terms
% * $v_t + a v_x = \nu v_{xx}, x \in (0, 1), t > 0$
% * $v(x, 0) = f(x), x \in [0, 1]$
% * $v(0, t) = 0, v(1, t) = 0, t \geq 0$
%% Numerical Solution
%%%
% Using $\nu = 1$
f_ = @(x)(sin(4 * pi * x));
a = 2; nu = 1;
%%%
% In domain $[0, 1] \times [0, 1]$, let $N = 250$, $M = 10$,
% a.k.a $\Delta t = 0.004$, $\Delta x = 0.1$
dx = 0.05; dt = 0.001;
T = round(1 / dt);
X = round(1 / dx);
x = 1; t = 1;
M = x * X; N = t * T;
%%%
% Let $r = \nu \frac{\Delta t}{\Delta x^2}$
r = nu * dt / dx / dx;
%%%
% Let $s = a \frac{\Delta t}{2 \Delta x}$
s = a * dt / dx / 2;
%%%
% Setting domain and boundary value
x_ = 0 : dx : 1;
t_ = 0 : dt : t;
[X_, T_] = meshgrid(x_, t_);
f = f_(x_);
u = zeros(size(X_));
u(1, :) = f;
u(:, 1) = 0; u(:, M + 1) = 0;

va = (-r + s) * ones(1, M - 2);
vc = (-r - s) * ones(1, M - 2);
vb = bsxfun(@plus, ones(1, M - 1), r * 2);

for jj = 1 : N - 1
    u(jj + 1, 2 : end - 1) = thomas_algo(va, vb, vc, u(jj, 2 : end - 1))';
end
figure; mesh(X_, T_, u);
xlabel('x'); ylabel('t'); zlabel('u(numerical)');
%% Value of Slices
%%%
% 
%   function plotSlice_2_6_3( t_slice, dt, u, x_ )
%       global u u_exact x_
%       t_slice_id = round(t_slice / dt) + 1;
%       u_slice_num = u(t_slice_id, :);
%       u_slice_num = u_slice_num(:);
%       figure; hold on;
%       plot(x_, u_slice_num, 'r-');
%       legend('numerical');
%       xlabel('x'); ylabel('u');
%   end
%
%%%
% $t = 0.06$
plotSlice_2_6_3( 0.06, dt, u, x_ );
%%%
% $t = 0.1$
plotSlice_2_6_3( 0.1, dt, u, x_ );
%%%
% $t = 0.9$
plotSlice_2_6_3( 0.9, dt, u, x_ );
u_slice_1 = u(round(0.9 / dt) + 1, :);
u_slice_1 = u_slice_1(:);
dt = 0.1; T = round(1 / dt); N = t * T;
r = nu * dt / dx / dx;
s = a * dt / dx / 2;
t_ = 0 : dt : t;
[X_, T_] = meshgrid(x_, t_);
f = f_(x_);
u = zeros(size(X_));
u(1, :) = f;
u(:, 1) = 0; u(:, M + 1) = 0;

va = (-r + s) * ones(1, M - 2);
vc = (-r - s) * ones(1, M - 2);
vb = bsxfun(@plus, ones(1, M - 1), r * 2);

for jj = 1 : N - 1
    d = (r + s) * u(jj, 1 : M - 1) + ...
        (1 - 2 * r) * u(jj, 2 : M) + ...
        (r - s) * u(jj, 3 : M + 1);
    d(1) = d(1) + (r + s) * 0;
    d(end) = d(end) + (r - s) * 0;
    u(jj + 1, 2 : end - 1) = thomas_algo(va, vb, vc, u(jj, 2 : end - 1))';
end
u_slice_2 = u(round(0.9 / dt) + 1, :);
u_slice_2 = u_slice_1(:);
figure; hold on;
title(sprintf('t = %.2f', 0.9));
plot(x_, u_slice_1, 'r-', x_, u_slice_2, 'b-');
legend('dt = 0.001', 'dt = 0.1');
xlabel('x'); ylabel('u');
fprintf('max error = %e\n', max(abs(u_slice_1  - u_slice_2)));