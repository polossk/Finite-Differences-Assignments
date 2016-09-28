%% HW 1.5.4
%% Eqution
% * $v_t = \nu v_{xx}, x \in (0, 1), t > 0$
% * $v(x, 0) = f(x), x \in [0, 1]$
% * $v(0, t) = a(t), v(1, t) = b(t), t \geq 0$
%% Numerical Solution
%%%
% Using $\nu = 1$
nu = 0.1;
%%%
% In domain $[0, 1] \times [0, 1]$, let $N = 250$, $M = 10$,
% a.k.a $\Delta t = 0.004$, $\Delta x = 0.1$
dx = 0.1; dt = 0.05;
T = round(1 / dt);
X = round(1 / dx);
x = 1; t = 2;
M = x * X; N = t * T;
%%%
% Let $r = \nu \frac{\Delta t}{\Delta x^2}$
r = nu * dt / dx / dx;
%%%
% Setting domain and boundary value
x_ = 0 : dx : 1;
t_ = 0 : dt : t;
[X_, T_] = meshgrid(x_, t_);
u = zeros(size(X_));
u(1, :) = 0;
u(:, 1) = sin(4 * pi * t_); u(:, M + 1) = 0;
for jj = 1 : N
    for kk = 2 : M
        u(jj + 1, kk) = (1 - 2 * r) * u(jj, kk) ...
            + r * (u(jj, kk - 1) + u(jj, kk + 1));
    end
end
figure; mesh(X_, T_, u);
xlabel('x'); ylabel('t'); zlabel('u(numerical)');
%% Value of Slices
%%%
%
% <include>plotSlice_1_5_1.m</include>
%
%%%
% $t = 0.06$
plotSlice_1_5_1( 0.1, dt, u, x_ );
%%%
% $t = 0.1$
plotSlice_1_5_1( 0.1, dt, u, x_ );
%%%
% $t = 0.9$
plotSlice_1_5_1( 2.0, dt, u, x_ );