%% HW 1.2.1
%% Initial-boundary-value Problem
% * $v_t = \nu v_{xx}, x \in (0, 1), t > 0$
% * $v(x, 0) = f(x), x \in [0, 1]$
% * $v(0, t) = a(t), v(1, t) = b(t), t \geq 0$
%%%
% where $f(0) = a(0)$ , and $f(1) = b(0)$
%% Numerical Solution
%%%
% Using $f(x) = \sin(2 \pi x)$, $a = b = 0$, and $\nu = 1 / 6$
f_ = @(x)(sin(2 * pi * x));
a = 0; b = 0;
nu = 1 / 6;
%%%
% In domain $[0, 1] \times [0, 1]$, let $N = 50$, $M = 10$,
% a.k.a $\Delta t = 0.02$, $\Delta x = 0.1$
dx =  0.1; dt = 0.02;
T = round(1 / dt);
X = round(1 / dx);
x = 1; t = 1;
M = x * X; N = t * T;
%%%
% Let $r = \nu \frac{\Delta t}{\Delta x^2}$
r = nu * dt / dx / dx;
%%%
% Setting domain and boundary value
x_ = 0 : dx : 1;
t_ = 0 : dt : t;
[X_, T_] = meshgrid(x_, t_);
f = f_(x_);
u = zeros(size(X_));
u(1, :) = f;
%%%
% It is acknowledged that using central difference method to solve this
% eqution, which indicating the solution like
%%%
% $$u^{n+1}_k = u^n_k + r(u^n_{k + 1} - 2 u^n_k + u^n_{k - 1})$$
%%%
% where notation $u^n_k$ indicating the numerical solution of $v$ at the
% point $(k\Delta x, n\Delta t)$
for jj = 1 : N
    for kk = 2 : M
        u(jj + 1, kk) = (1 - 2 * r) * u(jj, kk) ...
            + r * (u(jj, kk - 1) + u(jj, kk + 1));
    end
end
figure; mesh(X_, T_, u);
xlabel('x'); ylabel('t'); zlabel('u(numerical)');
%% Analytical Solution
%%%
% Notice that the analytical solutions following the form like
%%%
% $$v = \alpha \sin (\omega x + \phi) \exp \{ -(\beta t + \theta) \}$$
%%%
% By using the method of undetermined coefficients, the solution is as
% follows
%%%
% $$ v = \sin (2 \pi x) \exp \{ -4 \pi^2 \nu t \}$$
u_exact = sin(2 * pi * X_) .* exp( (-4 * pi * pi * nu) * T_);
figure; mesh(X_, T_, u_exact);
xlabel('x'); ylabel('t'); zlabel('u(exact)');
%% Error
delta = u_exact - u;
err = max(abs(delta(:)));
figure; mesh(X_, T_, u_exact - u);
xlabel('x'); ylabel('t'); zlabel('\Deltau');
fprintf('max error = %f\n', err);
%% Value of Slices
%%%
% 
%   function plotSlice( t_slice, dt, u, u_exact, x_ )
%       t_slice_id = round(t_slice / dt) + 1;
%       u_slice_num = u(t_slice_id, :);
%       u_slice_num = u_slice_num(:);
%       u_slice_ext = u_exact(t_slice_id, :);
%       u_slice_ext = u_slice_ext(:);
%       figure; hold on;
%       plot(x_, u_slice_num, 'r-', x_, u_slice_ext, 'g-');
%       legend('numerical', 'analytical');
%       xlabel('x'); ylabel('u');
%       u_slice_err = max(abs(u_slice_num  - u_slice_ext));
%       fprintf('when t = %.2f, max error = %f\n', t_slice, u_slice_err);
%   end
%
%%%
% $t = 0.06$
plotSlice( 0.06, dt, u, u_exact, x_ );
%%%
% $t = 0.1$
plotSlice( 0.1, dt, u, u_exact, x_ );
%%%
% $t = 0.9$
plotSlice( 0.9, dt, u, u_exact, x_ );
%% Slice of $t = 50$
% It's known that the formulation is stable if and only if $r \leq 1/2$,
% which indicating $\Delta t$ is defined under the constraint that
% $\Delta t \leq \frac{\Delta x^2}{6 \nu} = 0.01$.
t = 50; dt = 0.01;
N = round(t / dt);
t_ = 0 : dt : t;
r = nu * dt / dx / dx;
[X_, T_] = meshgrid(x_, t_);
u_1 = f; u_2 = f;
u_2(1) = 0; u_2(M + 1) = 0;
for jj = 1 : N
    for kk = 2 : M
        u_1(kk) = (1 - 2 * r) * u_2(kk) + r * (u_2(kk - 1) + u_2(kk + 1));
    end
    u_2 = u_1;
end
u_exact_2 = sin(2 * pi * x_) .* exp( (-4 * pi * pi * nu) * 50);
u_slice_err = max(abs(u_1  - u_exact_2));
fprintf('when t = %.2f, max error = %f\n', 50, u_slice_err);