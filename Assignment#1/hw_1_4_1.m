%% HW 1.4.1
%% Neumann-boundary-value Problem
% * $v_t = \nu v_{xx}, x \in (0, 1), t > 0$
% * $v(x, 0) = \cos(\pi x / 2), x \in [0, 1]$
% * $v(1, t) = 0, t \geq 0$
% * $v_x(0, t) = 0,  t \geq 0$
%% Numerical Solution
%%%
% Using $\nu = 1$
f_ = @(x)(cos(pi * x / 2));
a = 0; b = 0;
nu = 1;
%%%
% In domain $[0, 1] \times [0, 1]$, let $N = 250$, $M = 10$,
% a.k.a $\Delta t = 0.004$, $\Delta x = 0.1$
dx = 0.1; dt = 0.004;
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
u(:, M + 1) = b;
%%%
% It is acknowledged that using central difference method to solve this
% eqution, which indicating the solution like
%%%
% $$u^{n+1}_k = u^n_k + r(u^n_{k + 1} - 2 u^n_k + u^n_{k - 1})$$
%%%
% where notation $u^n_k$ indicating the numerical solution of $v$ at the
% point $(k\Delta x, n\Delta t)$
for jj = 1 : N
    u(jj + 1, 1) = (1 - 2 * r) * u(jj, 1) + 2 * r * u(jj, 2);
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
% $$v = \alpha \cos (\omega x + \phi) \exp \{ -(\beta t + \theta) \}$$
%%%
% By using the method of undetermined coefficients, the solution is as
% follows
%%%
% $$ v = \cos (\frac{\pi}{2} x) \exp \{ -\frac{\pi^2}{4} \nu t \}$$
u_exact = cos(pi * X_ / 2) .* exp( (-pi * pi * nu / 4) * T_);
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
%       global u u_exact x_
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