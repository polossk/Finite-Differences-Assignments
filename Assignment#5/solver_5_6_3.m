function [u, v, X_, T_] = solver_5_6_3( dx, dt, t )
    f_ = @(x)(sin(pi * x) .^ 40);
    v_ = @(X, T)(sin(pi * (X + T)) .^ 40);
    nu = -1;

    X = round(1 / dx);
    T = round(1 / dt);
    x = 1; Mx = x * X; N = t * T;

    rx = nu * dt / dx;

    x_ = 0 : dx : x;
    t_ = 0 : dt : t;
    [X_, T_] = meshgrid(x_, t_);
    u = zeros(size(X_));
    v = v_(X_, T_);
    u(1, :) = f_(x_);
    u_t = zeros(size(x_'));

    A_0 = (1 + rx) * speye(length(x_));
    A_1 = -rx * ones(1, length(x_) - 1);
    A = A_0 + sparse(diag(A_1, 1)); A(end, 1) = -rx;

    for kk = 1 : N
        u_t = (A * u(kk, :)')';
        u_t(end) = u_t(1);
        u(kk + 1, :) = u_t;
    end
end
