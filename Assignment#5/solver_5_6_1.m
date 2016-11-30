function [u, v, X_, T_] = solver_5_6_1( dx, dt )
    f_ = @(x)(sin(pi * x) .^ 80);
    v_ = @(X, T)(sin(pi * (X - T)) .^ 80);
    nu = 1;

    X = round(1 / dx);
    T = round(1 / dt);
    x = 1; t = 1;
    Mx = x * X; N = t * T;

    rx = nu * dt / dx;

    x_ = 0 : dx : x;
    t_ = 0 : dt : t;
    [X_, T_] = meshgrid(x_, t_);
    u = zeros(size(X_));
    v = v_(X_, T_);
    u(1, :) = f_(x_);
    u_t = zeros(size(x_'));

    A_0 = (1 - rx) * speye(length(x_));
    A_1 = rx * ones(1, length(x_) - 1);
    A = A_0 + sparse(diag(A_1, -1)); A(1, end) = rx;

    for kk = 1 : N
        u_t = (A * u(kk, :)')';
        u_t(1) = u_t(end);
        u(kk + 1, :) = u_t;
    end
end
