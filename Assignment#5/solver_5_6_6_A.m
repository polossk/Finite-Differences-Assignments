function [u, v, X_, T_] = solver_5_6_6_A( dx, dt, t )
    f_ = @(x) ((x >= 0) & (x <= 1)) * 1;
    g_ = @(x) abs(rem(x + 1, 3)) - 1;
    v_ = @(X, T) f_(g_(X - T)) * 1;
    nu = 1;

    X = round(1 / dx);
    T = round(1 / dt);
    x = 1; Mx = x * X; N = t * T;

    rx = nu * dt / dx;

    x_ = -1 : dx : 2;
    t_ = 0 : dt : t;
    [X_, T_] = meshgrid(x_, t_);
    u = zeros(size(X_));
    v = v_(X_, T_);
    u(1, :) = f_(x_);
    u_t = zeros(size(x_'));

    D_d = (1 + rx) / 2 * ones(1, length(x_) - 1);
    D_u = (1 - rx) / 2 * ones(1, length(x_) - 1);
    A = sparse(diag(D_d, -1)) + sparse(diag(D_u, 1));
    A(1, end) = (1 + rx) / 2; A(end, 1) = (1 - rx) / 2;

    for kk = 1 : N
        u_t = (A * u(kk, :)')';
        u_t(1) = u_t(end);
        u(kk + 1, :) = u_t;
    end
end
