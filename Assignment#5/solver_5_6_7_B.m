function [u, v, X_, T_] = solver_5_6_7_B( dx, dt, t )
    f_ = @(x) ((x >= 0.4) & (x <= 0.6)) * 1;
    g_ = @(x) abs(rem(x, 1));
    v_ = @(X, T) f_(g_(X + T)) * 1;
    nu = -1;

    X = round(1 / dx);
    T = round(1 / dt);
    x = 1; Mx = x * X; N = t * T;

    rx = nu * dt / dx;

    x_ = 0 : dx : 1;
    t_ = 0 : dt : t;
    [X_, T_] = meshgrid(x_, t_);
    u = zeros(size(X_));
    v = v_(X_, T_);
    u(1, :) = f_(x_);
    u_t = zeros(size(x_'));

    D_0 = (1 - rx * rx) * speye(length(x_));
    D_1 = ones(1, length(x_) - 1);
    A = rx * (rx + 1) * sparse(diag(D_1, -1)) + rx * (rx - 1) * sparse(diag(D_1, 1));
    A = 0.5 * A + D_0;
    A(1, end) = rx * (rx + 1) / 2; A(end, 1) = rx * (rx - 1) / 2;

    for kk = 1 : N
        u_t = (A * u(kk, :)')';
        u_t(end) = u_t(1);
        u(kk + 1, :) = u_t;
    end
end
