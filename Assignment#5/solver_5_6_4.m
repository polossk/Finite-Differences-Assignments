function [u, v, X_, T_] = solver_5_6_4( dx, dt, t )
    f_ = @(x)(sin(pi * x) .^ 40);
    v_ = @(X, T)(sin(pi * (X + T)) .^ 40);
    nu = -1;

    X = round(1 / dx);
    T = round(1 / dt);
    x = 1; Mx = x * X; N = t * T;

    qrx = nu * dt / dx / 4;

    x_ = 0 : dx : x;
    t_ = 0 : dt : t;
    [X_, T_] = meshgrid(x_, t_);
    u = zeros(size(X_));
    v = v_(X_, T_);
    u(1, :) = f_(x_);
    u_t = zeros(size(x_'));

    D_0 = speye(length(x_));
    D_1 = ones(1, length(x_) - 1);
    A = D_0 - qrx * sparse(diag(D_1, -1)) + qrx * sparse(diag(D_1, 1));
    A(1, end) = -qrx; A(end, 1) = qrx;
    B = D_0 + qrx * sparse(diag(D_1, -1)) - qrx * sparse(diag(D_1, 1));
    B(1, end) = qrx; B(end, 1) = -qrx;
    C = inv(A) * B;

    % w = zeros(length(x_), 1); w(1) = 1; w(end) = -1;
    % z = abs(w) .* qrx;
    % EE = D_0 - qrx * sparse(diag(D_1, -1)) + qrx * sparse(diag(D_1, 1));
    % EE(1, 1) = EE(1, 1) + qrx;
    % EE(end, end) = EE(end, end) - qrx;

    for kk = 1 : N
        u_t = (C * u(kk, :)')';
        u_t(end) = u_t(1);
        u(kk + 1, :) = u_t;
        % b = B * u(kk, :)';
        % y_1 = EE \ b;
        % y_2 = EE \ w;
        % temp = dot(z, y_1) / (1 - dot(z, y_2));
        % u_t = y_1 + temp * y_2;
    end
end
