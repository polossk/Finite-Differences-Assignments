function [ u, X_, T_ ] = solver_1_5_9( Fxt, a, b, f, M, dt, varargin )
    if (nargin > 6)
        t = varargin{1};
    else
        t = 2;
    end
    nu = 0.1; x = 1;
    dx = 1 / M;
    T = round(1 / dt);
    X = round(1 / dx);
    N = t * T; t_real = t * dt * T;
    r = nu * dt / dx / dx;
    s = nu * dt / dx;
    x_ = 0 : dx : x;
    t_ = 0 : dt : t_real;
    [X_, T_] = meshgrid(x_, t_);
    u = zeros(size(X_));
    u(1, :) = f(x_);
    u(:, M + 1) = b(t_);
    for jj = 1 : N
        u(jj + 1, 1) = (1 - 2 * r) * u(jj, 1) ...
            + 2 * r * u(jj, 2) - 2 * s * a(t_(jj));
        for kk = 2 : M
            u(jj + 1, kk) = (1 - 2 * r) * u(jj, kk) ...
                + r * (u(jj, kk - 1) + u(jj, kk + 1)) ...
                + dt * Fxt(x_(kk), t_(jj));
        end
    end
    figure; mesh(X_, T_, u);
    xlabel('x'); ylabel('t'); zlabel('u(numerical)');
    ylim([0 t_real]);
    plotSlice_1_5_1( 0.1, dt, u, x_ );
    plotSlice_1_5_1( 0.9, dt, u, x_ );
    plotSlice_1_5_1( 2.0, dt, u, x_ );
end