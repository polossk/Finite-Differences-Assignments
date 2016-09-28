function [ u, u_exact, X_, T_ ] = solver_1_4_1( M, dt )
    f_ = @(x)(cos(pi * x / 2));
    a = 0; b = 0;
    nu = 1;
    x = 1; t = 1;
    dx = 1 / M;
    T = round(1 / dt);
    X = round(1 / dx);
    N = t * T; t_real = t * dt * T;
    r = nu * dt / dx / dx;
    x_ = 0 : dx : x;
    t_ = 0 : dt : t_real;
    [X_, T_] = meshgrid(x_, t_);
    u = zeros(size(X_));
    u(1, :) = f_(x_);
    u(:, M + 1) = b;
    for jj = 1 : N
        u(jj + 1, 1) = (1 - 2 * r) * u(jj, 1) + 2 * r * u(jj, 2);
        for kk = 2 : M
            u(jj + 1, kk) = (1 - 2 * r) * u(jj, kk) ...
                + r * (u(jj, kk - 1) + u(jj, kk + 1));
        end
    end
    u_exact = cos(pi * X_ / 2) .* exp( (-pi * pi * nu / 4) * T_);
    figure; mesh(X_, T_, u);
    xlabel('x'); ylabel('t'); zlabel('u(numerical)');
    figure; mesh(X_, T_, u_exact);
    xlabel('x'); ylabel('t'); zlabel('u(exact)');
    delta = u_exact - u;
    err = max(abs(delta(:)));
    figure; mesh(X_, T_, u_exact - u);
    xlabel('x'); ylabel('t'); zlabel('\Deltau');
    fprintf('max error = %f\n', err);
end