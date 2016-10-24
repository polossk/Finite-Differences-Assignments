function u = solver_2_6_4( va, vb, vc, fvd, u_exact )
    f_ = @(x)(cos(pi * x / 2));
    a = 0; b = 0; nu = 1;
    dx =  0.1; dt = 0.004;
    T = round(1 / dt);
    X = round(1 / dx);
    x = 1; t = 1;
    M = x * X; N = t * T;
    r = nu * dt / dx / dx;
    x_ = 0 : dx : 1;
    t_ = 0 : dt : t;
    [X_, T_] = meshgrid(x_, t_);
    f = f_(x_);
    u = zeros(size(X_));
    u(1, :) = f; u(:, end) = b;

    for jj = 1 : N
        vd = fvd(u(jj, 1 : end - 1));
        A = diag(vb) + diag(va, -1) + diag(vc, 1);
        [u(jj + 1, 1 : end - 1), ~] = gmres(A, vd');
        % u(jj + 1, 1 : end - 1) = thomas_algo(va, vb, va, vd)';
    end

    figure; mesh(X_, T_, u);
    xlabel('x'); ylabel('t'); zlabel('u(numerical)');
    figure; mesh(X_, T_, u_exact);
    xlabel('x'); ylabel('t'); zlabel('u(exact)');
    delta = u_exact - u;
    err = max(abs(delta(:)));
    figure; mesh(X_, T_, u_exact - u);
    xlabel('x'); ylabel('t'); zlabel('\Deltau');
    fprintf('max error = %f\n', err);

    plotSlice( 0.06, dt, u, u_exact, x_ );
    plotSlice( 0.1, dt, u, u_exact, x_ );
    plotSlice( 0.9, dt, u, u_exact, x_ );
end