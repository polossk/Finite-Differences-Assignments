function [u, X_, Y_] = solver_4_2_3( dx, dy, dt, slices )
    f_ = @(X, Y)( zeros(size(X)) );
    F_ = @(X, Y, t)( sin(t) .* sin(2 * pi * X) .* sin(4 * pi * Y) );
    nu = 0.5;
    slices_id = round(slices / dt);

    X = round(1 / dx);
    Y = round(1 / dy);
    x = 1; y = 1;
    Mx = x * X; My = y * Y;
    
    rx = nu * dt / dx / dx;
    ry = nu * dt / dy / dy;

    x_ = 0 : dx : x;
    y_ = 0 : dy : y;
    [X_, Y_] = meshgrid(x_, y_);
    u_fin = f_(X_, Y_);

    x__ = x_(2 : end - 1); y__ = y_(2 : end - 1);
    [X__, Y__] = meshgrid(x__, y__);
    u = zeros(size(X__));
    u_t = zeros(size(X__));

    A_0 = ry * ones(1, length(y__));
    A_1 = ry * ones(1, length(y__) - 1);
    A = -2 * diag(A_0) + diag(A_1, 1) + diag(A_1, -1);

    B_0 = rx * ones(1, length(x__));
    B_1 = rx * ones(1, length(x__) - 1);
    B = -2 * diag(B_0) + diag(B_1, 1) + diag(B_1, -1);

    for id = 1 : max(slices_id)
        u_t = u + A * u + u * B + F_(X__, Y__, id * dt);
        if (any(id == slices_id))
            u_fin(2 : end - 1, 2 : end - 1) = u_t;
            plotSlice_4_2_3(u_fin, X_, Y_, dt, id);
        end
        u = u_t;
    end
end
