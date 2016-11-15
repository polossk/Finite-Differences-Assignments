function [u, X_, Y_] = solver_4_2_1( dx, dy, dt )
    f_ = @(X, Y)(sin(pi * X) .* sin(2 * pi * Y));
    nu = 1;

    X = round(1 / dx);
    Y = round(1 / dy);
    T = round(1 / dt);
    x = 1; y = 1; t = 1;
    Mx = x * X; My = y * Y; N = t * T;

    rx = nu * dt / dx / dx;
    ry = nu * dt / dy / dy;

    x_ = 0 : dx : x;
    y_ = 0 : dy : y;
    [X_, Y_] = meshgrid(x_, y_);
    f = f_(X_, Y_);
    u = zeros([size(X_), N + 1]);
    u(:, :, 1) = f;

    for kk = 1 : N
        for ii = 2 : My
            for jj = 2 : Mx
                u(ii, jj, kk + 1) = (1 - 2 * (rx + ry)) * u(ii, jj, kk) ...
                    + ry * (u(ii + 1, jj, kk) + u(ii - 1, jj, kk)) ...
                    + rx * (u(ii, jj + 1, kk) + u(ii, jj - 1, kk));
            end
        end
    end
end
