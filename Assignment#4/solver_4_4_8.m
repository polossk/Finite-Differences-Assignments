function [u, X_, Y_] = solver_4_4_8( dx, dy, dt )
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
    u = zeros([size(X_), N + 1]);
    u(:, :, 1) = f_(X_, Y_);

    u_t = zeros(size(u(2 : end - 1, 2 : end - 1, 1)));
    u_h = zeros(size(u_t));

    x__ = x_(2 : end - 1); y__ = y_(2 : end - 1);
    [X__, Y__] = meshgrid(x__, y__);

    A_0 = ry * speye(length(y__));
    A_1 = ry * ones(1, length(y__) - 1);
    A = -2 * A_0 + sparse(diag(A_1, 1)) + sparse(diag(A_1, -1));

    B_0 = rx * speye(length(x__));
    B_1 = rx * ones(1, length(x__) - 1);
    B = -2 * B_0 + sparse(diag(B_1, 1)) + sparse(diag(B_1, -1));

    coef_Ix = speye(length(x__));
    coef_Iy = speye(length(y__));

    coef_mrdx2 = coef_Ix - B;
    coef_prdy2 = coef_Iy + A;
    coef_mrdy2 = coef_Iy - A;
    coef_inv_mrdx2 = inv(coef_mrdx2);

    for kk = 1 : N
        u__ = u(:, :, kk);
        u_ = u__(2 : end - 1, 2 : end - 1);
        u_star = coef_prdy2 * u_ * coef_inv_mrdx2;
        coef_Right_Hand = u_star - A * u_;
        u_t = mldivide(coef_mrdy2, coef_Right_Hand);
        u(2 : end - 1, 2 : end - 1, kk + 1) = u_t;
    end
end
