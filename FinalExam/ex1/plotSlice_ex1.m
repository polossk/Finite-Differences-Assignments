function plotSlice_ex1( t_slices, dt, u, X_, Y_ )
    t_slices_id = round(t_slices / dt) + 1;
    for ii = 1 : length(t_slices)
        figure; hold on;
        v = u(:, :, :, t_slices_id(ii));
        subplot(1, 3, 1); mesh(X_, Y_, v(:, :, 1));
        view(3); xlabel('x'); ylabel('y'); zlabel('v_1(numerical)');
        subplot(1, 3, 2); mesh(X_, Y_, v(:, :, 2));
        view(3); xlabel('x'); ylabel('y'); zlabel('v_2(numerical)');
        subplot(1, 3, 3); mesh(X_, Y_, v(:, :, 3));
        view(3); xlabel('x'); ylabel('y'); zlabel('v_3(numerical)');
        suptitle(sprintf('t = %.2f', t_slices(ii)));
    end
end