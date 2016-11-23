function plotSlice_4_4_3( t_slices, dt, u, X_, Y_ )
    t_slices_id = round(t_slices / dt) + 1;
    for ii = 1 : length(t_slices)
        figure; hold on;
        title(sprintf('t = %.2f', t_slices(ii)));
        mesh(X_, Y_, u(:, :, t_slices_id(ii)));
        view(3); xlabel('x'); ylabel('y'); zlabel('u(numerical)');
    end
end
