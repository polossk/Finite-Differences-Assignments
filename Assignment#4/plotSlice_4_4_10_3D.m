function plotSlice_4_4_3( t_slices, dt, u, X_, Y_, dz )
    t_slices_id = round(t_slices / dt) + 1;
    z_slices = [0.1, 0.35, 0.65, 0.9];
    idx = round(z_slices / dz) + 1;
    for ii = 1 : length(t_slices)
        figure; hold on;
        suptitle(sprintf('t = %.3f', t_slices(ii)));
        for jj = 1 : 4
            subplot(2, 2, jj); hold on;
            [~, h] = contourf(X_(:, :, idx(jj)), Y_(:, :, idx(jj)), ...
                u(:, :, idx(jj), t_slices_id(ii)));
            h.ShowText = 'on';
            text(-0.019, 0.9, sprintf('z = %.2f', z_slices(jj)));
        end
        hold off;
    end
end
