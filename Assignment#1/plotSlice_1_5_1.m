function plotSlice_1_5_1( t_slice, dt, u, x_ )
    t_slice_id = round(t_slice / dt) + 1;
    u_slice_num = u(t_slice_id, :);
    u_slice_num = u_slice_num(:);
    figure; hold on;
    title(sprintf('t = %.2f', t_slice));
    plot(x_, u_slice_num, 'r-');
    legend('numerical');
    xlabel('x'); ylabel('u');
end