function plotSlice( t_slice, dt, u, u_exact, x_ )
    t_slice_id = round(t_slice / dt) + 1;
    u_slice_num = u(t_slice_id, :);
    u_slice_num = u_slice_num(:);
    u_slice_ext = u_exact(t_slice_id, :);
    u_slice_ext = u_slice_ext(:);
    figure; hold on;
    title(sprintf('t = %.2f', t_slice));
    plot(x_, u_slice_num, 'r-', x_, u_slice_ext, 'g-');
    legend('numerical', 'analytical');
    xlabel('x'); ylabel('u');
    u_slice_err = max(abs(u_slice_num  - u_slice_ext));
    fprintf('when t = %.2f, max error = %f\n', t_slice, u_slice_err);
end