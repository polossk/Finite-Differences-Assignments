function plotSlice_5_6_1( t_slices, dt, u, v, X_, T_ )
    figure; hold on;
    title('Analytic Solution');
    surf(X_, T_, v); shading flat; view(3);
    xlabel('x'); ylabel('t'); zlabel('v(analytic)');
    figure; hold on;
    title('Numerical Solution');
    surf(X_, T_, u); shading flat; view(3);
    xlabel('x'); ylabel('t'); zlabel('u(numerical)');
    t_slices_id = round(t_slices / dt) + 1;
    for ii = 1 : length(t_slices)
        figure; hold on;
        title(sprintf('t = %.2f', t_slices(ii)));
        plot(X_(1, :), u(t_slices_id(ii), :), 'b-');
        plot(X_(1, :), v(t_slices_id(ii), :), 'r--');
        xlabel('x'); ylabel('y');
        legend('numerical', 'analytic');
    end
end
