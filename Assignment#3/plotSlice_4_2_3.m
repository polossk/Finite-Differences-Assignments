function plotSlice_4_2_3(u, X_, Y_, dt, id)
    figure; hold on;
    title(sprintf('t = %.2f', id * dt));
    mesh(X_, Y_, u);
    view(3); xlabel('x'); ylabel('y'); zlabel('u(numerical)');
end
