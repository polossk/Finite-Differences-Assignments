 %% Demonstrate hyperbolic solvers

function hyper_demo()
    close all;
    global s;
    a = 0;
    b = 1;
    N = 80;
    x = linspace(a,b,N+1);
    xf = linspace(a,b,201);
    dx = (b-a)/N;
    s = 1;  % wave speed
    init_cond = @periodic_wave;
    cfl = 0.9;
    dt_stable = cfl*dx/s;
    Tfinal = 2;
    Nfinal = round(Tfinal/dt_stable) + 1;
    dt = Tfinal/Nfinal;
    method = 'upwind';    % 'ftcs','upwind','lf','lw'
    t0 = 0;
    un = init_cond(x,t0);
    plot(x,un,'bo');
    hold on;
    ue = init_cond(xf,t0);
    plot(xf,ue,'r');
    pause;
    hold off;
    unp1 = un;
    for n = 1:Nfinal
        tn = (n-1)*dt;
        tnp1 = (n)*dt;
        un_ext = [un(end-1) un un(2)];
        switch method
            case 'ftcs'
                unp1 = un_ext(2:end-1) - (s*dt/(2*dx))*(un_ext(3:end) - un_ext(1:end-2));
            case 'upwind'
                unp1 = un_ext(2:end-1) - (s*dt/dx)*(un_ext(2:end-1)-un_ext(1:end-2));
            case 'lf'
                unavg = (un_ext(3:end) + un_ext(1:end-2))/2;
                unp1 = unavg-(s*dt/(2*dx))*(un_ext(3:end) - un_ext(1:end-2));
            case 'lw'
                unp1 = un_ext(2:end-1) - (s*dt/(2*dx))*(un_ext(3:end)-un_ext(1:end-2)) ...
                    + (s^2*dt^2/(2*dx^2))*(un_ext(3:end) - 2*un_ext(2:end-1) + un_ext(1:end-2));
            otherwise
        end
        plot(x,unp1,'bo');
        hold on;
        plot(xf,init_cond(xf,tnp1),'r');
        % plot(x,0*x,'k-');
        title(sprintf('t = %12.4f',tnp1),'fontsize',18);
        axis([0 1 -1.5 1.5]);
        pause;
        hold off;
        un = unp1;
    end

    ue = init_cond(x,tnp1);
    err = norm(unp1-ue,1)*dx;
    fprintf('%3d %12.4e\n',N,err);
end

function u = square_wave(x,t)
    global s;
    d = mod(x - s*t,1);
    u = 0.3 <= d & d <= 0.7;
end

function u = gauss_wave(x,t)
    global s;
    d = mod(x-s*t,1);
    u = exp(-20*(d-0.5).^2);
end

function u = periodic_wave(x,t)
    global s;
    d = mod(x-s*t,1);
    % u = sin(4*pi*d);
    u = cos(8 * pi * d);
end