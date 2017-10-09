function wave_eqn()
    close all;
    global c f_init G_init make_periodic;
    % Initial conditions
    N = 320;
    Tfinal = 1;
    f_init = @zero_wave;
    G_init = @ramp;
    pause_on = false;
    make_periodic = true;
    % ------------------------------------------------
    % Start of code
    % ------------------------------------------------
    if (pause_on)
        fprintf('Hit enter to see each frame.\n');
    end
    % Physical parameters
    a = 0;
    b = 1;
    c = 1;  % advective speed
    % Numerical parameters
    dx = (b-a)/N;
    xe = linspace(a,b,N+1);
    xc = xe(1:end-1) + dx/2;
    method = 'lw';    % 'upwind', 'lw'
    % Time step size
    cfl = 0.9;
    dt_stable = cfl*dx/c;
    Nfinal = round(Tfinal/dt_stable + 1);
    dt = Tfinal/Nfinal;
    % For plotting
    xef = linspace(a,b,200);
    xcf = (xef(2:end) + xef(1:end-1))/2;
    xlim = [a b];
    ylim = [-2 2];
    % Initial conditions
    qn = soln(xc,0);
    % Plot initial conditions
    mq = 2;
    plot(xc,qn(mq,:),'bo');
    hold on;
    ue = soln(xcf,0);
    plot(xcf,ue(mq,:),'r');
    axis([xlim ylim]);
    if (pause_on)
        input('Hit enter to continue : ');
    end
    hold off;
    R = [1 1; -c c];
    Rinv = (1/(2*c))*[c -1; c 1];
    t_slice = [0, 0.02, 0.1, 0.6]; tt = 1;
    % Begin time stepping
    qnp1 = qn;
    for n = 1:Nfinal,
        tn = (n-1)*dt;
        tnp1 = n*dt;
        % periodic boundary conditions
        ql = soln(a-dx/2,tn);
        qr = soln(b+dx/2,tn);
        qn_bc = [ql qn qr];
        % Construct characteristic decomposition at each edge
        for i = 2:N+2
            delta = qn_bc(:,i) - qn_bc(:,i-1);
            alpha = Rinv*delta;
            % Waves and speeds
            wave(:,i,1) = alpha(1)*R(:,1);  % eigenvector 1
            s(i,1) = -c;    % eigenvalue 1
            wave(:,i,2) = alpha(2)*R(:,2);  % eigenvector 2
            s(i,2) = c;     % eigenvalue 2
            % wave speeds
            Amdq(:,i-1) = s(i,1)*wave(:,i,1);   % left going characterisic
            Apdq(:,i-1) = s(i,2)*wave(:,i,2);   % right going characteristic
            F_corr(:,i-1) = (c/2)*(1 - dt*c/dx)*(wave(:,i,1) + wave(:,i,2));
        end
        % Update solution
        upwind = Amdq(:,2:end) + Apdq(:,1:end-1);
        switch method
            case 'upwind'
                qnp1 = qn - (dt/dx)*upwind;
            case 'lw'
                % Add a second order correction term
                qnp1 = qn - (dt/dx)*upwind - ...
                    (dt/dx)*(F_corr(:,2:end)-F_corr(:,1:end-1));
            otherwise
        end
        if (tt <= 4 && tnp1 > t_slice(tt))
            tt = tt + 1;
            figure; plot(xc,qnp1(mq,:),'bo');
            hold on;
            ue = soln(xcf,tnp1);
            plot(xcf,ue(mq,:),'r');
            axis([xlim ylim]);
            title(sprintf('q(%d) : t = %12.4f',mq,tnp1),'fontsize',18);
            hold off;
        end
        qn = qnp1;
    end
    qe = soln(xc,tnp1);
    err = norm(qnp1(mq,:)-qe(mq,:),1)*dx;
    fprintf('%3d %12.4e\n',N,err);
end

function d = get_x(x)
    global make_periodic;
    if (make_periodic)
        d = mod(x,1);
    else
        d = x;
    end
end

function u = square_wave(x)
    d = get_x(x);
    u = 0.4 <= d & d <= 0.6;
end

function u = gauss_wave(x)
    d = get_x(x);
    u = exp(-80*(d-0.5).^2);
end

function u = periodic_wave(x)
    d = get_x(x);
    u = cos(8*pi*d);
end

function u = zero_wave(x)
    d = get_x(x);
    u = zeros(size(d));
end

function G = ramp(x)
    % G(x) = -integral of g(x) from to b;
    a = 0.4;
    b = 0.6;
    d = get_x(x);
    m1 = d < a;
    m2 = a <= d & d <= b;
    m3 = d > b;
    G = zeros(size(x));
    G(m1) = 0;
    G(m2) = 2;
    G(m3) = 0;
end

function G = sin_wave(x)
    % G(x) = -integral of g(x) from to b;
    G = cos(8*pi*x);
end


function q = soln(x,t)
    global c f_init G_init;
    % Exact solution goes here
    u = f_init(x);
    v= G_init(x);
    q = [u;v];
end