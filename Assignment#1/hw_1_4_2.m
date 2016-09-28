%% HW 1.4.2
%% Neumann-boundary-value Problem
% * $v_t = \nu v_{xx}, x \in (0, 1), t > 0$
% * $v(x, 0) = \cos(\pi x / 2), x \in [0, 1]$
% * $v(1, t) = 0, t \geq 0$
% * $v_x(0, t) = 0,  t \geq 0$
%% Numerical Solver
%
% <include>solver_1_4_1.m</include>
%

%% Numerical Solution of M = 20 and dt = 0.001
solver_1_4_1( 20, 0.001 );

%% Numerical Solution of M = 40 and dt = 0.0003
solver_1_4_1( 40, 0.0003 );