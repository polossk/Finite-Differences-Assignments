%% HW 1.5.5
%% Eqution
% * $v_t = \nu v_{xx}, x \in (0, 1), t > 0$
% * $v(x, 0) = f(x), x \in [0, 1]$
% * $v(0, t) = a(t), v(1, t) = b(t), t \geq 0$
%% Numerical Solver
%
% <include>solver_1_5_4.m</include>
%

%% Numerical Solution
% * $a(t) = \sin(4 \pi t)$
% * $b(t) = 0$
% * $f(x) = \sin(2 \pi x)$
% * $M = 10$, $dt = 0.05$
a_ = @(t)( sin(4 * pi * t) );
b_ = @(t)( 0 );
f_ = @(x)( sin(2 * pi * x) );
solver_1_5_4(a_, b_, f_, 10, 0.05);
