%% HW 1.5.9
%% Eqution
% * $v_t = \nu v_{xx} + F(x, t), x \in (0, 1), t > 0$
% * $v(x, 0) = f(x), x \in [0, 1]$
% * $v_x(0, t) = a(t), v(1, t) = b(t), t \geq 0$
%% Numerical Solver
%
% <include>solver_1_5_9.m</include>
%

%% Numerical Solution
% * $a(t) = 10 \sin(t)$
% * $b(t) = 4 \sin(6t)$
% * $f(x) = x(1 - x)$
% * $F(x, t) = \sin(2 \pi x) \sin(4 \pi t)$
% * $M = 20$, $dt = 0.01$
a_ = @(t)( 10 * sin(t) );
b_ = @(t)( 4 * sin(6 * t) );
f_ = @(x)( x .* (1 - x) );
F_ = @(x, t)( sin(2 * pi * x) * sin(4 * pi * t) );
solver_1_5_9(F_, a_, b_, f_, 10, 0.05, 2 * pi);