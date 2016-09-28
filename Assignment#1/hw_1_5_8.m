%% HW 1.5.8
%% Eqution
% * $v_t = \nu v_{xx} + F(x, t), x \in (0, 1), t > 0$
% * $v(x, 0) = f(x), x \in [0, 1]$
% * $v(0, t) = a(t), v(1, t) = b(t), t \geq 0$
%% Numerical Solver
%
% <include>solver_1_5_7.m</include>
%

%% Numerical Solution #1
% * $a(t) = 10 \sin(t)$
% * $b(t) = 4 \sin(6t)$
% * $f(x) = x(1 - x)$
% * $F(x, t) = \sin(2 \pi x) \sin(4 \pi t)$
% * $M = 10$, $dt = 0.05$
a_ = @(t)( 10 * sin(t) );
b_ = @(t)( 4 * sin(6 * t) );
f_ = @(x)( x .* (1 - x) );
F_ = @(x, t)( sin(2 * pi * x) * sin(4 * pi * t) );
solver_1_5_7(F_, a_, b_, f_, 10, 0.05, 2 * pi);

%% Numerical Solution #2
% * $M = 20$, $dt = 0.01$
solver_1_5_7(F_, a_, b_, f_, 20, 0.01, 2 * pi);