%% HW 1.5.7
%% Eqution
% * $v_t = \nu v_{xx} + F(x, t), x \in (0, 1), t > 0$
% * $v(x, 0) = f(x), x \in [0, 1]$
% * $v(0, t) = a(t), v(1, t) = b(t), t \geq 0$
%% Numerical Solver
%
% <include>solver_1_5_7.m</include>
%

%% Numerical Solution
% * $a(t) = 0$
% * $b(t) = 0$
% * $f(x) = 0$
% * $F(x, t) = \sin(2 \pi x) \sin(4 \pi t)$
% * $M = 10$, $dt = 0.05$
a_ = @(t)( 0 );
b_ = @(t)( 0 );
f_ = @(x)( 0 );
F_ = @(x, t)( sin(2 * pi * x) * sin(4 * pi * t) );
solver_1_5_7(F_, a_, b_, f_, 10, 0.05);

%% Numerical Solution #2
% * $M = 20$, $dt = 0.01$
solver_1_5_7(F_, a_, b_, f_, 20, 0.01);