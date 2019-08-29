syms t x A delta0 delta1 delta2 delta3 alpha beta gamma a b mu mu0 delta beta1 beta2 beta3 beta4

% x(1) = I
% x(2) = S_a
% x(3) = N
% x(4) = M

beta1 = 0;
beta2 = 0.001;
beta3 = 0.002;
beta4 = 0.003;

timespan = [0 400];
x0 = [242; 200; 200; 200];

sys1 = system(250, 0.2, 0.005, 0.002, 0.18, 0.0005, beta1, 0.008, 0.2, 0.0003, 0.005, 0.24);
[T1, sol1] = ode45(sys1, timespan, x0);

sys2 = system(250, 0.2, 0.005, 0.002, 0.18, 0.0005, beta2, 0.008, 0.2, 0.0003, 0.005, 0.24);
[T2, sol2] = ode45(sys2, timespan, x0);

sys3 = system(250, 0.2, 0.005, 0.002, 0.18, 0.0005, beta3, 0.008, 0.2, 0.0003, 0.005, 0.24);
[T3, sol3] = ode45(sys3, timespan, x0);

sys4 = system(250, 0.2, 0.005, 0.002, 0.18, 0.0005, beta4, 0.008, 0.2, 0.0003, 0.005, 0.24);
[T4, sol4] = ode45(sys4, timespan, x0);

plot(T1, sol1(:, 1));
hold on;
plot(T2, sol2(:, 1));
hold on;
plot(T3, sol3(:, 1));
hold on;
plot(T4, sol4(:, 1));
grid on;







function f = system(A, delta0, delta1, delta2, delta3, alpha, beta, gamma, a, b, mu, mu0)
    delta = delta0 + delta1 + delta2;
    f = @(t, x) [alpha*(x(3) - x(1) - x(2))*x(1) - delta*x(1) - a*x(1)/(1+b*x(1));
                 beta*(x(3) - x(1) - x(2))*x(4)/(1+gamma*x(4)) - (delta0 + delta3)*x(2);
                 A - delta0*x(3) - (delta1+delta2)*x(1) - a*x(1)/(1+b*x(1));
                 mu*x(1) - mu0*x(4)];
end



