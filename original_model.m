syms t x A delta0 delta1 delta2 delta3 alpha beta gamma a b mu mu0 delta

% x(1) = I
% x(2) = S_a
% x(3) = N
% x(4) = M

timespan = [0 100];
x0 = [200; 200; 200; 200];

sys = system(250, 0.2, 0.005, 0.002, 0.18, 0.0005, 0.0022, 0.008, 0.2, 0.0003, 0.005, 0.24);
[T, sol] = ode45(sys, timespan, x0);

plot(T, sol(:, 1));
hold on;
plot(T, sol(:, 2));
hold on;
plot(T, sol(:, 3));
hold on;
plot(T, sol(:, 4));
grid on;







function f = system(A, delta0, delta1, delta2, delta3, alpha, beta, gamma, a, b, mu, mu0)
    delta = delta0 + delta1 + delta2;
    f = @(t, x) [alpha*(x(3) - x(1) - x(2))*x(1) - delta*x(1) - a*x(1)/(1+b*x(1));
                 beta*(x(3) - x(1) - x(2))*x(4)/(1+gamma*x(4)) - (delta0 + delta3)*x(2);
                 A - delta0*x(3) - (delta1+delta2)*x(1) - a*x(1)/(1+b*x(1));
                 mu*x(1) - mu0*x(4)];
end



