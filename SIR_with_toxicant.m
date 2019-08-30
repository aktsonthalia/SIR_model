timespan = [0 300];

x0 = [100; 100; 100; 100; 100; 100; 100];

A = 250;
delta0 = 0.2;
delta1 = 0.005;
delta2 = 0.002;
delta3 = 0.18;
alpha = 0.0005;
beta = 0.054;
gamma = 0.008;
r = 0.0002;
a = 0.2;
b = 0.0003;
mu = 0.005;
mu0 = 0.24;
B = 150;
m = 0.2;
n = 0.05;
p = 0.2;

sys = system(A, delta0, delta1, delta2, delta3, alpha, beta, gamma, r, a, b, mu, mu0, B, m, n, p);

[T, sol] = ode45(sys, timespan, x0);

plot(T, sol(:, 1));
hold on;
plot(T, sol(:, 2));
hold on;
plot(T, sol(:, 3));
hold on;
plot(T, sol(:, 6));
hold on;
plot(T, sol(:, 7));
grid on;
plot(T, sol(:, 1) + sol(:, 2) + sol(:, 3) + sol(:, 5));

legend('Susceptible', 'Infected', 'Susceptible aware','Toxicant in Environment', 'Toxicant inside Organism', 'Total Population');

function f = system(A, delta0, delta1, delta2, delta3, alpha, beta, gamma, r, a, b, mu, mu0, B, m, n, p)
f = @(t, x) [A - delta0*x(1) - alpha*x(1)*x(2) - beta*x(1)*x(4)/(1+gamma*x(4)) + delta3*x(3) - r*x(7)*x(1);
             alpha*x(1)*x(2) - (delta0+delta1+delta2)*x(2) - a*x(2)/(1+b*x(2)) - r*x(7)*x(2);
             beta*x(1)*x(4)/(1+gamma*x(4)) - (delta0+delta3)*x(3) - r*x(7)*x(3);
             mu*x(2) - mu0*x(4);
             delta2*x(2) - delta0*x(5) + a*x(2)/(1+b*x(2)) - r*x(7)*x(5);
             B - m*x(6);
             n*x(6) - p*x(7);
            ];
end
