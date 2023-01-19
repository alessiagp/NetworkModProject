% Model without treatment (d=0)
% X(1) = A(t), X(2) = E(t)
tspan = [0;6000];
X0 = [1000,3.5e5];

[t,X] = ode45(@model,tspan,X0);

figure;
plot(t,X)
legend('Cancer Cells','Effector Cells')
xlabel('Time (h)')
ylabel('Number')

function dX = model(t,X)
r = 0.01;
K = 4e6;
mu_a = 2e-12;
mu_e = 4e-5;
p = 4e-14;
c = 1e2;
mu_ea = 4e-15;

dX = [r*X(1)*(1-(X(1)/K)) - mu_a*X(1)*X(2);
    -mu_e*X(2) + (p*X(1)*X(2))/(c+X(1)) - mu_ea*X(1)*X(2)];
end