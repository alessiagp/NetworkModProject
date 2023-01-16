clear
% Model one treatment
% X(1) = A(t), X(2) = E(t), X(3) = C(t)
d = 1.25 * 2.4e18; % for Cyt 62.5 mg/kg

Tf = 1000;
tspan = [0;Tf];
X0 = [1e10,1e6,d];
op1 = odeset('MaxStep',1);

[t,X] = ode45(@model,tspan,X0);

figure;
plot(t,X)
legend('Cancer Cells','Effector Cells','Drug Molecules')
xlabel('Time (h)')
ylabel('Number')

function dX = model(t,X)
r = 0.01;
K = 4e6;
mu_ac = 0.012; % for Cyt 62.5 mg/kg
mu_ca = mu_ac * 10;
d = 1.25 * 2.4e18; % for Cyt 62.5 mg/kg
mu_c = 0.231;
mu_a = 2e-12;
a = 2e3;
p = 4e-14;
mu_ec = 417;
c = 1e2;
b = 5e6;
mu_ea = 4e-15;
mu_e = 4e-5;
tau = 24;

A = r*X(1)*(1-(X(1)/K)) - mu_a*X(1)*X(2) - ((mu_ac*X(1)*X(3))/(a+X(3)));
E = -mu_e*X(2) + (p*X(1)*X(2))/(c+X(1)) - mu_ea*X(1)*X(2)  - (mu_ec*X(2)*X(3))/(b+X(3));
if rem(round(t(1),0),tau) == 0
    C = d - mu_c*X(3) - (mu_ca*X(3)*X(1))/(a+X(1)); 
else
    C = - mu_c*X(3) - (mu_ca*X(3)*X(1))/(a+X(1));
end
dX = [A;E;C];
end