clear
% Model one treatment
% X(1) = A(t), X(2) = E(t), X(3) = C(t)
d = 1.25 * 2.4e18; % for Cyt 62.5 mg/kg

Tf = 100;
tspan = [0;Tf];
X0 = d;
op1 = odeset('MaxStep',1);
[t,X] = ode45(@model,tspan,X0,op1);

figure;
plot(t,X,'+')
legend('Cancer Cells','Effector Cells','Drug Molecules')
xlabel('Time (h)')
ylabel('Number')

function dX = model(t,X)
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
nInject = round(100/24, 0);

if rem(round(t(1),0),tau) == 0
    dX = d - mu_c*X(1); 
else
    dX = - mu_c*X(1);
end
end