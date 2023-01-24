clear;
%%%%%%%%%%%%% Model with Cytarabine infusion, Cyt_low group

% X(1) = A(t), X(2) = E(t), X(3) = C(t)

Tf = 1440; % Final time
tspan = [0;Tf];
X0 = [5e4,2500,0]; % Starting conditions
op1 = odeset('MaxStep',1,'RelTol',1e-2,'AbsTol',1e-4);

[t,X] = ode45(@model,tspan,X0,op1);

figure;
plot(t,X)
legend('Cancer Cells','Effector Cells','Drug Molecules')
xlabel('Time (h)')
ylabel('Number')

function dX = model(t,X)
r = 0.01;   % A20 growth rate [h^-1]
K = 4e6;    % Max A20 number [cells/mouse]

%%% Cytotoxicity rate with drug [h^−1]

mu_ac1 = 0.0015; % Cyt 200 mg * 24h / m^2, 7 days of infusion
mu_ac2 = 0.021; % Cyt 1000 mg * 3h / m^2, 6 days with 3h long infusions every 12h

mu_ca1 = mu_ac1 * 10; % Deactivation of drug due to killing of A20 [h^-1]
mu_ca2 = mu_ac2 * 10;

%%% Dose [molecules] 

dt1 = 0.00198849634391836; % mean step size when the treatment is delivered
dt2 = 0.0019462106821853881;
d1 = (9.8 * 2.4e18 * dt1) / (7 * 24) ; % Cyt 200 mg/m^2 nr of molecules per time step
d2 = (2.3 * 2.4e18 * dt2) / 3 ; % Cyt 1000 mg/m^2 nr of molecules per time step

%%% Chemical deactivation rate [h^-1]

mu_c = 0.231; % Cyt

mu_a = 2e-12;   % Effectors-A20 interaction coefficient [h^-1]
a = 2e3;    % Drug amount producing 50% of max effect on A20 [molecules]
p = 4e-14; % Production rate of effectors stimulated by A20 [h^-1]
mu_ec = 417; % Mortality rate of drug on effector cells [h^-1]
c = 1e2;    % Num of A20 producing 50% of max immune activation [cells]
b = 5e6;    % Drug amount producing 50% max effect on healty cells [molecules]
mu_ea = 4e-15; % A20-effectors interaction coefficient [h^-1]
mu_e = 4e-5; % Death rate of effecors [h^-1]

e0 = 4.2; % Natural production of effector cells [cells/h]

E = e0 -mu_e*X(2) + (p*X(1)*X(2))/(c+X(1)) - mu_ea*X(1)*X(2)  - (mu_ec*X(2)*X(3))/(b+X(3));
if t < 336
    C = 0;
    A = r*X(1)*(1-(X(1)/K)) - mu_a*X(1)*X(2);
elseif t > 336 && t < 504
    C = d1 - mu_c*X(3) - (mu_ca1*X(3)*X(1))/(a+X(1));
    A = r*X(1)*(1-(X(1)/K)) - mu_a*X(1)*X(2) - ((mu_ac1*X(1)*X(3))/(a+X(3)));
elseif t > 504 && t < 672
    C = - mu_c*X(3) - (mu_ca1*X(3)*X(1))/(a+X(1));
    A = r*X(1)*(1-(X(1)/K)) - mu_a*X(1)*X(2) - ((mu_ac1*X(1)*X(3))/(a+X(3)));
elseif t > 672 && t < 675 || t > 684 && t < 687 || t > 696 && t < 699 || t > 708 && t < 711 || t > 720 && t < 723 || t(1) > 732 && t(1) < 735 || t(1) > 744 && t(1) < 747 || t(1) > 756 && t(1) < 759 || t(1) > 768 && t(1) < 771 || t(1) > 780 && t(1) < 783 || t > 792 && t < 795 || t > 804 && t < 807
    C = d2 - mu_c*X(3) - (mu_ca2*X(3)*X(1))/(a+X(1));
    A = r*X(1)*(1-(X(1)/K)) - mu_a*X(1)*X(2) - ((mu_ac2*X(1)*X(3))/(a+X(3)));
else
    C = - mu_c*X(3) - (mu_ca2*X(3)*X(1))/(a+X(1));
    A = r*X(1)*(1-(X(1)/K)) - mu_a*X(1)*X(2) - ((mu_ac2*X(1)*X(3))/(a+X(3)));
end
dX = [A;E;C];
end
