clear;
%%%%%%%%%%%%% Complete model

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

%mu_ac = 0.012; % Cyt 62.5 mg/kg 3 days
%mu_ac = 0.001; % Cyt 0.12 mg/kg 5 days
%mu_ac = 0.0041; % Ibr 9 mg/kg on days 1-5 and 8-10
%mu_ac = 0.0042; % Ibr 18 mg/kg on days 1-5 and 8-10
mu_ac = 0.0161; % Cyt 62.5 + Ibr 9 mg/kg days 1-5 and 8-10
%mu_ac = 0.0043; % Ibr 25 mg/kg on days 1-5 and 8-10

mu_ca = mu_ac * 10; % Deactivation of drug due to killing of A20 [h^-1]

%%% Dose [molecules/mouse]

%d = 1.25 * 2.4e18; % Cyt 62.5 mg/kg
%d = 2.4e-3 * 2.4e18; % Cyt 0.12 md/kg
d = 0.18 * 1.4e18; % Ibr 9 mg/kg
%d = 0.36 * 1.4e18; % Ibr 18 mg/kg
%d = 1.25 * 2.4e18 + 0.18 * 1.4e18; %Cyt 62.5 + Ibr 9
%d = 0.5 * 1.4e18; % Ibr 25 mg/kg

%%% Chemical deactivation rate [h^-1]

%mu_c = 0.231; % Cyt
mu_c = 0.116; % Ibr
%mu_c = 0.221; % Cyt 62.5 + Ibr 9

mu_a = 2e-12;   % Effectors-A20 interaction coefficient [h^-1]
a = 2e3;    % Drug amount producing 50% of max effect on A20 [molecules]
p = 4e-14; % Production rate of effectors stimulated by A20 [h^-1]
mu_ec = 417; % Mortality rate of drug on effector cells [h^-1]
c = 1e2;    % Num of A20 producing 50% of max immune activation [cells]
b = 5e6;    % Drug amount producing 50% max effect on healty cells [molecules]
mu_ea = 4e-15; % A20-effectors interaction coefficient [h^-1]
mu_e = 4e-5; % Death rate of effecors [h^-1]


e0 = 0; % Natural production of effector cells [cells/h]
%e0 = 4.2;
tau = 24;

A = r*X(1)*(1-(X(1)/K)) - mu_a*X(1)*X(2) - ((mu_ac*X(1)*X(3))/(a+X(3)));
E = e0 -mu_e*X(2) + (p*X(1)*X(2))/(c+X(1)) - mu_ea*X(1)*X(2)  - (mu_ec*X(2)*X(3))/(b+X(3));
if t(1) > 336 && rem(round(t(1),0),tau) == 0 && t(1) < 432 || rem(round(t(1),0),tau) == 0 && (t(1) > 480)&&(t(1) < 528)
%if t(1) > 336 && rem(round(t(1),0),tau) == 0 && t(1) < 432
%if t(1) > 336 && rem(round(t(1),0),tau) == 0 && t(1) < 384
%if t(1) < 0
    C = d - mu_c*X(3) - (mu_ca*X(3)*X(1))/(a+X(1)); 
else
    C = - mu_c*X(3) - (mu_ca*X(3)*X(1))/(a+X(1));
end
dX = [A;E;C];
end
