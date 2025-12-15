clear; clc;

% Given 
Tc = 270;          % in K
Th = 320;          
alpha = 0.02;
k = 50;
rp_min = 1.1;
rp_max = 6;
Sgen_max = 0.05;   % in kJ/K (entropy generation limit)

COP = @(rp) (Tc./(Th - Tc)) .* (1 - alpha*(rp - 1).^2 ./ rp);
Wc  = @(rp) k*(sqrt(rp) - 1);

% Heat delivered
Qh = @(rp) COP(rp).*Wc(rp);

% Entropy generation model 
% Entropy generation = Qh*(1/Tc - 1/Th)
Sgen = @(rp) Qh(rp).*(1/Tc - 1/Th);

% Optimization 
% Objective: maximize Qh/Wc => minimize negative
objfun = @(rp) -Qh(rp)./Wc(rp);

% Nonlinear entropy constraint
nonlcon = @(rp) deal(Sgen(rp) - Sgen_max, []);

rp0 = 2.5; % initial guess

% Solving optimization
rp_opt = fmincon(objfun, rp0, [], [], [], [], rp_min, rp_max, nonlcon);

% Parameter sweep for plotting 
rp = linspace(rp_min, rp_max, 400);
COP_vals  = COP(rp);
Sgen_vals = Sgen(rp);
QhW_vals  = Qh(rp)./Wc(rp);
feasible = (Sgen_vals <= Sgen_max);

% Plotting figure
figure;

subplot(3,1,1)
plot(rp, COP_vals, 'b','LineWidth',2); hold on;
xline(rp_opt,'r--','LineWidth',1.5);
xlabel('Pressure ratio r_p');
ylabel('COP');
title('Coefficient of Performance');
grid on;

subplot(3,1,2)
plot(rp, Sgen_vals, 'm','LineWidth',2); hold on;
yline(Sgen_max,'k--','LineWidth',1.5);
xline(rp_opt,'r--','LineWidth',1.5);
xlabel('Pressure ratio r_p');
ylabel('Entropy Generation (kJ/K)');
title('Entropy Generation Constraint');
grid on;

subplot(3,1,3)
plot(rp, QhW_vals, 'g','LineWidth',2); hold on;
plot(rp(feasible), QhW_vals(feasible), 'k','LineWidth',3);
xline(rp_opt,'r--','LineWidth',1.5);
xlabel('Pressure ratio r_p');
ylabel('Q_h / W_c');
title('Heat Delivered per Unit Work');
legend('All','Feasible region','Optimal r_p','Location','Best');
grid on;

% results 
fprintf('Optimal pressure ratio = %.3f\n', rp_opt);
fprintf('COP at optimum        = %.3f\n', COP(rp_opt));
fprintf('Entropy generation    = %.4f kJ/K\n', Sgen(rp_opt));
