clear all;

% Given 
R = 8.314e-3;        
c = 1.8;             
A_true = 1e5;        % true pre-exponential factor
E_true = 45;         % true activation energy (kJ/mol)
T0 = 300;            % initial temperature (in K)
tspan = linspace(0, 5000, 200);   % time grid (s)

% Generating synthetic data 
params_true = [A_true, E_true];
[t_true, T_true] = ode45(@(t,T) reactor_ode(t,T,params_true,c,R), ...
                         tspan, T0);

% Adding 1% noise
noise_level = 0.01;
T_noisy = T_true .* (1 + noise_level*randn(size(T_true)));

% Parameter estimation 
% Initial guesses
params_guess = [5e4, 30];   % [A, E]

% Lower and upper bounds
lb = [1e3, 10];
ub = [1e7, 100];

% Least-squares fitting
params_est = lsqcurvefit(@(p,t) model_wrapper(p,t,T0,c,R), params_guess, tspan, T_noisy, lb, ub);

A_est = params_est(1);
E_est = params_est(2);

% Model prediction with estimated params 
[t_est, T_est] = ode45(@(t,T) reactor_ode(t,T,params_est,c,R), tspan, T0);

% Plotting results 
figure;
plot(tspan, T_noisy, 'ko', 'MarkerSize',4); hold on;
plot(tspan, T_true, 'b-', 'LineWidth',2);
plot(tspan, T_est, 'r--', 'LineWidth',2);
xlabel('Time (s)');
ylabel('Temperature (K)');
legend('Noisy data','True model','Fitted model');
title('Reaction Heat Parameter Estimation');
grid on;

% results 
fprintf('True A = %.2e, Estimated A = %.2e\n', A_true, A_est);
fprintf('True E = %.2f kJ/mol, Estimated E = %.2f kJ/mol\n', E_true, E_est);


function dTdt = reactor_ode(t, T, params, c, R)
    A = params(1);
    E = params(2);
    rT = A * exp(-E/(R*T));          % reaction heat
    qext = 2000 * exp(-0.001*t);     % external heating  
    dTdt = (rT + qext) / c;
end


function Tmodel = model_wrapper(params, tspan, T0, c, R)
    [~, Tmodel] = ode45(@(t,T) reactor_ode(t,T,params,c,R),tspan, T0);
end
