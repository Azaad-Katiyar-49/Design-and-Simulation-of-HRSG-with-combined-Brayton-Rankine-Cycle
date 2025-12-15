clear all;

% Given 
T0 = 300;             
Tf = 900;             
Tb = 300;              % boundary temperature (K)
tf = 2000;             % total time (s)
c  = 2.5;              
Qtot = 5e5;            % total heat input (kJ)

% Time discretization 
N = 100;                       % number of time intervals
dt = tf/N;                     % time step
t  = linspace(0, tf, N+1);     % time vector

% Initial guess: uniform heating 
q0 = (Qtot/tf)*ones(N,1);      % kJ/s (uniform heat rate)

% Objective function (to minimize total entropy generation)
objfun = @(q) entropy_objective(q, T0, c, Tb, dt);

% Equality constraint (total heat fixed) 
Aeq = dt*ones(1,N);    % Integral of q(t)
beq = Qtot;

% Bounds on heat rate 
lb = zeros(N,1);       % No negative heating
ub = [];               % No upper bound

% Optimization options 
options = optimoptions('fmincon','Display','iter','Algorithm','sqp');

% Solving optimization problem 
q_opt = fmincon(objfun, q0, [], [], Aeq, beq, lb, ub, [], options);

% Temperature evolution 
T_opt = temperature_profile(q_opt, T0, c, dt);
T_uni = temperature_profile(q0,    T0, c, dt);

% Entropy generation comparison 
Sgen_opt = objfun(q_opt);
Sgen_uni = objfun(q0);

% Plotting results 
figure;

subplot(2,1,1)
plot(t(1:end-1), q0, '--k', 'LineWidth',1.5); hold on;
plot(t(1:end-1), q_opt, 'r', 'LineWidth',2);
xlabel('Time (s)');
ylabel('Heat rate q(t) (kJ/s)');
legend('Uniform heating','Optimal heating');
title('Heating Profiles');
grid on;

subplot(2,1,2)
plot(t, T_uni, '--k', 'LineWidth',1.5); hold on;
plot(t, T_opt, 'r', 'LineWidth',2);
xlabel('Time (s)');
ylabel('Temperature (K)');
legend('Uniform heating','Optimal heating');
title('Temperature Evolution');
grid on;

% results 
fprintf('Total entropy generation (uniform) = %.2f kJ/K\n', Sgen_uni);
fprintf('Total entropy generation (optimal) = %.2f kJ/K\n', Sgen_opt);

% defining entropy objective function 
function Sgen = entropy_objective(q, T0, c, Tb, dt)
    % Computes total entropy generation for a given heat profile q
    N = length(q);
    T = zeros(N+1,1);
    T(1) = T0;
    for i = 1:N
        T(i+1) = T(i) + (dt/c)*q(i);
    end
    sgen_inst = q./T(1:end-1) - q/Tb;
    Sgen = sum(sgen_inst)*dt;
end

function T = temperature_profile(q, T0, c, dt)
    % Computes temperature evolution
    N = length(q);
    T = zeros(N+1,1);
    T(1) = T0;
    for i = 1:N
        T(i+1) = T(i) + (dt/c)*q(i);
    end
end
