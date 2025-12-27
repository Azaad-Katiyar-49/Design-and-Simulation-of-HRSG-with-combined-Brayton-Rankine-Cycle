clc; clear all;

%% Question 2 – Chasing Power in the Thin Air (BRAYTON CYCLE)

%% TASK(a): IDEAL BRAYTON CYCLE (Assuming constant specific heats)

% Input parameters
rp = 10;                 % Compressor pressure ratio (P2/P1)
T1 = 300;                % Ambient inlet temperature (in K)
P1 = 100e3;              % Ambient inlet pressure (in Pa)
T3 = 1400;               % Turbine inlet temperature (K)

gamma = 1.4;             % Ratio of specific heats for air
cp = 1.005;              % in kJ/kg-K
R  = 0.287;              % in kJ/kg-K

% State 1: Compressor inlet
P2 = rp * P1; % since rp = P2/P1

% State 2: Compressor outlet (after isentropic compression)
T2 = T1 * rp^((gamma-1)/gamma); % since for isentropic processes, T2/T1 = rp^((gamma-1)/gamma)

% State 3: Turbine inlet (after constant pressure heat addition)
P3 = P2; % since constant pressure is maintained

% State 4: Turbine outlet (after isentropic expansion)
P4 = P1; % since when heat is removed from 4->1, pressure remains constant
T4 = T3 * (1/rp)^((gamma-1)/gamma); % since 3->4 is an isentropic process


fprintf('T1 = %.2f K, P1 = %.2f kPa\n', T1, P1/1000);
fprintf('T2 = %.2f K, P2 = %.2f kPa\n', T2, P2/1000);
fprintf('T3 = %.2f K, P3 = %.2f kPa\n', T3, P3/1000);
fprintf('T4 = %.2f K, P4 = %.2f kPa\n\n', T4, P4/1000);

%% TASK(b): REAL BRAYTON CYCLE (Non-Ideal Compressor and Turbine)

eta_c = 0.85;   % Compressor isentropic efficiency
eta_t = 0.88;   % Turbine isentropic efficiency

% since pressure ratio is fixed by design, P2 and P4 will be unchanged

% Actual compressor outlet temperature
T2_actual = T1 * (1 + (rp^((gamma-1)/gamma) - 1)/eta_c);

% Actual turbine outlet temperature
T4_actual = T3 * (1 - eta_t*(1 - (1/rp)^((gamma-1)/gamma)));


fprintf('T2_actual = %.2f K\n', T2_actual);
fprintf('T4_actual = %.2f K\n\n', T4_actual);

%% TASK(c): VARIABLE SPECIFIC HEAT MODEL

% Temperature dependent cp(T) model (assume a polynomial)
cp_T = @(T) (1.003 + 1.5e-4*(T-300) + 1e-8*(T-300).^2); % in kJ/kg-K

% Numerical integration for work and heat
W_comp_var = integral(@(T) cp_T(T), T1, T2_actual); % 1->2
W_turb_var = integral(@(T) cp_T(T), T4_actual, T3); % 3->4
Qin_var = integral(@(T) cp_T(T), T2_actual, T3);    % 2->3

Wnet_var = W_turb_var - W_comp_var;   % Net work output (kJ/kg)
eta_th_var = Wnet_var / Qin_var;                  % Thermal efficiency


fprintf('VARIABLE Cp MODEL RESULTS\n');
fprintf('Compressor Work (variable cp) = %.3f kJ/kg\n', W_comp_var);
fprintf('Turbine Work (variable cp) = %.3f kJ/kg\n', W_turb_var);
fprintf('Net Work Output (variable cp) = %.3f kJ/kg\n', Wnet_var);
fprintf('Heat Added (variable cp) = %.3f kJ/kg\n', Qin_var);
fprintf('Thermal Efficiency (variable cp) = %.3f\n\n', eta_th_var);

% Comparison (Variable and constant cp models):
% Variable cp generally predicts lower net work and efficiency because 
% as cp increases on increasing T, Qin increases faster, W_turb increases 
% less rapidly, and thus ratio Wnet/Qin decreases

%% TASK(d): PERFORMANCE CALCULATIONS

% Constant cp model 
Wc = cp * (T2_actual - T1);       % Compressor work (kJ/kg)
Wt = cp * (T3 - T4_actual);       % Turbine work (kJ/kg)
Wnet = Wt - Wc;                   % Net work output (kJ/kg)
Qin = cp * (T3 - T2_actual);      % Heat added (kJ/kg)
eta_th = Wnet / Qin;              % Thermal efficiency
work_ratio = Wnet / Wt;           % Work ratio

fprintf('CONSTANT Cp MODEL RESULTS\n');
fprintf('Compressor Work = %.3f kJ/kg\n', Wc);
fprintf('Turbine Work = %.3f kJ/kg\n', Wt);
fprintf('Net Work Output = %.3f kJ/kg\n', Wnet);
fprintf('Heat Added = %.3f kJ/kg\n', Qin);
fprintf('Thermal Efficiency = %.3f\n', eta_th);
fprintf('Work Ratio = %.3f\n\n', work_ratio);

% Relevance of work ratio in gas turbine operation :
% Work ratio is the ratio of net work and turbine work. It indicates the 
% fraction of turbine work available as useful output after driving the compressor.

%% TASK(e): T–s AND P–v DIAGRAMS

% Entropy calculation (reference at state 1)
s1 = 0;
s2 = cp*log(T2/T1) - R*log(P2/P1);
s3 = cp*log(T3/T2);
s4 = cp*log(T4/T3) - R*log(P4/P3);

s = [s1 s2 s3 s4 s1];
T = [T1 T2 T3 T4 T1];

% T–s plot
figure;
plot(s, T, '-o', 'LineWidth', 2);
xlabel('Entropy (kJ/kg·K)');
ylabel('Temperature (K)');
title('T–s Diagram of Brayton Cycle');
grid on;

% Specific volume calculation
v1 = R*T1/P1;
v2 = R*T2/P2;
v3 = R*T3/P3;
v4 = R*T4/P4;

v = [v1 v2 v3 v4 v1];
P = [P1 P2 P3 P4 P1]/1000;

% P–v plot
figure;
plot(v, P, '-o', 'LineWidth', 2);
xlabel('Specific Volume (m^3/kg)');
ylabel('Pressure (kPa)');
title('P–v Diagram of Brayton Cycle');
grid on;

%% TASK(f): PRESSURE RATIO SWEEP

rp_range = linspace(2, 25, 40);
eta_array = zeros(size(rp_range));
Wnet_array = zeros(size(rp_range));

for i = 1:length(rp_range)
    rp = rp_range(i);

    T2 = T1 * (1 + (rp^((gamma-1)/gamma) - 1)/eta_c);
    T4 = T3 * (1 - eta_t*(1 - (1/rp)^((gamma-1)/gamma)));

    Wc = cp * (T2 - T1);
    Wt = cp * (T3 - T4);

    Wnet_array(i) = Wt - Wc;
    eta_array(i) = (Wt - Wc)/(cp*(T3 - T2));
end

% Plotting net work
figure;
plot(rp_range, Wnet_array, 'LineWidth', 2);
xlabel('Pressure Ratio');
ylabel('Net Work Output (kJ/kg)');
title('Net Work Output vs Pressure Ratio');
grid on;

% Plotting thermal efficiency
figure;
plot(rp_range, eta_array, 'LineWidth', 2);
xlabel('Pressure Ratio');
ylabel('Thermal Efficiency');
title('Thermal Efficiency vs Pressure Ratio');
grid on;

% Identifying the pressure ratio that maximizes net work output 
[Wnet_max, idx] = max(Wnet_array);
rp_opt = rp_range(idx);
fprintf('Maximum net work occurs at pressure ratio rp = %.2f\n', rp_opt);

% Interpretation:
% Net work output peaks at an optimal pressure ratio where turbine work 
% gain balances compressor work increase. Efficiency continues to increase 
% but the increment is continuously diminishing.



%% ALTITUDE EFFECT :
% Reduction in inlet pressure due to altitude lowers ambient air density 
% which lowers mass flow rate.
% This also reduces pressure ratio effectiveness because when P1 decreases,
% to maintain the same rp, the compressor must produce very low exit 
% pressures.
% These leads to reducing both net work output and thermal efficiency of 
% the Brayton cycle.
