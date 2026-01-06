%% Final Project Submission: Optimization of a Combined Brayton–Rankine Cycle
% Our goal is to maximize efficiency by varying rp(Brayton cycle pressure ratio) and Ps( Rankine cycle pump outlet pressure) with the given constraints
clc; clear all;

%% Fixed Parameters 
T1 = 298;            % in K (compressor inlet temperature)
P1 = 1;                 % in bar
eta_c = 0.88;           % compressor isentropic efficiency
eta_gt = 0.90;          % gas turbine isentropic efficiency
eta_st = 0.92;          % steam turbine isentropic efficiency
eta_p = 0.60;           % pump isentropic efficiency
T_cond = 30 + 273.15;   % in K (condenser temp)
dT_pinch = 10;          % in K (pinch point difference) 
dT_super = 20;          % in K (T4 - T_steam_max) 
T4_limit = 650 + 273.15;% in K (Exhaust constraint since T4<650 as mentioned in the paper) 

% Thermodynamic constants (air/exhaust)
Cp_g = 1.005;           % in kJ/kgK
gamma = 1.4;
R = 0.287;              % in kJ/kgK

%% Optimization and Parametric Study 
T3_range = (900:100:1400) + 273.15; % Range of T3 is from 900C to 1400C
rp_range = 5:0.5:35;                % Brayton pressure ratio sweep
Ps_range = 10:2:100;                % Rankine pump pressure sweep (in bar)

parametric_results = []; % to store max efficiency for each temp T3
best_overall.eta = 0;

for T3 = T3_range
    best_eta_for_T3 = 0;
    for rp = rp_range
        % Brayton cycle calculations
        T2s = T1 * (rp)^((gamma-1)/gamma); % since for isentropic processes, T2/T1 = rp^((gamma-1)/gamma)
        T2 = T1 + (T2s - T1)/eta_c;
        T4s = T3 * (1/rp)^((gamma-1)/gamma); % isentropic processes => T2/T1 = rp^((gamma-1)/gamma)
        T4 = T3 - eta_gt * (T3 - T4s);
        
        % Constraint check: exhaust temperature
        if T4 > T4_limit, continue; end
        
        W_comp = Cp_g * (T2 - T1);
        W_gt = Cp_g * (T3 - T4);
        Q_in_fuel = Cp_g * (T3 - T2);
        
        for Ps = Ps_range
            % Rankine cycle coupling
            % Approximation of water saturation temp (saturation line)
            T_sat = 100 + 15 * log(Ps); % this is a logarithmic approximation of the saturated liquid-vapor curve for water
            
            % Checking pinch point constraint
            if T4 < T_sat + dT_pinch, continue; end
            
            % Enthalpy approximations (kJ/kg) for mw/mg calculation 
            T10 = T4 - dT_super; % Superheated steam temperature
            h10 = 2500 + 2.0 * (T10 - 273.15); % Approx steam enthalpy
            h8 = 4.18 * (T_sat - 273.15);     % Approx liquid enthalpy
            h7 = 4.18 * (T_cond - 273.15) + (Ps/10)/eta_p; % Pump outlet
            
            % Mass flow ratio (mw/mg)  
            mw_mg = (Cp_g * (T4 - (T_sat + dT_pinch))) / (h10 - h8); % corresponds to eq (8) from the paper
            
            % Work and efficiency 
            W_st = mw_mg * (h10 - 2300) * eta_st; % Approx steam turbine work
            W_net_total = (W_gt - W_comp) + W_st;
            eta_comb = W_net_total / Q_in_fuel;
            
            % Storing best for this T3 (Parametric analysis)
            if eta_comb > best_eta_for_T3
                best_eta_for_T3 = eta_comb;
            end
            
            % Storing best overall in the struct
            if eta_comb > best_overall.eta
                best_overall.eta = eta_comb;
                best_overall.T3 = T3;
                best_overall.rp = rp;
                best_overall.Ps = Ps;
                best_overall.states_B = [T1, T2, T3, T4];
                best_overall.P_B = [P1, rp, rp, 1];
                best_overall.states_R = [T_cond, T_sat, T10, T_cond];
            end
        end
    end
    parametric_results = [parametric_results; T3-273.15, best_eta_for_T3*100];
end

%% Displaying Results 
fprintf('Optimized Combined Cycle Results ->\n');
fprintf('Maximum Overall Efficiency: %.2f%%\n', best_overall.eta * 100);
fprintf('Optimal GT Inlet Temperature (T3): %.2f K (%.0f C)\n', best_overall.T3, best_overall.T3-273.15);
fprintf('Optimal Brayton Pressure Ratio (rp): %.2f\n', best_overall.rp);
fprintf('Optimal Rankine Pump Pressure (Ps): %.2f bar\n', best_overall.Ps);

%% Plotting the diagrams
figure('Color', 'w', 'Position', [50, 50, 1200, 400]);

% Subplot 1: Combined T-s Diagram 
subplot(1,3,1); hold on; grid on;
s_B = Cp_g * log(best_overall.states_B/T1) - R * log(best_overall.P_B/P1);
plot([s_B, s_B(1)], [best_overall.states_B, best_overall.states_B(1)], 'r-o', 'LineWidth', 1.5);
% Rankine plotting 
s_R = [0.5, 0.7, 1.8, 1.6]; % Representative entropy for water
plot([s_R, s_R(1)], [best_overall.states_R, best_overall.states_R(1)], 'b-s', 'LineWidth', 1.5);
title('Combined T-s diagram'); xlabel('s [kJ/kgK]'); ylabel('T [K]');
legend('Brayton cycle', 'Rankine cycle', 'Location', 'best');

% Subplot 2: Brayton P-v diagram
subplot(1,3,2); hold on; grid on;
v_B = (R * best_overall.states_B) ./ (best_overall.P_B * 100); 
plot([v_B, v_B(1)], [best_overall.P_B, best_overall.P_B(1)], 'k-d', 'LineWidth', 1.5);
title('Brayton p-v diagram'); xlabel('v [m^3/kg]'); ylabel('P [bar]');

% Subplot 3: Parametric analysis (Influence of T3 on efficiency) 
subplot(1,3,3); plot(parametric_results(:,1), parametric_results(:,2), 'm-p', 'LineWidth', 2);
grid on; title('Efficiency vs T_3');
xlabel('Gas turbine inlet temp T_3 [^\circC]'); ylabel('Max combined efficiency [%]');

sgtitle(['Optimization Results for T3 = ', num2str(best_overall.T3-273.15), ' ^\circC']);