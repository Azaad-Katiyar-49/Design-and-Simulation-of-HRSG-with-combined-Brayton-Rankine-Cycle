clc; clear all; 

%% Question 1 - The Steam That Refused to Go to Waste (RANKINE CYCLE)

%% TASK(a): IDEAL RANKINE CYCLE 

% Input parameters
P_boiler = 15e6;        % Boiler pressure (in Pa)
P_cond   = 10e3;        % Condenser pressure (in Pa)
T3       = 823;         % Turbine inlet temperature (in K)

% State 1: Condenser outlet 
% (since this is a saturated liquid state, values can be evaluated from steam table)
h1 = XSteam('hL_p', P_cond/1e5); % (in kJ/kg), P_cond/1e5 to adjust units
s1 = XSteam('sL_p', P_cond/1e5);
v1 = XSteam('vL_p', P_cond/1e5);
% T1 = Saturation temp at condenser pressure
T1 = XSteam('Tsat_p', P_cond/1e5);   % in C

% State 2: Pump outlet (after isentropic compression)
% since pressure increases - small work input
W_pump_ideal = v1 * (P_boiler - P_cond) / 1000;  % in kJ/kg
h2 = h1 + W_pump_ideal; %  since W_pump_ideal = h2-h1
s2 = s1; % isentropic - entropy remains same
T2 = T1; % since temperature rise across pump is very small

% State 3: Boiler outlet (after constant P heat add)
% (this state is superheated steam)
h3 = XSteam('h_pT', P_boiler/1e5, T3-273); % gives enthalpy at the specified superheated state
s3 = XSteam('s_pT', P_boiler/1e5, T3-273); % gives entropy
% T3​ = specified turbine inlet temperature
T3_C = T3 - 273;   % Kelvin to C

% State 4: Turbine outlet (after isentropic expansion)
s4 = s3; % since isentropic so entropy remains same
h4 = XSteam('h_ps', P_cond/1e5, s4);
% since expansion is isentropic so temperature is obtained using pressure and entropy
T4 = XSteam('T_ps', P_cond/1e5, s3);   % in C


fprintf('Enthalpy at state 1 = %.3f\n',h1);
fprintf('Entropy at state 1 = %.3f\n',s1);
fprintf('Temperature at state 1 = %.3f\n\n',T1);

fprintf('Enthalpy at state 2 = %.3f\n',h2);
fprintf('Entropy at state 2 = %.3f\n',s2);
fprintf('Temperature at state 2 = %.3f\n\n',T2);

fprintf('Enthalpy at state 3 = %.3f\n',h3);
fprintf('Entropy at state 3 = %.3f\n',s3);
fprintf('Temperature at state 3 = %.3f\n\n',T3_C);

fprintf('Enthalpy at state 4 = %.3f\n',h4);
fprintf('Entropy at state 4 = %.3f\n',s4);
fprintf('Temperature at state 4 = %.3f\n\n',T4);

%% TASK(b): REAL RANKINE CYCLE 

eta_t = 0.85;   % Turbine efficiency
eta_p = 0.75;   % Pump efficiency

% Turbine (actual)
h4_actual = h3 - eta_t * (h3 - h4); % since eta_t = (h3-h4_actual)/(h3-h4)

% Pump (actual)
W_pump_actual = W_pump_ideal / eta_p;
h2_actual = h1 + W_pump_actual; % since W_pump_actual = h2_actual - h1



% In the real Rankine cycle, irreversibilities are introduced only in the 
% pump and turbine through isentropic efficiencies. Therefore, only the 
% outlet states of the pump and turbine (States 2 and 4) differ between 
% ideal and real cycles.

s2_actual = XSteam('s_ph', P_boiler/1e5, h2_actual);

s4_actual = XSteam('s_ph', P_cond/1e5, h4_actual);
T4_actual = XSteam('T_ph', P_cond/1e5, h4_actual);


fprintf('Actual Enthalpy at state 2 = %.3f\n',h2_actual);
fprintf('Actual Entropy at state 2 = %.3f\n',s2_actual);
fprintf('Actual Temperature at state 2 = %.3f\n\n',T2); % T2 actual will be ~T1 (same as T2 ideal)

fprintf('Actual Enthalpy at state 4 = %.3f\n',h4_actual);
fprintf('Actual Entropy at state 4 = %.3f\n',s4_actual);
fprintf('Actual Temperature at state 4 = %.3f\n\n',T4_actual); 

%% TASK(c): PERFORMANCE CALCULATIONS 

Wt = h3 - h4_actual;             % Turbine work (in kJ/kg)
Wp = h2_actual - h1;             % Pump work (kJ/kg)
Wnet = Wt - Wp;                  % Net work (kJ/kg)
Qin = h3 - h2_actual;            % Heat added (kJ/kg)

eta_th = Wnet / Qin;             % Thermal efficiency of the cycle
BWR = Wp / Wt;                   % Back-work ratio

fprintf('Turbine Work = %.3f\n', Wt);
fprintf('Pump Work = %.3f\n', Wp);
fprintf('Net Work output = %.3f\n', Wnet);
fprintf('Heat added in the boiler = %.3f\n', Qin);
fprintf('Thermal Efficiency = %.3f\n', eta_th);
fprintf('Back Work Ratio = %.4f\n\n', BWR);

% Physical significance of Back-Work Ratio : 
% Back-work ratio is the ratio of pump work and the total turbine work. It 
% tells us how much turbine work is consumed by the pump. In Rankine 
% cycles, BWR is usually small which indicates efficient power production.

%% TASK(d): T-s and h-s DIAGRAMS 

% State vectors (Ideal / Actual as needed)
T = [XSteam('Tsat_p', P_cond/1e5), XSteam('Tsat_p', P_boiler/1e5), T3-273, XSteam('T_ps', P_cond/1e5, s3) ];

s = [s1 s2 s3 s4];
h = [h1 h2 h3 h4];

% T-s plot : 
figure;
plot(s, T, '-o', 'LineWidth', 2);
xlabel('Entropy (kJ/kg·K)');
ylabel('Temperature (C)');
title('T–s Diagram of Rankine Cycle');
grid on;

% h-s plot :
figure;
plot(s, h, '-o', 'LineWidth', 2);
xlabel('Entropy (kJ/kg·K)');
ylabel('Enthalpy (kJ/kg)');
title('h–s Diagram of Rankine Cycle');
grid on;


%% TASK(e): PARAMETRIC STUDY 

P_boiler_range = linspace(8e6, 20e6, 20);
eta_array = zeros(size(P_boiler_range));

for i = 1:length(P_boiler_range)
    P = P_boiler_range(i);

    h3 = XSteam('h_pT', P/1e5, T3-273);
    s3 = XSteam('s_pT', P/1e5, T3-273);
    h4 = XSteam('h_ps', P_cond/1e5, s3);
    h4a = h3 - eta_t*(h3 - h4);

    Wt = h3 - h4a;
    Wp = v1*(P - P_cond)/1000/eta_p;
    Qin = h3 - (h1 + Wp);

    eta_array(i) = (Wt - Wp) / Qin;
end

% plotting efficiency vs Boiler pressure
figure;
plot(P_boiler_range/1e6, eta_array, 'LineWidth', 2);
xlabel('Boiler Pressure (MPa)');
ylabel('Thermal Efficiency');
title('Effect of Boiler Pressure on Efficiency');
grid on;

% Explanation for the plot : 
% Figure shows the variation of thermal efficiency with boiler pressure 
% for the Rankine cycle at constant condenser pressure and turbine inlet 
% temperature. As boiler pressure increases, thermal efficiency increases 
% due to a rise in the average temperature of heat addition and an 
% increase in turbine work output. However, at higher pressures, the rate 
% of increase in efficiency diminishes because pump work increases and the 
% incremental rise in mean heat addition temperature becomes smaller. 
% This results in a gradually flattening efficiency curve.