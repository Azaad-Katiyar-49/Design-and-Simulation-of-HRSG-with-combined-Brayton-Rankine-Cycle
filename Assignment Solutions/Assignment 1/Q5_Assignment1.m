clear all;

% Given 
T1 = 350;        
T2 = 900;        
T0 = 298; % (dead-state temperature)
cp = @(T) 1200 + 0.4*T - 1.2e-4*T.^2;

% Entropy change (numerical integration)
integrand = @(T) cp(T)./T;
delta_s = integral(integrand, T1, T2)   

% Irreversibility levels
irreversibility = linspace(0, 0.10, 100);  % 0% to 10%

% Exergy destruction
% Xdest = T0*S_gen_dot
% S_gen_dot = (irreversibility)*delta_s
Xdest = T0 .* irreversibility .* delta_s;  

% Plotting exergy destruction vs irreversibility levels
figure;
% x axis -> irreversibility percentage and y axis -> Xdest in kJ/kg
plot(irreversibility*100, Xdest/1000, 'LineWidth', 2); 
xlabel('Irreversibility Level (%)');
ylabel('Exergy Destruction (kJ/kg)');
title('Exergy Destruction vs Irreversibility');
grid on;
