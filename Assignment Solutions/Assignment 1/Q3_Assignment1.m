clear all;

% Given 
P1 = 1;          
P2 = 20;         
n  = 1.28;
T1 = 300;        % in K (reference inlet temperature)
cp = 1.05;       % kJ/kg-K
R  = 0.287;      % kJ/kg-K

% Pressure discretization
P = linspace(P1, P2, 1000);   

% Temperature along polytropic path
% For ideal gas: T * P^{(1-n)/n} = constant
T = T1 * (P/P1).^((n-1)/n);

% Real-gas compressibility factor
Z = 1 + 0.0008*P - 120./T;

% Numerical derivatives
dT = gradient(T, P);
dP = gradient(P);

% Real-gas entropy change (numerical integration)
ds_real = cp * (dT ./ T) - R * Z .* (dP ./ P);
delta_s_real = trapz(P, ds_real);

% Ideal-gas entropy change
ds_ideal = cp * (dT ./ T) - R * (dP ./ P);
delta_s_ideal = trapz(P, ds_ideal);

% Percent deviation
percent_deviation = (delta_s_real - delta_s_ideal)/delta_s_ideal * 100;

% results
fprintf('Real-gas entropy change  = %.4f kJ/kg-K\n', delta_s_real);
fprintf('Ideal-gas entropy change = %.4f kJ/kg-K\n', delta_s_ideal);
fprintf('Percent deviation        = %.2f %%\n', percent_deviation);

% Plotting results
figure;
plot(P, T, 'LineWidth', 2)
xlabel('Pressure (bar)')
ylabel('Temperature (K)')
title('Temperature Variation Along Polytropic Compression')
grid on
