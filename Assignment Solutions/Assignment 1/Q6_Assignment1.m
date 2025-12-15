clear all;

% defining a time grid
t_start = 0;          
t_end   = 8000;       
N = 5000;             
t = linspace(t_start, t_end, N);

% Boundary temperatures(given) - 
Th = 900 - 300*exp(-0.0008*t);
Tc = 300 + 40*sin(0.002*t);

% Efficiency - 
eta = 1 - Tc./Th;

% Heat input - 
Qin = 20000*(1 + 0.3*sin(0.003*t));

% Instantaneous power output - 
P = eta .* Qin;

% Total work output - 
W_total = trapz(t, P);   % since P is in kW and t is in seconds

% rate of entropy generation 
Sgen = Qin./Th - Qin./Tc;

% Result : 
fprintf('Total work output = %.2f MJ\n', W_total/1000); % in MJ

 
% Plotting the results
figure;

subplot(3,1,3)
plot(t, Sgen, 'm', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('{S}_{gen} (kW/K)');
title('Entropy Generation Rate');
grid on;

subplot(3,1,1)
plot(t, Th, 'r', t, Tc, 'b', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Temperature (K)');
legend('T_h(t)', 'T_c(t)');
title('Reservoir Temperatures');
grid on;

subplot(3,1,2)
plot(t, eta, 'k', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Efficiency');
title('Time-Varying Efficiency');
grid on;


