clear all;

tspan = [0 4000];   % time span

T0 = 300;   % in K (assumed initial temperature)
u0 = 450 + 1.1*T0 + 0.0012*T0^2; 

% ODE solution for internal energy
[t, u] = ode45(@energy_balance, tspan, u0);

% Inverting u(T) to obtain temperature history
T = zeros(size(u));

for i = 1:length(u)
    T(i) = fsolve(@(Tvar) 450 + 1.1*Tvar + 0.0012*Tvar^2 - u(i), T0);
end

q_ext = 5000 * exp(-0.002 * t);
r_T = 1500 * (1 - exp(-0.01 * T));

% Finding time when reaction heat exceeds external heat
idx = find(r_T > q_ext, 1, 'first');
t_cross = t(idx);

% Plotting figures 
figure;
subplot(2,1,1)
plot(t, u, 'LineWidth', 2)
xlabel('Time (s)')
ylabel('Internal Energy (kJ/kg)')
title('Internal Energy Evolution')
grid on

subplot(2,1,2)
plot(t, T, 'r', 'LineWidth', 2)
xlabel('Time (s)')
ylabel('Temperature (K)')
title('Recovered Temperature History')
grid on

figure;
plot(t, q_ext, 'b', t, r_T, 'r--', 'LineWidth', 2)
xlabel('Time (s)')
ylabel('Heat Rate (kJ/s)')
legend('External Heating', 'Reaction Heat')
title('Comparison of Heat Contributions')
grid on

% result
fprintf('Reaction heat exceeds external heating at t = %.2f s\n', t_cross);


% Energy balance ODE function
function dudt = energy_balance(t, u)
    % estimating temperature from u using quadratic inversion
    % u = 450 + 1.1T + 0.0012T^2
    coeffs = [0.0012 1.1 (450 - u)];
    T = max(roots(coeffs));   % physical root
    q_ext = 5000 * exp(-0.002 * t);
    r_T = 1500 * (1 - exp(-0.01 * T));
    dudt = q_ext + r_T;
end
