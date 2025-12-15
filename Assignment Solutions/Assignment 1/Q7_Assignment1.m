clear all;

% Given 
PA = 1;            
PB = 10;           
TA = 300;          
m  = 1.25;         
R  = 0.287;        
u = @(T) 500 + 0.8*T + 1.5e-3*T.^2;

% final temperature TB using polytropic relation
% for an ideal gas: T*P^{(1-m)/m} = constant
polytropic_eq = @(TB) TB*(PB^((1-m)/m)) - TA*(PA^((1-m)/m));
TB = fsolve(polytropic_eq, 700);   % initial guess = 700 K


% Change in internal energy
delta_u = u(TB) - u(TA);   


% Numerical work calculation along the polytropic path :

% Expressing pressure as a function of temperature
P_T = @(T) PA*(T/TA).^(m/(m-1));   % in bar
% Converting pressure to kPa (since R is in kJ/kg-K)
P_T_kPa = @(T) 100*P_T(T);
% Specific volume from ideal gas law
v_T = @(T) R*T ./ P_T_kPa(T);   % in m^3/kg
% dv/dT (numerical derivative)
dv_dT = @(T) gradient(v_T(T), T);

% Work integral: W = integral of (P*dv)
integrand = @(T) P_T_kPa(T).*dv_dT(T);

T_path = linspace(TA, TB, 2000);
W_num = trapz(T_path, integrand(T_path));  


% Analytical polytropic work (ideal gas)
W_analytical = (R*(TB - TA))/(1 - m);   % in kJ/kg


% Heat transfer
Q = delta_u + W_num;   % from first Law


% Results
fprintf('Final temperature TB        = %.2f K\n', TB);
fprintf('Change in internal energy   = %.3f kJ/kg\n', delta_u);
fprintf('Numerical work              = %.3f kJ/kg\n', W_num);
fprintf('Analytical polytropic work  = %.3f kJ/kg\n', W_analytical);
fprintf('Heat transfer               = %.3f kJ/kg\n', Q);
