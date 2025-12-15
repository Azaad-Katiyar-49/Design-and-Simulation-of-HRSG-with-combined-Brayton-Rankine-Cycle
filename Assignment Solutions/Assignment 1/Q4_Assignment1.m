clear all;

% Given 
T1 = 310;        % in K
T2 = 670;        % in K
Tb = 300;        % in K

h = @(T) 300 + 2.5*T + 0.0007*T.^2;
s = @(T) 2.0*log(T) + 0.001*T;

% Enthalpy and entropy change
dh = h(T2) - h(T1);        
ds = s(T2) - s(T1);        

% Heat input and mass flow ranges
Qdot = linspace(20, 100, 400);     % kW
mdot = linspace(0.01, 2, 400);     % kg/s

% Creating a 2D grid using meshgrid
[Qdot_grid, mdot_grid] = meshgrid(Qdot, mdot);

% Entropy generation rate (using entropy balance relation)
Sgen = mdot_grid .* ds - Qdot_grid ./ Tb;

% Feasible region 
feasible = (Sgen >= 0);

% Plotting feasible operating region
figure;
contourf(Qdot_grid, mdot_grid, feasible, 1);
colormap([1 1 1; 0.2 0.6 1]);
colorbar;
xlabel('Heat Input Rate(kW)');
ylabel('Mass Flow Rate (kg/s)');
title('Feasible Operating Region (Sgen)');
grid on;

% Explanation :

% why a feasible region exists >
% larger Qdot means more entropy rejected 
% therefore, to satisfy second law, mdot(s2-s1)>=Qdot/Tb
% => higher heat input requires higher mass flow

% connection with internal energy rise >
% higher Qdot demands higher mdot
% entropy balance further restricts this relationship
% hence, only a subset of energy balanced states are thermodynamically
% stable