clear all;

% symbolic integration (analytical solution)
syms T                          % symbolic variable
cv = 700 + 0.35*T - 2e-4*T^2;
du_sym = int(cv,T,320,820);     % exact integration
du_sym = double(du_sym)


% numerical integration using integral()
cv_num = @(T) 700 + 0.35*T - 2e-4*T.^2;
cv_num_int = integral(cv_num,320,820)

% numerical integration using trapz()
N = 1000;
T_grid = linspace(320,820,N);
cv_vals = cv_num(T_grid);
cv_num_trapz = trapz(T_grid,cv_vals)

% error analysis and minimum grid resolution
N_vals = 50:50:5000;
error = zeros(size(N_vals));

for i = 1:length(N_vals)
    Tgrid = linspace(320,820,N_vals(i));
    du_num = trapz(Tgrid,cv_num(Tgrid));
    error(i) = abs((du_num-du_sym)/du_sym)*100;
end

N_required = N_vals(find(error<0.1,1)) % minimum N that satisfies the accuracy requirement

