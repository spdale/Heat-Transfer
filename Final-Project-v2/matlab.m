height = 200;
width = 500;

num_time_steps = 20000;

cylinder_diameter = 50;
cylinder_radius = cylinder_diameter / 2;
cylinder_center_x = height / 2;
cylinder_center_y = 100;

error_limit = 0.01;                  % 1% maximum change for convergence

U_inf = 2;                           % m/s uniform inflow
F = 1.9;                             % over-relaxation factor
free_lid = U_inf * (height / 2);     % free-lid streamfunction constant

Re_D = 200;                          % Given Reynolds number

T_surface = 400;                     % K
T_boundary = 300;                    % K
T_init = min(T_surface, T_boundary); % Bulk fluid initial temp

% Constants picked for air around room temp
alpha = 22.07 * 10^(-6);    % m^2/s     Thermal Diffusivity at 300K
% alpha = 0.1463 * 10**(-6)           % water at 300K
nu = 1.48 * 10^(-5);        % m^2/s     Kinematic Viscosity at 300K
% nu = 8.56 * 10**(-7)                % water at 300K


h_1 = (10 - 1) * nu / U_inf;
h_2 = (10 - 1) * alpha / U_inf;
h = min(h_1, h_2)       % grid spacing

dt = (h / U_inf) / 2


omega = zeros(width, height);
psi = zeros(width, height);
temps = zeros(width, height);


solid_rows = int16.empty;
solid_cols = int16.empty;
for i = 1:width
    for j = 1:height
        dist = sqrt((i - cylinder_center_x)^2 + (j - cylinder_center_y)^2);
        if (dist <= cylinder_radius)
            omega(i, j) = 1
        end
    end
end

solid_indices = 

omega(find(omega == 1)) = 2


