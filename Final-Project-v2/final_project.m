close all

height = 200;
width = 500;

num_time_steps = 5000;

cylinder_diameter = 50;
cylinder_radius = cylinder_diameter / 2;
cylinder_center_x = height / 2;
cylinder_center_y = height / 2;

error_limit = 0.01;                  % 1% maximum change for convergence

U_inf = 2;                        % m/s      uniform inflow
alpha = 22.07 * 10^(-6);             % m^2/s    Thermal Diffusivity at 300K
nu = 1.48 * 10^(-5);                 % m^2/s    Kinematic Viscosity at 300K
F = 1.9;                             %          over-relaxation factor
free_lid = U_inf * (height / 2);     %          free-lid streamfunction constant

Re_D = 200;                          % Given Reynolds number

T_surface = 400;                     % K
T_boundary = 300;                    % K
T_init = min(T_surface, T_boundary); % Bulk fluid initial temp

h_1 = (10 - 1) * nu / U_inf;
h_2 = (10 - 1) * alpha / U_inf;
h = min(h_1, h_2)       % grid spacing

U_max = 5*U_inf;

dt = (h / U_max) / 2



omega = zeros(width, height);
psi = zeros(width, height);
temps = zeros(width, height);

u = zeros(width, height);
v = zeros(width, height);


solid_points = zeros(width, height);
for i = 1:width
    for j = 1:height
        dist = sqrt((i - cylinder_center_x)^2 + (j - cylinder_center_y)^2);
        if dist <= cylinder_radius
            solid_points(i, j) = 1;
            temps(i, j) = T_surface;
        else
            temps(i, j) = T_boundary;
%             psi = U_inf * j - free_lid;
            psi(i, j) = (U_inf * j - free_lid) * h;
        end
        
        
        
    end
end

u(1,:) = U_inf; 


error_flag = true;
while error_flag
    psi_old = psi;
    
    for i = 2:(width - 1)
        for j = 2:(height - 1)
            if ~solid_points(i,j)
                psi(i, j) = psi(i, j) + (F / 4) * (psi(i - 1, j) + psi(i + 1, j) + psi(i, j - 1) + psi(i, j + 1) - 4 * psi(i, j));
            end
        end
    end
    
%     psi(0, :) = psi(3, :);
    
    error_array = abs(psi - psi_old) ./ psi_old;
    error_array(isnan(error_array)) = 0;
    
    error_term = max(error_array);
    
    if (error_term <= error_limit)
        error_flag = false;
    end
    
end



for time_step = 1:num_time_steps
    
    disp("Time step " + time_step + " of " + num_time_steps)
    
    psi_old = psi;
    omega_old = omega;
    temps_old = temps;
    
    for i = 2:(width - 1)
        for j = 2:(height - 1)
            if solid_points(i, j)
                omega(i, j) = (-2 / (h * h)) * (psi(i - 1, j) + psi(i + 1, j) + psi(i, j - 1) + psi(i, j + 1));
            end
        end
    end
    
    for i = 2:(width - 1)
        for j = 2:(height - 1)
            u(i, j) = (psi(i, j + 1) - psi(i, j - 1)) / (2 * h);
            v(i, j) = (psi(i - 1, j) - psi(i + 1, j)) / (2 * h);
        end
    end
    
    u(:,1)=U_inf;
    u(:,height)=U_inf;
    
    for i = 2:(width - 1)
        for j = 2:(height - 1)
            if ~solid_points(i, j)
                laplacian_vorticity = (omega_old(i - 1, j) + omega_old(i + 1, j) + omega_old(i, j - 1) + omega_old(i, j + 1) - 4 * omega_old(i, j)) / (h * h);

                delta_u_omega = 0;
                if (u(i, j) < 0)
                    delta_u_omega = u(i + 1, j) * omega_old(i + 1, j) - u(i, j) * omega_old(i, j);
                elseif (u(i, j) > 0)
                    delta_u_omega = u(i, j) * omega_old(i, j) - u(i - 1, j) * omega_old(i - 1, j);
                end

                delta_v_omega = 0;
                if (v(i, j) < 0)
                    delta_v_omega = v(i, j + 1) * omega_old(i, j + 1) - v(i, j) * omega_old(i, j);
                elseif (v(i, j) > 0)
                    delta_v_omega = v(i, j) * omega_old(i, j) - v(i, j - 1) * omega_old(i, j - 1);
                end

                omega(i, j) = omega_old(i, j) + dt * (-delta_u_omega / h - delta_v_omega / h + nu * laplacian_vorticity);
%                 if laplacian_vorticity > 5
%                     disp(laplacian_vorticity + " " + time_step)
%                 end
            end
        end
    end
    

    psi(width, :) = 2 * psi(width - 1, :) - psi(width - 2, :);
    omega(width, :) = omega(width - 1, :);
    
    
    
    error_flag = true;
    while error_flag
        psi_old = psi;

        for i = 2:(width - 1)
            for j = 2:(height - 1)
                if ~solid_points(i,j)
                    psi(i, j) = psi_old(i, j) + (F / 4) * (psi(i - 1, j) + psi(i + 1, j) + psi(i, j - 1) + psi(i, j + 1) + 4 * h * h * omega(i, j) - 4 * psi(i, j));
                end
            end
        end

    %     psi(0, :) = psi(3, :);

        error_array = abs(psi - psi_old) ./ psi_old;
        error_array(isnan(error_array)) = 0;

        error_term = max(error_array);

        if (error_term <= error_limit)
            error_flag = false;
        end
    end
    
    temps_old = temps;
    for i = 2:(width -1)
        for j = 2:(height - 1)
            if ~solid_points(i, j)
                laplacian_temps = (temps_old(i - 1, j) + temps_old(i + 1, j) + temps_old(i, j - 1) + temps_old(i, j + 1) - 4 * temps_old(i, j)) / (h * h);
        
                u_delta_T = 0;
                if u(i, j) < 0
                    u_delta_T = u(i, j) * (temps_old(i + 1, j) - temps_old(i, j));
                elseif u(i, j) > 0
                    u_delta_T = u(i, j) * (temps_old(i, j) - temps_old(i - 1, j));
                end
                
                v_delta_T = 0;
                if v(i, j) < 0
                    v_delta_T = v(i, j) * (temps_old(i, j + 1) - temps_old(i, j));
                elseif v(i, j) > 0
                    v_delta_T = v(i, j) * (temps_old(i, j) - temps_old(i, j - 1));
                end
%                 
                temps(i, j) = temps_old(i, j) + dt * (-u_delta_T / h - v_delta_T / h + alpha * laplacian_temps);
                
            end
        end
    end
    
    
end



















figure(1)
hold on
plot_data = flipud(rot90(psi));
s = pcolor(plot_data);
daspect([1 1 1]);
colormap(gray);
set(s, 'EdgeColor', 'none');
colorbar
contour(plot_data, 32, 'black');
title("Streamfunction");
hold off





figure(2)
hold on
plot_data = flipud(rot90(omega));
s = pcolor(plot_data);
daspect([1 1 1]);
colormap(gray);
set(s, 'EdgeColor', 'none');
colorbar
title("Vorticity");
hold off



figure(3)
hold on
plot_data = flipud(rot90(temps));
s = pcolor(plot_data);
daspect([1 1 1]);
colormap(jet);
set(s, 'EdgeColor', 'none');
colorbar
title("Temperature");
hold off
