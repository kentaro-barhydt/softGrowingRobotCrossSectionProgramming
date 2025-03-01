% Soft Growing Robotic Sling Cross-Section Programming
% Geometry Model - Numerical 
% 
% This script runs an alternative version of the model presented in the
% paper, "Mechanically Programming the Cross-Sectional Shape of Soft 
% Growing Robotic Structures for Patient Transfer", Section III.C. Unlike 
% the analytical model presented in the paper, which calculates the 
% fabrication parameters required to meet a set of design specifications, 
% this model calculates the cross-sectional geometry given only a set of 
% fabrication parameters. This model utilizes a numerical solver (fsolve) 
% to solve for transcendental equations that arise from this derivation. 

clear all; close all; clc;


%% INPUTS AND VARIABLES

% input design parameter values for magic sheet prototype 
s_c_design = 5.5
L_design = 3
s_s_design = 12



%% SETUP

% variables: 
syms s_c L th_c s_s th_s


% variable constraints:
th_c_range_min = 0;
th_c_range_max = 2*pi;
th_s_range_min = 0;
th_s_range_max = 2*pi;



%% AREA OF CENTER CHANNEL

% Define function for area of center channel A_c

r_c = s_c / th_c;

w_c = 2*r_c*sin(th_c/2);

A_c = 2*( (th_c)/(2*pi)*pi*r_c^2 - (1/2)*r_c^2*sin(th_c) ) ...
    + w_c * L;
A_c_p = diff(A_c,th_c);
A_c_p_fun = matlabFunction(A_c_p, 'Vars', [s_c, L, th_c]);


% Substitute in design parameter values for L and s_c:

A_c_p_fun_design = @(th_c) A_c_p_fun(s_c_design, L_design, th_c);


% Solve for th_c value that maximizes A_c

th_c_sol_local = [];
th_c_init_int = 0.1;
options = optimset('TolFun',1e-15,'Display','off');
for th_c_init = th_c_range_min:th_c_init_int:th_c_range_max
    if th_c_init == 0
        th_c_init = th_c_init_int * 0.001;
    end
    sol = fsolve(A_c_p_fun_design,th_c_init,options);
    th_c_sol_local = [th_c_sol_local; sol];
end
th_c_sol_local = unique(round(th_c_sol_local,4));

A_c_fun = matlabFunction(A_c, 'Vars', [s_c, L, th_c]);
A_c_fun_design = @(th_c) A_c_fun(s_c_design, L_design, th_c);
A_c_design_sol = A_c_fun_design(th_c_sol_local);

[A_c_design_sol_global, A_c_design_sol_global_i] = max(A_c_design_sol);
th_c_sol_global = th_c_sol_local(A_c_design_sol_global_i);



%% GEOMETRY OF CENTER CHANNEL 

% Calculate remaining geometric parameters for center channel 

r_c_fun = matlabFunction(r_c);
r_c_design = r_c_fun(s_c_design, th_c_sol_global);

w_c_fun = matlabFunction(w_c);
w_c_design = w_c_fun(s_c_design, th_c_sol_global);



%% GEOMETRY OF SIDE CHANNELS

% Solve for th_s value that satisfies straight segment length equality
% condition L_eq=0

r_s = s_s / (2*pi - th_s);

L_eq = L - ( 2 * r_s * sin(th_s/2) ); % Length equality condition
L_eq_fun = matlabFunction(L_eq, 'Vars', [s_s, th_s, L]);
L_eq_fun_design = @(th_s) L_eq_fun(s_s_design, th_s, L_design);

th_s_sol_local = [];
th_s_init_int = 0.1;
options = optimset('TolFun',1e-15,'Display','off');
for th_s_init = th_s_range_min:th_s_init_int:th_s_range_max
    if th_s_init == 0
        th_s_init = th_s_init_int * 0.001;
    end
    sol = fsolve(L_eq_fun_design,th_s_init,options);
    th_s_sol_local = [th_s_sol_local; sol];
end
th_s_sol_global = [];
th_s_sol_local = unique(round(th_s_sol_local,4));
for i = 1:length(th_s_sol_local)
    if abs(L_eq_fun_design(th_s_sol_local(i))) < 1e-4
        th_s_sol_global = [th_s_sol_global; th_s_sol_local(i)];
    end
end
if isempty(th_s_sol_global)
    error("ERROR: NO SOLUTION FOR th_s_sol_global")
elseif length(th_s_sol_global) > 1
    [L_eq_fun_design_min,L_eq_fun_design_min_i] = ...
        min(abs(th_s_sol_global), [], "all");
    th_s_sol_global = th_s_sol_global(L_eq_fun_design_min_i);
end


% Calculate radius of side channel 

r_s_fun = matlabFunction(r_s, 'Vars', [s_s, th_s]);
r_s_design = r_s_fun(s_s_design, th_s_sol_global);



%% PLOT FULL GEOMETRY 

figure('name','cross-section geometry')
hold on; axis equal
plot([-w_c_design/2, -w_c_design/2], [-L_design/2, L_design/2],'k')
plot([w_c_design/2, w_c_design/2], [-L_design/2, L_design/2],'k')

p_c_upper = [0, L_design/2 - r_c_design*cos(th_c_sol_global/2)];
p_c_lower = [0, r_c_design*cos(th_c_sol_global/2) - L_design/2];
p_s_left = [-w_c_design/2 - r_s_design*sin((pi-th_s_sol_global)/2), 0];
p_s_right = [w_c_design/2 + r_s_design*sin((pi-th_s_sol_global)/2), 0];
plot(p_c_upper(1),p_c_upper(2),'r.')
plot(p_c_lower(1),p_c_lower(2),'g.')
plot(p_s_left(1),p_s_left(2),'c.')
plot(p_s_right(1),p_s_right(2),'m.')

th_c_plot_upper = linspace(pi/2-th_c_sol_global/2, ...
    pi/2+th_c_sol_global/2, 100);
s_c_plot_upper_x = p_c_upper(1) + r_c_design.*cos(th_c_plot_upper);
s_c_plot_upper_y = p_c_upper(2) + r_c_design.*sin(th_c_plot_upper);
plot(s_c_plot_upper_x, s_c_plot_upper_y,'r')

th_c_plot_lower = linspace(3*pi/2-th_c_sol_global/2, ...
    3*pi/2+th_c_sol_global/2, 100);
s_c_plot_lower_x = p_c_lower(1) + r_c_design.*cos(th_c_plot_lower);
s_c_plot_lower_y = p_c_lower(2) + r_c_design.*sin(th_c_plot_lower);
plot(s_c_plot_lower_x, s_c_plot_lower_y,'g')

th_s_plot_left = linspace(th_s_sol_global/2, ...
    2*pi - th_s_sol_global/2, 100);
s_s_plot_left_x = p_s_left(1) + r_s_design.*cos(th_s_plot_left);
s_s_plot_left_y = p_s_left(2) + r_s_design.*sin(th_s_plot_left);
plot(s_s_plot_left_x, s_s_plot_left_y,'c')

th_s_plot_right = linspace(-pi + th_s_sol_global/2, ...
    pi - th_s_sol_global/2, 100);
s_s_plot_right_x = p_s_right(1) + r_s_design.*cos(th_s_plot_right);
s_s_plot_right_y = p_s_right(2) + r_s_design.*sin(th_s_plot_right);
plot(s_s_plot_right_x, s_s_plot_right_y,'m')

plot(0,0,'k.')

r_max = max([r_s_design r_c_design]);
axis(1.2 .* ...
    [p_s_left(1)-r_s_design ...
    p_s_right(1)+r_s_design ...
    -r_max ...
    r_max])




