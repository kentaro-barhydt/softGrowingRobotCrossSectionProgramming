% Soft Growing Robotic Sling Cross-Section Programming
% Geometry Model - Analytical
% 
% This script runs the analytical model presented in the paper, 
% "Mechanically Programming the Cross-Sectional Shape of Soft Growing 
% Robotic Structures for Patient Transfer", Section III.C. 

clear all; close all; clc;


%% INPUTS

% Input geometric parameter values:
H_c = 3.5
H_s = 4.0
w = 9.0


%% CHECK INPUTS 

% Check if input geometric parameter values are admissible
w_min = 2 * (0.5*H_s);
w_max = 2*H_s + H_c;
if w < w_min
    error('INVALID INPUT: w < w_min')
end
if w > w_max
    error('INVALID INPUT: w > w_max')
end
if H_c <= 0
    error('INVALID INPUT: H_c <=0')
end
if H_s <= 0
    error('INVALID INPUT: H_s <=0')
end
if w <= 0
    error('INVALID INPUT: w <=0')
end


%% ANALYTIC SOLUTIONS

s_c = H_c * asin( ( w^2 + H_c^2 - 2*w*H_s ) ...
    / ( 2*H_c*(w - H_s) ) )                                         % Eq.19

L = abs( 2*(w-H_s) )^-1 ...
    * sqrt( 4*H_s*( H_c^2*H_s - H_c^2*w - w^2*H_s + w^3 ) ...
    - ( w^2 - H_c^2 )^2 )                                           % Eq.20

s_s_p = H_s * asin( abs( 2*H_s*(w-H_s) )^-1 ...
    * sqrt( 4*H_s*( H_c^2*H_s - H_c^2*w - w^2*H_s + w^3 ) ...
    - ( w^2 - H_c^2 )^2 ) );                                        % Eq.22

if w-H_c*sin(s_c/H_c) > H_s
    s_s = pi*H_s - s_s_p
    disp('s_s_p <= 0.5*pi*H_s')                                     % Eq.21
elseif w-H_c*sin(s_c/H_c) <= H_s
    s_s = s_s_p
    disp('s_s_p > 0.5*pi*H_s')                                      % Eq.21
else
    error('INVALID INEQUALITY s_s_p')
end


%% CALCULATE REMAINING TERMS TO DEFINE FULL GEOMETRY

w_c = H_c*sin(s_c/H_c)                                              % Eq.15
r_c = 0.5*H_c                                                       % Eq.12
th_c = 2*s_c/H_c                                                    % Eq.13
r_s = 0.5*H_s                                                       % Eq.7
th_s = 2*pi - 2*s_s/H_s                                             % Eq.8


%% PLOT GEOMETRY

figure('name','cross-section geometry')
plot(0,0,'k.')
hold on; axis equal
H_max = max([H_s H_c]);
axis(1.2 .*[-w/2 w/2 -H_max/2 H_max/2])

plot([-w_c/2, -w_c/2], [-L/2, L/2],'k')
plot([w_c/2, w_c/2], [-L/2, L/2],'k')

p_c_upper = [0, L/2 - r_c*cos(th_c/2)]
p_c_lower = [0, r_c*cos(th_c/2) - L/2]
p_s_left = [-w_c/2 - r_s*sin((pi-th_s)/2), 0]
p_s_right = [w_c/2 + r_s*sin((pi-th_s)/2), 0]
plot(p_c_upper(1),p_c_upper(2),'k.')
plot(p_c_lower(1),p_c_lower(2),'k.')
plot(p_s_left(1),p_s_left(2),'k.')
plot(p_s_right(1),p_s_right(2),'k.')

th_c_plot_upper = linspace(pi/2-th_c/2, pi/2+th_c/2, 100);
s_c_plot_upper_x = p_c_upper(1) + r_c.*cos(th_c_plot_upper);
s_c_plot_upper_y = p_c_upper(2) + r_c.*sin(th_c_plot_upper);
plot(s_c_plot_upper_x, s_c_plot_upper_y,'k')

th_c_plot_lower = linspace(3*pi/2-th_c/2, 3*pi/2+th_c/2, 100);
s_c_plot_lower_x = p_c_lower(1) + r_c.*cos(th_c_plot_lower);
s_c_plot_lower_y = p_c_lower(2) + r_c.*sin(th_c_plot_lower);
plot(s_c_plot_lower_x, s_c_plot_lower_y,'k')

th_s_plot_left = linspace(th_s/2, 2*pi - th_s/2, 100);
s_s_plot_left_x = p_s_left(1) + r_s.*cos(th_s_plot_left);
s_s_plot_left_y = p_s_left(2) + r_s.*sin(th_s_plot_left);
plot(s_s_plot_left_x, s_s_plot_left_y,'k')

th_s_plot_right = linspace(-pi + th_s/2, pi - th_s/2, 100);
s_s_plot_right_x = p_s_right(1) + r_s.*cos(th_s_plot_right);
s_s_plot_right_y = p_s_right(2) + r_s.*sin(th_s_plot_right);
plot(s_s_plot_right_x, s_s_plot_right_y,'k')


%% CALCULATE ADDITIONAL RESULTING GEOMETRIC PARAMETERS 

A_cs_s = pi*r_s^2 - ( (pi*r_s^2)*(th_s/(2*pi)) - 0.5*r_s^2*sin(th_s) )
A_cs_c = 2*( (pi*r_c^2)*(th_c/(2*pi)) - 0.5*r_c^2*sin(th_c) ) + L*w_c
A_cs = 2*A_cs_s + A_cs_c

