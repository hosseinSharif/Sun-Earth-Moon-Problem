% In The Name of Allah
% In Persuit of The Moon Project
% Hossein Sharif 1401/11/28 14:00
% Location: Islamic Republic of Iran - Tehran

%% Sun-Earth Problem:
%    Solving The Differential Equation of Earth Displacement And Velocity Vector
%    With Runge-Kutta 4th Order Method
%%

%% Note!
%    In genral, we can say that any other system of two masses with mass of the one
%    of them very greater than the other, will go this way! %%

% using Plots

close all
clear all
clc

% Constants: (SI Units)
M_S = 1.989e30; % The Mass of The Sun (Kg)
M_E = 5.972e24; % The Mass of The Earth (Kg)
M_M = 7.34767309e22; % The Mass of The Moon (Kg)

G = 6.6743e-11; % Gravitational Constant (m^3kg^-1s^-2)
AU = 149597870700; % Mean Distance Between The Sun And The Earth (m)

SEMD = AU; % Sun-Earth Mean Distance = 149597870.700 Km
EMMD = SEMD / 391.1055443137255; % Earth-Moon Mean Distance = 382,500 Km
M_i = 5.15; % The Moon moves in an approximately elliptic orbit inclined
% at about five degrees to the plane of the ecliptic. (In Degrees)

%% Linear Velocity of The Earth In Orbit Around The Sun:
% This equation can be proved with the help of polar
% coordinates and two assumptions:
% 1- r is constant in orbital motion (r = AU)
% 2- Theta_Dot is constant in orbital motuon (Theta_Dot = cte.)
V_E0 = (G * M_S / SEMD)^0.5; % Earth Initial Velocity
V_M0 = (G * M_E / EMMD)^0.5 + V_E0; % Moon Initial Velocity

% Initialize System:
Day_Hours = 23 + (56 / 60) + (4 / 3600); % Hours of A Day
Year_Days = 365.2425; % Days of

Year_MIN = zeros(1, length(Year));
Year_MIN(1, :) = (Year_Days * Day_Hours * 60) .* Year(1, :);

Year_Seconds = Year_Days * Day_Hours * 3600;

dD = 0.0006944445; % Simulation Step In Days 0.0006944445 = 1 minute
Stop_Year = 2; % The year that the simulation will stop.

dY = dD / 365.2425; % Simulation  Step In Years
dt = dY * Year_Seconds; % In Seconds

Year = zeros(1, length(0:dY:Stop_Year));
Year(1, :) = 0:dY:Stop_Year; % In Years

% JPL Horizons Initial Velocities of Moon And Earth (1401-01-01: Solar New Year)
xM0 = -150614032638.7242; % -SEMD - (EMMD * cosd(M_i));
yM0 = 1002563524.411460; % 0;
zM0 = 41625742.77642026; % EMMD * sind(M_i);
VxM0 = -112.8832840254308; % 0;
VyM0 = -30808.89757527497; % -V_M0;
VzM0 = -84.00881299342267; % 0;

xE0 = -150289457012.5996; % -SEMD;
yE0 = 1186655852.197647; % 0;
zE0 = 28135493.33878030; % 0;
VxE0 = -652.1856667917787; % 0;
VyE0 = -29906.06102741717; % -V_E0;
VzE0 = 2.776852076703307; % 0;

x = zeros(12, length(Year));
x(:, 1) = [xM0 yM0 zM0 xE0 yE0 zE0 VxM0 VyM0 VzM0 VxE0 VyE0 VzE0]';

% Integrating With RK4:
for i = 1:(length(Year)-1)
    k1 = dt .* f(x(:, i));
    k2 = dt .* f(x(:, i) + k1 ./ 2);
    k3 = dt .* f(x(:, i) + k2 ./ 2);
    k4 = dt .* f(x(:, i) + k3);
    
    x(:, i+1) = x(:, i) + (k1 + 2 * k2 + 2 * k3 + k4) ./ 6;
    
end

%%
r_SM = x(1:3, :);
v_SM = x(7:9, :);

r_SE = x(4:6, :);
v_SE = x(10:12, :);

r_EM = r_SM - r_SE;
v_EM = v_SM - v_SE;

theta = zeros(1, length(Year));
PML = zeros(1, length(Year));

for i = 1:(length(Year))
    theta(1, i) = rad2deg(acos(dot(r_SM(1:3, i), r_EM(1:3, i)) / (norm(r_SM(1:3, i)) * norm(r_EM(1:3, i)))));
    PML(1, i) = 0.5 * (1 - cosd(theta(1, i)));
end

%% Plots:
figure(1);
hold on;
grid on;
axis equal;
plot(r_SE(1, :), r_SE(2, :));
plot(r_SM(1, :), r_SM(2, :),'r');

figure(2);
hold on;
grid on;
plot(Year_MIN(1, :), 100 * PML(1, :));