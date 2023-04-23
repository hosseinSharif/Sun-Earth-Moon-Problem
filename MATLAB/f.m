% State-Space Matrix Calculation Function:
function f_x = f(x_var)
% Constants: (SI Units)
M_S = 1.989e30; % The Mass of The Sun (Kg)
M_E = 5.972e24; % The Mass of The Earth (Kg)
M_M = 7.34767309e22; % The Mass of The Moon (Kg)

G = 6.6743e-11; % Gravitational Constant (m^3kg^-1s^-2)

D_A1_A2 = (((x_var(1) - x_var(4))^2 + (x_var(2) - x_var(5))^2 + (x_var(3) - x_var(6))^2)^(3 / 2));
D_B1 = ((x_var(1)^2 + x_var(2)^2 + x_var(3)^2)^(3 / 2)); % Inverse of The Distance of The Moon
D_B2 = ((x_var(4)^2 + x_var(5)^2 + x_var(6)^2)^(3 / 2)); % Inverse of The Distance of The Earth

A1 = (G * M_E) / D_A1_A2;
A2 = (G * M_M) / D_A1_A2;
B1 = (G * M_S) / D_B1;
B2 = (G * M_S) / D_B2;

Z1 = -B1 - A1;
Z2 = -B2 - A2;

% Stat-Space Matrix
A = [0 0 0 0 0 0 1 0 0 0 0 0
    0 0 0 0 0 0 0 1 0 0 0 0
    0 0 0 0 0 0 0 0 1 0 0 0
    0 0 0 0 0 0 0 0 0 1 0 0
    0 0 0 0 0 0 0 0 0 0 1 0
    0 0 0 0 0 0 0 0 0 0 0 1
    Z1 0 0 A1 0 0 0 0 0 0 0 0
    0 Z1 0 0 A1 0 0 0 0 0 0 0
    0 0 Z1 0 0 A1 0 0 0 0 0 0
    A2 0 0 Z2 0 0 0 0 0 0 0 0
    0 A2 0 0 Z2 0 0 0 0 0 0 0
    0 0 A2 0 0 Z2 0 0 0 0 0 0];

f_x = A * x_var;
end