# In The Name of Allah
# In Persuit of The Moon Project
# Hossein Sharif 1401/11/28 14:00
# Location: Islamic Republic of Iran - Tehran

#= Sun-Earth Problem:
    Solving The Differential Equation of Earth Displacement And Velocity Vector 
    With Runge-Kutta 4th Order Method
=#

#= Note!
    In genral, we can say that any other system of two masses with mass of the one
    of them very greater than the other, will go this way! =#

using Plots, Plots.Measures

# Constants: (SI Units)
const G = 6.6743e-11; # Gravitational Constant (m^3kg^-1s^-2)
const M_S = 1.989e30; # The Mass of The Sun (kg)
const AU = 149597870700; # Mean Distance Between The Sun And The Earth (m)

## Linear Velocity of The Earth In Orbit Around The Sun:
# This equation can be proved with the help of polar
# coordinates and two assumptions:
# 1- r is constant in orbital motion (r = AU)
# 2- Theta_Dot is constant in orbital motuon (Theta_Dot = cte.)
const V0 = (G * M_S / AU)^0.5;

# State-Space Matrix Calculation Function:
function f(x_var)
    Alpha = 1 / ((x_var[1]^2 + x_var[2]^2 + x_var[3]^2)^(3 / 2))
    Beta_var = -G * M_S * Alpha

    # Stat-Space Matrix
    A = [0 0 0 1 0 0
        0 0 0 0 1 0
        0 0 0 0 0 1
        Beta_var 0 0 0 0 0
        0 Beta_var 0 0 0 0
        0 0 Beta_var 0 0 0]

    B = A * x_var
    return B
end

# Initialize System:
Day_Hours = 23 + (56 / 60) + (4 / 3600); # Hours of A Day
Year_Days = 365.2425; # Days of 
Year_Seconds = Year_Days * Day_Hours * 3600;

dD = 5; # Simulation Step In Days
Stop_Year = 3; # The year that the simulation will stop.

dY = dD / 365.2425; # Simulation  Step In Years
Δt = dY * Year_Seconds; # In Seconds

Year = 0:dY:Stop_Year; # In Years
Δt = dY * Year_Seconds; # In Seconds

xE0 = -AU;
yE0 = 0;
zE0 = 0;

VxE0 = 0;
VyE0 = -V0;
VzE0 = 0;

x = zeros(6, length(Year));
x[:, 1] = [xE0 yE0 zE0 VxE0 VyE0 VzE0]';

# Integrating With RK4:
for i = 1:(length(Year)-1)
    k1 = Δt .* f(x[:, i])
    k2 = Δt .* f(x[:, i] + k1 ./ 2)
    k3 = Δt .* f(x[:, i] + k2 ./ 2)
    k4 = Δt .* f(x[:, i] + k3)

    x[:, i+1] = x[:, i] + (k1 .+ 2 * k2 .+ 2 * k3 .+ k4) ./ 6

end

plotlyjs()

plot(x[1, :], x[2, :])
