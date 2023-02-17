# In The Name of Allah

#= Sun-Earth Problem:
    Solving The Differential Equation of Earth Displacement And Velocity Vector 
    With Runge-Kutta 4th Order Method
=#

#= Note!
    In genral we can say that any other system of two masses with mass of the one
    of them very greater than the other, will go this way! =#

using Plots

# Variables And Prameters Zeroings


# Constants: (SI Units)
const G = 6.6743e-11; # m^3kg^-1s^-2
const M_S = 1.989e30; # kg
const AU = 149597870700; # 1.5e11
const V0 = (G * M_S / AU)^0.5; #SI

# State-Space Matrix Calculation Function:
function f(x_var)
    Alpha = 1 / ((x_var[1]^2 + x_var[2]^2 + x_var[3]^2)^(3 / 2))
    Beta_var = -G * M_S * Alpha

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
Stop_Year = 10; # The year that the simulation will stop.

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
