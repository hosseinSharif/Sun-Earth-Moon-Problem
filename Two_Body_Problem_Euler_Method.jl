closeall

using Plots

dt = 0.001
t = 0:dt:748.5
G = 0.01
M_S = 0.1

x0 = 1
y0 = 1
z0 = 0
xdot0 = -0.03
ydot0 = 0.01
zdot0 = 0

alpha = zeros(1,1)
beta = zeros(1,1)

x = zeros(6,length(t))
x[:,1] = [x0 y0 z0 xdot0 ydot0 zdot0]'

for i = 1:(length(t)-1)
    alpha = 1/((x[1,i]^2 + x[2,i]^2 + x[3,i]^2)^(3/2))
    beta = -G*M_S*alpha

    A = [0         0         0         1         0         0;
         0         0         0         0         1         0;
         0         0         0         0         0         1;
         beta      0         0         0         0         0;
         0         beta      0         0         0         0;
         0         0         beta      0         0         0]

    x[:,i+1] = x[:,i] + (dt.*A)*x[:,i]
end

gr()

plot(x[1,:],x[2,:])
