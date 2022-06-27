using OrdinaryDiffEq

monod_equation(r, R, Ks) = r*R/(Ks+R)
μ(i, R, r, K) = minimum( monod_equation.(r, R, K[:, i]) ) 

#K, C: k x n (resources x species)
function resource_consumption(j, C, R, r, K, N) 
    sum = 0.0 
    n = length(N)
    for i=1:n
        sum += C[j,i] * μ(i, R, r, K) * N[i]
    end 
    return sum 
end

@inbounds function resourcecompetition!(du, u, p, t)
    K, C, S, r, m, D = p
    k, n = size(C)
    N = u[1:n]
    R = u[n+1:n+k]
    for i=1:n 
        du[i] = N[i]*(μ(i, R, r, K) - m[i])
    end

    for j=1:k
        du[n+j] = D*(S[j] - R[j]) - resource_consumption(j, C, R, r, K, N)
    end
end



# Fig 1
## params

### 1 a,b, d
# S = [10, 10, 10]
# K = [1.0 0.75 0.25 0.7 0.2 0.65 0.68 0.38 0.46;
#      0.25 1.0 0.75 0.2 1.01 0.55 0.83 1.10 0.85;
#      0.75 0.25 1 1.10 0.7 0.95 0.6 0.5 0.77]

# C = [0.10 0.20 0.15 0.05 0.01 0.40 0.30 0.20 0.25;
#     0.15 0.10 0.20 0.15 0.30 0.35 0.25 0.02 0.35;
#      0.20 0.15 0.10 0.25 0.05 0.20 0.40 0.15 0.10]
# K = K[1:3, 1:3]
# C = C[1:3, 1:3]


### Fig 2 
S = [6, 10, 14, 4, 9]
K = [0.39 0.34 0.30 0.24 0.23 0.41 0.20 0.45 0.14 0.15 0.38 0.28;
     0.22 0.39 0.34 0.30 0.27 0.16 0.15 0.05 0.38 0.29 0.37 0.31;
     0.27 0.22 0.39 0.34 0.30 0.07 0.11 0.05 0.38 0.41 0.24 0.25;
     0.30 0.24 0.22 0.39 0.34 0.28 0.12 0.13 0.27 0.33 0.04 0.41;
     0.34 0.30 0.22 0.20 0.39 0.40 0.50 0.26 0.12 0.29 0.09 0.16]

C = [0.04 0.04 0.07 0.04 0.04 0.22 0.10 0.08 0.02 0.17 0.25 0.03;
     0.08 0.08 0.08 0.10 0.08 0.14 0.22 0.04 0.18 0.06 0.20 0.04;
     0.10 0.10 0.10 0.10 0.14 0.22 0.24 0.12 0.03 0.24 0.17 0.01;
     0.05 0.03 0.03 0.03 0.03 0.09 0.07 0.06 0.03 0.03 0.11 0.05;
     0.07 0.09 0.07 0.07 0.07 0.05 0.24 0.05 0.08 0.10 0.02 0.04]

K = K[:, 1:5]
C = C[:, 1:5]



k, n = size(C)
r = [1 for i=1:n] #1/day 
D = 0.25 #1/day 
m = [D for i=1:n]
p = [K, C, S, r, m, D]
# p = Dict(:K=>K, :C=>C, :S=>S, :r=>r, :m=>m, :D=>D)

## ics 
R0 = S
N0 = 0.1 .+ (1:n) ./ 100
u0 = [N0; R0]

## solve
# T = 200 #fig 1
T = 300 #fig 2
tspan = (0, T)
prob = ODEProblem(resourcecompetition!, u0, tspan, p)
sol = solve(prob, Tsit5(), saveat=0:0.1:T, abstol=1e-8, reltol=1e-8); 
t = sol.t; u = sol[:, :] #u is numvars x T

using GLMakie
fig = GLMakie.Figure()
ax1 = Axis(fig[1, 1])
for i=1:n lines!(ax1, t, u[i, :], color=["black", "red", "green", "cyan", "purple"][i]) end 
# ax.ylabel = 

ax2 = Axis3(fig[2, 1])
lines!(ax2, u[1,:], u[2,:], u[3,:])


fig = GLMakie.Figure()
ax1 = Axis(fig[1, 1:2])
for i=1:n lines!(ax1, t, u[i, :], color=["black", "red", "green", "cyan", "purple"][i]) end 
# ax.ylabel = 

ax2 = Axis3(fig[2, 1])
lines!(ax2, u[1,:], u[2,:], u[3,:])

ax3 = Axis(fig[2, 2])
total_biomass = sum(u[1:n, :], dims=1)[1,:]
lines!(ax3, t, total_biomass)