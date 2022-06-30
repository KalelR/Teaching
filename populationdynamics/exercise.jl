using OrdinaryDiffEq

monod_equation(r, R, Ks) = r*R/(Ks+R)
μ(R, r, K) = minimum( monod_equation.(r, R, K) ) 

#K, C: k x n (resources x species)
function resource_consumption(j, C, R, r, K, N) 
    sum = 0.0 
    n = length(N)
    for i=1:n
        sum += C[j,i] * μ(R, r[i], K[:,i]) * N[i]
    end 
    return sum 
end

@inbounds function resourcecompetition!(du, u, p, t)
    K, C, S, r, m, D = p
    k, n = size(C)
    N = u[1:n]
    R = u[n+1:n+k]
    for i=1:n 
        du[i] = N[i]*(μ(R, r[i], K[:,i]) - m[i])
    end

    for j=1:k
        du[n+j] = D*(S[j] - R[j]) - resource_consumption(j, C, R, r, K, N)
    end
end



# Fig 1
## params

### 1 a,b, d
S = [10, 10, 10]
K = [1.0 0.75 0.25 0.7 0.2 0.65 0.68 0.38 0.46;
     0.25 1.0 0.75 0.2 1.01 0.55 0.83 1.10 0.85;
     0.75 0.25 1 1.10 0.7 0.95 0.6 0.5 0.77]

C = [0.10 0.20 0.15 0.05 0.01 0.40 0.30 0.20 0.25;
    0.15 0.10 0.20 0.15 0.30 0.35 0.25 0.02 0.35;
     0.20 0.15 0.10 0.25 0.05 0.20 0.40 0.15 0.10]
K = K[1:3, 1:3]
C = C[1:3, 1:3]


### Fig 2 
# S = [6, 10, 14, 4, 9]
# K = [0.39 0.34 0.30 0.24 0.23 0.41 0.20 0.45 0.14 0.15 0.38 0.28;
#      0.22 0.39 0.34 0.30 0.27 0.16 0.15 0.05 0.38 0.29 0.37 0.31;
#      0.27 0.22 0.39 0.34 0.30 0.07 0.11 0.05 0.38 0.41 0.24 0.25;
#      0.30 0.24 0.22 0.39 0.34 0.28 0.12 0.13 0.27 0.33 0.04 0.41;
#      0.34 0.30 0.22 0.20 0.39 0.40 0.50 0.26 0.12 0.29 0.09 0.16]

# C = [0.04 0.04 0.07 0.04 0.04 0.22 0.10 0.08 0.02 0.17 0.25 0.03;
#      0.08 0.08 0.08 0.10 0.08 0.14 0.22 0.04 0.18 0.06 0.20 0.04;
#      0.10 0.10 0.10 0.10 0.14 0.22 0.24 0.12 0.03 0.24 0.17 0.01;
#      0.05 0.03 0.03 0.03 0.03 0.09 0.07 0.06 0.03 0.03 0.11 0.05;
#      0.07 0.09 0.07 0.07 0.07 0.05 0.24 0.05 0.08 0.10 0.02 0.04]
# K = K[:, 1:5]
# C = C[:, 1:5]



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
T = 5000.0 #fig 2
tspan = (0.0, T)
prob = ODEProblem(resourcecompetition!, u0, tspan, p)
sol = solve(prob, Tsit5(), saveat=0:0.1:T, abstol=1e-8, reltol=1e-8); 
t = sol.t; u = sol[:, :] #u is numvars x T

using GLMakie
# fig = GLMakie.Figure()
# ax1 = Axis(fig[1, 1], ylabel="Ni", xlabel="t")
# for i=1:n lines!(ax1, t, u[i, :], color=["black", "red", "green", "cyan", "purple"][i]) end 
# # ax.ylabel = 

# ax2 = Axis3(fig[2, 1])
# lines!(ax2, u[1,:], u[2,:], u[3,:], color=t)


# Tplot=200 #fig 1
Tplot=300 #fig 2
fig = Figure()
ax1 = Axis(fig[1, 1], ylabel="Ni", xlabel="t")
for i=1:n lines!(ax1, t, u[i, :], color=["black", "red", "green", "cyan", "purple"][i]) end 
ax1.xticks=0:25:Tplot
xlims!(0, Tplot)

ax2 = Axis3(fig[1:3, 2])
l = lines!(ax2, u[1,:], u[2,:], u[3,:], color=t)
Colorbar(fig[2,3], l, label="t")

ax3 = Axis(fig[3, 1], ylabel="Total biomass", xlabel="t")
total_biomass = sum(u[1:n, :], dims=1)[1,:]
lines!(ax3, t, total_biomass, color=:black)
ax3.xticks=0:25:Tplot
xlims!(0, Tplot)

ax4 = Axis(fig[2, 1], ylabel="Rj", xlabel="t")
for i=1:k lines!(ax4, t, u[n+i, :], color=["black", "red", "green", "cyan", "purple"][i]) end 
ax4.xticks=0:25:Tplot
xlims!(0, Tplot)
save("populationdynamics-allinfo-n_$(n)-k_$(k).png", fig, px_per_unit=3)

#=
#Lyapunovs 
using DynamicalSystems
ds = ContinuousDynamicalSystem(prob)
Q = zeros(Float64, (n+k, n+k))
for i=1:n+k Q[i,i] = 1.0 end
λs = lyapunovspectrum(ds, 20000, Q; Ttr=5000)


#Bifurcation diagram
T = 4000.0 #fig 2
Ttr = 750.0
tspan = (0.0, T)
prob = ODEProblem(resourcecompetition!, u0, tspan, p)
K41s = 0.1:0.001:0.5
all_maxs = [Float64[] for i ∈ K41s];
all_mins = [Float64[] for i ∈ K41s];
for (idx, K41) ∈ enumerate(K41s)
    K[4,1] = K41
    p = [K, C, S, r, m, D]
    sol = solve(prob, Tsit5(), saveat=Ttr:0.1:T, abstol=1e-8, reltol=1e-8); 
    t = sol.t; u = sol[:, :] #u is numvars x T
    idx_maxs, maxs = findmaxima(u[1,:])
    idx_mins, mins = findminima(u[1,:])
    all_maxs[idx] = maxs
    all_mins[idx] = mins
end
fig = Figure()
ax = Axis(fig[1, 1], ylabel="Ni", xlabel="K41")
for (idx, K41) ∈ enumerate(K41s)
    scatter!(ax, [K41 for i=1:length(all_maxs[idx])], all_maxs[idx], color=(:black, 0.5), markersize=3)
end

save("bifurcationdiagram-full-n_$(n)-k_$(k).png", fig, px_per_unit=3)

=#








### test 
# lines(t, u[1,:])
# using Peaks
# pks, vals = findmaxima(u[1,:])
# scatter!(t[pks], vals, color=:red)
# pks, vals = findminima(u[1,:])
# scatter!(t[pks], vals, color=:green)