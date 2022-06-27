using DynamicalSystems
σ = 10.0; ρ = 28.0; β = 8/3
u0 = rand(3)
lor = Systems.lorenz(u0; σ, ρ, β)


# Dissipation rate
λs = lyapunovspectrum(lor, 100000; Ttr=10000)
##ans:    0.9046781537145369,   2.6143953086660426e-5,  -14.57128897720743
mean_diss_rate = sum(λs)
divF = -(σ + 1 + β) #since divF is constant, its average should be equal to itself. The average is also equal to mean_diss_rate, as we can verify!

# Chaotic attractor 
T = 1000
Δt = 0.001
u0 = [1.0, 1.0, 1.0]
tr = trajectory(lor, T, u0; Δt)
t = 0:Δt:T
using GLMakie 

fig = Figure(resolution=(1920, 1080))
ax = Axis3(fig[1, 1], azimuth=2.1, elevation=π/8)
l = lines!(ax, tr[:,1], tr[:,2], tr[:,3], color=t)
Colorbar(fig[1,2], l, label="t")
save("lorenz-attractor-colorastime.png", fig, px_per_unit=3)

## plot with switching 
function indices_switching(tr)
xs = tr[:,1]
idxs = zeros(Int64, length(xs))
for i = 2:length(xs) 
    if xs[i] * xs[i-1] < 0 #changed signs
        idxs[i] = idxs[i-1] + 1
    else 
        idxs[i] = idxs[i-1]
    end
end
return idxs
end

idxs = indices_switching(tr)
fig = Figure(resolution=(1920, 1080))
ax = Axis3(fig[1, 1], azimuth=2.1, elevation=π/8)
# lines!(ax, tr[:,1], tr[:,2], tr[:,3], color=idxs, colormap:tab20)
lines!(ax, tr[:,1], tr[:,2], tr[:,3], color=idxs, colormap=:gist_rainbow)
scatter!([u0[1]], [u0[2]], [u0[3]], color=:black, markersize=6000)
ax = Axis(fig[2,1])
lines!(ax, t, tr[:,1], color=idxs, colormap=:gist_rainbow)
ax.xlabel="t"; ax.ylabel="x"
save("lorenz-attractor-switching.png", fig, px_per_unit=3)


## compare ics 
u02 = deepcopy(u0)
u02[1] += 1e-8 
tr2 = trajectory(lor, T, u02; Δt)

fig = Figure(resolution=(1920, 1080))
# ax = Axis3(fig[1, 1:2])
ax = Axis3(fig[1:2, 2], azimuth=2.1, elevation=π/8)
alpha = 0.5
lines!(ax, tr[:,1], tr[:,2], tr[:,3], color=(:black, 0.3))
lines!(ax, tr2[:,1], tr2[:,2], tr2[:,3], color=(:red, 0.8))
scatter!([u0[1]], [u0[2]], [u0[3]], color=(:black, 0.5), markersize=6000)
scatter!([u02[1]], [u02[2]], [u02[3]], color=(:red, 0.8), markersize=5000)
ax.xlabel="x"; ax.ylabel="y"; ax.zlabel="z"

ax = Axis(fig[1,1])
lines!(ax, t, tr[:,1], color=:black)
lines!(ax, t, tr2[:,1], color=:red)
xlims!(0, 200)
ax.xlabel="t"; ax.ylabel="x"

ax = Axis(fig[2,1])
norm(x, y) = sqrt.(sum((x .- y) .^ 2, dims=2)[:,1])
distance = norm(Matrix(tr), Matrix(tr2))
lines!(ax, t, distance, color=:blue)
xlims!(0, 200)
ax.xlabel="t"; ax.ylabel="euclidian distance"

save("lorenz-attractor-twoics.png", fig, px_per_unit=3)


# Eigenvalues of C- and C+ 
using LinearAlgebra
# ρs = range(1,2, length=100)
ρs = range(1,30, length=1000)
all_λs = zeros(ComplexF64, (length(ρs), 3));
for (idx, ρ) ∈ enumerate(ρs)
    C₊ = [sqrt(β*(ρ-1)), sqrt(β*(ρ-1)), ρ-1]
    C₋ = [-sqrt(β*(ρ-1)), -sqrt(β*(ρ-1)), ρ-1]
    J = Matrix(lor.jacobian(C₊, [σ, ρ, β], 0))
    λs = eigvals(J)
    all_λs[idx,:] = λs
end

fig = Figure() 
for i=1:3
ax = Axis(fig[i, 1], ylabel="Re(λ_$(i))", xlabel="ρ")
lines!(ax, ρs, real.(all_λs[:,i]))
hlines!(ax, 0, color=:red, linestyle="--")
vlines!(ax, 1.346, color=:purple, linestyle="--")
vlines!(ax, 24.74, color=:purple, linestyle="--")
xlims!(ρs[1], ρs[end])
ax = Axis(fig[i, 2], ylabel="Im(λ_$(i))", xlabel="ρ")
lines!(ax, ρs, imag.(all_λs[:,i]))
vlines!(ax, 1.346, color=:purple, linestyle="--")
vlines!(ax, 24.74, color=:purple, linestyle="--")
xlims!(ρs[1], ρs[end])
end
save("lorenz-eigenvalues-twootherfixedpoints-near-r_1.png", fig, px_per_unit=3)
save("lorenz-eigenvalues-twootherfixedpoints-near-until-r_30-hopf.png", fig, px_per_unit=3)

## Interactive 
using InteractiveDynamics
using DynamicalSystems, GLMakie
using OrdinaryDiffEq

diffeq = (alg = Tsit5(), adaptive = false, dt = 0.01)
ps = Dict(
    1 => 1:0.1:30,
    2 => 10:0.1:50,
    3 => 1:0.01:10.0,
)
pnames = Dict(1 => "σ", 2 => "ρ", 3 => "β")

lims = (
    (-30, 30),
    (-30, 30),
    (0, 100),
)

ds = Systems.lorenz()

u1 = [10,20,40.0]
u3 = [20,10,40.0]
u0s = [u1, u3]

idxs = (1, 2, 3)
diffeq = (alg = Tsit5(), dt = 0.01, adaptive = false)

figure, obs, step, paramvals = interactive_evolution(
    ds, u0s; ps, idxs, tail = 1000, diffeq, pnames, lims
)

# Use the `slidervals` observable to plot fixed points
lorenzfp(ρ,β) = [
    Point3f(sqrt(β*(ρ-1)), sqrt(β*(ρ-1)), ρ-1),
    Point3f(-sqrt(β*(ρ-1)), -sqrt(β*(ρ-1)), ρ-1),
]

fpobs = lift(lorenzfp, slidervals[2], slidervals[3])
ax = content(figure[1,1][1,1])
scatter!(ax, fpobs; markersize = 5000, marker = :diamond, color = :black)