using NamedArrays
using Pkg          
Pkg.instantiate()
Pkg.activate(".")

using BenchmarkTools
#using Integrate
using CSV, DataFrames, Interpolations
using Plots
using OrdinaryDiffEq
using DifferentialEquations
using Printf, DelimitedFiles
using ForwardDiff
using LaTeXStrings

#include("source/GW.jl")
#include("source/install.jl")
#include("source/Integrate_r.jl")
include("source/Integrate.jl")
#include("source/make_poly_tables.jl")
include("source/math.jl")
#include("source/Pipeline.jl")
include("source/polytrope.jl")
#include("source/save.jl")
##
n = 4
Krange = [1e9, 5e9, 1e10, 5e10, 1e11]
M_DMrange = [0.01,0.05,0.1,0.5,1] #solar masses
M_NS1 = 0

radii = zeros(Float64,5,5)
rho_0s = zeros(Float64,5,5)
profiles = Matrix{Tuple{Vector{Float64},Vector{Float64}}}(undef,length(Krange),length(M_DMrange))

for Ki in range(1,length(Krange))
    for Mi in range(1,length(M_DMrange))
        compute = Polytrope.compute_polytrope(n); #calculate polytrope
        apply = Polytrope.apply_polytrope(compute,M_DMrange[Mi],M_NS1,Krange[Ki],n)
        radii[Ki,Mi] = apply.R_DM
        rho_0s[Ki,Mi] = apply.rho_r[1]
        profiles[Ki,Mi] = ((apply.rs), (apply.rho_r))
    end
end
rho_0s = Math.cgs_density(rho_0s)

println("Maximum radius [km]")
radii_named = NamedArray(
    radii,
    (string.(Krange), string.(M_DMrange)),
    ("K", "M_DM")
)
display(radii_named)

println("")
println("Central densities [cgs]")
rho_0_named = NamedArray(
    rho_0s,
    (string.(Krange), string.(M_DMrange)),
    ("K", "M_DM")
)
display(rho_0_named)

rs, rho_r = profiles[1,1]
plot(rs,Math.cgs_density(rho_r), xlabel = "Radius [km]", ylabel = L"Density [M_\odot / km^3]")