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
n=3
rho_0 = 1*10^-8 #Assume a central density
titlestring = L"\rho_0 = 1E-8, n = 4"

output = Polytrope.compute_polytrope(n); #calculate polytrope
# Plot theta vs xi

xi1 = output.root
xis = output.sol.t
thetas = output.sol[1,:];
dthetas = output.sol[2,:];
p1 = plot(xis,thetas, label="Polytrope solution")

thetas_cut = Polytrope.polytrope_profile.(Ref(output),xis)
# println(typeof(thetas_cut))
# println(thetas_cut[1,:])
thetas_cut = reduce(hcat,thetas_cut)
plot!(xis,thetas_cut[1,:],label = "Cut Solution")
plot!(xlabel = L"\xi", ylabel = L"\theta",title = titlestring)
display(p1)

#  Plot M(dimensionless) vs xi
mass_r_nondim = output.mass
#println(mass_r_nondim)
p2 = plot(range(0,xi1,length=1000),mass_r_nondim(range(0,xi1,length=1000)))
plot!(xlabel = L"\xi", ylabel = L"M_{\rm enclosed}", title = titlestring)
display(p2)

# Non dimensionalize
M_DM = 2 #M_solar
M_nondim = mass_r_nondim(xi1)

K = (1/(4*pi*rho_0)*M_DM/M_nondim)^(2/3) *4*pi*Math.G*rho_0^(1-(1/n)) / (n+1) #calculate polytrope constant
alpha = ((n+1)*K/(4*pi*Math.G*rho_0^(1-(1/n))))^(1/2)   #alpha is nondim. constant, easier to use
 
#distance
rs = xis.*alpha #an array of distances in km
R_max = xi1*alpha  #outer edge of halo in km
println(R_max)

#mass/density
rho_r = rho_0.*(thetas_cut).^n   #the density as a function of distance
mass_r = 4*pi* alpha^3 * rho_0 * mass_r_nondim(range(0,xi1,length=1000))

p3 = plot(rs,rho_r[1,:], xlabel = "distance [km]", ylabel = L"\rho(r) [M_\odot/km^3]", title = titlestring)
display(p3)

p4 = plot(range(0,R_max,length=1000),mass_r,
xlabel = "distance [km]", ylabel = L"M_{\rm enclosed}\;[M_\odot]", title  = titlestring)
display(p4)
##
n=3
rho_0 = 1*10^-8 #Assume a central density
M_DM = 2

output = Polytrope.compute_polytrope(n); #calculate polytrope
output2 = Polytrope.apply_polytrope(output,M_DM,rho_0,n)

p5 = plot(output2.rs,output2.rho_r)
display(p5)
p6 = plot(output2.rs,output2.mass_r)
display(p6)