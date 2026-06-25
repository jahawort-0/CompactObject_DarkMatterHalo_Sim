include("dependencies.jl")
##
#---- This script takes in rho_0 and M_DM, and finds the associated K and R_DM
n = 5
rhorange = [1e11,5e11,1e12,5e12,1e13]
M_DMrange = [0.02,0.04,0.06,0.08,0.1,0.15,0.2,0.25,0.3] #solar masses
M_NS1 = 0

radii = zeros(Float64,length(rhorange),length(M_DMrange))   #Outer Radius
Ks = zeros(Float64,length(rhorange),length(M_DMrange))      #K value
radius_Mh = zeros(Float64,length(rhorange),length(M_DMrange))       #Half mass radius
profiles = Matrix{Tuple{Vector{Float64},Vector{Float64}}}(undef,length(Krange),length(M_DMrange))   #saves density profile for plotting

for rhoi in range(1,length(rhorange))
    for Mi in range(1,length(M_DMrange))
        compute = Polytrope.compute_polytrope(n); #calculate polytrope
        apply = Polytrope.apply_polytrope_rho(compute,M_DMrange[Mi],M_NS1,Math.cgs_density_inv(rhorange[rhoi]),n)   #redim/apply polytrope

        half_mass_interp = linear_interpolation(apply.mass_r,apply.rs)  #create interpolation for half mass radius

        radius_Mh[rhoi,Mi] = half_mass_interp(M_DMrange[Mi]*0.5)    #half mass radius
        radii[rhoi,Mi] = apply.R_DM     #outer radius
        Ks[rhoi,Mi] = apply.K       #K value
        profiles[rhoi,Mi] = ((apply.rs), (apply.rho_r))
    end
end

println("")
println("Maximum radius [km]")
radii_named = NamedArray(
    radii,
    (string.(rhorange), string.(M_DMrange)),
    ("rho_0[cgs]", "M_DM[Msun]")
)
display(radii_named)

println("")
println("Half Mass radius [km]")
radius_Mh_named = NamedArray(
    radius_Mh,
    (string.(rhorange), string.(M_DMrange)),
    ("rho_0[cgs]", "M_DM[Msun]")
)
display(radius_Mh_named)

println("")
println("K values")
K_named = NamedArray(
    Ks,
    (string.(rhorange), string.(M_DMrange)),
    ("rho_0[cgs]", "M_DM[Msun]")
)
display(K_named)

#rs, rho_r = profiles[1,1]
#plot(rs,Math.cgs_density(rho_r), xlabel = "Radius [km]", ylabel = L"Density [M_\odot / km^3]")