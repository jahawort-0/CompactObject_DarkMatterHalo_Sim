include("dependencies.jl")
##
# ----- This script take in K and M_DM, returns rho_0 and radius
n = 5
Krange = [2e9,4e9,6e9,8e9,1e10,1.2e10,1.4e10,1.6e10,1.8e10,2e10,4e10,6e10]
M_DMrange = [0.02,0.04,0.06,0.08,0.1,0.15,0.2,0.25,0.3] #solar masses
M_NS1 = 0

radii = zeros(Float64,length(Krange),length(M_DMrange))
rho_0s = zeros(Float64,length(Krange),length(M_DMrange))
radius_Mh = zeros(Float64,length(Krange),length(M_DMrange))       #Half mass radius
profiles = Matrix{Tuple{Vector{Float64},Vector{Float64}}}(undef,length(Krange),length(M_DMrange))

for Ki in range(1,length(Krange))
    for Mi in range(1,length(M_DMrange))
        compute = Polytrope.compute_polytrope(n); #calculate polytrope
        apply = Polytrope.apply_polytrope(compute,M_DMrange[Mi],M_NS1,Krange[Ki],n)

        half_mass_interp = linear_interpolation(apply.mass_r,apply.rs)  #create interpolation for half mass radius

        radius_Mh[Ki,Mi] = half_mass_interp(M_DMrange[Mi]*0.5)    #half mass radius
        radii[Ki,Mi] = apply.R_DM
        rho_0s[Ki,Mi] = apply.rho_r[1]
        profiles[Ki,Mi] = ((apply.rs), (apply.rho_r))
    end
end
rho_0s = Math.cgs_density(rho_0s)

println("")
println("Maximum radius [km]")
radii_named = NamedArray(
    radii,
    (string.(Krange), string.(M_DMrange)),
    ("K", "M_DM[Msun]")
)
display(radii_named)

println("")
println("Half Mass radius [km]")
radius_Mh_named = NamedArray(
    radius_Mh,
    (string.(Krange), string.(M_DMrange)),
    ("K", "M_DM[Msun]")
)
display(radius_Mh_named)

println("")
println("Central densities [cgs]")
rho_0_named = NamedArray(
    rho_0s,
    (string.(Krange), string.(M_DMrange)),
    ("K", "M_DM[Msun]")
)
display(rho_0_named)

#rs, rho_r = profiles[1,1]
#plot(rs,Math.cgs_density(rho_r), xlabel = "Radius [km]", ylabel = L"Density [M_\odot / km^3]")