include("dependencies.jl")
##
n = 5
Krange = [2e9,4e9,6e9,8e9,1e10,1.2e10,1.4e10,1.6e10,1.8e10,2e10,4e10,6e10]
M_DMrange = [0.02,0.04,0.06,0.08,0.1,0.15,0.2,0.25,0.3] #solar masses
M_NS1 = 0

radii = zeros(Float64,length(Krange),length(M_DMrange))
rho_0s = zeros(Float64,length(Krange),length(M_DMrange))
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

#rs, rho_r = profiles[1,1]
#plot(rs,Math.cgs_density(rho_r), xlabel = "Radius [km]", ylabel = L"Density [M_\odot / km^3]")