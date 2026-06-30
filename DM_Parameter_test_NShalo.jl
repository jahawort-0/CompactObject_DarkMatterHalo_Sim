include("dependencies.jl")
##
# ----- This script is for a NS-centered halo, takes in K and rho_0, returns M_DM
n = 5
#Krange = [1e11,2e11,4e11,6e11,8e11, 1e12]
Krange = [1e9, 1e10, 1e11, 1e12, 1e13, 1e14,1e15]
#rho_range_0 = [2e14,4e14,6e14,8e14,1e15] #g/cm3    Nuclear density : 2.8e14
rho_range_0 = [1e9,1e10,1e11,1e12,1e13,1e14,1e15] #g/cm3    Nuclear density : 2.8e14
rho_range = Math.cgs_density_inv(rho_range_0)
M_NS1 = 1.4
R_NS1 = 12
rend =1000

M_DMs = zeros(Float64,length(Krange),length(rho_range))
#radius_Mh = zeros(Float64,length(Krange),length(rho_range))       #Half mass radius
profiles = Matrix{Tuple{Vector{Float64},Vector{Float64},Vector{Float64}}}(undef,length(Krange),length(rho_range))

for Ki in range(1,length(Krange))
    for rhoi in range(1,length(rho_range))
        output = Polytrope.solve_halo(n,Krange[Ki],M_NS1,R_NS1,rho_range[rhoi],rend)

        #half_mass_interp = linear_interpolation(output.mass_r,apply.rs)  #create interpolation for half mass radius

        #radius_Mh[Ki,Mi] = half_mass_interp(M_DMrange[Mi]*0.5)    #half mass radius
        #radii[Ki,Mi] = apply.R_DM
        rho_bool = BitArray(output[2].>0)
        lastrho = ((output[2])[rho_bool])[end]
        ind = findfirst(==(lastrho), output[2])
        if output[1][ind]<R_NS1
            M_DMs[Ki,rhoi] = NaN
        elseif output[2][end]>(0.1/(4*pi*1000^2))
            M_DMs[Ki,rhoi] = NaN
        else
        M_DMs[Ki,rhoi] = output[4][end]
        end
        profiles[Ki,rhoi] = ((output[1]), (output[2]), (output[4]))
    end
end
#rho_0s = Math.cgs_density(rho_0s)

println("")
println("Enclosed Mass (w/in 1000km) [Msun]")
radii_named = NamedArray(
    M_DMs,
    (string.(Krange), string.(rho_range_0)),
    ("K", "rho_0 [g/cm3]")
)
display(radii_named)

# rs, rho_r, mass_r = profiles[2,5]
# p1 = plot(rs,Math.cgs_density(rho_r), xlabel = "Radius [km]", ylabel = L"Density [g/cm^3]")
# plot!(xlim=(0,1000))
# display(p1)
# p2 = plot(rs,mass_r, xlabel = "Radius [km]", ylabel = "Enclosed Mass [Msun]")
# display(p2)