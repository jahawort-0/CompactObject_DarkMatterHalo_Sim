include("dependencies.jl")
##
# ----- This script is for a NS-centered halo, takes in K and rho_0, returns M_DM
n = 5
#rho_range_0 = [5e13, 6e13, 7e13, 8e13, 9e13, 1e14] #g/cm3  Nuclear density : 2.8e14
#Krange = [2.5e10, 2.75e10, 3e10, 3.25e10, 3.5e10] 
rho_range_0 = 10 .^(range(12,16,length = 50))
Krange = 10 .^(range(9.8,12.3,length = 50))    
#rho_range_0 = 10 .^ [13.5, 13.6, 13.7, 13.8, 13.9]
#Krange = 10 .^ [11.1, 11.2, 11.3, 11.4, 11.5]
rho_range = Math.cgs_density_inv(rho_range_0)
M_NS1 = 1.4
R_NS1 = 12
rend =2000

M_DMs = zeros(Float64,length(Krange),length(rho_range))
R_DMs = zeros(Float64,length(Krange),length(rho_range))
valid_bit = zeros(Bool,length(Krange),length(rho_range))
#radius_Mh = zeros(Float64,length(Krange),length(rho_range))       #Half mass radius
profiles = Matrix{Tuple{Vector{Float64},Vector{Float64},Vector{Float64}}}(undef,length(Krange),length(rho_range))

for Ki in range(1,length(Krange))
    for rhoi in range(1,length(rho_range))
        output = Polytrope.solve_halo(n,Krange[Ki],M_NS1,R_NS1,rho_range[rhoi],rend)

        #half_mass_interp = linear_interpolation(output.mass_r,apply.rs)  #create interpolation for half mass radius

        #radius_Mh[Ki,Mi] = half_mass_interp(M_DMrange[Mi]*0.5)    #half mass radius
        #radii[Ki,Mi] = apply.R_DM
        rho_bool = BitArray(output.rho_r.>0)
        lastrho = ((output.rho_r)[rho_bool])[end]
        ind = findfirst(==(lastrho), output.rho_r)

        if output.rs[ind]<R_NS1
            valid_bit[Ki,rhoi] = 0
        elseif output.R_DM>100
            valid_bit[Ki,rhoi] = 0
        else
            valid_bit[Ki,rhoi] = 1
        end
        #M_DMs[Ki,rhoi] = output.mass_r[end]
        M_DMs[Ki,rhoi] = output.M_DM
        R_DMs[Ki,rhoi] = output.R_DM

        profiles[Ki,rhoi] = ((output.rs), (output.rho_r), (output.mass_r))
    end
end
#rho_0s = Math.cgs_density(rho_0s)

## Named Array
println("")
println("Enclosed Mass (w/in 2000km) [Msun]")
radii_named = NamedArray(
    M_DMs,
    (string.(Krange), string.(rho_range_0)),
    ("K", "rho_0 [g/cm3]")
)
display(radii_named)

##Heatmap
M_DM_filter = similar(M_DMs)
M_DM_filter[valid_bit] .= M_DMs[valid_bit]
M_DM_filter[.!valid_bit] .= NaN
phm = heatmap(log10.(rho_range_0),log10.(Krange),M_DM_filter,
    ylabel = ("log10 K value"), xlabel = "log10 central density [g/cm^3]", 
    colorbar_title = "enclosed mass", dpi = 500, clim = (0,1)
    )
#heatmap(log10.(rho_range_0), log10.(Krange), R_DMs, ylabel = ("log10 K"), xlabel = ("log10 density"), colorbar_title = "DM radius",dpi = 500)
## plot 1 profile
rs, rho_r, mass_r = profiles[3,5]
p1 = plot(rs,Math.cgs_density(rho_r), xlabel = "Radius [km]", ylabel = L"Density [g/cm^3]")
plot!(p1, xlim=(0,2000))
display(p1)
p2 = plot(rs,mass_r, xlabel = "Radius [km]", ylabel = "Enclosed Mass [Msun]",xlim = (0,2000))
display(p2)