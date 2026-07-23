include("dependencies.jl")
##
# ----- This script is for a NS-centered halo, takes in K and rho_0, returns M_DM
n = 5
#rho_range_0 = [5e13, 6e13, 7e13, 8e13, 9e13, 1e14] #g/cm3  Nuclear density : 2.8e14
#Krange = [2.5e10, 2.75e10, 3e10, 3.25e10, 3.5e10] 
rho_range_0 = 10 .^(range(13.6,14.5,length = 200))
Krange = 10 .^(range(9.9,10.6,length = 200))    
rho_range = Math.cgs_density_inv(rho_range_0)
M_NS1 = 1.4
R_NS1 = 12
rend =2000

M_DMs = zeros(Float64,length(Krange),length(rho_range))
R_DMs = zeros(Float64,length(Krange),length(rho_range))
halo_bit = trues(size(M_DMs))
bound_bit = trues(size(M_DMs))
#profiles = Matrix{Tuple{Vector{Float64},Vector{Float64},Vector{Float64}}}(undef,length(Krange),length(rho_range))

let iteration = 0
for Ki in range(1,length(Krange))
    for rhoi in range(1,length(rho_range))
        iteration += 1
        output = Polytrope.solve_halo(n,Krange[Ki],M_NS1,R_NS1,rho_range[rhoi],rend)
    
        # if output.R_DM<=R_NS1       #Check for DM not escaping NS surface
        #     halo_bit[Ki,rhoi] = false
        # elseif output.R_DM>2000          #Check for DM exceeding 500km
        #     bound_bit[Ki,rhoi] = false
        # end
        
        M_DMs[Ki,rhoi] = output.M_DM
        R_DMs[Ki,rhoi] = output.R_DM

        #profiles[Ki,rhoi] = ((output.rs), (output.rho_r), (output.mass_r))
        print(iteration," ")
    end
end
end
CSV.write("heatmap_Rs.csv",DataFrame(R_DMs, :auto))
CSV.write("heatmap_Ms.csv",DataFrame(M_DMs, :auto))


##Heatmap
halo_bit = trues(size(M_DMs))
halo_bit = (R_DMs .> R_NS1)

bound_bit = trues(size(M_DMs))
bound_bit = (R_DMs .<1000.)

M_DM_filter = M_DMs;     R_DM_filter = R_DMs
M_DM_filter[.!BitArray(halo_bit)] .= NaN
M_DM_filter[.!BitArray(bound_bit)] .= NaN
R_DM_filter[.!BitArray(halo_bit)] .= NaN

sampleK = log10(3.25e10)
samplerho = log10(9e13)

mass_hm = heatmap(log10.(rho_range_0),log10.(Krange),M_DM_filter,
    ylabel = (L"\mathrm{log}_{10} K"), xlabel = L"\mathrm{log}_{10} \rho _0\ \mathrm{[g/cm^3]}", xlim=(13.6,14.5),
    colorbar_title = L"\mathrm{halo\ mass}\ \mathrm{[M_\odot]}", dpi = 500,clim=(0,0.501),
    guidefontsize=16,tickfontsize=11,legendfontsize=16, colorbar_titlefontsize=16,
    right_margin=5px,left_margin=5px,bottom_margin=5px)   #,size=(800,400)
contour!(mass_hm,log10.(rho_range_0),log10.(Krange),M_DM_filter,levels=[0.1], color=:white, linewidth=1, linestyle=:solid, colorbar_entry = false,label = "0.1M")
contour!(mass_hm,log10.(rho_range_0),log10.(Krange),M_DM_filter,levels=[0.2], color=:white, linewidth=1, linestyle=:dot, colorbar_entry = false,label = "0.2M")
# contour!(mass_hm,log10.(rho_range_0),log10.(Krange),R_DM_filter,levels=[50], color=:black, linewidth=1, linestyle=:solid, colorbar_entry = false)
# contour!(mass_hm,log10.(rho_range_0),log10.(Krange),R_DM_filter,levels=[100], color=:black, linewidth=1, linestyle=:solid, colorbar_entry = false)
annotate!(mass_hm,13.85,9.95,"DM halo contained within NS",10)
annotate!(mass_hm,14.1,10.58,"DM halo radius exceeds 1000km",10)
annotate!(mass_hm,14.45,10.12,(L"0.1\ \mathrm{M_\odot}",10,:white))
annotate!(mass_hm,14.45,10.32,(L"0.2\ \mathrm{M_\odot}",10,:white))
plot!(mass_hm,[samplerho],[sampleK],color = :green,marker = :o,markersize = 5, markerstrokewidth=3, legend = false)
annotate!(mass_hm,samplerho+0.046,sampleK+0.025,("|sample",8,:black))

rad_hm = heatmap(log10.(rho_range_0),log10.(Krange),R_DM_filter,
    ylabel = ("log10 K value"), xlabel = "log10 central density [g/cm^3]", 
    colorbar_title = "99% mass radius [km]", dpi = 500, clim = (0,1000)
    )
plot!(rad_hm,[samplerho],[sampleK],color = :green2,marker = :circle,markersize = 5, legend = false)


display(mass_hm)
#display(rad_hm)
# ## Named Array
# println("")
# println("Enclosed Mass (w/in 2000km) [Msun]")
# radii_named = NamedArray(
#     M_DMs,
#     (string.(Krange), string.(rho_range_0)),
#     ("K", "rho_0 [g/cm3]")
# )
# display(radii_named)

## plot one profile
# rs, rho_r, mass_r = profiles[15,10]
# p1 = plot(rs,Math.cgs_density(rho_r), xlabel = "Radius [km]", ylabel = L"Density [g/cm^3]")
# display(p1)
# p2 = plot(rs,mass_r, xlabel = "Radius [km]", ylabel = "Enclosed Mass [Msun]")
# display(p2)