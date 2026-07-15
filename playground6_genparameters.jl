include("dependencies.jl")
##
n = 5 
rho_range_0 = 10 .^(range(13.6,14.3,length = 10))   # cgs units, Nuclear density : 2.8e14
rho_range = Math.cgs_density_inv(rho_range_0)
M_NS1 = 1.4
R_NS1 = 12.
rend =2000

Ktarget = 2.7e10 #10^(10.45)   
M_DM_target = 0.2

M_DMs = zeros(length(rho_range))
R_DMs = zeros(length(rho_range))

for rhoi in range(1,length(rho_range))
    output = Polytrope.solve_halo(n,Ktarget,M_NS1,R_NS1,rho_range[rhoi],rend)

    rho_bool = BitArray(output.rho_r.>0)
    lastrho = ((output.rho_r)[rho_bool])[end]
    ind = findfirst(==(lastrho), output.rho_r)

    if output.rs[ind]<R_NS1
        valid_bit[rhoi] = 0
    elseif output.R_DM>500
        valid_bit[rhoi] = 0
    else
        valid_bit[rhoi] = 1
    end
    M_DMs[rhoi] = output.M_DM
    R_DMs[rhoi] = output.R_DM
end

M_interp = linear_interpolation(M_DMs,rho_range_0)  #pass in goal mass, return central density in cgs
R_interp = linear_interpolation(rho_range_0,R_DMs)

rho_result = M_interp(M_DM_target)
R_result = R_interp(rho_result)
println("K = ",Ktarget)
println("M halo = ", M_DM_target, " Msun")
println("central density = ", rho_result, " g/cm3")
println("Truncated radius ~= ", R_result, " km")
println("")
## plot halo to check
output = Polytrope.solve_halo(n,Ktarget,M_NS1,R_NS1,Math.cgs_density_inv(rho_result),rend)
p1 = plot(output.rs,output.mass_r,xlabel = "radius [km]", ylabel = "Mass enclosed [Msun]")
p2 = plot(output.rs,output.rho_r,xlabel = "radius [km]", ylabel = "density")
display(p1)
display(p2)