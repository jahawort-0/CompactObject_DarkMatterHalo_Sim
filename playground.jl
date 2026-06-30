include("dependencies.jl")


##
n=5
K = 1e10 #Assume polytropic constant
titlestring = "n = 5"

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

#  Plot M(dimensionless) vs xi
mass_r_nondim = output.mass
#println(mass_r_nondim)
p2 = plot(range(0,xi1,length=1000),mass_r_nondim(range(0,xi1,length=1000)))
plot!(xlabel = L"\xi", ylabel = L"M_{\rm enclosed}", title = titlestring)

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

p4 = plot(range(0,R_max,length=1000),mass_r,
xlabel = "distance [km]", ylabel = L"M_{\rm enclosed}\;[M_\odot]", title  = titlestring)

display(p1)
display(p2)
#display(p3)
#display(p4)



##
"""This does all of the above more compactly"""
n=5
#rho_0 = 1*10^-8 #Assume a central density
n = 5
K = 1.5e9
M_DM = 0.1
M_NS1 = 1.35

output = Polytrope.compute_polytrope(n); #calculate polytrope
output2 = Polytrope.apply_polytrope(output,M_DM,M_NS1,K,n)

#p5 = plot(output2.rs,output2.rho_r, xlabel = "Radius [km]", ylabel = L"Density [M_\odot / km^3]")
#display(p5)
mass_interp = output2.mass_interp
p6 = plot(output2.rs,mass_interp(output2.rs), xlabel = "Radius [km]", ylabel = L"Mass\ enclosed [M_\odot]",dpi =200)
display(p6)
#print(output2.rho_r[1])


## Enclosed mass

n=5
K = 1e10
M_DM = 0.1
M_NS1 = 0

output = Polytrope.compute_polytrope(n); #calculate polytrope
output2 = Polytrope.apply_polytrope(output,M_DM,M_NS1,K,n)

M_enclosed = output2.mass_interp(1)

# using analytical solution
analytical_output = Polytrope.n5apply_polytrope(M_DM,K)
radius_range = range(0,5,100)

p1 = plot(radius_range,output2.mass_interp.(radius_range),label="numerical solver")
plot!(radius_range,analytical_output.mass_interp.(radius_range),label="analytical solver")
plot!(xlabel = "Radius [km]", ylabel = L"Mass\ enclosed [M_\odot]")
display(p1)

#residuals
res = output2.mass_interp.(radius_range) .- analytical_output.mass_interp.(radius_range)
p2 = plot(radius_range, res,xlabel = "Radius [km]", ylabel = L"Mass\ enclosed [M_\odot]")
plot!(radius_range,zeros(length(radius_range)))


## 06/26 New polytrope solution
#old first
n = 5
K = 1.5e9
M_DM = 0.1

output = Polytrope.compute_polytrope(n); #calculate polytrope
output2 = Polytrope.apply_polytrope(output,M_DM,0,K,n)

p1 = plot(output2.rs,output2.rho_r, xlabel = "Radius [km]", ylabel = L"Density [M_\odot / km^3]",dpi = 200, label = "no NS")
p2 = plot(output2.rs,output2.mass_r, xlabel = "Radius [km]", ylabel = L"Mass [M_\odot]",dpi = 200, label = "no NS")


#Then new
rho_0 = 1e-3#4.7796e-6
M_NS1 = 1.4
R_NS1 = 12
NS_density = M_NS1 / (4/3*pi*R_NS1^3)
K_NS = 7e10

rend = 1000
# p = [n,K,M_NS1,R_NS1,rho_0]
# u = [NS_density + rho_0, 0]
# du = [Float64(0),Float64(0)]
# Polytrope.halo_polytrope_problem!(du,u,p,1e-8)
# println(u)
# println(du)

output3 = Polytrope.solve_halo(n,K_NS,M_NS1,R_NS1,rho_0,rend)

plot!(p1,output3[1,:],output3[2,:],xlim = (0,rend), ylim = (0,5e-6), label = "w/ NS")
plot!(p2,output3[1,:],output3[4,:],xlim = (0,rend), ylim = (0,0.1), label = "w/ NS")
display(p2)
display(p1)
