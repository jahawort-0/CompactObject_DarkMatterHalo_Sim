include("dependencies.jl")

## Set parameters
n = 5
K = 1.5e9
M_DM = 0.1
M_NS1 = 1.35
M_NS2 = 1.5

#Initialize Polytrope
poly = Polytrope.compute_polytrope(n)   #integrate nondimensional polytrope
realpoly = Polytrope.apply_polytrope(poly,M_DM,M_NS1,K,n)   # redimensionalize polytrope
circ_orbit=Integrate.initial_circular_orbit(M_NS1, M_NS2, M_DM, realpoly);  #Solve for Roche sepeartion, ciruclar orbit
# initial_circular_orbit returns [a, M_NS1, M_NS2, M_DM, R_DM, J, 0.]

ICs = [circ_orbit[4],circ_orbit[1]];    #[M_DM, a]
parameters = [circ_orbit[2],circ_orbit[3],realpoly.mass_interp];    #[M_NS1, M_NS2, mass(r)]

M1 = M_DM + M_NS1
η = (M1*M_NS2)/((M1+M_NS2)^2)
t_decay = 5/256*((Math.c)^5) * (circ_orbit[1]^4) /(Math.G*(M1+M_NS2))^3/ η
## Integration
tend = 1.65729454e5
times=range(0,stop=tend,length=1000)

integration_sol=Integrate.integrate_halo(ICs,parameters,(0,tend)); #The -1 tells the function that we are decreasing in separation

sol_interpolated=integration_sol(times)

M_DM_t = sol_interpolated[1,:]
a_t = sol_interpolated[2,:]

p1 = plot(times,a_t, xlabel = "Time [s]", ylabel = "seperation [km]",dpi = 200,ylim=(0,circ_orbit[1]*1.1))
p2 = plot(times,M_DM_t, xlabel = "Time [s]", ylabel = "Mass [Msun]",dpi = 200,ylim=(0,M_DM*1.1))
display(p2)
display(p1)