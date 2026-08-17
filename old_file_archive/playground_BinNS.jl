include("dependencies.jl")

## Set parameters
n = 5
K = 1.5e9
M_DM = 0.101
M_NS1 = 1.3
M_NS2 = 1.4
a_start = 20.

#Initialize Polytrope
poly = Polytrope.compute_polytrope(n)   #integrate nondimensional polytrope
realpoly = Polytrope.apply_polytrope(poly,M_DM,M_NS1,K,n)   # redimensionalize polytrope
circ_orbit=Integrate.initial_circular_orbit(M_DM, M_NS1, M_NS2, a_start, realpoly);  #Solve for Roche sepeartion, ciruclar orbit
# initial_circular_orbit returns [a, M_DM, M_NS1, M_NS2, R_DM, J, 0.]

ICs = [circ_orbit[1],circ_orbit[2],circ_orbit[3],circ_orbit[4],0,circ_orbit[6]];    #[a, M_DM, M_NS1, M_NS2, phase, J]
parameters = [realpoly.mass_interp];    #[M_NS1, M_NS2, mass(r)]

# extra, calculate approx decay time
M1 = M_DM + M_NS1
η = (M1*M_NS2)/((M1+M_NS2)^2)
t_decay = 5/256*((Math.c)^5) * (circ_orbit[1]^4) /(Math.G*(M1+M_NS2))^3/ η
## Integration
tend = 1e-3
times=range(0,stop=tend,length=10000)
#times=range(tend-2.5,stop=tend,length=1000)

# du = zeros(Float64, 6)
# Integrate.halo_problem!(du,ICs,parameters,0)
# println(u)
# println(du)
#integration_sol = solve(ODEProblem(Integrate.halo_problem!, ICs, (0,tend), parameters))
integration_sol=Integrate.integrate_halo(ICs,parameters,(0,tend)); #The -1 tells the function that we are decreasing in separation

sol_interpolated=integration_sol(times)

# find point of a_min and cut all data after that point
a_min = 12
a_diff = sol_interpolated[1,:] .- a_min
keep = (a_diff.>1e-6)
sol_interpolated = sol_interpolated[:,keep]
times = times[keep]

M_DM_t = sol_interpolated[2,:]
a_t = sol_interpolated[1,:]

p1 = plot(times,a_t, xlabel = "Time [s]", ylabel = "seperation [km]",dpi = 200,ylim=(0,circ_orbit[1]*1.1))
p2 = plot(times,M_DM_t, xlabel = "Time [s]", ylabel = "Mass [Msun]",dpi = 200,ylim=(0,M_DM*1.1))
display(p2)
display(p1)
## GW
full_out = zeros(Float64,length(times),8)
full_out[:,1:6] .= sol_interpolated'
full_out[:,7].= ones(length(times))
full_out[:,8] = times
waveform_out=Pipeline.package_quadrupole_frequency_new(full_out)
df_waveform=DataFrame(waveform_out, ["t", "f", "ddI_p", "ddI_c", "a", "M_DM", "M_NS1", "M_NS2", "Phase"])

## plot

freq(t,f1) = f1 / ((1-t/t_decay)^(3/8))

p1 = plot(df_waveform.t,df_waveform.f,xlabel = "time [s]", ylabel = "frequency [Hz]", dpi = 200, ylim = (0,500),xlim = (t_decay - 1e1,t_decay+1), label = "simulated")
plot!(p1,df_waveform.t,freq.(df_waveform.t,df_waveform.f[1]), label = "expected", ls=:dash)
vline!(p1,[t_decay],label = "decay time")
display(p1)
# p2 = plot(df_waveform.Phase,df_waveform.ddI_p,xlabel = "phase", ylabel = "ddI", label = "+", dpi = 300)
# plot!(p2, df_waveform.Phase,df_waveform.ddI_c,xlabel = "phase", ylabel = "ddI", label = "x", dpi = 300)
#display(p2)
