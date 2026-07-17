include("dependencies.jl")

## Set parameters
n = 5
K = 3.25e10
rho_0 = Math.cgs_density_inv(9e13)
#M_DM = 0.1
M_NS1 = 1.4 #Solar masses
R_NS1 = 12  #km
M_NS2 = 1.4 #Solar masses
a_start = 200. #km

#Initialize Polytrope
realpoly = Polytrope.solve_halo(n,K,M_NS1,R_NS1,rho_0,2000) #realpoly = [rs,R_DM,rho_r,M_DM_r,rho_interp,mass_interp,M_DM]

#check plot of polytrope
p1 = plot(realpoly.rs, realpoly.rho_r, xlabel = "Radius [km]", ylabel = "Density [Msun/km^3]", xlim = (0,realpoly.R_DM))
display(p1)

M_DM = realpoly.M_DM
println("Mass of DM halo: ",M_DM," Solar Masses")
println("Mass of halo outside NS: ",M_DM - realpoly.mass_interp(R_NS1), " Solar Masses")

circ_orbit = Integrate.initial_circular_orbit(M_DM, M_NS1, M_NS2, a_start, realpoly);  #Solve for Roche seperation, ciruclar orbit   circ_orbit = [a_start, M_DM, M_NS1, M_NS2, R_DM, J, a_RL, 0.]

ICs = [circ_orbit[7],circ_orbit[2],circ_orbit[3],circ_orbit[4],0,circ_orbit[6]];    #ICs = [a, M_DM, M_NS1, M_NS2, phase, J]
parameters = [realpoly.mass_interp];    #[mass(r)]
println("Dark matter outflow beginning at r = ", circ_orbit[7], " km")

# Not necessary, calculate approx decay time for binary nuetron stars
M1 = M_DM + M_NS1
η = (M1*M_NS2)/((M1+M_NS2)^2)
t_decay = 5/256*((Math.c)^5) * (circ_orbit[7]^4) /(Math.G*(M1+M_NS2))^3/ η
println("Decay time: ", t_decay, " seconds")
## Integration
tend = t_decay * 2
integration_sol=Integrate.integrate_halo(ICs,parameters,(0,tend));      #Run time integration of system

# Intepolate the integration solution at regular timesteps for plotting
t_a_min = integration_sol.t[end]    #solver quits once minimum seperation is reached, extract final time
times = range(0,stop=t_a_min,length=1000)
sol_interpolated=integration_sol(times)

#plotting seperation and mass
M_DM_t = sol_interpolated[2,:]
a_t = sol_interpolated[1,:]

p2 = plot(times,a_t, xlabel = "Time [s]", ylabel = "seperation [km]",dpi = 200,ylim=(0,circ_orbit[1]*1.1))
p3 = plot(times,M_DM_t, xlabel = "Time [s]", ylabel = "Mass [Msun]",dpi = 200,ylim=(0,M_DM*1.1))
display(p2)
display(p3)
## GW
#Initialize the variables to be passed into quadrupole package
full_out = zeros(Float64,length(times),8)
full_out[:,1:6] .= sol_interpolated'
full_out[:,7].= ones(length(times)) #leftover from R_WD
full_out[:,8] = times

#Find quadrupole evolution, calculated waveform for interpolated time values
waveform_out=Pipeline.package_quadrupole_frequency_new(full_out)
df_waveform=DataFrame(waveform_out, ["t", "f", "ddI_p", "ddI_c", "a", "M_DM", "M_NS1", "M_NS2", "Phase"])

## plot
#freq vs time
p4 = plot(df_waveform.t,df_waveform.f,xlabel = "time [s]", ylabel = "frequency [Hz]", dpi = 200, label = "simulated")
vline!(p4,[t_decay],label = "decay time")
# frequency_t(t,f1) = f1./(1 .-(t./t_decay))^(3/8)
# plot!(p4, df_waveform.t,frequency_t(df_waveform.f[1],df_waveform.t),label = "analytical")
display(p4)

# waveform vs phase
p5 = plot(df_waveform.Phase,df_waveform.ddI_p,xlabel = "phase", ylabel = "ddI", label = "+", dpi = 300)
plot!(p5, df_waveform.Phase,df_waveform.ddI_c,xlabel = "phase", ylabel = "ddI", label = "x", dpi = 300)
#display(p5)

## Save to CSV
CSV.write("Merger_wDM.csv", df_waveform)