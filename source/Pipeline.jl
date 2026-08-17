module Pipeline
    using Distributed
    using OrdinaryDiffEq
    using FLoops
    using DistributedArrays
    using CSV, DataFrames, Interpolations
    using Printf
    using FFTW
    include("math.jl")
    include("save.jl")
    include("Integrate.jl")
    include("Integrate_r.jl")
    include("polytrope.jl")

    #New Function for DM
    function package_quadrupole_frequency_new(full_out_test)
        #full_out_test [a, M_DM,M_NS1,M_NS2, phase, J, R_DM, t]
        out=zeros(length(full_out_test[:,1]),9)
        out[:,1].=full_out_test[:,8]
        out[:,5:8].=full_out_test[:,1:4]  # a, Masses  
        out[:,9].=full_out_test[:,5]    #phase
        for i in 1:length(full_out_test[:,1])
            out[i,2]=2/Math.period(full_out_test[i,1],full_out_test[i,2]+full_out_test[i,3]+full_out_test[i,4]) #frequency
            zero_theta=zeros(7)
            #theta = [a, M_DM, M_NS1, M_NS2, R_WD, J, phase]
            zero_theta[1]=full_out_test[i,1]    #a
            zero_theta[2]=full_out_test[i,2]    #M_DM
            zero_theta[3]=full_out_test[i,3]    #M_NS1
            zero_theta[4]=full_out_test[i,4]    #M_NS2
            zero_theta[5]=full_out_test[i,7]    #R_WD??
            zero_theta[6]=full_out_test[i,6]    #J
            zero_theta[7]=full_out_test[i,5]    #phase
            t_ddI=Math.ddI_dddI_from_theta(zero_theta; ddI_only=true)
            out[i,3:4].=t_ddI[2][1:2]
        end
        return(out)
        #out["t", "f", "ddI_p", "ddI_c", "a", "M_DM", "M_NS1", "M_NS2, "Phase"]
    end

    function evolve_halo(n::Float64 ,K::Float64 ,rho_0::Float64 ,R_NS1::Float64 ,M_NS1::Float64 ,M_NS2::Float64 ,
        filename::String)

        #Initialize Polytrope
        println("Initializing dark matter halo...")
        realpoly = Polytrope.solve_halo(n,K,M_NS1,R_NS1,rho_0,2000) #realpoly = [rs,R_DM,rho_r,M_DM_r,rho_interp,mass_interp,M_DM]
        M_DM = realpoly.M_DM

        #Solve for Roche seperation, ciruclar orbit
        println("Setting up the system...")
        circ_orbit = Integrate.initial_circular_orbit(M_DM, M_NS1, M_NS2, realpoly);  #circ_orbit = [a_RL, M_DM, M_NS1, M_NS2, R_DM, J, 0.]

        #Setup Initial Conditions for integration
        ICs = [circ_orbit[1],circ_orbit[2],circ_orbit[3],circ_orbit[4],0,circ_orbit[6]];    #ICs = [a, M_DM, M_NS1, M_NS2, phase, J]
        parameters = [realpoly.mass_interp];    #[mass(r)]

        #Print Status
        println("Total dark matter halo mass: M_DM = ", M_DM, " Solar masses")
        println("Dark matter outflow beginning at seperation r = ", circ_orbit[1], " km")

        #Find time interval to integrate over
        M1 = M_DM + M_NS1;      η = (M1*M_NS2)/((M1+M_NS2)^2)
        t_decay = 5/256*((Math.c)^5) * (circ_orbit[1]^4) /(Math.G*(M1+M_NS2))^3/ η  #decay time for a system w/o mass transfer
        tend = t_decay * 1.2    #system will decay slower with mass transfer, tend must extend beyond true end of system behavior
        println("Max simulation time: ",tend," seconds")

        #Run time integration of system
        println("Evolving the system over time...")
        integration_sol=Integrate.integrate_halo(ICs,parameters,(0,tend));

        # Intepolate the integration solution at regular timesteps for plotting
        t_a_min = integration_sol.t[end]    #solver quits once minimum seperation is reached, extract final time
        #times = vcat(range(0,stop=(t_a_min*0.95),length=800),range((t_a_min*0.95),t_a_min,length=201)[2:end])
        #times = range(0,t_a_min,length=Int(floor(100*2*t_a_min)))        
        #times = range(0,t_a_min,length=10000)

        times = Float64[]   #Variable timestep interpolation times, important for 
            push!(times,0.0)
            P = Math.period(circ_orbit[1], M1+M_NS2)
            while times[end]<t_a_min
                dt = (P/4.139579320)*(1-(times[end]/t_a_min))^(3/8)
                push!(times,times[end]+dt)
            end
            push!(times,t_a_min)

        sol_interpolated=integration_sol(times)

        #Initialize the variables to be passed into quadrupole package
        full_out = zeros(Float64,length(times),8)
        full_out[:,1:6] .= sol_interpolated'
        full_out[:,7].= ones(length(times)) #leftover from R_WD
        full_out[:,8] = times

        #Find quadrupole evolution, calculated waveform for interpolated time values
        println("Finding gravitational waves...")
        waveform_out=package_quadrupole_frequency_new(full_out)
        df_waveform=DataFrame(waveform_out, ["t", "f", "ddI_p", "ddI_c", "a", "M_DM", "M_NS1", "M_NS2", "Phase"])
        CSV.write(filename, df_waveform)

        #Print end result of integration
        println(M_DM - df_waveform.M_DM[end]," Solar masses of dark matter were removed in ", df_waveform.t[end], " seconds.")
    end

    #Calculate GW strain 
    function calc_strain(ddI,DL)
        return (2*Math.G/(Math.c^4*DL)).*ddI
    end

    #FFT GW strain
    # function FT_strain(h_t,times)
    #     FT_h = fft(h_t);              #FFT of strain
    #     dt = times[2] - times[1];          # time step
    #     fs = 1/dt;                    # sampling frequency
    #     N = length(h_t);

    #     freqs = (0:N-1) * (fs/N);
    #     half = 1:div(N,2)
    #     freqs = freqs[half]
    #     FT_h = FT_h[half]
    #     FT_h = real.(FT_h)
    #     return(freqs,FT_h)
    # end
    function FT_strain(h_t,phase_t, freq,times)     #calculates h tilde (without complex phase)
        A = abs.(h_t./(-cos.(2 .*phase_t)))
        fdot = similar(freq)
        fdot[2:end-1] = (freq[3:end] .- freq[1:end-2]) ./ (times[3:end] .- times[1:end-2])
        fdot[1] = fdot[2]
        fdot[end] = fdot[end-1]

        return(A,A .* sqrt.((2*pi)./abs.(fdot)))
    end

    function htilde_phase(phase,freq,times)
        return(phase .- ((2*pi) .* freq .* times .+ phase[1]))
    end

    function def_integral_df(h1,h2,freq)
        h = conj(h1) .* h2
        sum=0
        for i = 1:length(freq)-1
            df = freq[i+1]-freq[i]
            area = df * ((h[i]+h[i+1])/2)
            sum+=area
        end
        return sum
    end


end