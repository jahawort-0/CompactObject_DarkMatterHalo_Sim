module Integrate
    using Distributed
    using OrdinaryDiffEq
    using FLoops
    using DistributedArrays
    include("math.jl")
    include("save.jl")
    
    """
        Returns initial conditions for the beginning of Roche Lobe overflow as a vector.
        Vector output is [a, M_WD, M_NS, R_WD, J, θ], where:
        a = total orbital separation of the binary.
        M_WD = mass of the white dwarf.
        M_NS = mass of the neutron star.
        R_WD = radius of the white dwarf.
        J = total angular momentum of the binary.
        θ = angular parameter.
        ```julia
        initial_circular_orbit(M_WD::Float64;M_NS::Float64=1.35,mu_e::Float64=2.)
        ```
        mu_e = mean molecular mass per electron (default 2 for white dwarfs of all but the most exotic chemical compositions (e.g., ONeMG)).
    """


    #const M_ch=1.44 #Chandresekhar mass in M⊙
    const optional_default = (mu_e = 2.0, equation = 1, A = 10.0, poly_shell = false, a_ring_frac = 1.0, mode = :isotropic, eta_acc = 0.1, period_calc = true)
    #`optional` contains optional parameters. In order, these optional parameters are:
    #mu_e=2.0:                     the chemical composition of the white dwarf (for CO/ONeMg, mu_e=2.0/2.35, respectively)
    #equation=1:                   choice of empirical equation for white dwarf mass vs. radius relation (1=Nauenberg/2=Eggleton
    #A=10.0:                       the rate at which mass flows off the WD surface (characteristic number of orbits for mass to leave)
    #poly_shell=true:              if true/false, uses polytropic/uniform density distribution to model mass in process of being stripped from WD
    #a_ring_frac=1.0:              radius of circumbinary ring formed by matter decreted by WD
    #mode=:isotropic: mode of mass loss assumed in calculating angular momentum loss
    #eta_acc=0.1:                  accretion efficiency of the NS
    #period_calc=true:             if true/false, calculates angular velocity as a function of radius/angular momentum

    #Of these, I expect A, poly_shell, mode, and period_calc to have considerable effect on the evolution of the system

    function initial_circular_orbit(M_NS1::Float64, M_NS2::Float64, M_DM::Float64, RealPoly::Any; optional=optional_default)::Vector{Float64}
    
        #returns keplerian orbital parameters, theta, at initial Roche Lobe overflow

        a=Math.RL_contact(M_NS1,M_NS2,M_DM,RealPoly;optional=optional)
        J=Math.circular_J(a,(M_NS1 + M_DM),M_NS2)
        R_DM=RealPoly.R_DM
        return([a, M_NS1, M_NS2, M_DM, R_DM, J, 0.])
    end
    
    #"Below is the heart of the integration process"
    """Integrates theta forward one step with time interval either given explicitly (dt!=-1) or as a fraction of the orbital period determined by dt_period. Returns the absolute time step taken and the updated theta parameters."""

    


    #NEW CODE--given optional parameters, returns an integrator which can then be reused at each new timestep
    function make_one_step(optional=optional_default; parallel::Bool=true, uneven::Bool=false)
        #First we set up the integration problem once which will then be used every subsequent time when the temporary function (ftn!) is called
        prob = nothing

        t_acceleration! = nothing
    
        # There's a warning from SciMLBase here:
        #   Warning: Using arrays or dicts to store parameters of different types can hurt performance.
        #   Consider using tuples instead.
        #   @ SciMLBase ~/.julia/packages/SciMLBase/MzuF2/src/performance_warnings.jl:33
        # but we can't use tuples because ForwardDiff uses only arrays. So we suppress all warnings temporarily here.
        redirect_stderr(devnull) do 
             if uneven #define the acceleration problem based on parallelization choices
                t_acceleration! = Math.acceleration_uneven!
            else
                parallel ?
                    t_acceleration! = (du, θ, p, t) -> Math.acceleration_even!(du, θ, p, t; optional=optional) :
                    t_acceleration! = (du, θ, p, t) -> Math.acceleration_serial!(du, θ, p, t; optional=optional)
            end
        end

        #integ = init(prob, alg=Vern9()) #We build the integrator once and will reuse it every time
        
        function ftn!(theta::Vector{Float64}; dt_period::Float64=1e-2, dt::Float64=-1.)::Tuple{Float64, Vector{Float64}}
            @assert length(theta)==6
            #integrates forward one step with time interval either given explicitly (dt!=-1)
            #or as a fraction of the orbital period determined by dt_period
            if dt==-1
                @inbounds begin
                    a = theta[1]
                    M_WD = theta[2]
                    M_NS = theta[3]
                    dt=dt_period*Math.period(a, M_WD + M_NS)
                end
            end

            #resetting initial values for prob to current values
            # Construct a tiny problem
            prob  = ODEProblem(t_acceleration!, theta, (0., dt), Any[])
            integ = init(prob, Vern9(); save_everystep=false, save_start=false, save_end=false, tstops=Float64[])
    
            step!(integ)
            theta .= integ.u        
            
            return integ.t, copy(theta)
        end
    
        return ftn!
        
    end

    #NEWER code--tries to eliminate the problem of making a new integrator every step to avoid unnecessary memory allocations
   function make_one_step_new(theta, optional=optional_default; parallel::Bool=true, uneven::Bool=false, time_limit=3.16e7)
        #First we set up the integration problem once which will then be used every subsequent time when the temporary function (ftn!) is called
        prob = nothing

        t_acceleration! = nothing
    
        # There's a warning from SciMLBase here:
        #   Warning: Using arrays or dicts to store parameters of different types can hurt performance.
        #   Consider using tuples instead.
        #   @ SciMLBase ~/.julia/packages/SciMLBase/MzuF2/src/performance_warnings.jl:33
        # but we can't use tuples because ForwardDiff uses only arrays. So we suppress all warnings temporarily here.
        redirect_stderr(devnull) do 
             if uneven #define the acceleration problem based on parallelization choices
                t_acceleration! = Math.acceleration_uneven!
            else
                parallel ?
                    t_acceleration! = (du, θ, p, t) -> Math.acceleration_even!(du, θ, p, t; optional=optional) :
                    t_acceleration! = (du, θ, p, t) -> Math.acceleration_serial!(du, θ, p, t; optional=optional)
            end
        end

        #integ = init(prob, alg=Vern9()) #We build the integrator once and will reuse it every time


        @assert length(theta)==6
            #integrates forward one step with time interval either given explicitly (dt!=-1)
            #or as a fraction of the orbital period determined by dt_period
            if dt==-1
                @inbounds begin
                    a = theta[1]
                    M_WD = theta[2]
                    M_NS = theta[3]
                    dt=dt_period*Math.period(a, M_WD + M_NS)
                end
            end

            #resetting initial values for prob to current values
            # Construct a tiny problem
            prob  = ODEProblem(t_acceleration!, theta, (0., time_limit), Any[])
            integ = init(prob, Vern9(); save_everystep=false, save_start=false, save_end=false, tstops=Float64[])
        
        function ftn!(theta::Vector{Float64}; dt_period::Float64=1e-2, dt::Float64=-1.)::Tuple{Float64, Vector{Float64}}
            integ.dt = 0.01 
        
            step!(integ)
            theta .= integ.u        
            
            return integ.t, copy(theta)
        end
    
        return ftn!
        
    end 


    #NEWEST code--resets the initial conditions of the integrator for each time step while controlling size of the time step

    function make_one_step_reuse(optional=optional_default;
                             parallel::Bool=true,
                             uneven::Bool=false,
                             alg=Vern9())

        # There's a warning from SciMLBase here:
        #   Warning: Using arrays or dicts to store parameters of different types can hurt performance.
        #   Consider using tuples instead.
        #   @ SciMLBase ~/.julia/packages/SciMLBase/MzuF2/src/performance_warnings.jl:33
        # but we can't use tuples because ForwardDiff uses only arrays. So we suppress all warnings temporarily here.
        t_acceleration! = nothing
        redirect_stderr(devnull) do
            if uneven
                t_acceleration! = Math.acceleration_uneven!
            else
                parallel ?
                    t_acceleration! = (du, θ, p, t) ->
                        Math.acceleration_even!(du, θ, p, t; optional=optional) :
                    t_acceleration! = (du, θ, p, t) ->
                        Math.acceleration_serial!(du, θ, p, t; optional=optional)
            end
        end
    
        # Dummy initial state (will be overwritten by reinit!)
        θ0 = zeros(6)
    
        prob = ODEProblem(t_acceleration!, θ0, (0.0, Inf), Any[])
        integ = init(prob, alg;
                     save_everystep=false,
                     save_start=false,
                     save_end=false)
    
        function ftn!(theta::Vector{Float64};
                      dt_period::Float64=1e-2,
                      dt::Float64=-1.0)::Tuple{Float64,Vector{Float64}}
    
            @assert length(theta) == 6
    
            # Determines time step as a fraction of the physical step
            if dt == -1.0
                @inbounds begin
                    a     = theta[1]
                    M_WD  = theta[2]
                    M_NS  = theta[3]
                    dt = dt_period * Math.period(a, M_WD + M_NS)
                end
            end
    
            # Reset integrator state 
            #reinit!(integ; u0=theta, t0=0.0) 
            reinit!(integ, theta; t0=0.0, reset_dt=true)
            #no explicit time-dependent variable here 
            #(i.e., time translation invariant) 
            #so it is okay to start over at t=0 
    
            # Force integration to exactly dt
            empty!(integ.opts.tstops)
            add_tstop!(integ, dt)
    
            advance_to_tstop!(integ)
    
            theta .= integ.u
            return dt, copy(theta)
        end
    
        return ftn!
    end


    #EXTRA MOST NEWESTEST 
    #set up a problem and integrate it for a
    #set number of steps with set size
    function solve_one_orbit(theta::Vector{Float64};
                             t_init=0,
                             optional=optional_default,
                             parallel::Bool=true,
                             uneven::Bool=false,
                             alg=Vern9(),
                             steps=100)

        # There's a warning from SciMLBase here:
        #   Warning: Using arrays or dicts to store parameters of different types can hurt performance.
        #   Consider using tuples instead.
        #   @ SciMLBase ~/.julia/packages/SciMLBase/MzuF2/src/performance_warnings.jl:33
        # but we can't use tuples because ForwardDiff uses only arrays. So we suppress all warnings temporarily here.
        t_acceleration! = nothing
        redirect_stderr(devnull) do
            if uneven
                t_acceleration! = Math.acceleration_uneven!
            else
                parallel ?
                    t_acceleration! = (du, θ, p, t) ->
                        Math.acceleration_even!(du, θ, p, t; optional=optional) :
                    t_acceleration! = (du, θ, p, t) ->
                        Math.acceleration_serial!(du, θ, p, t; optional=optional)
                end
            end

        a = theta[1]
        M_WD = theta[2]
        M_NS = theta[3]
        P=Math.period(a, M_WD + M_NS)
        dt=P/steps
        
    
        prob = ODEProblem(t_acceleration!, theta, (t_init, t_init+P), Any[])
        sol=solve(prob,alg,dt=dt,adaptive=false,save_everystep=true)
    
        return(sol.t,sol.u)
    end

 function integrate_chunk(theta;
                            orbit_limit::Float64=1e5, step_limit::Int=10^7,
                                mass_limit::Float64=1e-2, dt_period=1e-2,
                                parallel::Bool=true, uneven::Bool=false,
                                write_file::Bool=true, filename::Union{String,Int}=-1,
                                optional=optional_default,
                                checkpoint_every::Int=500,
                                verbose::Bool=false,
                                time_limit=3.16e7,
                                theta_init::Vector{Float64}=[-1.0],
                                steps=100,
                                alg=Vern9())

        t=0
        i=1

        out=zeros(step_limit,7)

    
        while theta[6] < 2*π*orbit_limit && i < step_limit-steps && theta[2] > mass_limit && t<time_limit && theta[1]> 6 * Math.G * theta[3]/Math.c^2
            new_t,new_thetas=solve_one_orbit(theta; 
                             t_init=t,
                             optional=optional,
                             parallel=parallel,
                             uneven=uneven,
                             alg=alg,
                             steps=steps)
            println(new_t[end],new_thetas[end])
            for j in 1:steps 
                    #save the output--
                    #overwrite the last step of the last integration since it is the same as the first in this one 
                    out[i , 1:6] .= new_thetas[j]
                    i+=1 
                    out[i , 7] = new_t[j]
            end

            t=new_t[end]
            theta=out[i,:]
        
        end

        return(out[1:i,:])
    end

    
    """Given an initial set of masses, `integrate` solves for the equations of motion from initial Roche Lobe overflow until the eventual complete dissipation of the white dwarf or the integration finishes according to user defined cutoff (either a limit on the number of orbits completed or on the integration steps taken). `integrate` stores and returns the second derivative of the quadrupole moment as a function of time for the orbit. The function also saves the run to a file named by the date and initial conditions if `write_file`. """
    function integrate(M_WD::Float64=1.0; M_NS::Float64=1.35, theta_init::Vector{Float64}=[-1.0], mu_e::Float64=2., 
            orbit_limit::Float64=1e5, step_limit::Int=10^7, mass_limit::Float64=1e-3, dt_period=1e-2,
            parallel::Bool=true, uneven::Bool=false, write_file::Bool=true, verbose=false, filename::Union{String,Int}=-1,
        optional=optional_default,time_limit=3.16e7)::AbstractMatrix{<:Real}
        #returns time series second derivatives of the quadrupole moment
        @assert M_WD > mass_limit
        @assert M_WD < Math.M_ch # Chandrasekhar mass
        #@assert M_NS > 1.1 # very low lower limit, slightly below mass of the secondary neutron star in PSR J0453+1559
        #@assert M_NS < 2.9 # TOV limit (a generous estimate) 
        #Disabling asserts so that primordial/stellar mass black holes can be tested
        @assert orbit_limit > 0
        @assert step_limit > 1
        @assert mass_limit > 0
        @assert dt_period > 0

        one_step! = make_one_step(optional; parallel=parallel, uneven=uneven) # NEW 
    
        t=0.
        if theta_init[1]==-1.0
            theta=initial_circular_orbit(M_WD,M_NS=M_NS,optional=optional)
        else
            @assert length(theta_init)==6
            theta=theta_init
        end
        ddI=Math.ddI_from_theta(theta,optional=optional)
        datumses=zeros(step_limit+1,9)
        datumses[1,:].=[t, ddI[1,1], ddI[1,2], theta... ]
        i=1


        
    
        @inbounds begin
            phase = theta[6]
            M_WD = theta[2]
            while phase<2*π*orbit_limit && i<step_limit && M_WD>mass_limit
                #integrate until time or step limit are elapsed or WD completely dissipates
                
                #dt,theta=one_step!(theta; dt_period=dt_period, parallel=parallel, uneven=uneven, optional=optional) OLD
                dt,theta=one_step!(theta; dt_period=dt_period) #NEW
            
                ddI=Math.ddI_from_theta(theta; optional=optional)
                i+=1
                t+=dt
                datumses[i,:].=[t, ddI[1,1], ddI[1,2], theta... ]
                if mod(i,10)==0 && verbose
                    println(i," ",  theta)
                end
            end
        end
    
        datumses=datumses[1:i,:]
    
        if write_file #saves datumses to a csv in directory "Output"
            Save.save_as_csv(datumses,M_WD,M_NS,optional.A,optional.mode,filename=filename)
        end
        
        return(datumses)
    end


    """Same as above but with checkpoints """
    function integrate_checkpointed(M_WD::Float64=1.0; M_NS::Float64=1.35, mu_e::Float64=2.,
                                orbit_limit::Float64=1e5, step_limit::Int=10^7,
                                mass_limit::Float64=1e-2, dt_period=1e-2,
                                parallel::Bool=true, uneven::Bool=false,
                                write_file::Bool=true, filename::Union{String,Int}=-1,
                                optional=optional_default,
                                checkpoint_every::Int=500,
                                verbose::Bool=false,
                                time_limit=3.16e7,
                                theta_init::Vector{Float64}=[-1.0])

        if verbose
            println("I am using the new version.")
        end
        @assert M_WD > mass_limit
        @assert M_WD < Math.M_ch
    
        #one_step! = make_one_step(optional; parallel=parallel, uneven=uneven) # NEW  # new -- defines the problem once here 


        if theta_init[1]==-1.0
            theta=initial_circular_orbit(M_WD,M_NS=M_NS,optional=optional)
        else
            @assert length(theta_init)==6
            theta=theta_init
        end
    
        if filename==-1
            filename=Save.format_filename(M_WD,M_NS,optional.A,optional.mode)
        end
        
        t = 0.
        theta = initial_circular_orbit(M_WD, M_NS=M_NS, optional=optional)
        ddI = Math.ddI_from_theta(theta, optional=optional)

        one_step! = make_one_step_reuse(optional; parallel=parallel, uneven=uneven) 
        # NEWER--implements improved memory allocation technique
    
        # buffer for checkpointing
        buffer = zeros(checkpoint_every, 9)
        buffer_count = 1
        buffer[1,:] .= [t, ddI[1,1], ddI[1,2], theta...]
        
        i = 1
       
        rss_mb=0
        
        while theta[6] < 2*π*orbit_limit && i < step_limit && theta[2] > mass_limit && rss_mb < 30000 && t<time_limit
            
            dt, theta = one_step!(theta; dt_period=dt_period) #old
        
            #t_dt= dt_period*Math.period(theta[1], theta[2] + theta[3])#NEWER
        
            #dt,theta=one_step!(theta, dt_period=dt_period) #new--calls the function defined internally
            #dt,theta=one_step!(theta, dt_period=dt_period) #newer--calls the function defined internally
            
            ddI = Math.ddI_from_theta(theta, optional=optional)
            t += dt
            i += 1
    
            buffer_count += 1
            buffer[buffer_count, :] .= [t, ddI[1,1], ddI[1,2], theta...]
            
    
            if buffer_count >= checkpoint_every
                Save.save_as_csv_append(buffer, M_WD, M_NS, optional.A, optional.mode; filename=filename)
                buffer .= 0.0          # clear buffer
                buffer_count = 0       # reset counter
                #GC.gc()                # release unused memory
            end
        
            if mod(i, 100) == 0
                rss_mb = parse(Int, read(`ps -o rss= -p $(getpid())`, String))/1024
                if verbose
                    println("Step $i, θ = $theta, Memory usage: $(rss_mb) MB")
                end
            end
    
            
        end
    
        # Save any leftover rows in the buffer
        if buffer_count > 0 && write_file
            Save.save_as_csv_append(buffer[1:buffer_count, :], M_WD, M_NS, optional.A, optional.mode; filename=filename)
        end
    
        return nothing  
    end

    """Performs many WD-NS merger simulations at once. (This is another place where we take advantage of parallelization). Takes in sim_conditions (NamedTuple), which has:
    - M_WDs (array/tuple of white dwarf mass(es))
    - M_NSs (array/tuple of neutron star mass(es))
    and optional_conditions (NamedTuple), which has:
    - mu_es (array/tuple of mu_e(s), the mean molecular mass per electron)
    - WD_mass_radius_equations ((1), (2), or (1,2); specifies whether to use the Nauenberg WD mass-radius relation (1), Eggleton relation (2), or both (1,2))
    - As (array/tuple of A value(s), a constant related to mass loss of the white dwarf)
    - use_polytropic_density ((true), (false), or (true,false); specifies whether to use the polytropic density distribution to model mass in process of being stripped from WD (true), the uniform density distribution (false), or both (true,false))
    - a_ring_fracs (array/tuple of a_ring_frac(s), the radius of circumbinary ring formed by matter decreted by WD. Only relevant if 'modes' includes "Circumbinary Ring")
    - mass_loss_modes (array/tuple of mass loss mode specifications for calculating angular momentum loss; possible values include :jeans, :isotropic, and :circumbinary. :jeans = Jeans mode, :isotropic = isotropic re-emission, :circumbinary = circumbinary ring)
    - eta_accs (array/tuple of eta_acc(s), the accretion efficiency of the neutron star)
    - use_period_calc ((true), (false), or (true,false); specifies whether to calculate the angular velocity of white dwarf as a function of the period (true), the angular momentum (false), or both (true,false))
    many_integrate then iterates through all possible combinations of these user-specified parameters, with one integrate run dedicated to each combination. Any parameters not specified in optional_conditions will fall back to the default values, which are:
    - mu_es = (2.0)
    - WD_mass_radius_equations = (1)
    - As = (10.0)
    - use_polytropic_density = (true)
    - a_ring_fracs = (1.0)
    - mass_loss_modes = (:isotropic)
    - eta_accs = (0.1)
    - use_period_calc = (true)
    The ranges of M_WD and M_NS must be physical. Some of the checks we have are:
    @assert M_WD < Math.M_ch # Chandrasekhar mass
    @assert M_NS > M_WD
    @assert M_NS > 1.1 # very low lower limit, slightly below mass of the secondary neutron star in PSR J0453+1559
    @assert M_NS < 2.9 # TOV limit (a generous estimate)"""
    function many_integrate(sim_conditions::NamedTuple; mu_e::Float64=2., 
            orbit_limit::Float64=1e5, step_limit::Int=10^7, mass_limit::Float64=1e-3, dt_period=1e-3,
            parallel::Bool=true, distributed::Bool=true, uneven::Bool=false, write_file::Bool=true, filename::Union{String,Int}=-1, optional_conditions::NamedTuple=NamedTuple())
        pairs = initialize_pairs(sim_conditions.M_WDs, sim_conditions.M_NSs, optional_conditions)
        if parallel
            if distributed
                pairs_distr = distribute(pairs)
                @distributed (+) for pair in pairs_distr
                    M_WD, M_NS, mu_e, equation, A, poly_shell, a_ring_frac, mode, eta_acc, period_calc = pair
                    optional = (mu_e=mu_e, equation=equation, A=A, poly_shell=poly_shell, a_ring_frac=a_ring_frac, mode=mode, eta_acc=eta_acc, period_calc=period_calc)
                    integrate(M_WD; M_NS, mu_e, 
                        orbit_limit, step_limit, mass_limit, dt_period,
                        parallel, uneven, write_file, filename, optional)
                    GC.gc() # Free up memory
                    0 # trivial return value for reduction, needed to time distributed memory calls
                end    
                # Need to be careful here because apparently having each process read and write to disk can worsen performance if not done carefully
            else
                @floop ThreadedEx(basesize = 8) for pair in pairs
                    M_WD, M_NS, mu_e, equation, A, poly_shell, a_ring_frac, mode, eta_acc, period_calc = pair
                    optional = (mu_e=mu_e, equation=equation, A=A, poly_shell=poly_shell, a_ring_frac=a_ring_frac, mode=mode, eta_acc=eta_acc, period_calc=period_calc)
                    integrate(M_WD; M_NS, mu_e, 
                        orbit_limit, step_limit, mass_limit, dt_period,
                        parallel, uneven, write_file, filename, optional)
                    GC.gc()
                end
            end
        else
            for pair in pairs
                M_WD, M_NS, mu_e, equation, A, poly_shell, a_ring_frac, mode, eta_acc, period_calc = pair
                optional = (mu_e=mu_e, equation=equation, A=A, poly_shell=poly_shell, a_ring_frac=a_ring_frac, mode=mode, eta_acc=eta_acc, period_calc=period_calc)
                integrate(M_WD; M_NS, mu_e, 
                    orbit_limit, step_limit, mass_limit, dt_period,
                    parallel, uneven, write_file, filename, optional)
                GC.gc() # Free up memory
            end
        end
    end

    """
        Helper function for distributed memory parallelization mode. Makes distributed pairs for all inputted parameter arguments specified by the user, including if the user specified to iterate over the optional parameter arguments.
    """
    function initialize_pairs(M_WDs, M_NSs, optional_conditions)
        mu_es = hasproperty(optional_conditions, :mu_es) ? optional_conditions.mu_es : (optional_default[1],)
        equations = hasproperty(optional_conditions, :WD_mass_radius_equations) ? optional_conditions.WD_mass_radius_equations : (optional_default[2],)
        As = hasproperty(optional_conditions, :As) ? optional_conditions.As : (optional_default[3],)
        poly_shells = hasproperty(optional_conditions, :use_polytropic_density) ? optional_conditions.use_polytropic_density : (optional_default[4],)
        a_ring_fracs = hasproperty(optional_conditions, :a_ring_fracs) ? optional_conditions.a_ring_fracs : (optional_default[5],)
        modes = hasproperty(optional_conditions, :mass_loss_modes) ? optional_conditions.mass_loss_modes : (optional_default[6],)
        eta_accs = hasproperty(optional_conditions, :eta_accs) ? optional_conditions.eta_accs : (optional_default[7],)
        period_calcs = hasproperty(optional_conditions, :use_period_calc) ? optional_conditions.use_period_calc : (optional_default[8],)

        # Apologies for the syntax below, it isn't elegant but it makes the code run
        pairs = [(M_WD, M_NS, mu_e, equation, A, poly_shell, a_ring_frac, mode, eta_acc, period_calc) for M_WD in M_WDs for M_NS in M_NSs for mu_e in mu_es for equation in equations for A in As for poly_shell in poly_shells for a_ring_frac in a_ring_fracs for mode in modes for eta_acc in eta_accs for period_calc in period_calcs]
        return pairs
    end

    ## New DM halo functions

    function update_mass_r(mass_r)

    end

    function halo_problem!(du,u,p,t)
        #--------Setup inputs ---------

        #du = [dM_DM_dt, da_dt]
        # u = [M_DM, a]
        # p = [M_NS1, M_NS2, mass_r]
        M_DM = u[1]     
        a = u[2]        

        M_NS1 = p[1]    
        M_NS2 = p[2]    
        mass_r = p[3]

        #--------Setup the model of the system, makes many assumptions ---------

        R_RL = Math.Roche_Limit(a, (M_NS1 + M_DM), M_NS2)

        Mtot = M_NS1 + M_NS2 + M_DM
        mu = (M_NS1 + M_DM)*M_NS1/Mtot  #reduced mass

        P = Math.period(a, Mtot)

        #dM_DM_dt = -1e-200
        dM_DM_dt = -10/P * (M_DM - mass_r(R_RL))  #will need A = 10
        #This works under the assumption that the change in mass is monotonic
        #For more complicated mass transfer (second halo forming) we need a more complex treatment
        
        da_rr_dt = -64 * Math.G^3 * mu * Mtot^2 / (5*Math.c^5*a^3)

        dM_NS2_dt = 0 #Simple assumption that removed mass will leave system

        beta = abs(dM_NS2_dt/dM_DM_dt)
        gamma = (M_NS1+M_DM)/M_NS2  #isotropic reemission

        da_mt_dt = -2*a * dM_DM_dt/(M_NS1+M_DM) * (1- beta*(M_DM+M_NS1)/M_NS2 -
            (1-beta)*(gamma+0.5) * (M_NS1+M_DM)/Mtot)

        da_dt = da_rr_dt + da_mt_dt

        #dJ_dt = 0  #No change in angular momentum, not a very good assumption

        #---------Define outputs ----------
        du[1] = dM_DM_dt
        du[2] = da_dt
    end

    function integrate_halo(u0,p,tspan)
        #u0: Initial conditions
        #p: system parameters
        #tspan: time over which we integrate
    
        prob = ODEProblem(halo_problem!, u0, tspan, p)
        sol = solve(prob)

        return(sol)
    end

end