module Pipeline
    using Distributed
    using OrdinaryDiffEq
    using FLoops
    using DistributedArrays
    using CSV, DataFrames, Interpolations
    using Printf
    include("math.jl")
    include("save.jl")
    include("Integrate.jl")
    include("Integrate_r.jl")



    """The following function is meant to integrate from start to finish a system with given M_WD and M_NS. It is not always as flexible as it needs to be, and it does not actually take optional keywords, but rather these are hardcoded for time being. In particular, this pipeline works well for Jeans mode, but struggles with Isotropic Re-emission (which is the default option). This is why there are more flexible methods included below for specific purposes."""
    function integrate_r_full(M_WD::Float64, M_NS::Float64;
        optional=Integrate_r.optional_default)
    
        base=zeros(7)
        base[1:6].=Integrate.initial_circular_orbit(M_WD,M_NS=M_NS)
        out_minus=Integrate_r.integrate_r_sneaky(base,-1,res=1e-7)
    
        new_base=base.+out_minus[end-5,:]
        out_turnaround=turn_around_dtheta_dt(new_base[1:6],out_minus[end,7]-out_minus[end-5,7])
       
    
        newer_base=new_base.+out_turnaround[end,:]
        println("NEWER BASE", newer_base)
        newer_base[7]=0
        out_plus=Integrate_r.integrate_r_sneaky(newer_base,+1,res=0.1)
        #out_plus=Integrate_r.integrate_r_sneaky(new_base,+1,res=0.1)
        
        #=flag=true
        while flag
            out_plus=Integrate_r.integrate_r_sneaky(new_base,+1,res=0.1)
            if length(out_plus[:,7])<=1 || out_plus[end,7]==0
                new_base[1]=new_base[1]+0.01
            else
                flag=false
            end
        end=#
    
        
        waveform_out=zeros(length(out_minus[:,7])+length(out_plus[:,7]),8)
        
        waveform_out[1:length(out_minus[:,7]),1].=out_minus[:,7]
        waveform_out[length(out_minus[:,7])+1:end,1].=out_plus[:,7].+out_minus[end,7]
        
        for i in 1:length(out_minus[:,7])
            t_val=base[1:end-1].+out_minus[i,1:end-1]
            t_val[4]=Math.R0_WD(t_val[2])
            t_val[6]=0
            t_ddI=Math.ddI_from_theta(t_val)
            waveform_out[i,2]=abs(t_ddI[1])
            waveform_out[i,3]=abs(t_ddI[2])
            waveform_out[i,4]=4*pi/Math.period(t_val[1],t_val[2]+t_val[3])
            waveform_out[i,5:7].=out_minus[i,1:3].+base[1:3]
            waveform_out[i,8]=out_minus[i,6]
        end
    
        for i in 1:length(out_plus[:,7])
            t_val=new_base[1:end-1].+out_plus[i,1:end-1]
            t_val[4]=Math.R0_WD(t_val[2])
            t_val[6]=0
            t_ddI=Math.ddI_from_theta(t_val)
            waveform_out[length(out_minus[:,7])+i,2]=abs(t_ddI[1])
            waveform_out[length(out_minus[:,7])+i,3]=abs(t_ddI[2])
            waveform_out[length(out_minus[:,7])+i,4]=4*pi/Math.period(t_val[1],t_val[2]+t_val[3])
            waveform_out[length(out_minus[:,7])+i,5:7].=out_plus[i,1:3].+new_base[1:3]
            waveform_out[length(out_minus[:,7])+i,8]=out_plus[i,6]+new_base[6]
        end
        
        return(waveform_out)
    
    end

    function dtheta_dt_sneaky(theta_init::Vector{Float64})
            function ftn!(du::Vector{Float64}, delta_theta::Vector{Float64},p::Vector{Any},t::Float64; 
            parallel::Bool=true, uneven::Bool=false, optional=optional_default,evolve_R_WD=false)
                current_theta=theta_init.+delta_theta
                if evolve_R_WD
                    #println(current_theta[4])
                     Math.acceleration!(du, current_theta,p,t; parallel=parallel, uneven=uneven, optional=optional)
                else
                    #println(current_theta[2])
                    current_theta[4]=Math.R0_WD(current_theta[2])
                    Math.acceleration!(du, current_theta,p,t; parallel=parallel, uneven=uneven, optional=optional)
                    du[4]=0
                end
                return(du)
            end
            return ftn!
        end

    """ Normal parallelization implementation of in-place acceleration! function""" 
    function dtheta_dt_even_sneaky!(du::Vector{Float64}, theta_init::Vector{Float64}, delta_theta::Vector{Float64},p::Vector{Any},t::Float64; optional=optional_default,evolve_R_WD=evolve_R_WD)
        (dtheta_dt_sneaky(theta_init))(du, delta_theta, p, t; optional=optional,evolve_R_WD=evolve_R_WD)
    end
    
        
    """ Serial implementation of in-place acceleration! function"""
    function dtheta_dt_serial_sneaky!(du::Vector{Float64}, theta_init::Vector{Float64}, delta_theta::Vector{Float64},p::Vector{Any},t::Float64; optional=optional_default,evolve_R_WD=false)
        (dtheta_dt_sneaky(theta_init))(du, delta_theta, p, t; parallel=false, optional=optional,evolve_R_WD=evolve_R_WD)
    end
    
    function integration_dt_sneaky(theta_init,
                                dt::Float64, t_end::Float64;
                                delta=zeros(7),
                                optional=Integrate_r.optional_default,
                                 parallel::Bool=true,
                                 uneven::Bool=false,
                                 alg=Vern9(),
                                evolve_R_WD=false)
            #println(theta_init)
            t_acceleration! = nothing
            redirect_stderr(devnull) do
                if uneven
                    throw("We don't do that here.") 
                else
                    parallel ?
                        t_acceleration! = (du, θ, p, t) ->
                            dtheta_dt_even_sneaky!(du, theta_init, θ, p, t; optional=optional,evolve_R_WD=evolve_R_WD) :
                        t_acceleration! = (du, θ, p, t) ->
                            dtheta_dt_serial_sneaky!(du, theta_init, θ, p, t; optional=optional,evolve_R_WD=evolve_R_WD)
                end
            end
    
            prob = ODEProblem(t_acceleration!, zeros(6), (0, t_end), Any[]) 
            #defines the problem based on starting position (delta) and step count
            sol=solve(prob,alg,dt=dt,adaptive=false,save_everystep=true)
    
            return(sol.t,sol.u)
    end

    """Carries out the turnaround stage in  time basis given a starting condition turn_around"""
    function turn_around_dtheta_dt(turn_around,test_t_end)
        local ts, thetas
        try
            ts,thetas=integration_dt_sneaky(turn_around,test_t_end/1000,5*test_t_end;evolve_R_WD=true)
            #e.g., ts is a vector with length ~1001 and thetas has ~1001 rows each of length 6
        catch
            #println("EXCEPTION")
            ts,thetas=integration_dt_sneaky(turn_around,test_t_end/100,5*test_t_end;evolve_R_WD=false)
        end
        thetas = reduce(hcat, thetas)'
        
        if thetas[end,1]>thetas[end-1,1]
            #println("No Recursion")
            out=zeros(length(ts),7)
            out[:,1:6].=thetas
            out[:,7].=ts
            return(out)
        end
    
    
        #println("Recursion")
        thetas_2=turn_around_dtheta_dt(turn_around.+thetas[end,:],test_t_end)
    
        out=zeros(length(ts)+length(thetas_2[:,1])-1,7)
        out[1:length(ts),1:6].=thetas
        out[1:length(ts),7].=ts
    
        out[length(ts)+1:end,1:6].=(thetas_2[2:end,1:6].+thetas[end,:]')
        out[length(ts)+1:end,7].=(thetas_2[2:end,7].+ts[end])
        return(out)
    end

    """Function more specifically optimized for isotropic re-emission"""
    function integrate_r_full_new(M_WD::Float64, M_NS::Float64;
        optional=Integrate_r.optional_default)
    
        base=zeros(7)
        base[1:6].=Integrate.initial_circular_orbit(M_WD,M_NS=M_NS)
        out_minus=(Integrate_r.integrate_r_sneaky(base,-1,res=1e-7).+reshape(base, 1, :))
        
    
        new_base=out_minus[end-5,:]
        new_base[4]=Math.R0_WD(new_base[2])
        println(new_base)
        #return(out_minus)
        
        T=out_minus[end,7]-out_minus[end-5,7]
        println("Second Step")
        out_turnaround=(turn_around_dtheta_dt(new_base[1:6],T)).+reshape(new_base, 1, :)
    
        newer_base=out_turnaround[end,:]
        newer_base[4]=Math.R0_WD(newer_base[2])
        println("Third Step")
        println(newer_base)
        out_plus=turn_around_and_out(newer_base,T)
    
        newest_base=out_plus[end,:]
        newest_base[4]=Math.R0_WD(newest_base[2])
        dt=(out_plus[end,7]-out_plus[end-1,7])/10
        println("Final Step")
        out_end=final_stage_integrate(newest_base[1:6],dt)
        out_end[:,7].+=newest_base[7]
    
        out_full=vcat(out_minus[1:end-6,:],out_turnaround)
        out_full=vcat(out_full[1:end-1,:],out_plus)
        out_full=vcat(out_full[1:end-1,:],out_end)
        return(out_full)
    end


    """Forces integration until R0_WD(0.05) for UCXB-like systems under Isotropic Re-emission assumption """
    function integrate_r_full_new_terminal(M_WD::Float64, M_NS::Float64;
        optional=Integrate_r.optional_default)
    
        base=zeros(7)
        base[1:6].=Integrate.initial_circular_orbit(M_WD,M_NS=M_NS)
        out_minus=(Integrate_r.integrate_r_sneaky(base,-1,res=1e-7).+reshape(base, 1, :))
    
        if Math.R0_WD(out_minus[end,2])>out_minus[end,1] || merger_check(out_minus) || out_minus[end,2]<0.05
            return(out_minus)
        end
    
        println(size(out_minus))
    
        new_base=out_minus[end-5,:]
        new_base[4]=Math.R0_WD(new_base[2])
        println(new_base)
        
        
        
        T=out_minus[end,7]-out_minus[end-5,7]
        println("Second Step")
        out_turnaround=(turn_around_dtheta_dt(new_base[1:6],T)).+reshape(new_base, 1, :)
    
        out_full=vcat(out_minus[1:end-6,:],out_turnaround)
        if Math.R0_WD(out_full[end,2])>out_full[end,1] || merger_check(out_full) || out_full[end,2]<0.05
            return(out_full)
        end
        
        println(size(out_full))
        
        newer_base=out_turnaround[end,:]
        newer_base[4]=Math.R0_WD(newer_base[2])
        println("Third Step")
        println(newer_base)
        out_plus=turn_around_and_out2(newer_base,T)
    
        
    
        println("Saved out_plus")
        println(size(out_plus))
        out_full=vcat(out_full[1:end-1,:],out_plus)
        print("Concatenated out_plus")
        if Math.R0_WD(out_full[end,2])>out_full[end,1] || merger_check(out_full) || out_full[end,2]<0.05
            return(out_full)
        end
    
        println(size(out_full))
        
        newest_base=out_plus[end,:]
        newest_base[4]=Math.R0_WD(newest_base[2])
        dt=(out_plus[end,7]-out_plus[end-1,7])
        println("Final Step")
        out_end=final_stage_integrate_terminal(newest_base[1:6],100*dt)
        out_end[:,7].+=newest_base[7]
        out_full=vcat(out_full[1:end-1,:],out_end)
    
    
        if Math.R0_WD(out_full[end,2])>out_full[end,1] || merger_check(out_full) || out_full[end,2]<0.05
            return(out_full)
        end
    
        try
            out_epilogue=turn_around_and_out_terminate(out_end[end,:],1.0e5,dr=1)
            out_full = vcat(out_full[1:end-1,:], out_epilogue)
            
            return(out_full)
        catch
            return(out_full)
        end
    end

    """Given the output of one of the above pipelines, formats the data in a CSV for human consumption. """
    function package_quadrupole_frequency(full_out_test)
        out=zeros(length(full_out_test[:,1]),8)
        out[:,1].=full_out_test[:,7]
        out[:,5:7].=full_out_test[:,1:3]
        out[:,8].=full_out_test[:,6]
        for i in 1:length(full_out_test[:,1])
            out[i,2]=2/Math.period(full_out_test[i,1],full_out_test[i,2]+full_out_test[i,3])
            zero_theta=zeros(6)
            zero_theta[1:5].=full_out_test[i,1:5]
            t_ddI=Math.ddI_dddI_from_theta(zero_theta; ddI_only=true)
            out[i,3:4].=t_ddI[2][1:2]
        end
        return(out)
    end

    """Saves output as a .csv file. NOTE THE SPECIFIC FORMATTING USED--you may want to change this based on your system."""
    function save_thing(full_out_test3,M_WD,M_NS)
        df=DataFrame(package_quadrupole_frequency(full_out_test3), ["t", "f", "ddI_p", "ddI_c", "a", "M_WD", "M_NS", "Phase"])
        label=@sprintf("Outputs_6/R_int_%.2f_%.2f_iso_re.csv",M_WD,M_NS )
        CSV.write(label, df)#return here
        println(label)
    end

    function turn_around_and_out2(newer_base,test_t_end;loop_limit=10,dr=0.1)

        i = 0
        flag = true
    
        t_out_plus = nothing
        turn_around_collection = nothing
    
        while i < loop_limit && flag
    
            i += 1
    
            t_out_plus =
                Integrate_r.integrate_r_sneaky_end(newer_base,+1;dr=dr)
    
            if size(t_out_plus,1) <= 1 || t_out_plus[end] == 0
    
                turn_around_int =
                    turn_around_dtheta_dt(newer_base[1:6],test_t_end)
    
                t_shifted_out =
                    turn_around_int .+ reshape(newer_base,1,:)
    
                if isnothing(turn_around_collection)
                    turn_around_collection = t_shifted_out
                else
                    turn_around_collection =
                        vcat(
                            turn_around_collection[1:end-1,:],
                            t_shifted_out
                        )
                end
    
                newer_base .+= turn_around_int[end,:]
    
                test_t_end *= 2
    
            else
                flag = false
            end
        end
    
        if isnothing(t_out_plus)
            error("Loop never executed; t_out_plus was never assigned.")
        end
    
        shifted_out =
            Array(t_out_plus) .+ reshape(newer_base,1,:)
    
        if i > 1
            return vcat(
                turn_around_collection[1:end-1,:],
                shifted_out
            )
        end
    
        return shifted_out
    end



    function turn_around_and_out(newer_base,test_t_end;loop_limit=10,dr=0.1)
        i=0
        flag=true
        turn_around_collection=0
        
        while i<loop_limit && flag
            i+=1
            println(i)
            t_out_plus=Integrate_r.integrate_r_sneaky_end(newer_base,+1;dr=dr)
            #println(t_out_plus)
            if length(t_out_plus[:,1])<=1 || t_out_plus[end]==0
                turn_around_int=turn_around_dtheta_dt(newer_base[1:6],test_t_end)
                #print(turn_around_int)
                t_shifted_out=turn_around_int.+ reshape(newer_base, 1, :)
                try
                    turn_around_collection=vcat(turn_around_collection[1:end-1,:],t_shifted_out)
                catch
                    turn_around_collection=t_shifted_out
                    #@show typeof(reshape(newer_base, 1, :))
                    #@show typeof(turn_around_int)
                    #@show typeof(t_shifted_out)
                    #@show typeof(turn_around_collection)
                    #println(turn_around_collection[1:10,:])
                    #println(turn_around_int[1:10,:])
                    #println((turn_around_int.+newer_base')[1:10,:])
                end
                newer_base.+=turn_around_int[end,:]
                println(newer_base)
                test_t_end*=2
                
            else
                flag = false;dr
                if i>1
                    #println(turn_around_collection)
                    println(size(turn_around_collection))
                    println(size(t_out_plus))
                    println(size(t_out_plus.+newer_base'))
                    shifted_out = Array(t_out_plus) .+ reshape(newer_base, 1, :)
                    #@show typeof(t_out_plus)
                    #@show typeof(shifted_out)
                    #@show typeof(turn_around_collection)
                    return (vcat(turn_around_collection[1:end-1,:],shifted_out))
                end
                return Array(t_out_plus) .+ reshape(newer_base, 1, :)
            end
        end
    end

    function turn_around_and_out_terminate(newer_base,test_t_end;loop_limit=10,dr=0.1)
        i=0
        flag=true
        turn_around_collection=0
    
        R=Math.R0_WD(0.05)-newer_base[1]
        
        t_out_plus=Integrate_r.integrate_r_specific(newer_base,+1, R ;dr=dr)     
                
        return Array(t_out_plus)'.+ reshape(newer_base, :)'
        
    end








end