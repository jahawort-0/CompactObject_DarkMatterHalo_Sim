module Polytrope
    using OrdinaryDiffEq
    using QuadGK
    using Interpolations
    using CSV, DataFrames
    using Memoization
    using ChainRulesCore

    const G::Float64 = 1.327e11   # km^3 / M⊙ / s^2

    
    """ Returns a function f(xi, theta) where theta = [θ, θ′]
        #Computes Lane–Emden 2nd-order equation in first-order form."""
    function polytrope_a(n)
    	function f!(dθ, θ, p, xi)
    		# θ = [theta, dtheta/dxi]
    		# dθ = [d theta/dxi, d² theta/dxi²]
    		# Lane-Emden: (1/xi^2) d/dxi(xi^2 θ′) = -θ^n
    		θ_0=max(θ[1],0) #avoid negative exponentiation
    		a = -2 * θ[2] / xi - θ_0^n
    		dθ[1] = θ[2]
    		dθ[2] = a
    	end
    	return f!
    end
    """continuous approximation for dimensionless polytrope density"""
    function lane_emden_int(xi, sol)
    	if xi > sol.t[end]
    		return 0.0
    	end
    	θ = sol(xi)       # sol(x) gives interpolated vector [θ, θ′]
    	return max(θ[1],0)
    end

    """Analytical Solution to n=5 Lane–Emden eqn"""
    function LEn5_sol(ξ)
        θ = (1+(ξ^2)/3)^(-1/2)
        dθ = (-ξ/3)*(1+(ξ^2)/3)^(-3/2)
        return[θ,dθ]
    end

    """ dm/dxi = xi^2 * θ^n  (standard Lane–Emden mass integrand)"""
    function dm_dxi(xi, sol, n)
    	θ = lane_emden_int(xi, sol)
    	return xi^2 * θ^n
    end

    function mass_enclosed(sol, root, n; N=2000)
        rs = range(0, root; length=N)   #sample points to surface
        f = [dm_dxi(r, sol, n) for r in rs] #evalute integrand at each r
        ms = zeros(length(rs))  #initialize enclosed mass

        for i in 2:length(rs)   #ms[1]=M(r=0) = 0, start at ms[2]
            dr = rs[i] - rs[i-1]    #radius step
            ms[i] = ms[i-1] + (f[i] + f[i-1])/2 * dr    #evaluate enclosed mass
        end
        # Continuous interpolation
        M_enclosed = interpolate((rs,), ms, Gridded(Linear()))
        return M_enclosed
    end


    """Used to find xi1, where the polytrope reaches zero value, corresponds to surface of DM halo """
    # may be a possible issue with oscillating solutions (n=1), finding wrong root, nonphysical radius
    function find_roots(ftn;low=0,high=20,iterations=20)
    	soln=(low+high)/2   #average between high and low
    	if iterations==0    #kills the loop after n iterations
    		return(soln)
    	end
    	if ftn(soln)>0      #moves lower bound
    		low=soln
    	elseif ftn(soln)<0  #moves upper bound
    		high=soln
    	else
    		return soln     #if zero is exactly found
    	end
    	return(find_roots(ftn;low=low,high=high,iterations=iterations-1))   #loops function for n iterations
    end

    #stores the solution of polytropic equations of state so that they do not 
    #need to be solved each time that they are used
    struct CachedPolytrope
        sol   :: Any      #Differential equation solution table
        root  :: Float64  #scale factor if the surface of the white dwarf (for normalization)
        tnorm :: Float64  #solution to the mass integral (for normalization)
        mass  :: Any  #mass enclosed as a function of radius
    end

    #stores all CashedPolytrope types previously calculated
    const polytrope_cache = Dict{Float64, CachedPolytrope}()

    """Calculates the solution to a polytropic differential equation and its normalized radius and mass """
    function compute_polytrope(n)
        f! = polytrope_a(n)

        if n == 5
            θ0 = [1.0, 0.0]
            tspan = (0.001, 100)
            prob = ODEProblem(f!, θ0, tspan, nothing)
            sol = solve(prob; reltol=1e-9, abstol=1e-9, dtmax=1e-3)
    
            # The n=5 polytrope never goes to zero, however its area does converge to sqrt(3)
            # Taking xi=100 as the edge of the polytrop encapsulates 99.95% of the total mass
            root = 100
    
            # normalization integral ∫ xi^2 θ^n dxi
            tnorm = quadgk(x -> dm_dxi(x, sol, n), 0, root)[1]

            mass = mass_enclosed(sol,root,n);
    
            return CachedPolytrope(sol, root, tnorm, mass)

        else
            θ0 = [1.0, 0.0]
            tspan = (0.001, 20)
            prob = ODEProblem(f!, θ0, tspan, nothing)
            sol = solve(prob; reltol=1e-9, abstol=1e-9, dtmax=1e-3)
    
            # root-finding for θ=0
            root = find_roots(xi -> sol(xi)[1])
    
            # normalization integral ∫ xi^2 θ^n dxi
            tnorm = quadgk(x -> dm_dxi(x, sol, n), 0, root)[1]

            mass = mass_enclosed(sol,root,n);
    
            return CachedPolytrope(sol, root, tnorm, mass)
        end
    end

    function polytrope_profile(output,xi)
        xi1 = output.root
        if xi > xi1
            return [0,0]
        else
            return output.sol(xi)
        end
    end

    struct RealPolytrope
        rs   :: Vector{Float64}
        R_DM  :: Float64 
        rho_r :: Vector{Float64}  
        mass_r  :: Vector{Float64}  
        rho_interp :: Any
        mass_interp :: Any
        M_DM :: Float64
    end

    """Nondimensionalize Polytrope solution"""
    function apply_polytrope(CachedPoly, M_DM,M_NS1,K,n)
        xi1 = CachedPoly.root       #outer radius of DM halo
        xi_range = range(0,xi1,length = 1000)     #array of dimensionless radii
        M_tot_nondim = CachedPoly.mass(xi1)                         #total dimensionless mass
        if n==5
            M_tot_nondim = sqrt(3)
        end
        thetas = Polytrope.polytrope_profile.(Ref(CachedPoly),xi_range)
        thetas = reduce(hcat,thetas)[1,:]

        rho_0 = (K*(n+1)/(4*pi*G)*((4*pi*M_tot_nondim/M_DM)^(2/3)))^(1/((1/3)-(1/n)))
        alpha = ((n+1)*K/(4*pi*G*rho_0^(1-(1/n))))^(1/2)   #alpha is nondim. constant, easier to use
        
        #redimensionalize
        rs = xi_range.*alpha #an array of distances in km
        R_DM = xi1*alpha  #outer edge of halo in km
        rho_r = rho_0.*((thetas).^n)  #density profile
        mass_r = 4*pi* alpha^3 * rho_0 * CachedPoly.mass(xi_range)  #mass profile
        #mass_r = mass_r .+ M_NS1    #add in central point mass

        rho_interp = linear_interpolation(rs, rho_r, extrapolation_bc=Flat())
        mass_interp = linear_interpolation(rs, mass_r, extrapolation_bc=Flat())

        df = DataFrame([rs,rho_r,mass_r], ["rs", "rho_r", "mass_r"])

        CSV.write("Polytrope_sol.csv", df)

        return RealPolytrope(rs,R_DM,rho_r,mass_r,rho_interp,mass_interp,M_DM)
    end

    function n5apply_polytrope(M_DM::Float64,K::Float64)   #temporary
        M_tot_nondim = sqrt(3)
        n = 5

        rho_0 = (K*(n+1)/(4*pi*G)*((4*pi*M_tot_nondim/M_DM)^(2/3)))^(1/((1/3)-(1/n)))
        alpha = ((n+1)*K/(4*pi*G*rho_0^(1-(1/n))))^(1/2)   #alpha is nondim. constant, easier to use
        
        rho_interp(r) = rho_0 * (LEn5_sol(r/alpha)[1])^n
        mass_interp(r) = 4*pi*alpha^3 * rho_0 * (r/alpha)^2 * abs(LEn5_sol(r/alpha)[2])

        return RealPolytrope([],NaN,[],[],rho_interp,mass_interp,M_DM)
    end

    struct RealPolytrope_wK
        rs   :: Vector{Float64}
        R_DM  :: Float64 
        rho_r :: Vector{Float64}  
        mass_r  :: Vector{Float64}  
        K :: Float64
        rho_interp :: Any
        mass_interp :: Any
        M_DM :: Float64
    end

    #For a given central density
    function apply_polytrope_rho(CachedPoly, M_DM,M_NS1,rho_0,n)
        xi1 = CachedPoly.root       #outer radius of DM halo
        xi_range = range(0,xi1,length = 1000)     #array of dimensionless radii
        M_tot_nondim = CachedPoly.mass(xi1)                         #total dimensionless mass
        if n==5
            M_tot_nondim = sqrt(3)
        end
        thetas = Polytrope.polytrope_profile.(Ref(CachedPoly),xi_range)
        thetas = reduce(hcat,thetas)[1,:]

        K = (1/(4*pi*rho_0)*(M_DM/M_tot_nondim))^(2/3) *4*pi*G*rho_0^(1-(1/n)) / (n+1) #calculate polytrope constant
        alpha = ((n+1)*K/(4*pi*G*rho_0^(1-(1/n))))^(1/2)   #alpha is nondim. constant, easier to use
        
        #redimensionalize
        rs = xi_range.*alpha #an array of distances in km
        R_DM = xi1*alpha  #outer edge of halo in km
        rho_r = rho_0.*((thetas).^n)  #density profile
        mass_r = 4*pi* alpha^3 * rho_0 * CachedPoly.mass(xi_range)  #mass profile
        #mass_r = mass_r .+ M_NS1    #add in central point mass

        rho_interp = linear_interpolation(rs, rho_r, extrapolation_bc=Flat())
        mass_interp = linear_interpolation(rs, mass_r, extrapolation_bc=Flat())

        df = DataFrame([rs,rho_r,mass_r], ["rs", "rho_r", "mass_r"])

        CSV.write("Polytrope_sol.csv", df)

        return RealPolytrope_wK(rs,R_DM,rho_r,mass_r,K,rho_interp,mass_interp,M_DM)

    end

        #Numerical solution to Polytrope with Neutron star at center
    function halo_polytrope_problem_inside!(du,u,p,r)
        # u = [ρDM, M, M_DM]
        #du = [dρ, dM,dM_DM]
        # p = [n,K,M_NS1,R_NS1,ρDM_0]
        # r = radius, variable of integration
        #setup inputs
        ρDM = u[1]
        M = u[2]
        M_DM = u[3]
        n = p[1]
        K = p[2]
        M_NS1 = p[3]
        R_NS1 = p[4]

        #Make useful terms
        Γ = 1 + (1/n)
        ρNS = M_NS1 / (4/3*pi*R_NS1^3)
        ρ = ρDM + ρNS

        if ρDM<= 0  #"ending" the integration if density drops to zero
            u[1] = 0    #set DM density to zero
            du[1] = 0   #no change in DM density
            du[2] = 4*pi * r^2 * ρNS    #total enclosed mass comes only from NS
            du[3] = 0   #no change in DM mass
        else    #normal integration
        #The second order ODE
        dρ = -G/(K*Γ) * M / r^2 * ρDM / ρDM^(Γ-1)
        dM = 4*pi * r^2 * ρ
        dM_DM = 4*pi * r^2 * ρDM

        #outputs
        du[1] = dρ
        du[2] = dM
        du[3] = dM_DM
        end
    end

    function halo_polytrope_problem_outside!(du,u,p,r)
        # u = [ρDM, M]
        #du = [dρ, dM]
        # p = [n,K,M_NS1,R_NS1,ρDM_0]
        # r = radius, variable of integration
        #setup inputs
        ρDM = u[1]
        M = u[2]
        n = p[1]
        K = p[2]
        M_NS1 = p[3]
        R_NS1 = p[4]

        #Make useful terms
        Γ = 1 + (1/n)
        ρNS = 0
        ρ = ρDM + ρNS

        if ρDM<= 0  #"ending" the integration if density drops to zero
            u[1] = 0    #set DM density to zero
            du[1] = 0   #no change in DM density
            du[2] = 0   #total enclosed mass 
            du[3] = 0   #no change in DM mass
        else    #normal integration

        #The second order ODE
        dρ = -G/(K*Γ) * M / r^2 * ρDM / ρDM^(Γ-1)
        dM = 4*pi * r^2 * ρ

        #outputs
        du[1] = dρ
        du[2] = dM
        du[3] = dM
        end
    end

    function halo(p,rspan)
        # p = [n,K,M_NS1,R_NS1,ρDM_0]
        #rspan: radius over which we integrate

        rspan_in = (rspan[1],p[4])
        rspan_out = (p[4],rspan[end])

        M_NS1 = p[3]
        R_NS1 = p[4]
        ρDM_0 = p[5]
        NS_density = M_NS1 / (4/3*pi*R_NS1^3)
        M0 = 4π/3 * rspan[1]^3 * (NS_density + ρDM_0)
        M0_DM = 4π/3 * rspan[1]^3 * ρDM_0


        #u0 = [ρ0, dρ0 = 0] Setup Initial conditions
        u0 = [ρDM_0, M0,M0_DM]
    
        prob1 = ODEProblem(halo_polytrope_problem_inside!, u0, rspan_in, p)
        sol1 = solve(prob1, abstol = 1e-12, reltol = 1e-12)

        u0_out = sol1.u[end]
        prob2 = ODEProblem(halo_polytrope_problem_outside!,u0_out,rspan_out,p)
        sol2 = solve(prob2, abstol = 1e-12, reltol = 1e-12, dtmax = 0.01)

        return(sol1, sol2)
    end

    function solve_halo(n,K,M_NS1,R_NS1,ρDM_0,rend)
        p = [n,K,M_NS1,R_NS1,ρDM_0]
        rspan = (1e-6,rend)

        sol_in,sol_out = halo(p,rspan)

        rs_in = LinRange(0,R_NS1, 1000)
        rs_out = vcat(LinRange(R_NS1,R_NS1+100, 500),LinRange(R_NS1+100,2000, 1000)[2:end])
        rs = vcat(rs_in,rs_out[2:end])

        sol_in_interp = sol_in(rs_in)
        sol_out_interp = sol_out(rs_out)
        
        rho_r_in = sol_in_interp[1,:]
        rho_r_out = sol_out_interp[1,:]
        rho_r = vcat(rho_r_in,rho_r_out[2:end])

        M_r_in = sol_in_interp[2,:]
        M_r_out = sol_out_interp[2,:]
        M_r = vcat(M_r_in,M_r_out[2:end])

        M_DM_r_in = sol_in_interp[3,:]
        M_DM_r_out = sol_out_interp[3,:]
        M_DM_r = vcat(M_DM_r_in,M_DM_r_out[2:end])

        rho_interp = linear_interpolation(rs, rho_r, extrapolation_bc=Flat())
        mass_interp = linear_interpolation(rs, M_DM_r, extrapolation_bc=Flat())

        #solve for total mass   #not quite working yet
            # M_in = M_DM_r_in[end]
            # Γ = 1 + (1/n)
            # rho_surf = rho_interp(R_NS1)
            #     u0_back = [rho_r_out[1],M_r_out[1],M_DM_r_out[1]]
            #     rspan_back = (R_NS1,0)
            #     prob3 = ODEProblem(halo_polytrope_problem_outside!,u0_back,rspan_back,p)
            #     sol3 = solve(prob3, abstol = 1e-12, reltol = 1e-12)
            # rho0 = maximum(sol3(0))
            # α = ((n+1)*K/(4*pi*G*rho0^(1-(1/n))))^(1/2)
            # rho_prime = (rho_interp(R_NS1+0.00001)-rho_interp(R_NS1)/0.00001)
            # M_out = 4*pi*rho0^((93*Γ-4)/2) * (K*Γ/(4*pi*G*(Γ-1)))^(3/2) * (sqrt(3) - ((R_NS1/α)^2 *abs((Γ-1) * (rho_surf/rho0)^(Γ-2) * α / rho0 * rho_prime)))
            # Mtot = M_in + M_out

        #M_DM = NaN #Mtot

        frac = 0.99
        # examine change in enclose mass
        dM = diff(M_DM_r_out)
        plateau = findfirst(abs.(dM) .< 1e-12)  #will return spot where enclosed mass plateaus

        if plateau !== nothing  #if plateau is in the range we are searching
            R_DM = rs_out[plateau]
        else        #there may not be a platueau within the region we search, take the 95% mass boundary
            M_target = frac * M_DM_r_out[end]
            radius_mass_interp = linear_interpolation(M_DM_r_out, rs_out)
            R_DM = radius_mass_interp(M_target)
        end

        M_DM = mass_interp(R_DM)

        ind = findfirst(==(1),rs.>R_DM)  #clear halo outside 'max radius'
        rho_r = rho_r[1:ind];       rho_r[end] = 0
        M_DM_r = M_DM_r[1:ind];     M_DM_r[end] = M_DM
        rs = rs[1:ind]

        rho_r = rho_r[1:end-1]
        M_DM_r = M_DM_r[1:end-1]
        rs = rs[1:end-1]

        rho_interp = linear_interpolation(rs, rho_r, extrapolation_bc=Flat())
        mass_interp = linear_interpolation(rs, M_DM_r, extrapolation_bc=Flat())

        df = DataFrame([rs,rho_r,M_DM_r], ["rs", "rho_r", "mass_r"])

        CSV.write("Polytrope_sol.csv", df)

        return RealPolytrope(rs,R_DM,rho_r,M_DM_r,rho_interp,mass_interp,M_DM)
    end

    """Recalls polytrope if previously calculated, otherwise calculates and cashes the polytrope before returning it"""
    function get_polytrope(n)
        if haskey(polytrope_cache, n) #returns if previously calculated
            return polytrope_cache[n]
        end
        #calculates and cashes if necessary
        P = compute_polytrope(n) 
        polytrope_cache[n] = P
        return P
    end

    """Calculates the mass of dark matter shell outside of its Roche Lobe """    
    function mass_outside_radius(R_DM, R_RL, M_DM, n; ROUND=2)
        if ROUND!=-1
            n=round(n,digits=ROUND) 
            #so that the expensive calculation of a polytrope solution is not needed each time
            #in case n changes between time steps
        end
        P = get_polytrope(n)
        sol = P.sol
        root = P.root
        t_norm = P.tnorm
        m_table = P.mass
    
        #scaled_R_RL = R_RL * root / R_DM #fraction of radius within white dwarf


        # #integrates the fraction of the mass outside of the roche lobe and scales to the given mass
        # if scaled_R_RL<0
        #     val=t_norm
        # elseif scaled_R_RL>root
        #     val = 0
        # else
        #     val = m_table(scaled_R_RL/root)
        # end
    
        # return M_DM * val / t_norm

        M_outside = M_DM - mass_r(R_RL)
        return M_outside
    
        
        
    end

end