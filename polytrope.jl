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
    end

    """Nondimensionalize Polytrope solution"""
    function apply_polytrope(CachedPoly, M_DM,rho_0,n)
        xi1 = CachedPoly.root       #outer radius of DM halo
        xi_range = range(0,xi1,length = 1000)     #array of dimensionless radii
        M_tot_nondim = CachedPoly.mass(xi1)                         #total dimensionless mass
        thetas = Polytrope.polytrope_profile.(Ref(CachedPoly),xi_range)
        thetas = reduce(hcat,thetas)[1,:]

        K = (1/(4*pi*rho_0)*(M_DM/M_tot_nondim))^(2/3) *4*pi*G*rho_0^(1-(1/n)) / (n+1) #calculate polytrope constant
        alpha = ((n+1)*K/(4*pi*G*rho_0^(1-(1/n))))^(1/2)   #alpha is nondim. constant, easier to use
        
        #redimensionalize
        rs = xi_range.*alpha #an array of distances in km
        R_DM = xi1*alpha  #outer edge of halo in km
        rho_r = rho_0.*((thetas).^n)  #density profile
        mass_r = 4*pi* alpha^3 * rho_0 * CachedPoly.mass(xi_range)  #mass profile

        return RealPolytrope(rs,R_DM,rho_r,mass_r)

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
    
        scaled_R_RL = R_RL * root / R_DM #fraction of radius within white dwarf


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