module Math
	using ForwardDiff #package used for chain rules of large derivatives
    using LinearAlgebra
	using Memoization
    using Base.Threads

    # Constants
    const c::Float64 = 2.998e5    # km/s
    const G::Float64 = 1.327e11   # km^3 / M⊙ / s^2
    const period_factor::Float64 = sqrt(4.0 * π^2 / G)
    const dJ_rr_factor::Float64 = -2.0 * G / (5.0 * c^5)
    const Mpc_to_km_factor = 2.0*G*3.2e-20/c^4. #for conversion from Mpc to km
    const M_ch=1.44 #Chandresekhar mass in M⊙

    include("polytrope.jl")

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
    # A and mode are physically meaningful parameters while poly_shell and period_calc are simply alternate ways to calculate the same values
    #Default values of `true` on these latter parameters are the preferred means of calculation

    const optional_default = (mu_e = 2.0, equation = 1, A = 10.0, poly_shell = false, a_ring_frac = 1.0, mode = :isotropic, eta_acc = 0.1, period_calc = true)

    """
        Returns the natural radius of a white dwarf.
        ```julia
        R0_WD(M_WD)
        ```
        M_WD = mass of the white dwarf.
    """
    function R0_WD(M_WD;optional=optional_default) 
        #white dwarf radius in km as ftn of mass in M_odot and mean molecular weight/electron=2
    
        #assert comments can be uncommented for debugging, but they prove cumbersome in our integration.
        #for now, all `@assert` statements are commented out.
    	#@assert M_WD > 0 
    	#@assert mu_e > 0 
    
        
        x=M_WD/M_ch
        y=M_WD/0.00057
        #see Even (2009; https://arxiv.org/pdf/0908.2116) for summary
        if optional.equation==1 #Nauenberg equation
        	return(15584.0/optional.mu_e * x^(-1.0/3.0) * (1-x^(4.0/3.0))^(1.0/2.0))
        else #Eggleton equation
            return(7931.0*sqrt(x^(-2.0/3.0)-x^(2.0/3.0)) * (1.0+3.5*y^(-2.0/3.0)+1/y)^(-2.0/3.0)) 
        end
    end
    
    """
        Returns the dynamical timescale of the white dwarf.
        ```julia
        dynamical_timescale(M_WD,R_WD)
        ```
        M_WD = mass of the white dwarf.
        R_WD = radius of the white dwarf.
    """
    function dynamical_timescale(M_WD,R_WD)
        #returns dynamical timescale for perturbed white dwarf
    	#with radius R and mass M
    
    	#@assert M_WD > 0 
    	#@assert R_WD > 0 
    
    	return((R_WD^3.0/(G*M_WD))^(1.0/2.0))
    end
    
    """
        Returns the Roche limit of the binary, using the Eggleton approximation.
        ```julia
        Roche_Limit(a,M_WD,M_NS)
        ```
        a = total orbital separation of the binary.
        M_WD = mass of the white dwarf.
        M_NS = mass of the neutron star.
    """
    function Roche_Limit(a,M_1,M_NS2)
        #returns Roche limit of M_1 = M_NS1 + M_DM for companion of mass M_NS2 
    	#at a separation of distance a 
        
    	#@assert M_WD < M_NS 
    	#@assert a > 0
    	
    	q=M_1/M_NS2
        q_2_3rds = q^(2.0/3.0)
    	return(a * 0.49*q_2_3rds/(0.6*q_2_3rds+log(1.0+cbrt(q))))
    end
    
    
    """
        Returns the orbital period of the binary.
        ```julia
        period(a, M_tot)
        ```
        a = total orbital separation of the binary.
        M_tot = total mass of the binary (white dwarf + neutron star).
    """
    function period(a, M_tot)
        #@assert a > 0 
    	#@assert M_tot > 0 
    
        #equivalent to sqrt(4.0 * π^2 * a^3 / (G * Mtot))
        return period_factor * sqrt(a^3.0 / M_tot)
    end
    
    
    """
        Returns the decretion (mass loss) rate of the white dwarf.
        ```julia
        decretion_rate(a, R_WD, M_WD, M_NS)
        ```
        a = total orbital separation of the binary.
        R_WD = radius of the white dwarf.
        M_WD = mass of the white dwarf.
        M_NS = mass of the neutron star.
        If Roche lobe is greater than radius of white dwarf, returns 0.
    """
    function decretion_rate(a, R_WD,
                            M_WD, M_NS;optional=optional_default) 
    
    	#@assert R_WD > 0 
    
        RL = Roche_Limit(a, M_WD, M_NS)
        if RL > R_WD
            return 0.0
        end
    
        if optional.poly_shell 
        #use the mass in a shell assuming polytropic equation of state n=1.5
            shell_mass=Polytrope.mass_outside_radius(R_WD,RL,M_WD,1.5)
        else #simplification of a uniform density white dwarf
            shell_mass=M_WD * ((R_WD - RL) / R_WD)^3.0
        end
            
        return -optional.A * shell_mass / period(a, M_WD + M_NS) * ((R_WD - RL) / R_WD)^3.0
    end
    
    """
        Returns the specific angular momentum of material ejected from binary.
        ```julia
        gamma(M_WD, M_NS) 
        ```
        M_WD = mass of the white dwarf.
        M_NS = mass of the neutron star.
        mode = mass loss mode. One of :jeans, :isotropic, or :circumbinary.")
    
        Jeans mode - material leaves the donor star directly, carrying that star’s specific angular momentum.
        Isotropic re-emission - material first accretes onto the neutron star, then is ejected isotropically, carrying the accretor’s angular momentum.
        Circumbinary ring - material forms a ring around both stars at some radius a_ring = a x a_ring_frac, where a is the total binary separation.
    """
    function gamma(M_WD, M_NS;optional=optional_default) 
    
    	#@assert M_WD > 0
    	#@assert M_NS > 0
        #@assert optional.a_ring_frac > 0 
    
        if optional.mode == :jeans
            return M_NS / M_WD
        elseif optional.mode == :isotropic
            return M_WD / M_NS
        elseif optional.mode == :circumbinary
            return (M_WD + M_NS)^2.0 / (M_WD * M_NS) * sqrt(optional.a_ring_frac)
        else
            error("$mode not defined. Use :jeans, :isotropic, or :circumbinary")
        end
    end
    
    
    """
        Returns the Eddington accretion rate onto an object.
        ```julia
        eddington_rate(M)
        ```
        M = mass of the accretor.
        Default values apply to Eddington accretion onto the neutron star from the white dwarf.
    """
    function eddington_rate(M;optional=optional_default) 
    
        #@assert M > 0 
    	#@assert 0 < optional.eta_acc <= 1 
    
        # Equivalent to 4πG*m_p / (σ_T * eta_acc * c) * M
    
        return 7.05e-17 * M / optional.eta_acc  # M in solar masses, 7.05e-17/s
    end
    
    
    """
        Returns the (fraction of) mass accreted onto the neutron star over the mass lost from the white dwarf per unit time.
    """
    function beta(a, R_WD, M_WD, M_NS;optional=optional_default)         
        M_WD_dot = -decretion_rate(a, R_WD, M_WD, M_NS; optional=optional)
        lam_edd = eddington_rate(M_NS; optional=optional)
        if M_WD_dot < lam_edd
            return 1.0
        end
        return lam_edd / M_WD_dot
    end
    
    
    """
        Returns the rate of total angular momentum loss of the binary due to non-conservative mass transfer.
    """
    function J_dot(J, a, R_WD, M_WD, M_NS;optional=optional_default)         
        g = gamma(M_WD, M_NS; optional=optional)
        b = beta(a, R_WD, M_WD, M_NS; optional=optional)
    
    	#@assert b <= 1.0 
    
        M_WD_dot = decretion_rate(a, R_WD, M_WD, M_NS; optional=optional)
        return J * g * (1.0 - b) * M_WD_dot / (M_WD + M_NS)
    end

    """
        Returns the rate of change for total orbital separation.
    """
    function a_dot(a, R_WD, M_WD, M_NS;optional=optional_default) 
        g = gamma(M_WD, M_NS; optional=optional)
        b = beta(a, R_WD, M_WD, M_NS; optional=optional)
    
        # @assert b <= 1.0
    
        M_WD_dot = decretion_rate(a, R_WD, M_WD, M_NS; optional=optional)
        one_minus_b = 1.0 - b
        return -2.0 * a * M_WD_dot / M_WD * (1 - b * M_WD / M_NS - one_minus_b * (g + 0.5) * M_WD / (M_WD + M_NS))
    end
    
    
    """
        Returns rate of change for white dwarf radius.
        ```julia
        dot_R_WD(M_WD, R_WD)
        ```
        M_WD = mass of white dwarf.
        R_WD = radius of white dwarf.
    """
    function dot_R_WD(M_WD, R_WD; optional=optional_default) 
        natural_R = R0_WD(M_WD; optional=optional)
        tau = dynamical_timescale(M_WD, R_WD)
        return (natural_R - R_WD) / tau
    end
    
    
    """
        Returns reduced mass.
        ```julia
        mu(M_WD, M_NS)
        ```
        M_WD = mass of white dwarf.
        M_NS = mass of neutron star.
    """
    function mu(M_WD, M_NS)
    
    	#@assert M_WD > 0 
    	#@assert M_NS > 0 
    
        return M_WD * M_NS / (M_WD + M_NS)
    end
    
    """
        Returns the rate of total angular momentum loss of the binary due to the radiation reaction (see math notes).
        ```julia
        dJ_rr(I2::Matrix{Float64}, I3::Matrix{Float64})::Float64
        ```
        I2 = second order time derivative of quadrupole moment (quadrupole acceleration tensor).
        I3 = third order time derivative of quadrupole moment (quadrupole jerk tensor).
    """
    function dJ_rr(I2::Matrix{Float64}, I3::Matrix{Float64})::Float64
        # Corresponds to: -2G/(5c^5) * (dot(I2[0],I3[1]) - dot(I2[1],I3[0]))
    
        term::Float64 = dot(I2[:, 1], I3[:, 2]) - dot(I2[:, 2], I3[:, 1])
        return dJ_rr_factor * term
    end
    
    """
        Converts polar coordinates to Cartesian coordinates.
    """
    function x(r::T, phase::T)::T where {T<:Real}
    
        #@assert r > 0 
    
        return r * cos(phase)
    end
    
    """
        Converts polar coordinates to Cartesian coordinates.
    """
    function y(r::T, phase::T)::T where {T<:Real}
    
        #@assert r > 0
    
        return r * sin(phase)
    end
    
    """ 
        Returns rate of change of neutron star mass.
    """
    function dM_NS(a, R_WD, M_WD, M_NS ;optional=optional_default) 
        t_dM_WD=decretion_rate(a,R_WD,M_WD,M_NS; optional=optional)
        t_edd=eddington_rate(M_NS; optional=optional)
        return(min(-t_dM_WD,t_edd))
    end
    
    """ 
        Returns rate of change of orbital phase of the white dwarf.
    """
    function dphase(J,a,M_WD, M_NS;optional=optional_default) 
        if optional.period_calc #Calculates change in period under the adiabatic approximation
            return (2.0*π/period(a, M_WD+M_NS))
        end
        return(J/(a^2.0*mu(M_WD,M_NS)))
        
    end  
    
    
    """
        Returns I_{xx}, where I is the quadrupole moment tensor.
        ```julia
        I_xx(M_WD::T, M_NS::T, r::T, phase::T)
        ```
    """
    function I_xx(M_WD::T, M_NS::T, r::T, phase::T)::T where T<:Real
        return (mu(M_WD, M_NS) * (x(r,phase))^2.0)
    end           
    
    """Given polar orbital parameters, calculates I_{xy}, where I is the quadrupole moment tensor."""
    function I_xy(M_WD::T, M_NS::T, r::T, phase::T)::T where T<:Real
        return mu(M_WD, M_NS) * x(r,phase) * y(r,phase)
    end                
    
    """Fixes the higher derivatives of the function so that it is formatted in a way that ForwardDiff is comfortable with."""
    function fix_higher_derivatives(ftn::Function,higher_derivatives::AbstractVector{T})::Function where {T<:Real} 
    	function inner_ftn(gen_theta::AbstractVector{T}) where {T<:Real}
    		return(ftn(gen_theta,higher_derivatives))
    	end
    	return ftn 
    end
    
    """Given a function, its vectorized parameters, and the time derivatives of those parameters, returns the time derivative of the function as calculated via the chain rule."""
    function time_derivative(ftn::Function, ind_vars_and_derivatives::AbstractVector{T}; 
        optional=optional_default,
        ind_vars=zeros(T,length(ind_vars_and_derivatives)-6)) where {T<:Real}
        #ftn is the function to be taken a time derivative of
        #ind_vars is an array with the values of the arguments of ftn
        #derivatives is an array with the known time derivatives of ind_vars
        
        ind_vars.= ind_vars_and_derivatives[1:end-6] #we want a real array (as opposed to a view) to pass through gradient
        derivatives=@view ind_vars_and_derivatives[7:length(ind_vars_and_derivatives)]

        #defining t_ftn to be the function with optional parameters fixed to input optional
        function t_ftn(vars)
            return(ftn(vars; optional=optional))
        end
    
    	# ftn_at_ind_vars is the function evaluated at the parameters in ind_vars
        ftn_at_ind_vars = DiffResults.GradientResult(ind_vars)
    	ForwardDiff.gradient!(ftn_at_ind_vars,t_ftn,ind_vars)

        return(sum(DiffResults.gradient(ftn_at_ind_vars) .* derivatives))#chain rule
    end
    
    """Given a function, defines its next derivative based on the function's parameters and their time derivatives. `n` determines how many time derivatives of \vec{phase} (i.e., the parameters to be integrated) are parameters of the function itself. """
    function next_derivative(in_ftn::Function; n::Integer=1)::Function
        #returns a function which is the next higher derivative of the in_ftn
        #n is the total derivatives past the given analytic form
        @assert n>0
        function ftn(gen_theta_and_derivatives::AbstractVector{T}; optional=optional_default) where {T<:Real}
            
            @assert length(gen_theta_and_derivatives)==6*(n+1) 
    		
            return(time_derivative(in_ftn, gen_theta_and_derivatives; optional=optional))
        end
        return ftn
    end
    
    #The following functions wrap the original derivatives in vectors to make life easier when taking derivatives due to mass transfer only. In particular, `theta` is of the form `[a, M_WD, M_NS, R_WD, J, phase]`.
    """Wrapper function for time derivative of total orbital separation."""
    function da_vector(theta::AbstractVector{T}; optional=optional_default)::T where {T<:Real}    
        # promote the numeric keyword args to match the type of the inputs
        optional = (; optional..., A=T(optional.A), eta_acc=T(optional.eta_acc), a_ring_frac=T(optional.a_ring_frac))
    
        #unpacking theta
        a, M_WD, M_NS, R_WD, J, phase=theta
        
        return(a_dot(a,R_WD, M_WD, M_NS; optional = optional))
    end
    
    """Wrapper function for time derivative of white dwarf mass."""
    function dM_WD_vector(theta::AbstractVector{T}; optional=optional_default)::T where {T<:Real}        
        optional= (; optional..., A=T(optional.A))
        
        a, M_WD, M_NS, R_WD, J, phase=theta
        return(decretion_rate(a, R_WD, M_WD, M_NS; optional=optional))
    end
    
    """Wrapper function for time derivative of neutron star mass."""
    function dM_NS_vector(theta::AbstractVector{T}; optional=optional_default)::T where {T<:Real}    
        optional= (; optional..., A=T(optional.A), eta_acc=T(optional.eta_acc))
        
        a, M_WD, M_NS, R_WD, J, phase=theta
        return(dM_NS(a, R_WD, M_WD, M_NS; optional=optional))
    end
    
    """Wrapper function for time derivative of white dwarf radius."""
    function dR_WD_vector(theta::AbstractVector{T}; optional=optional_default)::T where {T<:Real}        
        optional = (; optional..., mu_e=T(optional.mu_e))
        
        a, M_WD, M_NS, R_WD, J, phase=theta
        return( dot_R_WD(M_WD, R_WD; optional=optional))
    end
    
    """Wrapper function for time derivative of total orbital momentum of the binary."""
    function dJ_vector(theta::AbstractVector{T}; optional=optional_default)::T where {T<:Real}        
        optional = (; optional..., A=T(optional.A), eta_acc=T(optional.eta_acc), a_ring_frac=T(optional.a_ring_frac))
        
        a, M_WD, M_NS, R_WD, J, phase=theta
        return(J_dot(J, a, R_WD, M_WD, M_NS; optional=optional))
    end
    
    """Wrapper function for time derivative of orbital phase of white dwarf."""
    function dphase_vector(theta::AbstractVector{T}; optional=optional_default)::T where {T<:Real} 
        a, M_WD, M_NS, R_WD, J, phase=theta
        return(dphase(J,a,M_WD, M_NS; optional=optional))
    end
    
    """Wrapper function for quadrupole moment tensor's xx component."""
    function Ixx_vector(theta::AbstractVector{T}; optional=optional_default)::T where {T<:Real} 
        a, M_WD, M_NS, R_WD, J, phase=theta
        return(I_xx(M_WD, M_NS, a, phase))
    end
    
    """Wrapper function for quadrupole moment tensor's xy component."""
    function Ixy_vector(theta::AbstractVector{T}; optional=optional_default)::T  where {T<:Real} 
        a, M_WD, M_NS, R_WD, J, phase=theta
        return(I_xy(M_WD, M_NS, a, phase))
    end
    
    #"The below block of code uses the vectorized derivatives above to output functions which take higher time derivatives of $\vec{\theta}$ and $I_{ij}$ as a function of phase. "
    dda_vector=next_derivative(da_vector)
    ddM_WD_vector=next_derivative(dM_WD_vector)
    ddM_NS_vector=next_derivative(dM_NS_vector)
    ddR_WD_vector=next_derivative(dR_WD_vector)
    ddJ_vector=next_derivative(dJ_vector)
    ddphase_vector=next_derivative(dphase_vector)
    
    ddda_vector=next_derivative(dda_vector,n=2)
    dddM_WD_vector=next_derivative(ddM_WD_vector,n=2)
    dddM_NS_vector=next_derivative(ddM_NS_vector,n=2)
    dddR_WD_vector=next_derivative(ddR_WD_vector,n=2)
    dddJ_vector=next_derivative(ddJ_vector,n=2)
    dddphase_vector=next_derivative(ddphase_vector,n=2)
    
    dI_xx_vector=next_derivative(Ixx_vector)
    ddI_xx_vector=next_derivative(dI_xx_vector,n=2)
    dddI_xx_vector=next_derivative(ddI_xx_vector,n=3)
    
    dI_xy_vector=next_derivative(Ixy_vector)
    ddI_xy_vector=next_derivative(dI_xy_vector,n=2)
    dddI_xy_vector=next_derivative(ddI_xy_vector,n=3)
    
    #=bundling derivatives together which will prove convenient later=#
    first_derivatives=[da_vector, dM_WD_vector, dM_NS_vector, dR_WD_vector, dJ_vector, dphase_vector]
    second_derivatives=[dda_vector, ddM_WD_vector, ddM_NS_vector, ddR_WD_vector, ddJ_vector, ddphase_vector]
    third_derivatives=[ddda_vector, dddM_WD_vector, dddM_NS_vector, dddR_WD_vector, dddJ_vector, dddphase_vector]

    """Constructs a symmetric 2x2 gravitational-wave strain tensor corresponding to given “plus” (+) and “cross” (×) polarization amplitudes."""
    function tensor_plus_cross(plus::Float64, cross::Float64)::Matrix{Float64}
        return([plus cross; cross -plus        
        ])
    end
    
    #"The below block of code is central to the interface between the mass transfer and radiation reaction effects. It also sits in the innermost loop of our integrator. This is a good place to peer review."
    #=Central Engine of the integrator=#
    """Calculates the second derivative of the quadrupole moment `t_ddI` due to kinematics of keplerian motion and mass transfer `dtheta`. Also returns higher derivatives of theta useful for calculating the third derivative of the quadrupole moment."""
    function ddI_helper(theta::AbstractVector{T}; parallel::Bool=true, optional=optional_default,
        dtheta=nothing,theta_dtheta=nothing, ddtheta=nothing, theta_ddtheta=nothing #pre-allocating memory
        )::Tuple{Vector{T},Vector{T},Vector{T},Array{T, 2}} where {T}
        #calulates a list of values which will be helpful in the next two functions
        # Allocate temporaries if not provided
        dtheta = dtheta === nothing ? zeros(T,6) : dtheta
        theta_dtheta = theta_dtheta === nothing ? zeros(T,12) : theta_dtheta
        ddtheta = ddtheta === nothing ? zeros(T,6) : ddtheta
        theta_ddtheta = theta_ddtheta === nothing ? zeros(T,18) : theta_ddtheta
        
    
    #calculating first time derivatives due to MT only

        if parallel
            Threads.@threads for i in 1:6
               @inbounds dtheta[i]=first_derivatives[i](theta; optional=optional)
            end
        else
            for i in 1:6
               @inbounds dtheta[i]=first_derivatives[i](theta; optional=optional)
            end
        end
    
        theta_dtheta[1:6] .= theta
        theta_dtheta[7:12] .= dtheta
        #ddtheta=zeros(6)
        #calculating second time derivatives due to MT only
        if parallel
            Threads.@threads for i in 1:6
                @inbounds ddtheta[i]=second_derivatives[i](theta_dtheta; optional=optional)
            end
        else
            for i in 1:6
                @inbounds ddtheta[i]=second_derivatives[i](theta_dtheta; optional=optional)
            end
        end
    
        theta_ddtheta[1:12] .= theta_dtheta
        theta_ddtheta[13:18] .= ddtheta
    
        #using orbital elements and first and second time derivatives 
        #to calculate second derivative of quadrupole moment
        #(for use in calculating RR acceleration)
        t_ddI=tensor_plus_cross(ddI_xx_vector(theta_ddtheta),
            ddI_xy_vector(theta_ddtheta)
        )
        return(dtheta, ddtheta, theta_ddtheta, t_ddI)
    end


    """The first derivative is much costlier than the following derivatives (about as costly as all the others put together), so this function parallelizes using only two workers by having the first worker calculate the first derivative while the second worker calculates the second through sixth derivatives. EDIT: It turns out that this parallelization is not great, so I am not updating it with the optional keyword."""
    function uneven_parallelization_helper(theta::AbstractVector{T}, derivatives)::Vector{T} where {T}
        dtheta::Vector{Float64} = zeros(6)

        # Task 1: expensive separation derivative (separation--many chain rules) 
        task_1 = Threads.@spawn begin
            dtheta[1] = derivatives[1](theta)
        end
    
        # Task 2: lighter subsequent derivatives
        task_2 = Threads.@spawn begin
            @inbounds for j in 2:6
                dtheta[j] = derivatives[j](theta)
            end
        end
    
        # Wait for completion
        fetch(task_1)
        fetch(task_2)
    
        return dtheta
    end
        

    """Same as above function, but attempting to parallelize by dedicating one worker to more intensive calculation of derivatives of separation (a) and one worker to calculating all other derivatives. EDIT: It turns out that this parallelization is not great, so I am not updating it with the optional keyword. """
    function ddI_helper_uneven(theta::AbstractVector{T})::Tuple{Vector{T},Vector{T},Vector{T},Array{T, 2}} where {T}
        #calulates a list of values which will be helpful in the next two functions
        theta_ddtheta::Vector{Float64} = zeros(18)
        theta_ddtheta[1:6] = @view theta[1:end]
        dtheta::Vector{Float64} = uneven_parallelization_helper(theta, first_derivatives)
        theta_ddtheta[7:12] = @view dtheta[1:end]
        #calculating first time derivatives due to MT only
    
        theta_dtheta::Vector{Float64} = @view theta_ddtheta[1:12]
        ddtheta::Vector{Float64} = uneven_parallelization_helper(theta_dtheta, second_derivatives)
        theta_ddtheta[13:18]= @view ddtheta[1:end]
        #calculating second time derivatives due to MT only
    
        #using orbital elements and first and second time derivatives 
        #to calculate second derivative of quadrupole moment
        #(for use in calculating RR acceleration)
        t_ddI::Array{Float64, 2}=tensor_plus_cross(ddI_xx_vector(theta_ddtheta),
            ddI_xy_vector(theta_ddtheta)
        )
        return(dtheta, ddtheta, theta_ddtheta, t_ddI)
    end


    """Returns the second derivative of the quadrupole moment as calculated in `ddI_Helper`. """
    function ddI_from_theta(theta::Vector{Float64}; parallel::Bool=true, optional=optional_default)::Array{Float64,2}
        return(ddI_helper(theta; parallel=parallel, optional=optional)[end])
    end
    	    
    """Calculates the third derivative of the quadrupole moment and returns it along with the time derivative of theta and the second derivative of the quadrupole moment as calculated in `ddI_Helper`."""
    function ddI_dddI_from_theta(theta::Vector{T}; ddI_only::Bool=false, parallel::Bool=true, optional=optional_default,
        dddtheta = nothing, theta_dddtheta = nothing #preallocating
        )::Tuple{Vector{T},Array{T,2},Array{T,2}} where {T}
        #returns time derivatives of theta and second and third derivatives of 
        #the quadrupole moment due to MT effects
        
        dtheta, ddtheta, theta_ddtheta, t_ddI = ddI_helper(theta; parallel=parallel, optional=optional)
        #using orbital elements and first through third time derivatives  
        #to calculate third derivative of quadrupole moment 
        #(for use in calculating RR acceleration)
        
        #dddtheta=zeros(6)
        dddtheta = dddtheta === nothing ? zeros(T,6) : dddtheta
        theta_dddtheta = theta_dddtheta === nothing ? zeros(T,24) : theta_dddtheta
    
        #calculating second time derivatives due to MT only
        if parallel
            Threads.@threads for i in 1:6
                @inbounds dddtheta[i]=third_derivatives[i](theta_ddtheta, optional=optional)
            end
        else
            for i in 1:6
                @inbounds dddtheta[i]=third_derivatives[i](theta_ddtheta, optional=optional)
            end
        end
    
        theta_dddtheta[1:18] .= theta_ddtheta
        theta_dddtheta[19:24] .= dddtheta
    
        t_dddI=tensor_plus_cross(dddI_xx_vector(theta_dddtheta),
            dddI_xy_vector(theta_dddtheta)
        )
        return(dtheta,t_ddI,t_dddI)
    end


"""Uneven parallelization to calculate third derivatives. EDIT: It turns out that this parallelization is not great, so I am not updating it with the optional keyword. """
    function ddI_dddI_from_theta_uneven(theta::Vector{Float64};ddI_only::Bool=false)::Tuple{Vector{Float64},Array{Float64,2},Array{Float64,2}}
        #returns time derivatives of theta and second and third derivatives of 
        #the quadrupole moment due to MT effects
        
        dtheta, ddtheta, theta_ddtheta, t_ddI = ddI_helper_uneven(theta)
        #using orbital elements and first through third time derivatives  
        #to calculate third derivative of quadrupole moment 
        #(for use in calculating RR acceleration)

        theta_dddtheta::Vector{Float64} = zeros(24)
        theta_dddtheta[1:18] = @view theta_ddtheta[1:end]
        dddtheta=uneven_parallelization_helper(theta_ddtheta, third_derivatives)
        theta_dddtheta[19:24] = @view dddtheta[1:end]
        #calculating second time derivatives due to MT only
        
        t_dddI::Array{Float64,2}=tensor_plus_cross(dddI_xx_vector(theta_dddtheta),
            dddI_xy_vector(theta_dddtheta)
        )
        return(dtheta,t_ddI,t_dddI)
    end

    
    """ An "in place" implementation of the full time derivative of theta including effects from both mass transfer and radiation reaction.""" 
    function acceleration!(du::Vector{Float64}, theta::Vector{Float64},p::Vector{Any},t::Float64; parallel::Bool=true, uneven::Bool=false, optional=optional_default,verbose=false)
        if parallel && uneven
            #keyword argument dictates the style of parallelization to be used 
            #(mostly for benchmarking purposes)
            dtheta,t_ddI,t_dddI=ddI_dddI_from_theta_uneven(theta)
        else
            dtheta,t_ddI,t_dddI=ddI_dddI_from_theta(theta; parallel=parallel, optional=optional)
        end
    
        #forces due to RR
        t_dJ_rr::Float64=dJ_rr(t_ddI,t_dddI)
        @inbounds begin
            a::Float64 = theta[1]
            J::Float64 = theta[5]
            t_da_rr::Float64 = 2.0 * a * t_dJ_rr / J

            #println(t_ddI,t_dddI)
            #println(t_dJ_rr)
            #println(t_da_rr)# For debugging--DELETE LATER
            #println()
        
            #adding RR forces to MT forces
            dtheta[1]+=t_da_rr
            dtheta[5]+=t_dJ_rr
        end

        if verbose
            println(dtheta)
        end
    
        #mutate du
        du .= dtheta
    end

    """ Normal parallelization implementation of in-place acceleration! function"""    
    function acceleration_even!(du::Vector{Float64}, theta::Vector{Float64},p::Vector{Any},t::Float64; optional=optional_default)
        acceleration!(du, theta, p, t; optional=optional)
    end

    """ Uneven parallelization implementation of in-place acceleration! function"""    
    function acceleration_uneven!(du::Vector{Float64}, theta::Vector{Float64},p::Vector{Any},t::Float64)
        acceleration!(du, theta, p, t, uneven=true)
    end
        
    """ Serial implementation of in-place acceleration! function"""
    function acceleration_serial!(du::Vector{Float64}, theta::Vector{Float64},p::Vector{Any},t::Float64; optional=optional_default)
        acceleration!(du, theta, p, t; parallel=false, optional=optional)
    end
    
    #The next several blocks define functions intended to initiate the orbit beginning at mass transfer for initial masses
    #Next, in Integrate.jl, we will define functions to take an integration step forward calculating waveform and finally to integrate until the White Dwarf dissipates completely.
    """Calculates in km the orbital separation for the white dwarf to overflow its Roche Lobe and begin mass transfer"""
    function RL_contact(M_NS1::Float64,M_NS2::Float64, M_DM::Float64, RealPoly::Any; optional=optional_default)::Float64
        #returns separation at which the donor DM halo overflows its Roche Lobe
        R_DM::Float64 = RealPoly.R_DM
        RL_a1::Float64 = Roche_Limit(1.,(M_NS1 + M_DM),M_NS2) #roche limit at separation of 1 km
        return(R_DM / RL_a1)#initial separation in km
    end
    
    """Calculates the angular momentum of the system at initial Roche Lobe overflow """	
    function circular_J(a::Float64,M_1::Float64,M_NS2::Float64)::Float64
        #returns angular momentum for a circular orbit
        return(( a * G * M_1 ^ 2.0 * M_NS2 ^ 2.0 / (M_1 + M_NS2))^(1.0/2.0))
    end

    function cgs_density(density)
        return (density * 1.988416e30 / 1e12)
    end

    function cgs_density_inv(density)
        return (density / 1.988416e30 * 1e12)
    end
end