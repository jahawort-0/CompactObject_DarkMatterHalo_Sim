include("dependencies.jl")

##
n = 5
K = 3e9
M_DM = 0.05
M_NS1 = 1.35
M_NS2 = 2.5

poly = Polytrope.compute_polytrope(n)   #integrate nondimensional polytrope
realpoly = Polytrope.apply_polytrope(poly,M_DM,M_NS1,K,n)   # redimensionalize polytrope
initial_conditions=Integrate.initial_circular_orbit(M_NS1, M_NS2, M_DM, realpoly);  #Solve for Roche sepeartion, ciruclar orbit

initial_conditions_w_time=zeros(8);
initial_conditions_w_time[1:7]=initial_conditions;

## Integration
phase_1_output=Integrate_r.integrate_r_sneaky(initial_conditions_w_time,-1); #The -1 tells the function that we are decreasing in separation
phase_1_output_total=phase_1_output.+initial_conditions_w_time'; #Need to add changes to isotropic