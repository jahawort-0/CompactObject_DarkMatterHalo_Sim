include("dependencies.jl")
## Set parameters
#have been using K = 3.25e10, rho_0 = Math.cgs_density_inv(9e13)
#Merger5 is for poster plotting, K = 2.7e10, rho_0 = Math.cgs_density_inv(1.46666e14), 99% halo
n = 5.
#K = 10^10.45
K = 3.25e10
#rho_0 = Math.cgs_density_inv(10^14)rho_0 = Math.cgs_density_inv(9e13)
#M_DM = 0.1
M_NS1 = 1.4 #Solar masses
R_NS1 = 12.  #km
M_NS2 = 3. #Solar masses
a_start = 200. #km
filename = "Merger_10per_jeans2.csv"

Pipeline.evolve_halo(n,K,rho_0,R_NS1,M_NS1,M_NS2,a_start,filename)
##
# df = CSV.File("Merger3_10per_iso.csv")
# p1 = plot(df.t, df.a, xlabel = "time [s]", ylabel = "semi major axis [km]")
# p2 = plot(df.t, df.f, xlabel = "time [s]", ylabel  = "frequency [Hz]")

# display(p1)
# display(p2)