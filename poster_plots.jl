include("dependencies.jl")

##
file1 = CSV.File("Merger_noDM.csv")    #no DM mass transfer
    t1 = file1.t
    freq1 = file1.f
    ddI_p1 = file1.ddI_p
    phase1 = file1.Phase
    a1 = file1.a
    M_DM1 = file1.M_DM

file2 = CSV.File("Merger_10per_iso.csv")     #10% DM accreted, isotropic    
    t2 = file2.t 
    freq2 = file2.f
    ddI_p2 = file2.ddI_p
    phase2 = file2.Phase
    a2 = file2.a
    M_DM2 = file2.M_DM

file3 = CSV.File("Merger_10per_jeans2.csv")    #10% DM accreted, jeans
    t3 = file3.t 
    freq3 = file3.f
    ddI_p3 = file3.ddI_p
    ddI_c3 = file3.ddI_c
    phase3 = file3.Phase
    a3 = file3.a
    M_DM3 = file3.M_DM

file4 = CSV.File("Merger_100per.csv")    #100% DM accreted
    t4 = file4.t 
    freq4 = file4.f
    ddI_p4 = file4.ddI_p
    phase4 = file4.Phase

phase1interp = linear_interpolation(sort(freq1),sort(phase1),extrapolation_bc=Flat())
dphase2 = phase2 .- phase1interp(freq2)
dphase3 = phase3 .- phase1interp(freq3)
dphase4 = phase4 .- phase1interp(freq4)

#phase vs freq
p1 = plot(freq1, phase1, label = "no DM", xlabel = "GW Frequency [Hz]", ylabel = "Phase [rad]", dpi = 300, xlim = (0,1100))
plot!(p1, freq2, phase2, label = "10%, isotropic")
plot!(p1, freq3, phase3, label = "10%, Jeans")
plot!(p1, freq4, phase4, label = "100% accretion")

#dphase vs freq
p2 = plot(freq2, dphase2, label = L"\mathrm{Isotropic}", xlabel = L"\mathrm{GW} f\ \mathrm{[Hz]}", ylabel = L"\Delta \phi\ \mathrm{[rad]}", dpi = 500, xlim = (150,1100),color = cb_palette[2],
titlefontsize=16,guidefontsize=16,tickfontsize=14,legendfontsize=10,linestyle=:dash, linewidth=3, size=(800,400),
bottom_margin=15px,left_margin=15px)
plot!(p2, freq3, dphase3, label = L"\mathrm{Jeans}", color = cb_palette[1],linestyle=:solid, linewidth=3)
plot!(p2, freq4, dphase4, label = L"100\%\ \mathrm{accretion}", color = cb_palette[3],linestyle=:dot, linewidth=3)
annotate!(p2, 280, 6.7, L"K=3.25\times 10^{10}",10)
annotate!(p2, 305, 5.9, L"\rho _0=9\times 10^{13}\ \mathrm{g/cm^3}",10)

#ddI vs time
p3 = plot(t1, ddI_p1, label = "no DM", xlabel = "Time [s]", ylabel = "ddI_p", dpi = 300)
plot!(p3, t2, ddI_p2, label = "10%, isotropic")
plot!(p3, t3, ddI_p3, label = "10%, Jeans")
plot!(p3, t4, ddI_p4, label = "100% accretion")

#freq vs time
p4 = plot(t1, freq1, label = "no DM", xlabel = "Time [s]", ylabel = "GW Frequency [Hz]", dpi = 300, ylim = (0,1100))
plot!(p4, t2, freq2, label = "10%, isotropic")
plot!(p4, t3, freq3, label = "10%, Jeans")
plot!(p4, t4, freq4, label = "100% accretion")

display(p1) #phase vs freq
display(p2) #dphase vs freq
display(p3) #ddI vs time
display(p4) #freq vs time
##

p5 = plot(t3,a3, xlabel = L"t\ \mathrm{[s]}", ylabel = L"a\ \mathrm{[km]}",dpi = 500,label="with DM",color = cb_palette[1],
guidefontsize=20,tickfontsize=14,legendfontsize=10,linestyle=:solid, linewidth=3,ylim=(0,110),legend=:bottomleft,
size=(800,300),bottom_margin=25px, top_margin=0px, left_margin=20px)
plot!(p5, t1, a1, linewidth = 3, linestyle=:dot, label = "without DM", color=cb_palette[5])

p6 = plot(t3,M_DM3, xlabel = L"t\ \mathrm{[s]}", ylabel = L"M_{DM}\ \mathrm{[M_\odot]}",dpi = 500,legend = false,color = cb_palette[1],
guidefontsize=20,tickfontsize=14,legendfontsize=10,linestyle=:solid, linewidth=3, ylim = (0,0.2),
size=(800,300),bottom_margin=25px, top_margin=0px, left_margin=20px)

p7 = plot(t3,freq3, xlabel = L"t\ \mathrm{[s]}",ylabel = L"f\ \mathrm{[Hz]}",dpi = 500,legend = false,color = cb_palette[1],
guidefontsize=20,tickfontsize=14,legendfontsize=10,linestyle=:solid, linewidth=3, ylim = (0,1100),
size=(800,300),bottom_margin=25px, top_margin=0px, left_margin=20px)

p8 = plot(p5,p6,p7,layout = (3,1),link = :x,margins=2Plots.mm, size=(800,600))

display(p5)
display(p6)
display(p7)
#display(p8)

## polytrope
n = 5
K = 3.25e10
rho_0 = Math.cgs_density_inv(9e13)
M_NS1 = 1.4 #Solar masses
R_NS1 = 12  #km
M_NS2 = 1.4 #Solar masses
a_start = 200. #kmßß

#Initialize Polytrope
realpoly = Polytrope.solve_halo(n,K,M_NS1,R_NS1,rho_0,2000) #realpoly = [rs,R_DM,rho_r,M_DM_r,rho_interp,mass_interp,M_DM]

#check plot of polytrope
p10 = vline([12],color=cb_palette[2],linewidth=3,linestyle=:dash)
plot!(p10,realpoly.rs, Math.cgs_density.(realpoly.rho_r)./1e13, xlabel = L"r\ \mathrm{[km]}", ylabel = L"\rho\ \mathrm{[g/cm^3]}", xlim = (0,40), ylim=(0,10),
dpi=500,guidefontsize=30,tickfontsize=16,legendfontsize=10,linestyle=:solid, linewidth=5, top_margin = 40px, legend=false,size = (800, 600), color = cb_palette[1],
bottom_margin = 5px, left_margin = 10px, right_margin=5px)
annotate!(p10,1,10.5,L"\times 10^{13}",16)
annotate!(p10,30,8,L"K=3.25\times 10^{10}",20)
annotate!(p10,30,7.0,L"\rho _0=9\times 10^{13}\ \mathrm{g/cm^3}",20)
annotate!(p10,15.5,9.7,"NS surface")
display(p10)

p11 = vline([12],color=cb_palette[2],linewidth=3,linestyle=:dash)
plot!(p11,realpoly.rs, realpoly.mass_r, xlabel = L"r\ \mathrm{[km]}", ylabel = L"M_{\mathrm{enclosed}}\ \mathrm{[M_\odot]}", xlim = (0,realpoly.R_DM), ylim=(0,0.205),
dpi=500,guidefontsize=30,tickfontsize=16,legendfontsize=10,linestyle=:solid, linewidth=5, legend=false, size=(800,600), bottom_margin=5px, color = cb_palette[1],
top_margin=0px, left_margin=10px)
annotate!(p11,35,0.195,"NS surface")

display(p11)

## GW strain
DL = 1 #Mpc
hp_t = Pipeline.calc_strain(ddI_p3,(DL*Math.Mpc_to_km))

plot(t3,hp_t)

freqs,FT_h = Pipeline.FT_strain(hp_t,t3)
plot(freqs,FT_h,xlabel = "frequency [Hz]", ylabel = "Amplitude",dpi=200)