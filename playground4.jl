#merger1->1.4,1.4, K=3.25
#merger2-> 1.4,3.0, K =3.25
#merger3->1.4,1.4, K = 2.7
#merger4-> 1.4,3, K =2.7
file1 = CSV.File("Merger4_noDM.csv")    #no DM mass transfer
t1 = file1.t 
freq1 = file1.f
ddI_p1 = file1.ddI_p
phase1 = file1.Phase

file2 = CSV.File("Merger4_10per_iso.csv")     #10% DM accreted, isotropic    
t2 = file2.t 
freq2 = file2.f
ddI_p2 = file2.ddI_p
phase2 = file2.Phase

file3 = CSV.File("Merger4_10per_jeans.csv")    #10% DM accreted, jeans
t3 = file3.t 
freq3 = file3.f
ddI_p3 = file3.ddI_p
phase3 = file3.Phase

file4 = CSV.File("Merger4_100per.csv")    #100% DM accreted
t4 = file4.t 
freq4 = file4.f
ddI_p4 = file4.ddI_p
phase4 = file4.Phase

phase1interp = linear_interpolation(freq1,phase1,extrapolation_bc=Flat())
dphase2 = phase2 .- phase1interp(freq2)
dphase3 = phase3 .- phase1interp(freq3)
dphase4 = phase4 .- phase1interp(freq4)

#phase vs freq
p1 = plot(freq1, phase1, label = "no DM", xlabel = "GW Frequency [Hz]", ylabel = "Phase [rad]", dpi = 300, xlim = (0,2000))
plot!(p1, freq2, phase2, label = "10%, isotropic")
plot!(p1, freq3, phase3, label = "10%, Jeans")
plot!(p1, freq4, phase4, label = "100% accretion")

#dphase vs freq
p2 = plot(freq2, dphase2, label = "10%, isotropic", xlabel = "GW Frequency [Hz]", ylabel = "Δ phase [rad]", dpi = 300, xlim = (0,2000),color = cb_palette[1], title = "Extra Accumulated Phase",
titlefontsize=16,guidefontsize=16,tickfontsize=14,legendfontsize=12)
plot!(p2, freq3, dphase3, label = "10%, Jeans", color = cb_palette[2])
plot!(p2, freq4, dphase4, label = "100% accretion", color = cb_palette[3])
annotate!(p2, 142, 5.9, "K=3.25e10",12)
annotate!(p2, 198, 4.9, "ρ₀=9e13 g/cm^3",12)

#ddI vs time
p3 = plot(t1, ddI_p1, label = "no DM", xlabel = "Time [s]", ylabel = "ddI_p", dpi = 300)
plot!(p3, t2, ddI_p2, label = "10%, isotropic")
plot!(p3, t3, ddI_p3, label = "10%, Jeans")
plot!(p3, t4, ddI_p4, label = "100% accretion")

#freq vs time
p4 = plot(t1, freq1, label = "no DM", xlabel = "Time [s]", ylabel = "GW Frequency [Hz]", dpi = 300, ylim = (0,2000))
plot!(p4, t2, freq2, label = "10%, isotropic")
plot!(p4, t3, freq3, label = "10%, Jeans")
plot!(p4, t4, freq4, label = "100% accretion")

#display(p1) #phase vs freq
display(p2) #dphase vs freq
#display(p3) #ddI vs time
display(p4) #freq vs time