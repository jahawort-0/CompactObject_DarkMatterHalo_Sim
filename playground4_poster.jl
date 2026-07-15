file1 = CSV.File("Merger_noDM.csv")    #no DM mass transfer
    t1 = file1.t 
    freq1 = file1.f
    ddI_p1 = file1.ddI_p
    phase1 = file1.Phase

file2 = CSV.File("Merger_10per_iso.csv")     #10% DM accreted, isotropic    
    t2 = file2.t 
    freq2 = file2.f
    ddI_p2 = file2.ddI_p
    phase2 = file2.Phase

file3 = CSV.File("Merger_10per_jeans.csv")    #10% DM accreted, jeans
    t3 = file3.t 
    freq3 = file3.f
    ddI_p3 = file3.ddI_p
    phase3 = file3.Phase

file4 = CSV.File("Merger_100per.csv")    #100% DM accreted
    t4 = file4.t 
    freq4 = file4.f
    ddI_p4 = file4.ddI_p
    phase4 = file4.Phase

phase1interp = linear_interpolation(freq1,phase1,extrapolation_bc=Flat())
dphase2 = phase2 .- phase1interp(freq2)
dphase3 = phase3 .- phase1interp(freq3)
dphase4 = phase4 .- phase1interp(freq4)

#phase vs freq
p1 = plot(freq1, phase1, label = "no DM", xlabel = "GW Frequency [Hz]", ylabel = "Phase [rad]", dpi = 300, xlim = (0,1100))
plot!(p1, freq2, phase2, label = "10%, isotropic")
plot!(p1, freq3, phase3, label = "10%, Jeans")
plot!(p1, freq4, phase4, label = "100% accretion")

#dphase vs freq
p2 = plot(freq2, dphase2, label = L"10\%, isotropic", xlabel = L"\mathrm{GW} f\ \mathrm{[Hz]}", ylabel = L"\Delta \phi\ \mathrm{[rad]}", dpi = 500, xlim = (150,1100),color = cb_palette[1],
titlefontsize=16,guidefontsize=16,tickfontsize=14,legendfontsize=10,linestyle=:solid, linewidth=3)
plot!(p2, freq3, dphase3, label = L"10\%, Jeans", color = cb_palette[2],linestyle=:dash, linewidth=3)
plot!(p2, freq4, dphase4, label = L"100\%\ accretion", color = cb_palette[3],linestyle=:dot, linewidth=3)
annotate!(p2, 280, 6.7, L"K=3.25\times 10^{10}",10)
annotate!(p2, 305, 5.9, L"\rho _0=9\times 10^{13}\ g/cm^3",10)

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

#display(p1) #phase vs freq
display(p2) #dphase vs freq
#display(p3) #ddI vs time
#display(p4) #freq vs time