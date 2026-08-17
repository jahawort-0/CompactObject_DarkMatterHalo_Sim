include("dependencies.jl")
## Import
file = CSV.File("Merger8_10per_jeans.csv");    #10% DM accreted, jeans
    time = file.t;
    freq = file.f;
    ddI_p = file.ddI_p;
    ddI_c = file.ddI_c;
    phase = file.Phase;
    a = file.a;
    M_DM = file.M_DM;

file2 = CSV.File("Merger8_noDM.csv");    #10% DM accreted, jeans
    time2 = file2.t;
    freq2 = file2.f;
    ddI_p2 = file2.ddI_p;
    ddI_c2 = file2.ddI_c;
    phase2 = file2.Phase;
    a2 = file2.a;
    M_DM2 = file2.M_DM;

## Plot strain vs time
DL = 1 #Mpc
hp_t = Pipeline.calc_strain(ddI_p,(DL*Math.Mpc_to_km))
hp_t2 = Pipeline.calc_strain(ddI_p2,(DL*Math.Mpc_to_km))

plot(time,hp_t*2e21, xlim = (time[end]-1,time[end]))
#plot!(time,-cos.(2 .*phase), xlim=(time[end]-1,time[end]))

## FFT
A_h,h_f_tilde = Pipeline.FT_strain(hp_t,phase,freq,time)
A_h2, h_f_tilde2 = Pipeline.FT_strain(hp_t2,phase2,freq2,time2)
#plot(freq,A_h,xlabel = "frequency [Hz]", ylabel = "Strain Amplitude",dpi=200,xlim=(0,1000), ylim = (0,3.5e-21))

## h tilde
p1 = plot(log10.(freq),log10.(h_f_tilde),xlabel = "log frequency [Hz]", ylabel = "log ~h(f)",dpi=200,xlim=(1,log10(1100)),ylim=(-23.5,log10(4e-21)),
label = "w/ DM")
plot!(p1,log10.(freq2),log10.(h_f_tilde2),linestyle=:dash, label = "no DM")
display(p1)

##phase of h tilde
hf_phase = Pipeline.htilde_phase(phase,freq,time)
h_f_tilde_full = h_f_tilde .* exp.(im .* hf_phase)

h_f_tilde_full_fiducial = h_f_tilde2 .* exp.(im .* Pipeline.htilde_phase(phase2,freq2,time2))
##  Plot complex FT of GW
plot(freq,real(h_f_tilde_full),xlabel = "frequency", ylabel = "amplitude", label = "real component", dpi=300, xlim = (0,1200))
plot!(freq,imag(h_f_tilde_full),linestyle =:dot, label = "imaginary component")

## Take inner products
IP_fid = Pipeline.def_integral_df(h_f_tilde_full_fiducial,h_f_tilde_full_fiducial,freq2)
IP_DM = Pipeline.def_integral_df(h_f_tilde_full,h_f_tilde_full,freq)
#IP_cross