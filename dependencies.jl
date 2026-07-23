using Pkg          
Pkg.instantiate()
Pkg.activate(".")

using BenchmarkTools
#using Integrate
using CSV, DataFrames, Interpolations
using Plots
using Plots.PlotMeasures
using OrdinaryDiffEq
using DifferentialEquations
using Printf, DelimitedFiles
using ForwardDiff
using LaTeXStrings
using NamedArrays
using ColorSchemes
using FFTW


#include("source/GW.jl")
#include("source/install.jl")
#include("source/Integrate_r.jl")
include("source/Integrate.jl")
#include("source/make_poly_tables.jl")
include("source/math.jl")
include("source/Pipeline.jl")
include("source/polytrope.jl")
#include("source/save.jl")

cb_palette = [
    "#0072B2", # blue
    "#D55E00", # vermillion
    "#009E73", # bluish green
    "#CC79A7", # reddish purple
    "#E69F00", # orange
    "#56B4E9", # sky blue
    "#F0E442"  # yellow
];