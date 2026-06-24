using Pkg          
Pkg.instantiate()
Pkg.activate(".")

using BenchmarkTools
#using Integrate
using CSV, DataFrames, Interpolations
using Plots
using OrdinaryDiffEq
using DifferentialEquations
using Printf, DelimitedFiles
using ForwardDiff
using LaTeXStrings
using NamedArrays

#include("source/GW.jl")
#include("source/install.jl")
#include("source/Integrate_r.jl")
include("source/Integrate.jl")
#include("source/make_poly_tables.jl")
include("source/math.jl")
#include("source/Pipeline.jl")
include("source/polytrope.jl")
#include("source/save.jl")