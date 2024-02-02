using Pkg
Pkg.activate(@__DIR__) # @__DIR__ = directory this script is in
Pkg.instantiate() # If a Manifest.toml file exist in the current project, download all the packages declared in that manifest. Else, resolve a set of feasible packages from the Project.toml files and install them.

using CSV
using XLSX
using DataFrames
using YAML
using Printf
using JuMP
using Gurobi
using Cbc
using DelimitedFiles

cd("C://Users//KENISM//OneDrive - VITO//Documents//_Research_2022a_OffShoreBiddingZone//_Models") #may have to change

a = Model(Gurobi.Optimizer)

include("data.jl")

include("function_nodalclearing.jl")

# include("function_zonalclearing_exact_HC.jl")

# include("function_zonalclearing_GSK_SHC.jl")
# include("function_zonalclearing_GSK_AHC.jl")

include("redispatch.jl")

include("process_results.jl")

