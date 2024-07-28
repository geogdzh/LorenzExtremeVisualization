include("analyze_util.jl")

# create directory for figures
isdir(pwd() * "/figs") ? nothing : mkdir(pwd() * "/figs")

## 
include("comparing_climates.jl")


##
include("holding_times.jl")


##
include("bayesian_windows.jl")