# download input files from: https://gist.github.com/isaacsas/7648fe065f1884307a2906a0a48ea12b
# then add ParseRxNetwork package (unregistered)
#] 
# add https://github.com/isaacsas/ParseRxNetwork.jl.git

using DiffEqBase, DiffEqBiological, Plots, OrdinaryDiffEq, Sundials, DataFrames, CSV
using ParseRxNetwork
using TimerOutputs

# parameters
doplot = true
networkname = "testbcrbng"
tf = 10000.

# BNG simulation data
fname    = "../data/bcr/bcr.net"
cdatfile = "../data/bcr/bcr.cdat"
gdatfile = "../data/bcr/bcr.gdat"
cdatdf = CSV.File(cdatfile, delim=" ", ignorerepeated=true) |> DataFrame
gdatdf = CSV.File(gdatfile, delim=" ", ignorerepeated=true) |> DataFrame

const to = TimerOutput()
reset_timer!(to)

@timeit to "bionetgen" prnbng = loadrxnetwork(BNGNetwork(),string(networkname,"bng"), fname); 
rnbng = prnbng.rn; u₀ = prnbng.u₀; p = prnbng.p; 

u0 = u₀   #convert.(Float64, prn.u₀)

# get the bng reaction network
@timeit to "baddodes" addodes!(rnbng; build_jac=false, build_symfuncs=false)
@timeit to "bODEProb" boprob = ODEProblem(rnbng, u0, (0.,tf), p)
show(to)

# BNG simulation data testing
asykgroups = prnbng.groupstoids[:Activated_Syk]

# note solvers run _much_ faster the second time 
# reset_timer!(to); @timeit to "BNG_CVODE_BDF" begin bsol = solve(boprob, CVODE_BDF(),dense=false, saveat=1., abstol=1e-8, reltol=1e-8); end; show(to)

# basyk = sum(bsol[prnbng.groupstoids[:Activated_Syk],:], dims=1)
# vars = sum(bsol[prnbng.groupstoids[:Ig_alpha_P],:],dims=1)

# if doplot
#     pyplot()
#     plot(bsol.t[2:end], basyk'[2:end], label=:Activated_Syk, xscale=:log10)
#     plot!(bsol.t[2:end], vars'[2:end], label=:Ig_alpha_P, xscale=:log10)
# end