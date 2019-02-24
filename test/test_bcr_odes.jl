# download input files from: https://gist.github.com/isaacsas/7648fe065f1884307a2906a0a48ea12b
# then add ParseRxNetwork package (unregistered)
#] 
# add https://github.com/isaacsas/ParseRxNetwork.jl.git

using DiffEqBase, DiffEqBiological,Plots,OrdinaryDiffEq,Sundials
using ParseRxNetwork
using TimerOutputs

# parameters
doplot = true
networkname = "testbcrbng"
tf = 10000.

# input files, specify path 
#speciesf = "SOME-PATH-TO/BCR_pop.txt"
#rxsf = "SOME-PATH-TO/BCR_rxn.txt"

const trssa = TimerOutput()
const to = TimerOutput()
reset_timer!(trssa)
reset_timer!(to)

#initialpop = convert.(Float64, prn.u₀)
# @timeit trssa "netgen" prn = loadrxnetwork(RSSANetwork(), networkname, speciesf, rxsf; printrxs = false)
# rn = prn.rn 
@timeit to "bionetgen" prnbng = loadrxnetwork(BNGNetwork(),string(networkname,"bng"), bngfname); 
rnbng = prnbng.rn; u₀ = prnbng.u₀; p = prnbng.p; 

u0 = u₀   #convert.(Float64, prn.u₀)

# @timeit trssa "raddodes" addodes!(rn; build_jac=false, build_symfuncs=false)
# @timeit trssa "rODEProb" roprob = ODEProblem(rn, u0, (0.,tf))
# show(trssa)

# get the bng reaction network
@timeit to "baddodes" addodes!(rnbng; build_jac=false, build_symfuncs=false)
@timeit to "bODEProb" boprob = ODEProblem(rnbng, u0, (0.,tf), p)
show(to)


# note solvers run _much_ faster the second time 
#reset_timer!(trssa); @timeit trssa "RSSA_CVODE_BDF" begin rsol = solve(roprob, CVODE_BDF(linear_solver=:GMRES),dense=false, saveat=tf/1000); end; show(trssa)
reset_timer!(to); @timeit to "BNG_CVODE_BDF" begin bsol = solve(boprob, CVODE_BDF(),dense=false, saveat=1., abstol=1e-8, reltol=1e-8); end; show(to)

#rasyk = sum(rsol[prnbng.groupstoids[:Activated_Syk],:], dims=1)
basyk = sum(bsol[prnbng.groupstoids[:Activated_Syk],:], dims=1)
vars = sum(bsol[prnbng.groupstoids[:Ig_alpha_P],:],dims=1)
#plot(rsol.t, rasyk', label=:rsol)

if doplot
    pyplot()
    plot(bsol.t[2:end], basyk'[2:end], label=:Activated_Syk, xscale=:log10)
    plot!(bsol.t[2:end], vars'[2:end], label=:Ig_alpha_P, xscale=:log10)
end