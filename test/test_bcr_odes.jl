# download input files from: https://gist.github.com/isaacsas/7648fe065f1884307a2906a0a48ea12b
# then add ParseRxNetwork package (unregistered)
#] 
# add https://github.com/isaacsas/ParseRxNetwork.jl.git

using DiffEqBase, DiffEqBiological,Plots,OrdinaryDiffEq, Sundials
using ParseRxNetwork
using TimerOutputs

# parameters
networkname = "testbcrbng"
tf = 10.

# input files, specify path 
#speciesf = "SOME-PATH-TO/BCR_pop.txt"
#rxsf = "SOME-PATH-TO/BCR_rxn.txt"

const to = TimerOutput()
reset_timer!(to)

# get the reaction network
@timeit to "bionetgen" prnbng = loadrxnetwork(BNGNetwork(),string(networkname,"bng"), bngfname); 
rnbng = prnbng.rn; u₀ = prnbng.u₀; p = prnbng.p; shortsymstosyms = prnbng.symstonames;
@timeit to "addodes" addodes!(rn; build_jac=false, build_symfuncs=false)
@timeit to "ODEProb" oprob = ODEProblem(rn, u₀, (0.,tf), p)
show(to)

# note solvers run _much_ faster the second time 
reset_timer!(to); @timeit to "CVODE_BDF" begin sol = solve(oprob, CVODE_BDF(),dense=false); end; show(to)

asyk = sum(sol[prnbng.groupstoids[:Activated_Syk],:], dims=1)
plot(sol.t, asyk, label=:Activated_Syk)


