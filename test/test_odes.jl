#] 
# add https://github.com/isaacsas/ParseRxNetwork.jl.git

using DiffEqBase, DiffEqBiological,Plots,OrdinaryDiffEq, Sundials
using ParseRxNetwork
using TimerOutputs


# parameters
networkname = "tester"
tf = 10.

# input files
#speciesf = "SOME-PATH-TO/BCR_pop.txt"
#rxsf = "SOME-PATH-TO/BCR_rxn.txt"

# Create a TimerOutput, this is the main type that keeps track of everything.
const to = TimerOutput()
reset_timer!(to)

# get the reaction network
@timeit to "netgen" rn,initialpop = get_rxnetwork_simple(networkname, speciesf, rxsf; printrxs = false)
@timeit to "addodes" addodes!(rn; build_symjac=false, build_symfuncs=false)
@timeit to "ODEProb" oprob = ODEProblem(rn,convert.(Float64,initialpop),(0.,tf))
show(to)

# note solvers run _much_ faster the second time 
reset_timer!(to); @timeit to "CVODE_BDF" begin sol = solve(oprob,CVODE_BDF(),dense=false); end; show(to)
reset_timer!(to); @timeit to "CVODE_BDF" begin sol = solve(oprob,CVODE_BDF(),dense=false); end; show(to)

reset_timer!(to); @timeit to "Rodas4" begin sol = solve(oprob,Rodas4(autodiff=false),dense=false); end;
reset_timer!(to); @timeit to "Rodas4" begin sol = solve(oprob,Rodas4(autodiff=false),dense=false); end;
reset_timer!(to); @timeit to "Rodas4" begin sol = solve(oprob,Rodas4(autodiff=false),dense=false,calck=false); end; show(to)