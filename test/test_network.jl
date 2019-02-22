using ParseRxNetwork
using DiffEqBase, DiffEqBiological, DiffEqJump, OrdinaryDiffEq, Plots

using TimerOutputs

# Create a TimerOutput, this is the main type that keeps track of everything.
const to = TimerOutput()

# parameters
method = RSSA()
nsteps = 1e5                  # time to simulate to
networkname = "tester"
tf = 10.
#speciesf = path to species initial population file
#rxsf = path to reaction network file


# get the reaction network
reset_timer!(to)
@timeit to "netgen" prn = loadrxnetwork(RSSANetwork(), networkname, speciesf, rxsf; printrxs = false)
rn = prn.rn 
initialpop = prn.u₀
println("network parsed")

# one simulation
@timeit to "addjumps" addjumps!(rn,build_regular_jumps=false, minimal_jumps=true)
println("added jumps")
@timeit to "dprob" prob = DiscreteProblem(rn, initialpop, (0.,tf))
@timeit to "jprob" jump_prob = JumpProblem(prob, method, rn; save_positions=(false,false))
# integrator = init(jump_prob, SSAStepper(), saveat=tf/1000)
# print("solving...")
# @timeit to "stepping" begin
#     for i = 1:nsteps
#         step!(integrator)
#     end
# end
# sol = integrator.sol
# println("done")
@timeit to "solve" sol = solve(jump_prob, SSAStepper(), saveat=tf/1000)
show(to)


# bng file
reset_timer!(to)
@timeit to "bionetgen" prnbng = loadrxnetwork(BNGNetwork(),string(networkname,"bng"), bngfname) 
rnbng = prnbng.rn; u₀ = prnbng.u₀; p = prnbng.p; shortsymstosyms = prnbng.symstonames;
@timeit to "addjumps" addjumps!(rnbng,build_regular_jumps=false, minimal_jumps=true)
@timeit to "dprob" prob = DiscreteProblem(rnbng, initialpop, (0.,tf), p)
@timeit to "jprob" jump_prob = JumpProblem(prob, method, rnbng; save_positions=(false,false))
@timeit to "solve" sol = solve(jump_prob, SSAStepper(), saveat=tf/1000)
show(to)

# plot
#labs = [String(sym) for i=1:1, sym in rn.syms];
#plot(sol.t, sol[:,:]', linetype=:steppost, labels=labs)
#plot(sol.t, sol[:,:]', labels=labs)
