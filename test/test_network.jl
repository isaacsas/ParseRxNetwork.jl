using ParseRxNetwork, DiffEqBase, DiffEqBiological, DiffEqJump, OrdinaryDiffEq, Plots

# parameters
method = Direct()       
tf = .1                  # time to simulate to
networkname = "tester"
#speciesf = path to species initial population file
#rxsf = path to reaction network file

# get the reaction network
rn,initialpop = get_rxnetwork_simple(networkname, speciesf, rxsf; printrxs = true)

# one simulation
prob = DiscreteProblem(initialpop, (0.,tf))
jump_prob = JumpProblem(prob, method, rn; save_positions=(false,false))
sol = solve(jump_prob, SSAStepper(), saveat=tf/1000)

# plot
labs = [String(sym) for i=1:1, sym in rn.syms]
#plot(sol.t, sol[:,:]', linetype=:steppost, labels=labs)
#plot(sol.t, sol[:,:]', labels=labs)