module ParseRxNetwork

using DataStructures, DiffEqBiological

abstract type NetworkFileFormat end

struct RSSANetwork <: NetworkFileFormat end
struct BNGNetwork <: NetworkFileFormat end

include("parsing_routines_rssafiles.jl")
include("parsing_routines_bngnetworkfiles.jl")

export RSSANetwork, BNGNetwork
export get_rxnetwork_simple

end # module
