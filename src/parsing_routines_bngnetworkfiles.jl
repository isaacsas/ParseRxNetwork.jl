"""
Parsing routines for BioNetGen .net files. Note, this only handles a subset
of the BioNetGen .net file format. In particular, any functions in the file
that use non-Julia expressions (i.e. like conditional logic) or time-dependent
features, will definitely not work. 
"""


# function build_rxnetwork(ft::BNGNetwork, networkname, rxstrs, rxrates; printrxs=false, kwargs...)
    
#     # string representing the network
#     rnstr = "@min_reaction_network $(networkname) begin\n"
#     for (i,rxstr) in enumerate(rxstrs)
#         rate = rxrates[i]
#         rnstr *= "\t $rate, $rxstr \n"
#     end
#     rnstr *= "end\n"

#     if printrxs
#         print(rnstr)
#     end

#     # build the network using DiffEqBiological
#     rn = eval( Meta.parse(rnstr) )

#     rn
# end

"""
Seek to a line matching `start_string`.
"""
function seek_to_block(lines, idx, start_string)
    while lines[idx] != start_string
        idx += 1        
        (idx > length(lines)) && error("Block: ", start_string, " was never found.")
    end
    idx += 1
end

const PARAM_BLOCK_START = "begin parameters"
const PARAM_BLOCK_END = "end parameters"
function parse_params(ft::BNGNetwork, lines, idx)

    idx = seek_to_block(lines, idx, PARAM_BLOCK_START)

    # parse params
    pvals = []
    ptoids = Dict{Symbol,Int}()
    while lines[idx] != PARAM_BLOCK_END
        vals = split(lines[idx])
        pidx = parse(Int,vals[1])
        psym = Symbol(vals[2])
        ptoids[psym] = pidx    

        # value could be an expression
        pval = Meta.parse(vals[3])    
        push!(pvals,pval)

        idx += 1
        (idx > length(lines)) && error("Block: ", PARAM_BLOCK_END, " was never found.")
    end

    ptoids,pvals,idx
end

const SPECIES_BLOCK_START = "begin species"
const SPECIES_BLOCK_END = "end species"
function parse_species(ft::BNGNetwork, lines, idx, ptoids, pvals)

    idx = seek_to_block(lines, idx, SPECIES_BLOCK_START)
    println("idx = ", idx, ", line is:", lines[idx])

    # parse species
    symstoids = Dict{Symbol,Int}()
    u0exprs   = Vector{Any}()
    while lines[idx] != SPECIES_BLOCK_END
        vals = split(lines[idx])
        sidx = parse(Int,vals[1])
        ssym = Symbol(vals[2])
        symstoids[ssym] = sidx
        push!(u0exprs,  Meta.parse(vals[3]))
        idx += 1
        (idx > length(lines)) && error("Block: ", SPECIES_BLOCK_END, " was never found.")
    end

    u0exprs,symstoids,idx
end

const REACTIONS_BLOCK_START = "begin reactions"
const REACTIONS_BLOCK_END = "end reactions"
function parse_reactions(ft, lines, idx, ptoids, pvals, symstoids)

    # reverse Dicts
    idstosyms = Vector{Symbol}(undef,length(symstoids))
    for (k,v) in symstoids
        idstosyms[v] = k
    end

    idstops = Vector{Symbol}(undef,length(ptoids))
    for (k,v) in ptoids
        idstops[v] = k
    end

    idx = seek_to_block(lines, idx, REACTIONS_BLOCK_START)

    idtosymstr = id -> (id==0) ? ∅ : string(idstosyms[id])

    rxstrs = "begin\n"
    while lines[idx] != REACTIONS_BLOCK_END
        vals        = split(lines[idx])        
        reactantids = (parse(Int,rsym) for rsym in split(vals[2],","))
        productids  = (parse(Int,psym) for psym in split(vals[3],","))
        rateconst   = Meta.parse(vals[4])
        
        pstr = "" # TODO
        rstr = join((idtosymstr(rid) for rid in reactantids), "+")
        pstr = join((idtosymstr(pid) for pid in productids), "+")
        
        # create string for this reaction
        push!(rxstrs, string(", "))

        idx += 1
        (idx > length(lines)) && error("Block: ", REACTIONS_BLOCK_END, " was never found.")
    end
end


# for parsing a subset of the BioNetGen .net file format
function loadrxnetwork(ft::BNGNetwork, networkname, rxfilename; kwargs...)

    file  = open(rxfilename, "r");
    lines = readlines(file)
    idx   = 1
    ptoids,pvals,idx      = parse_params(ft, lines, idx)
    u0exprs,symstoids,idx = parse_species(ft, lines, idx, ptoids, pvals)
    rxstrs,idx            = parse_reactions(ft, lines, idx, ptoids, pvals, symstoids)
    close(file)
    
    # build the DiffEqBiological representation of the network
    #rn = build_rxnetwork(ft, networkname, rxstrs, rxrates; kwargs...)

    #rn, initialpop
end