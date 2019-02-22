"""
Parsing routines for BioNetGen .net files. Note, this only handles a subset
of the BioNetGen .net file format. In particular, any functions in the file
that use non-Julia expressions (i.e. like conditional logic) or time-dependent
features, will definitely not work. 
"""

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

#Recursively traverses an expression and replaces a symbol with another symbol or expression
# here psyms stores the allowed symbols in the final expression
function recursive_replace!(expr::Any, ptoids, pvals, psyms)
    if typeof(expr) == Symbol
        if expr ∉ psyms
            haskey(ptoids,expr) && return pvals[ptoids[expr]]
        end
    elseif typeof(expr) == Expr
        for i = 1:length(expr.args)
            expr.args[i] = recursive_replace!(expr.args[i], ptoids, pvals, psyms)
        end
    end
    return expr
end

# recursively traverses an expression and replaces all Exprs and Symbols with numerical values
function recursive_replace!(expr::Any, ptoids, pvals)
    if typeof(expr) == Symbol
        haskey(ptoids,expr) && return recursive_replace!(pvals[ptoids[expr]], ptoids, pvals)
    elseif typeof(expr) == Expr
        for i = 1:length(expr.args)
            expr.args[i] = recursive_replace!(expr.args[i], ptoids, pvals, psyms)
        end
    end
    @assert typeof(expr) <: Number string("Error, expr type was not a number.")
    return expr
end


const PARAM_BLOCK_START = "begin parameters"
const PARAM_BLOCK_END = "end parameters"
function parse_params(ft::BNGNetwork, lines, idx)

    idx = seek_to_block(lines, idx, PARAM_BLOCK_START)

    # parse params
    pvals = []
    ptoids = Dict{Symbol,Int}()
    psyms = OrderedSet{Symbol}()
    while lines[idx] != PARAM_BLOCK_END
        vals = split(lines[idx])
        pidx = parse(Int,vals[1])
        psym = Symbol(vals[2])
        ptoids[psym] = pidx    

        # value could be an expression
        pval = Meta.parse(vals[3])
        
        # save numeric parameters
        if vals[end] == "Constant"
            push!(psyms, psym)
        else   # replace non-numeric Symbols or Exprs
            if (typeof(pval) <: Symbol) 
                (pval ∉ psyms) && (pval = pvals[ptoids[pval]])
            else
                @assert typeof(pval) <: Expr string("Expected expression for parameter value but got: ", typeof(pval), ", line idx = ", idx)
                pval = recursive_replace!(pval, ptoids, pvals, psyms)
            end
        end

        push!(pvals,pval)

        idx += 1
        (idx > length(lines)) && error("Block: ", PARAM_BLOCK_END, " was never found.")
    end

    ptoids,pvals,psyms,idx
end

stripchars = (s, r) -> replace(s, Regex("[$r]") => "")
const CHARS_TO_STRIP = ",~!()."
const SPECIES_BLOCK_START = "begin species"
const SPECIES_BLOCK_END = "end species"
function parse_species(ft::BNGNetwork, lines, idx, ptoids, pvals)

    idx = seek_to_block(lines, idx, SPECIES_BLOCK_START)

    # parse species
    symstoids = Dict{Symbol,Int}()
    u0exprs   = Vector{Any}()
    while lines[idx] != SPECIES_BLOCK_END
        vals = split(lines[idx])
        sidx = parse(Int,vals[1])
        ssym = Symbol(stripchars(vals[2], CHARS_TO_STRIP))
        symstoids[ssym] = sidx
        push!(u0exprs,  Meta.parse(vals[3]))
        idx += 1
        (idx > length(lines)) && error("Block: ", SPECIES_BLOCK_END, " was never found.")
    end

    u0exprs,symstoids,idx
end

const REACTIONS_BLOCK_START = "begin reactions"
const REACTIONS_BLOCK_END = "end reactions"
function parse_reactions!(rxiobuf, ft, lines, idx, ptoids, pvals, psyms, symstoids)

    # reverse Dicts
    idstosyms = Vector{Symbol}(undef,length(symstoids))
    for (k,v) in symstoids
        idstosyms[v] = k
    end

    idtosymstr = id -> (id==0) ? ∅ : string(idstosyms[id])
    idx = seek_to_block(lines, idx, REACTIONS_BLOCK_START)
    write(rxiobuf, "begin\n")
    while lines[idx] != REACTIONS_BLOCK_END
        vals        = split(lines[idx])        
        reactantids = (parse(Int,rsym) for rsym in split(vals[2],","))
        productids  = (parse(Int,psym) for psym in split(vals[3],","))
        rateconst   = Meta.parse(vals[4])
                
        rctype = typeof(rateconst)
        if rctype <: Number
            pstr = string(rateconst)
        elseif rctype <: Symbol
            # if true constant use symbol
            if rateconst ∈ psyms
                pstr = string(rateconst)
            else
                pstr = string(pvals[ptoids[rateconst]])
            end
        elseif rctype <: Expr
            rateconst = recursive_replace!(rateconst, ptoids, pvals, psyms)
        else
            error(string("Rate constant in reaction is not a Number, Symbol or Expr, at reaction: ", vals[1]))
        end

        reactstr = join((idtosymstr(rid) for rid in reactantids), " + ")
        productstr = join((idtosymstr(pid) for pid in productids), " + ")
        
        # create string for this reaction
        write(rxiobuf, "  ", pstr, ", ", reactstr, " --> ", productstr, "\n")

        idx += 1
        (idx > length(lines)) && error("Block: ", REACTIONS_BLOCK_END, " was never found.")
    end
    write(rxiobuf, "end")

    idx
end


# for parsing a subset of the BioNetGen .net file format
function loadrxnetwork(ft::BNGNetwork, networkname, rxfilename; kwargs...)

    file  = open(rxfilename, "r");
    lines = readlines(file)
    idx   = 1
    println("Parsing parameters...")
    ptoids,pvals,psyms,idx = parse_params(ft, lines, idx)
    println("done")
    println("Parsing species...")
    u0exprs,symstoids,idx = parse_species(ft, lines, idx, ptoids, pvals)
    println("done")
    println("Parsing reactions...")
    rxiobuf = IOBuffer()
    write(rxiobuf, string(networkname, " = @min_reaction_network "))
    idx = parse_reactions!(rxiobuf, ft, lines, idx, ptoids, pvals, psyms, symstoids)
    println("done")
    close(file)
    foreach(psym -> write(rxiobuf, " ", string(psym)), psyms)
    write(rxiobuf, "\n")
    rxstrs = String(take!(rxiobuf))

    # get the initial condition
    u₀ = [Float64(recursive_replace!(u0expr, ptoids, pvals)) for u0expr in u0exprs]

    # build the DiffEqBiological representation of the network
    #rn = rxstrs
    rn = eval(Meta.parse(rxstrs))    

    # parameters values for constant params
    p = [pvals[ptoids[psym]] for psym in psyms]

    rn,u₀,p
end