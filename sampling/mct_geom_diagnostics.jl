"""Additional edge-geography split diagnostics, including zero-population pieces.

Uses positive-area node/geography incidence derived from Census-block polygons.
Each unit touching nodes in >=2 districts is counted as geometrically split.
This is NOT the population-based unit-split statistic; the official plans must
be scored on the same geometry before direct official comparisons.
"""
function build_mct_geometry_split_diagnostics(partition, incidence_file::AbstractString)
    isfile(incidence_file) || throw(ArgumentError("Missing incidence file: $incidence_file"))
    attrs = partition.graph.node_attributes
    nodes = Dict{String,Int}()
    for (i,a) in enumerate(attrs)
        name=String(a["NAME"])
        haskey(nodes,name) && throw(ArgumentError("Duplicate graph NAME: $name"))
        nodes[name]=i
    end
    town_to_index=Dict{String,Int}()
    place_to_index=Dict{String,Int}()
    town_members=Vector{Tuple{Int,Int}}()
    place_members=Vector{Tuple{Int,Int}}()
    seen=Set{Tuple{Int,String,String}}()
    open(incidence_file,"r") do io
        eof(io) && throw(ArgumentError("Empty geometric incidence"))
        split(chomp(readline(io)), '\t') == ["node_name","geography","geoid"] ||
            throw(ArgumentError("Bad incidence TSV header"))
        for line in eachline(io)
            pieces=split(chomp(line),'\t';keepempty=true)
            length(pieces)==3 || throw(ArgumentError("Bad incidence row"))
            name,kind,id=pieces
            haskey(nodes,name) || throw(ArgumentError("Unknown graph node $name"))
            key=(nodes[name],kind,id)
            key in seen && throw(ArgumentError("Duplicate incidence $key"))
            push!(seen,key)
            if kind=="township"
                push!(town_members,(nodes[name],get!(town_to_index,id,length(town_to_index)+1)))
            elseif kind=="place"
                push!(place_members,(nodes[name],get!(place_to_index,id,length(place_to_index)+1)))
            else
                throw(ArgumentError("Bad geography kind $kind"))
            end
        end
    end
    function count_splits(pairs::Vector{Tuple{Int,Int}}, n_units::Int, assignments)
        # Number of DIFFERENT districts touching each geometry (not count edges).
        first_seen=zeros(Int,n_units)
        split= falses(n_units)
        for (node,unit) in pairs
            d=Int(assignments[node])
            1<=d<=partition.num_dists || throw(ArgumentError("Bad district $d"))
            if first_seen[unit]==0
                first_seen[unit]=d
            elseif first_seen[unit]!=d
                split[unit]=true
            end
        end
        return count(identity,split)
    end
    function township_splits(p,districts=nothing;update=nothing)
        assignment=update===nothing ? p.node_to_dist : p.node_to_dist_update
        return count_splits(town_members,length(town_to_index),assignment)
    end
    function cross_place_splits(p,districts=nothing;update=nothing)
        assignment=update===nothing ? p.node_to_dist : p.node_to_dist_update
        return count_splits(place_members,length(place_to_index),assignment)
    end
    return (township_splits=township_splits,
            cross_place_splits=cross_place_splits,
            n_townships=length(town_to_index),
            n_cross_places=length(place_to_index),
            n_memberships=length(seen))
end
