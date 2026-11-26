"""
    build_mct_energy_scores(partition, atoms_path; township_weight, cross_place_weight)

Build a township-first, geographically disjoint MCT energy from a validated
population-atom TSV produced by scripts/mct/build_mct_atoms.py.

The two penalized unit collections are:
  * qualifying township COUSUBs (CLASSFP T*, LSAD 44)
  * incorporated cities/villages not fully contained in one township

Nested incorporated places are *not* independently penalized. Their ousted
population remains available as a diagnostic to expose the tradeoff when a
parent township is split.

The TSV atoms assign each resident to exactly one (node, township, place)
category and populations are reconciled to every graph node. This does not
count municipal ousting twice for overlapping municipalities; however the
*additive* township and cross-municipality energies can penalize the same
resident under both geographic objectives when applicable. This is intended.
"""
function build_mct_energy_scores(
    partition,
    atoms_path::AbstractString;
    township_weight::Real,
    cross_place_weight::Real,
)
    for (label, weight) in (("township_weight", township_weight),
                            ("cross_place_weight", cross_place_weight))
        (isfinite(weight) && weight >= 0) ||
            throw(ArgumentError("$label must be finite and nonnegative"))
    end
    isfile(atoms_path) || throw(ArgumentError("MCT atoms TSV missing: $atoms_path"))

    nodes = partition.graph.node_attributes
    name_to_index = Dict{String, Int}()
    node_pop = zeros(Int, length(nodes))
    for (i, node) in enumerate(nodes)
        name = String(node["NAME"])
        haskey(name_to_index, name) && throw(ArgumentError("Duplicate graph NAME: $name"))
        name_to_index[name] = i
        node_pop[i] = Int(node["POP20"])
        node_pop[i] >= 0 || throw(ArgumentError("Negative graph node population at $name"))
    end

    township_ids = Dict{String, Int}()
    cross_place_ids = Dict{String, Int}()
    nested_place_ids = Dict{String, Int}()
    # Sparse, aggregated tuples in graph-node order; no full county-subdivision
    # matrix is needed in the TSV. Empty ID is represented by integer 0.
    atom_node = Int[]
    atom_township = Int[]
    atom_cross = Int[]
    atom_nested = Int[]
    atom_population = Int[]
    observed_node_pop = zeros(Int, length(nodes))
    row_keys = Set{Tuple{Int, String, String}}()

    open(atoms_path, "r") do io
        eof(io) && throw(ArgumentError("Empty MCT atoms TSV"))
        expected = ["node_name", "township_geoid", "place_geoid", "place_class", "population"]
        split(chomp(readline(io)), '\t'; keepempty=true) == expected ||
            throw(ArgumentError("Unexpected MCT atoms TSV header"))
        for (line_number, line) in enumerate(eachline(io))
            cols = split(chomp(line), '\t'; keepempty=true)
            length(cols) == 5 || throw(ArgumentError("MCT atoms row $(line_number+1): expected five columns"))
            node_name, t_geoid, p_geoid, place_class, pop_text = cols
            haskey(name_to_index, node_name) ||
                throw(ArgumentError("Unknown MCT graph NAME $node_name"))
            i = name_to_index[node_name]
            place_class in ("none", "nested", "cross") ||
                throw(ArgumentError("Unknown place_class $place_class"))
            (isempty(p_geoid) == (place_class == "none")) ||
                throw(ArgumentError("Place identifier and place_class disagree at node $node_name"))
            isempty(t_geoid) || occursin(r"^39\d{8}$", t_geoid) ||
                throw(ArgumentError("Invalid Ohio COUSUB GEOID $t_geoid"))
            isempty(p_geoid) || occursin(r"^39\d{5}$", p_geoid) ||
                throw(ArgumentError("Invalid Ohio PLACE GEOID $p_geoid"))
            pop = parse(Int, pop_text)
            pop > 0 || throw(ArgumentError("MCT atom populations must be positive"))
            rowkey = (i, t_geoid, p_geoid)
            rowkey in row_keys && throw(ArgumentError("Repeated graph-node/unit intersection $rowkey"))
            push!(row_keys, rowkey)
            observed_node_pop[i] += pop
            push!(atom_node, i)
            push!(atom_township, isempty(t_geoid) ? 0 : get!(township_ids, t_geoid, length(township_ids)+1))
            push!(atom_cross, place_class == "cross" ?
                get!(cross_place_ids, p_geoid, length(cross_place_ids)+1) : 0)
            push!(atom_nested, place_class == "nested" ?
                get!(nested_place_ids, p_geoid, length(nested_place_ids)+1) : 0)
            push!(atom_population, pop)
        end
    end
    isempty(atom_node) && throw(ArgumentError("MCT atoms TSV has no populated rows"))
    for (i, pop) in enumerate(node_pop)
        observed_node_pop[i] == pop ||
            throw(ArgumentError("MCT atoms population mismatch at $(nodes[i]["NAME"]): " *
                                "atoms=$(observed_node_pop[i]), graph=$pop"))
    end
    total_pop = sum(node_pop)
    total_pop == Int(partition.graph.total_pop) ||
        throw(ArgumentError("Graph and partition total populations disagree"))
    total_pop > 0 || throw(ArgumentError("MCT population is zero"))

    # A unit's largest district can be anywhere in the state; recompute the
    # full rows rather than restricting the maximization to moved districts.
    nd = partition.num_dists
    township_counts = zeros(Int, nd, length(township_ids))
    cross_counts = zeros(Int, nd, length(cross_place_ids))
    nested_counts = zeros(Int, nd, length(nested_place_ids))
    township_pop = sum(atom_population[k] for k in eachindex(atom_population)
                       if atom_township[k] != 0)
    cross_pop = sum(atom_population[k] for k in eachindex(atom_population)
                    if atom_cross[k] != 0)
    nested_pop = sum(atom_population[k] for k in eachindex(atom_population)
                     if atom_nested[k] != 0)

    # Count a geographic unit as cut when its represented positive population
    # occupies at least two districts. Do not count the number of boundary
    # segments, which is a different diagnostic.
    function count_split_units(counts)
        nsplit = 0
        for unit in axes(counts, 2)
            occupied = 0
            for district in axes(counts, 1)
                if counts[district, unit] > 0
                    occupied += 1
                    occupied > 1 && break
                end
            end
            nsplit += occupied > 1
        end
        return nsplit
    end

    function calculate_scores(p, update; with_splits::Bool=false)
        assignment = update === nothing ? p.node_to_dist : p.node_to_dist_update
        length(assignment) == length(nodes) ||
            throw(ArgumentError("Partition and MCT graph lengths disagree"))
        fill!(township_counts, 0)
        fill!(cross_counts, 0)
        fill!(nested_counts, 0)
        @inbounds for k in eachindex(atom_node)
            d = Int(assignment[atom_node[k]])
            1 <= d <= nd || throw(ArgumentError("Invalid assigned district $d"))
            pop = atom_population[k]
            t = atom_township[k]
            c = atom_cross[k]
            n = atom_nested[k]
            t == 0 || (township_counts[d, t] += pop)
            c == 0 || (cross_counts[d, c] += pop)
            n == 0 || (nested_counts[d, n] += pop)
        end
        township_ousted = township_pop - sum(maximum, eachcol(township_counts); init=0)
        cross_ousted = cross_pop - sum(maximum, eachcol(cross_counts); init=0)
        nested_ousted = nested_pop - sum(maximum, eachcol(nested_counts); init=0)
        scale = 100.0 / total_pop
        return (scale * township_ousted, scale * cross_ousted, scale * nested_ousted,
                with_splits ? count_split_units(township_counts) : 0,
                with_splits ? count_split_units(cross_counts) : 0,
                with_splits ? count_split_units(nested_counts) : 0)
    end

    # The MH kernel receives one combined callback: the two pieces share a
    # single pass through the atoms. The writer can inspect components at
    # output checkpoints without changing proposal acceptance behavior.
    function combined(p, districts=nothing; update=nothing)
        t, c, _, _, _, _ = calculate_scores(p, update)
        return township_weight * t + cross_place_weight * c
    end
    function township(p, districts=nothing; update=nothing)
        t, _, _, _, _, _ = calculate_scores(p, update)
        return t
    end
    function cross_place(p, districts=nothing; update=nothing)
        _, c, _, _, _, _ = calculate_scores(p, update)
        return c
    end
    function nested_place(p, districts=nothing; update=nothing)
        _, _, n, _, _, _ = calculate_scores(p, update)
        return n
    end
    function township_splits(p, districts=nothing; update=nothing)
        _, _, _, t, _, _ = calculate_scores(p, update; with_splits=true)
        return t
    end
    function cross_place_splits(p, districts=nothing; update=nothing)
        _, _, _, _, c, _ = calculate_scores(p, update; with_splits=true)
        return c
    end
    function nested_place_splits(p, districts=nothing; update=nothing)
        _, _, _, _, _, n = calculate_scores(p, update; with_splits=true)
        return n
    end
    return (
        combined=combined,
        township=township,
        cross_place=cross_place,
        nested_place=nested_place,
        township_splits=township_splits,
        cross_place_splits=cross_place_splits,
        nested_place_splits=nested_place_splits,
        n_atoms=length(atom_node),
        n_townships=length(township_ids),
        n_cross_places=length(cross_place_ids),
        n_nested_places=length(nested_place_ids),
        accounted_population=total_pop,
    )
end
