# import MultiScaleMapSampler

default_variables = Dict{String, Tuple{Union{Type, Nothing}, Any}}(
                        "base_seed" => (Int64, 1029472),
                        "mul_seed" => (Int64, 512339),
                        "steps" => (Int64, 4*10^7),
                        "gamma" => (Float64, 0.0),
                        "output_freq" => (Int64, 500),
                        "thread" => (Int64, 1),
                        "chamber" => (String, nothing),
                        "place_w" => (Float64, 1.0),
                        "cnty_w" => (Float64, 1.0),
                        "iso_weight" => (Float64, 0.0),
                        "popdev" => (Float64, 0.05),
                        "const_weight" => (Float64, 5.0),
                        "cousub_ousted_weight" => (Float64, 0.0),
                        "place_ousted_weight" => (Float64, 0.0),
                        "place_geoids_file" => (String, ""),
                        "target_maps" => (Int64, 0),
                        "max_districts_per_county" => (Int64, 2),
                    )

function parse_arg(
    args::Vector{String},
    arg::String
)
    if "--"*arg in args
        type = default_variables[arg][1]
        ind = findall(x->x=="--"*arg, args)[end]+1
        if type <: Real
            return parse(type, args[ind])
        elseif type <: String
            return args[ind]
        else
            println(arg, " is of unrecognized type ", type)
            @assert false
        end
    else
        if default_variables[arg][2] !== nothing
            return default_variables[arg][2]
        else
            println(arg, " must be specified")
            @assert false
        end
    end
end


function parse_districts(args::Vector{String})
    if "--num_dists" in args
        ind = findall(x->x=="--num_dists", args)[end]+1
        return parse(Int32, args[ind])
    end

    chamber = parse_arg(args, "chamber")

    if chamber == "lower"
        return 99
    elseif chamber == "upper"
        return 33
    elseif chamber == "congressional"
        return 15
    else
        println(" unrecognized chamber ", chamber)
        @assert false
    end

    println("Unspecified chamber or districts")
    @assert false
end


function parse_pop_dev(args::Vector{String})
    # Accept both spellings for compatibility with existing job files.
    popdev_flags = [flag for flag in ("--popdev", "--pop_dev") if flag in args]
    if !isempty(popdev_flags)
        flag = popdev_flags[end]
        ind = findall(x->x==flag, args)[end]+1
        return parse(Float64, args[ind])
    end

    chamber = parse_arg(args, "chamber")

    if chamber == "lower" || chamber == "upper"
        return 0.05
    elseif chamber == "congressional"
        return 0.02
    else
        println(" unrecognized chamber ", chamber)
        @assert false
    end

    println("Unspecified chamber or districts")
    @assert false
end


function parse_outfile(args::Vector{String})
    if "--outfile" in args
        ind = findall(x->x=="--outfile", args)[end]+1
        return args[ind]
    end

    thread_id = parse_arg(args, "thread")
    gamma = parse_arg(args, "gamma")
    popdev = parse_pop_dev(args)
    place_weight_scale = parse_arg(args, "place_w")
    cnty_weight_scale = parse_arg(args, "cnty_w")
    chamber = parse_arg(args, "chamber")
    outName = "run"*string(thread_id)
    base_seed = parse_arg(args, "base_seed")
    mul_seed = parse_arg(args, "mul_seed")
    iso_weight = parse_arg(ARGS, "iso_weight")
    const_weight = parse_arg(ARGS, "const_weight")

    if (chamber == "lower" || chamber == "upper") && popdev != 0.05
        outName *= "_popdev"*string(popdev)
    elseif chamber == "congressional" && popdev != 0.02
        outName *= "_popdev"*string(popdev)
    end

    if cnty_weight_scale != 1.0
        outName *= "_cntyw"*string(cnty_weight_scale)
    end
    if place_weight_scale != 1.0
        outName *= "_placew"*string(place_weight_scale)
    end
    if base_seed != 1029472
        outName *= "_baseseed"*string(base_seed)
    end
    if mul_seed != 512339
        outName *= "_mulseed"*string(mul_seed)
    end
    if gamma != 0.0
        outName *= "_gamma"*string(gamma)
    end
    if iso_weight != 0.0
        outName *= "_isoWeight"*string(iso_weight)
    end
    if const_weight != 0.0
        outName *= "_constWeight"*string(const_weight)
    end
    cousub_ousted_weight = parse_arg(args, "cousub_ousted_weight")
    if cousub_ousted_weight != 0.0
        outName *= "_cousubOustedWeight"*string(cousub_ousted_weight)
    end
    place_ousted_weight = parse_arg(args, "place_ousted_weight")
    if place_ousted_weight != 0.0
        outName *= "_placeOustedWeight"*string(place_ousted_weight)
    end
    if "--enforce_b7" in args
        outName *= "_B7"
    end

    outName *= ".jsonl";
    return outName
end


function parse_outpath(args::Vector{String})
    outName = parse_outfile(args)

    # Honor the explicit command-line argument if provided
    if "--outpath" in args
        ind = findall(x->x=="--outpath", args)[end]+1
        return joinpath(args[ind], outName)
    end

    # Default to the relative path architecture in repo
    chamber = parse_arg(args, "chamber")
    return joinpath("outputs", "chains", "soft_constraint", chamber, outName)
end


"""
Build a CycleWalk energy callback for the statewide share of population outside
the largest district piece of its Census county subdivision, in percentage points.
Each node contributes its full POP20 to its scalar COUSUB ID.
"""
function build_cousub_ousted_score(partition::LinkCutPartition)
    attrs = partition.graph.node_attributes
    ids = Dict{String, Int}()
    units = map(attrs) do node
        id = node["COUSUB"]
        id isa AbstractString && !isempty(id) ||
            throw(ArgumentError("Every graph node must have a scalar COUSUB ID"))
        get!(ids, id, length(ids) + 1)
    end
    pops = [Int(node["POP20"]) for node in attrs]
    total_pop = sum(pops)
    counts = zeros(Int, partition.num_dists, length(ids))

    # A unit's largest piece can be in a district outside the proposed move.
    function score(p, districts=nothing; update=nothing)
        assignment = update === nothing ? p.node_to_dist : p.node_to_dist_update
        fill!(counts, 0)
        @inbounds for i in eachindex(assignment, units, pops)
            counts[assignment[i], units[i]] += pops[i]
        end
        return 100 * (total_pop - sum(maximum, eachcol(counts))) / total_pop
    end

    return score
end


"""Restrict the Place dictionaries used by both the edge weights and the energy."""
function restrict_place_populations!(graph::BaseGraph, path::String)
    geoids = Set(strip(line) for line in readlines(path)
                 if !isempty(strip(line)) && !startswith(strip(line), "#"))
    !isempty(geoids) && all(id -> occursin(r"^39\d{5}$", id), geoids) ||
        throw(ArgumentError("Expected Ohio Place GEOIDs in $path"))
    present = Set{String}()
    for node in graph.node_attributes
        places = Dict(id => pop for (id, pop) in node["PLACE"] if id in geoids)
        union!(present, keys(places))
        node["PLACE"] = places
    end
    present == geoids || throw(ArgumentError("Place GEOIDs missing from graph: $(sort!(collect(setdiff(geoids, present))))"))
    return sort!(collect(geoids))
end


"""Score population outside each selected Place's largest district piece, as a percentage of statewide POP20."""
function build_place_ousted_score(partition::LinkCutPartition)
    ids = Dict{String, Int}()
    pieces = Tuple{Int, Int, Int}[]
    for (node, attrs) in enumerate(partition.graph.node_attributes)
        for (id, pop) in attrs["PLACE"]
            pop > 0 || continue
            push!(pieces, (node, get!(ids, id, length(ids) + 1), Int(pop)))
        end
    end
    isempty(pieces) && throw(ArgumentError("The selected Places have no population"))
    place_pop = sum(last, pieces)
    total_pop = partition.graph.total_pop
    counts = zeros(Int, partition.num_dists, length(ids))

    function score(p, districts=nothing; update=nothing)
        assignment = update === nothing ? p.node_to_dist : p.node_to_dist_update
        fill!(counts, 0)
        @inbounds for (node, place, pop) in pieces
            counts[assignment[node], place] += pop
        end
        return 100 * (place_pop - sum(maximum, eachcol(counts))) / total_pop
    end

    return score
end


"""Write the contracted cities' district labels under their original precinct keys."""
function build_precinct_output(graph::BaseGraph, expected_cities)
    city_nodes = filter(node -> haskey(node, "PROTECTED_CITY"), graph.node_attributes)
    length(city_nodes) == length(expected_cities) &&
        Set(node["PROTECTED_CITY"] for node in city_nodes) == Set(expected_cities) ||
        throw(ArgumentError("Expected one contracted node each for $(join(expected_cities, ", "))"))
    city_members = [((String(node["NAME"]),), [(String(name),) for name in node["MEMBER_NAMES"]])
                    for node in city_nodes]
    names = Set((String(node["NAME"]),) for node in graph.node_attributes)
    for (city, members) in city_members
        isempty(members) && throw(ArgumentError("No precinct members for $city"))
        for member in members
            member in names && throw(ArgumentError("Duplicate precinct key $member in $city"))
            push!(names, member)
        end
    end

    function write_map(partition, _measure, step, count, writer)
        map = CycleWalk.build_output_map(writer, partition, "step$(step - count)")
        for (city, members) in city_members
            district = pop!(map.districting, city)
            for member in members
                map.districting[member] = district
            end
        end
        CycleWalk.addMap(writer.atlas.io, map)
    end

    return write_map
end


function edge_weight_func(
    base_graph::BaseGraph,
    n1::Int, n2::Int;
    cnty_scale::Float64=1.0,
    place_scale::Float64=1.0,
    cnty_field::Union{String, Nothing}=nothing,
    place_field::Union{String, Nothing}=nothing
)
    @assert !(place_scale != 1.0 && place_field === nothing)
    @assert !(cnty_scale != 1.0 && cnty_field === nothing)
    @assert place_scale >= 1.0

    cur_weight = base_graph.simple_graph.weights[n1, n2]
    if cnty_field !== nothing && cnty_scale != 1.0
        @assert cnty_scale > 1
        # County weighting is orthogonal to Place weighting. Do not return here:
        # a Census Place can cross a county boundary.
        cnty1 = base_graph.node_attributes[n1][cnty_field]
        cnty2 = base_graph.node_attributes[n2][cnty_field]
        if cnty1 != cnty2
            cur_weight /= cnty_scale
        end
    end

    if place_field === nothing || place_scale == 1.0
        return cur_weight
    end

    # Restore the surviving South Carolina construction using Ohio's PLACE
    # population dictionaries. Each dictionary maps a Census Place GEOID to the
    # population of that Place inside the precinct.
    pop1 = Float64(base_graph.node_attributes[n1][base_graph.pop_col])
    pop2 = Float64(base_graph.node_attributes[n2][base_graph.pop_col])
    if pop1 <= 0 || pop2 <= 0
        return cur_weight
    end

    places1 = base_graph.node_attributes[n1][place_field]
    places2 = base_graph.node_attributes[n2][place_field]
    common_places = intersect(Set(keys(places1)), Set(keys(places2)))
    isempty(common_places) && return cur_weight

    shared_place_score = 0.0
    for place_geoid in common_places
        shared_place_score += (
            Float64(places1[place_geoid]) * Float64(places2[place_geoid])
        )
    end
    shared_place_score /= pop1 * pop2

    @assert -1e-12 <= shared_place_score <= 1.0 + 1e-12
    shared_place_score = clamp(shared_place_score, 0.0, 1.0)

    # Faithfully restore the surviving example's linear interpolation:
    # score=0 gives no Place multiplier; score=1 gives place_scale.
    place_multiplier = (
        place_scale * shared_place_score + 1.0 - shared_place_score
    )
    return cur_weight * place_multiplier
end


"""
    fixed_region_district_cap(graph, region, cap; num_dists)

Create a `CapRegionDistricts` check whose cap is the same for every region.
The constraint is intended for screening output and need not be added to the
sampler's constraints.
"""
function fixed_region_district_cap(
    graph,
    region::String,
    cap::Integer;
    num_dists::Integer,
)
    cap >= 1 || throw(ArgumentError("region district cap must be positive"))
    constraint = CapRegionDistricts(
        graph,
        region;
        num_dists=Int(num_dists),
        excess_split=0,
    )
    for region_name in keys(constraint.region_to_dist_cap)
        constraint.region_to_dist_cap[region_name] = Int(cap)
    end
    return constraint
end


"""
    build_original_county_output_screen(sampling_graph, assignment_graph_path;
                                        num_dists, max_districts_per_county)

Build a saved-state-only screen for Ohio Article XIX, Section 2(B)(7) and the
existing county district-count cap. It uses county identities and 2020 county
populations from the *original precinct graph*, never the scalar COUNTY field
on contracted city nodes. A protected-city node can represent precincts from
multiple counties; all of those counties must receive its district assignment.

A district pair violates B7 when it shares *at least two* counties with
population <= 400,000. Counties with populations > 400,000 are disregarded
for this particular rule. Return :pass, :county_cap, or :b7.

This is deliberately NOT a CycleWalk proposal constraint: the original MCMC
kernel remains free to visit states that will not be written to AtlasIO.
"""
function build_original_county_output_screen(
    sampling_graph::BaseGraph,
    assignment_graph_path::String;
    num_dists::Integer,
    max_districts_per_county::Integer,
)
    nd = Int(num_dists)
    1 <= nd <= 64 || throw(ArgumentError("County screen supports 1–64 districts"))
    cap = Int(max_districts_per_county)
    1 <= cap <= nd || throw(ArgumentError("Invalid county district-count cap"))

    # Load the original 7,928-node precinct graph only once at initialization.
    # A contracted city's scalar COUNTY attribute cannot safely be used here.
    original = BaseGraph(
        assignment_graph_path, "POP20";
        inc_node_data=Set(["NAME", "COUNTY", "POP20", "area", "border_length"]),
        area_col="area",
        node_border_col="border_length",
        edge_perimeter_col="length",
    )

    original_by_name = Dict{String, Tuple{String, Int}}()
    county_population = Dict{String, Int}()
    for node in original.node_attributes
        name = String(node["NAME"])
        county = String(node["COUNTY"])
        pop = Int(node["POP20"])
        pop >= 0 || throw(ArgumentError("Negative POP20 on original precinct $name"))
        haskey(original_by_name, name) &&
            throw(ArgumentError("Duplicate NAME in original graph: $name"))
        original_by_name[name] = (county, pop)
        county_population[county] = get(county_population, county, 0) + pop
    end
    isempty(original_by_name) && throw(ArgumentError("Original assignment graph is empty"))
    length(county_population) == 88 ||
        throw(ArgumentError("Expected 88 Ohio counties; got $(length(county_population))"))

    counties = sort!(collect(keys(county_population)))
    county_index = Dict(county => i for (i, county) in enumerate(counties))
    # For each sampling node, enumerate all original counties it intersects.
    # This representation accurately handles hard-contracted multi-county cities.
    county_ids_by_node = Vector{Vector{Int}}(undef, length(sampling_graph.node_attributes))
    seen_original_names = Set{String}()
    sampling_pop = 0

    for (i, node) in enumerate(sampling_graph.node_attributes)
        node_name = String(node["NAME"])
        names = if haskey(node, "PROTECTED_CITY")
            haskey(node, "MEMBER_NAMES") ||
                throw(ArgumentError("Protected city node $node_name lacks MEMBER_NAMES"))
            String.(node["MEMBER_NAMES"])
        else
            [String(node["NAME"])]
        end
        isempty(names) && throw(ArgumentError("Sampling node $i has no original members"))
        ids = Set{Int}()
        member_pop = 0
        for name in names
            haskey(original_by_name, name) ||
                throw(ArgumentError("Sampling member $name absent from original graph"))
            name in seen_original_names &&
                throw(ArgumentError("Original precinct $name appears in multiple sampling nodes"))
            push!(seen_original_names, name)
            county, pop = original_by_name[name]
            push!(ids, county_index[county])
            member_pop += pop
        end
        # Ensures that all county populations represented by a contracted node
        # really are included in the same district during screening.
        node_pop = Int(node["POP20"])
        node_pop == member_pop ||
            throw(ArgumentError("Population mismatch in sampling node $node_name: " *
                                "$node_pop vs original members $member_pop"))
        sampling_pop += member_pop
        county_ids_by_node[i] = sort!(collect(ids))
    end

    length(seen_original_names) == length(original_by_name) ||
        throw(ArgumentError("Sampling graph omits $(length(original_by_name)-length(seen_original_names)) original precincts"))
    sampling_pop == sum(values(county_population)) ||
        throw(ArgumentError("Sampling graph and original graph differ in total POP20"))

    # 2020 Census POP20, strictly greater than 400,000 is exempt from B7.
    nonexempt_counties = [i for (i, county) in enumerate(counties)
                          if county_population[county] <= 400_000]
    exempt_counties = [county for county in counties if county_population[county] > 400_000]
    expected_exempt_counties = Set(["CUYAHOGA", "FRANKLIN", "HAMILTON", "LUCAS",
                                   "MONTGOMERY", "SUMMIT"])
    Set(exempt_counties) == expected_exempt_counties ||
        throw(ArgumentError("Unexpected 2020 B7-exempt counties in original graph: $exempt_counties"))
    @info "Ohio B7 output screen initialized" original_nodes=length(original_by_name) sampling_nodes=length(county_ids_by_node) counties=length(counties) exempt_counties max_districts_per_county=cap

    county_masks = zeros(UInt64, length(counties))
    pair_shared_nonexempt = zeros(UInt8, nd, nd)
    function screen(p)
        assignment = p.node_to_dist
        length(assignment) == length(county_ids_by_node) ||
            throw(ArgumentError("Partition assignment length does not match sampling nodes"))
        fill!(county_masks, 0)
        @inbounds for node in eachindex(county_ids_by_node)
            district = Int(assignment[node])
            1 <= district <= nd ||
                throw(ArgumentError("Unexpected district label $district"))
            mask = UInt64(1) << (district - 1)
            for county in county_ids_by_node[node]
                county_masks[county] |= mask
            end
        end
        # First apply the existing output-only county district-count cap,
        # now evaluated against all constituent original precincts.
        @inbounds for mask in county_masks
            count_ones(mask) <= cap || return :county_cap
        end

        fill!(pair_shared_nonexempt, 0)
        # At most one nonexempt county may be shared by each district pair.
        @inbounds for county in nonexempt_counties
            mask = county_masks[county]
            for d1 in 1:(nd - 1)
                (mask & (UInt64(1) << (d1 - 1))) == 0 && continue
                for d2 in (d1 + 1):nd
                    (mask & (UInt64(1) << (d2 - 1))) == 0 && continue
                    pair_shared_nonexempt[d1, d2] += 1
                    pair_shared_nonexempt[d1, d2] <= 1 || return :b7
                end
            end
        end
        return :pass
    end
    return screen
end


"""
    run_until_compliant_quota!(partition, proposal, measure, rng,
                                compliance_constraint, writer;
                                output_freq, target_maps,
                                write_map=CycleWalk.output, screen_check=nothing)

Run the same CycleWalk Metropolis-Hastings transition kernel used by
`run_metropolis_hastings!`. At each ordinary output checkpoint, write the current
plan only when it satisfies `compliance_constraint`, or when `screen_check`
returns `:pass` if supplied. Noncompliant states remain in
the Markov chain. Stop after exactly `target_maps` maps have been written.

The initial plan is treated as the step-zero output checkpoint, matching
CycleWalk's standard writer behavior.
"""
function run_until_compliant_quota!(
    partition,
    proposal,
    measure,
    rng,
    compliance_constraint,
    writer;
    output_freq::Integer,
    target_maps::Integer,
    write_map=CycleWalk.output,
    screen_check=nothing,
)
    output_freq > 0 || throw(ArgumentError("output_freq must be positive"))
    target_maps > 0 || throw(ArgumentError("target_maps must be positive"))

    output_freq = Int(output_freq)
    target_maps = Int(target_maps)
    written_maps = 0
    screened_maps = 0
    rejected_maps = 0
    rejected_county_cap = 0
    rejected_b7 = 0
    total_steps = 0

    function screen_and_maybe_write!(step::Int)
        screened_maps += 1
        verdict = if screen_check === nothing
            satisfies_constraint(partition, compliance_constraint) ? :pass : :county_cap
        else
            screen_check(partition)
        end
        verdict in (:pass, :county_cap, :b7) ||
            throw(ArgumentError("Unknown output-screen verdict: $verdict"))
        if verdict !== :pass
            rejected_maps += 1
            if verdict === :county_cap
                rejected_county_cap += 1
            else
                rejected_b7 += 1
            end
            return false
        end

        write_map(partition, measure, step, 0, writer)
        written_maps += 1
        if written_maps == 1 || written_maps % 50 == 0 || written_maps == target_maps
            @info "Compliant Atlas progress" written_maps target_maps screened_maps rejected_maps rejected_county_cap rejected_b7 step
        end
        return true
    end

    screen_and_maybe_write!(0)

    while written_maps < target_maps
        # Run without a writer so all proposals and accept/reject decisions are
        # unchanged. The current state is screened only after the normal output
        # interval has elapsed.
        run_metropolis_hastings!(
            partition,
            proposal,
            measure,
            output_freq,
            rng;
            writer=nothing,
            output_freq=output_freq,
        )
        total_steps += output_freq
        screen_and_maybe_write!(total_steps)
    end

    return (
        written_maps=written_maps,
        screened_maps=screened_maps,
        rejected_maps=rejected_maps,
        rejected_county_cap=rejected_county_cap,
        rejected_b7=rejected_b7,
        total_steps=total_steps,
    )
end
