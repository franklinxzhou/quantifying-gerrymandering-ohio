using RandomNumbers
using CycleWalk
include(joinpath(@__DIR__, "aux_place.jl"))

@show ARGS

cycle_walk_frac = 0.1

### Parse variables ###
base_seed = parse_arg(ARGS, "base_seed")
mul_seed = parse_arg(ARGS, "mul_seed")
steps = parse_arg(ARGS, "steps")
gamma = parse_arg(ARGS, "gamma")
output_freq = parse_arg(ARGS, "output_freq")
thread_id = parse_arg(ARGS, "thread")
chamber = parse_arg(ARGS, "chamber")
place_w_scale = parse_arg(ARGS, "place_w")
cnty_w_scale = parse_arg(ARGS, "cnty_w")
iso_weight = parse_arg(ARGS, "iso_weight")
num_dists = parse_districts(ARGS) 
popdev = parse_pop_dev(ARGS) 
outPath = parse_outpath(ARGS)
const_weight = parse_arg(ARGS, "const_weight")
cousub_ousted_weight = parse_arg(ARGS, "cousub_ousted_weight")
place_ousted_weight = parse_arg(ARGS, "place_ousted_weight")
place_geoids_file = parse_arg(ARGS, "place_geoids_file")
preserve_five_cities = "--preserve_five_cities" in ARGS
preserve_cities = preserve_five_cities || "--preserve_cities" in ARGS
target_maps = parse_arg(ARGS, "target_maps")
max_districts_per_county = parse_arg(ARGS, "max_districts_per_county")
# Opt-in to the additional Article XIX §2(B)(7) saved-map filter.
# This does not restrict proposals, only which checkpoints enter AtlasIO.
enforce_b7 = "--enforce_b7" in ARGS
# Opt-in: split PLACE populations by COUNTY without changing the five hard-
# contracted city nodes or the saved-map B7 output filter.
place_county_units = "--place_county_units" in ARGS
# Dedicated Columbus-only pilot. Require opt-in and prevent accidental use
# with an obsolete full-incorporated-place configuration.
columbus_only = "--columbus_only" in ARGS
columbus_only || throw(ArgumentError("This dedicated runner requires --columbus_only"))
place_county_units && preserve_five_cities && enforce_b7 ||
    throw(ArgumentError("Columbus-only runner requires --place_county_units --preserve_five_cities --enforce_b7"))
selected_geoids = Set(strip(split(line, '#', limit=2)[1])
    for line in readlines(place_geoids_file) if !isempty(strip(split(line, '#', limit=2)[1])))
selected_geoids == Set(["3918000"]) ||
    throw(ArgumentError("Columbus-only runner requires a GEOID file containing only 3918000"))
if place_county_units
    preserve_five_cities || throw(ArgumentError(
        "--place_county_units requires --preserve_five_cities"))
    !isempty(place_geoids_file) || throw(ArgumentError(
        "--place_county_units requires --place_geoids_file for the intended incorporated-place ensemble"))
end

target_maps >= 0 || throw(ArgumentError("target_maps must be nonnegative"))
if enforce_b7
    chamber == "congressional" || throw(ArgumentError("--enforce_b7 requires --chamber congressional"))
    target_maps > 0 || throw(ArgumentError("--enforce_b7 requires a positive --target_maps"))
end
output_freq > 0 || throw(ArgumentError("output_freq must be positive"))
isfinite(cousub_ousted_weight) && cousub_ousted_weight >= 0 ||
    throw(ArgumentError("cousub_ousted_weight must be finite and nonnegative"))
isfinite(place_ousted_weight) && place_ousted_weight >= 0 ||
    throw(ArgumentError("place_ousted_weight must be finite and nonnegative"))
if preserve_cities
    target_maps > 0 || throw(ArgumentError("City preservation requires a compliant-map quota"))
    cousub_ousted_weight == 0 || throw(ArgumentError("The contracted graph has mixed COUSUB nodes. Use Place energy instead."))
end
if target_maps > 0
    chamber == "congressional" ||
        throw(ArgumentError("compliant-map quotas require chamber congressional"))
    1 <= max_districts_per_county <= num_dists ||
        throw(ArgumentError("max_districts_per_county must be between 1 and num_dists"))
    ispath(outPath) && error("Refusing to overwrite existing ensemble: $outPath")
    mkpath(dirname(outPath))
end

steps = Int(steps/cycle_walk_frac)
output_freq = Int(output_freq/cycle_walk_frac)
@show steps, output_freq
@show iso_weight

ad_param = Dict{String, Any}(
    "place_w_scale"=>place_w_scale,
    "place_weight_model"=>"linear_population_overlap",
    "place_population_field" => place_county_units ? "PLACE_COUNTY" : "PLACE",
    "place_unit_definition" => place_county_units ? "PLACE_GEOID|COUNTY_NAME" : "PLACE_GEOID",
    "place_county_units" => place_county_units,
    "columbus_only" => true,
    "edge_place_field" => "PLACE_COUNTY",
    "energy_place_field" => "PLACE_COUNTY",
    "energy_unit_definition" => "Columbus GEOID|COUNTY_NAME; each county component has its own largest district piece",
    "energy_score_key" => "place_ousted_pp",
    "edge_place_units" => ["3918000|FRANKLIN", "3918000|DELAWARE", "3918000|FAIRFIELD"],
    "energy_place_units" => ["3918000|FRANKLIN", "3918000|DELAWARE", "3918000|FAIRFIELD"],
    "cnty_w_scale"=>cnty_w_scale,
    "chamber" => chamber,
    "popdev" => popdev
)
if target_maps > 0
    merge!(ad_param, Dict{String, Any}(
        "nominal_steps" => steps,
        "target_compliant_maps" => target_maps,
        "max_districts_per_county" => max_districts_per_county,
        "output_filter" => enforce_b7 ?
            "original precinct counties <= $(max_districts_per_county) districts; Article XIX 2(B)(7) district-pair shared nonexempt counties <= 1" :
            "each county touches at most $(max_districts_per_county) districts",
        "enforce_b7" => enforce_b7,
        "b7_population_threshold" => 400_000,
        "b7_population_basis" => "2020 Census POP20 in original assignment_graph",
        "b7_screening_scope" => "output checkpoints only; MCMC proposals unchanged",
        "output_freq" => output_freq,
        "cycle_walk_frac" => cycle_walk_frac,
        "base_seed" => base_seed,
        "mul_seed" => mul_seed,
        "thread" => thread_id,
        "rng_seed" => base_seed + mul_seed * thread_id,
    ))
end



# Check the Columbus graph's mathematical invariants BEFORE initializing MCMC.
function validate_columbus_only_graph!(g::BaseGraph)
    expected_units = Dict("3918000|FRANKLIN" => 880830,
                          "3918000|DELAWARE" => 14570,
                          "3918000|FAIRFIELD" => 10849)
    units = Dict{String,Int}()
    protected = Set{String}()
    city_total = 0
    statewide_total = 0
    length(g.node_attributes) == 6975 || throw(ArgumentError("Wrong number of five-city nodes"))
    for node in g.node_attributes
        if haskey(node, "PROTECTED_CITY")
            push!(protected, String(node["PROTECTED_CITY"]))
        end
        pop = Int(node["POP20"])
        statewide_total += pop
        pp = node["PLACE"]
        pc = node["PLACE_COUNTY"]
        label = String(node["NAME"])
        all(k -> String(k)=="3918000", keys(pp)) ||
            throw(ArgumentError("Non-Columbus PLACE present at " * label))
        all(k -> String(k) in keys(expected_units), keys(pc)) ||
            throw(ArgumentError("Non-Columbus PLACE_COUNTY present at " * label))
        whole = sum((Int(v) for v in values(pp)); init=0)
        parts = sum((Int(v) for v in values(pc)); init=0)
        whole == parts && whole <= pop ||
            throw(ArgumentError("PLACE/PLACE_COUNTY population mismatch at " * label))
        city_total += whole
        for (key,val) in pc
            units[String(key)] = get(units,String(key),0) + Int(val)
        end
    end
    protected == Set(["Akron","Cincinnati","Cleveland","Dayton","Toledo"]) ||
        throw(ArgumentError("Five protected cities changed"))
    units == expected_units || throw(ArgumentError("Columbus county totals differ: $units"))
    city_total == 906249 && statewide_total == 11799448 ||
        throw(ArgumentError("Wrong population totals for Columbus/Ohio"))
    @info "Columbus-only graph PASS" statewide_total city_total units protected
    return nothing
end

### Set graph ###
assignment_graph = joinpath("data", "graph", "OH_pct20.json")
pctGraphPath = joinpath("data", "graph", "OH_pct20_preserved_five_cities_columbus_only.json")
protected_cities = ["Akron", "Cincinnati", "Cleveland", "Dayton", "Toledo"]
place_field = "PLACE_COUNTY"
nodeData = Set(["NAME", "COUNTY", "POP20", "PLACE", "area", "border_length"]);
push!(nodeData, "PLACE_COUNTY")
if cousub_ousted_weight > 0
    push!(nodeData, "COUSUB")
end
union!(nodeData, ["PROTECTED_CITY", "MEMBER_NAMES"])
ad_param["sampling_graph"] = pctGraphPath
ad_param["assignment_graph"] = assignment_graph
ad_param["protected_cities"] = protected_cities
base_graph = BaseGraph(pctGraphPath, "POP20", inc_node_data=nodeData, 
                       area_col="area", node_border_col="border_length", 
                       edge_perimeter_col="length");
if place_county_units
    # Validate original per-PLACE populations BEFORE filtering incorporated GEOIDs.
    # A county-specific representation must sum back to the preserved PLACE field
    # on every graph node, including all five hard-contracted city supernodes.
    # The legacy validator requires 1,341 PLACE_COUNTY units and therefore
    # cannot validate a deliberately filtered Columbus-only graph.
    validate_columbus_only_graph!(base_graph)
end
if !isempty(place_geoids_file)
    # The generated graph is ALREADY Columbus-only in both PLACE and
    # PLACE_COUNTY; the GEOID file was validated at startup. Do not mutate
    # node_attributes here or rely on a version-specific restrict helper.
    ad_param["weighted_place_geoids"] = sort!(collect(selected_geoids))
end
write_map = build_precinct_output(base_graph, protected_cities)
weight_func(g, n1, n2) = edge_weight_func(g, n1, n2, cnty_scale=cnty_w_scale,
                                          place_scale=place_w_scale,
                                          cnty_field="COUNTY",
                                          place_field=place_field)
modify_edge_weights!(base_graph, weight_func);
graph = MultiLevelGraph(base_graph, ["NAME"]);

rng = PCG.PCGStateOneseq(UInt64, base_seed + mul_seed*thread_id)
constraints = initialize_constraints()
add_constraint!(constraints, PopulationConstraint(graph, num_dists, popdev))
partition = LinkCutPartition(
    graph,
    constraints,
    num_dists;
    rng=rng,
    verbose=true,
)

cycle_walk = build_lifted_tree_cycle_walk(constraints)
internal_walk = build_internal_forest_walk(constraints)
proposal = [(cycle_walk_frac, cycle_walk), 
            (1.0-cycle_walk_frac, internal_walk)]
measure = Measure()
if gamma != 0
    push_energy!(measure, get_log_spanning_forests, gamma)
end
if iso_weight > 0
    push_energy!(measure, get_isoperimetric_score, iso_weight)
end
if chamber == "upper"
    push_energy!(measure, get_ohio_senate_county_violations, const_weight)
end
if cousub_ousted_weight > 0
    cousub_ousted_score = build_cousub_ousted_score(partition)
    push_energy!(measure, cousub_ousted_score, cousub_ousted_weight; desc="cousub_ousted_pp")
end
# Separate the three Columbus-in-county units when computing ousted pop.
# IMPORTANT: this is NOT the statewide largest piece of Columbus.
# U = 100 / N_OH * sum_c [N_Columbus,c - max_d N_Columbus,c,d].
# We define this helper in the dedicated runner instead of depending on the
# Grid-local aux_place.jl version supporting a place_field keyword.
function build_columbus_county_ousted_score(partition::LinkCutPartition)
    expected = Dict("3918000|FRANKLIN" => 880830,
                    "3918000|DELAWARE" => 14570,
                    "3918000|FAIRFIELD" => 10849)
    unit_ids = Dict{String,Int}()
    unit_totals = Dict{String,Int}()
    pieces = Tuple{Int,Int,Int}[]
    for (node_index, attrs) in enumerate(partition.graph.node_attributes)
        for (raw_key, raw_population) in attrs["PLACE_COUNTY"]
            key = String(raw_key)
            haskey(expected, key) ||
                throw(ArgumentError("Unexpected Columbus county component: $key"))
            pop = Int(raw_population)
            pop > 0 || throw(ArgumentError("Nonpositive Columbus county population"))
            unit_index = get!(unit_ids, key, length(unit_ids) + 1)
            push!(pieces, (node_index, unit_index, pop))
            unit_totals[key] = get(unit_totals, key, 0) + pop
        end
    end
    unit_totals == expected ||
        throw(ArgumentError("Incorrect Columbus county component populations: $unit_totals"))
    total_columbus = sum(values(unit_totals))
    total_columbus == 906249 || throw(ArgumentError("Wrong Columbus total"))
    total_ohio = Int(partition.graph.total_pop)
    total_ohio == 11799448 || throw(ArgumentError("Wrong statewide population"))
    counts = zeros(Int, partition.num_dists, length(unit_ids))

    function score(p, districts=nothing; update=nothing)
        assignment = update === nothing ? p.node_to_dist : p.node_to_dist_update
        fill!(counts, 0)
        @inbounds for (node_index, unit_index, pop) in pieces
            counts[assignment[node_index], unit_index] += pop
        end
        ousted = total_columbus - sum(maximum, eachcol(counts))
        return 100.0 * ousted / total_ohio
    end
    return score
end

place_ousted_score = if place_ousted_weight > 0 || !isempty(place_geoids_file)
    ad_param["place_ousted_population_basis"] = "statewide Ohio POP20 percentage points"
    ad_param["place_ousted_unit_definition"] = "Columbus-in-county PLACE_COUNTY units; each has its own majority district"
    build_columbus_county_ousted_score(partition)
end

if place_ousted_weight > 0
    push_energy!(measure, place_ousted_score, place_ousted_weight; desc="place_ousted_pp")
end

writer = Writer(measure, constraints, partition, outPath; 
                additional_parameters=ad_param)
push_writer!(writer, get_isoperimetric_scores)
if cousub_ousted_weight > 0
    push_writer!(writer, cousub_ousted_score; desc="cousub_ousted_pp")
end
if place_ousted_score !== nothing
    push_writer!(writer, place_ousted_score; desc="place_ousted_pp")
end

@show "running mcmc"
try
    if target_maps > 0
        # Same MCMC kernel as before: both screens act only on output checkpoints.
        # For B7, reconstruct county incidence from the *original precinct graph*,
        # since a protected-city sampling node can span multiple counties.
        county_b7_screen = enforce_b7 ? build_original_county_output_screen(
            base_graph, assignment_graph;
            num_dists=num_dists,
            max_districts_per_county=max_districts_per_county,
        ) : nothing
        county_output_filter = enforce_b7 ? nothing : fixed_region_district_cap(
            graph, "COUNTY", max_districts_per_county; num_dists=num_dists,
        )
        @time stats = run_until_compliant_quota!(
            partition, proposal, measure, rng, county_output_filter, writer;
            output_freq=output_freq,
            target_maps=target_maps,
            write_map=write_map,
            screen_check=county_b7_screen,
        )
        @show stats
    else
        @time run_metropolis_hastings!(partition, proposal, measure, steps, rng,
                                       writer=writer, output_freq=output_freq)
    end
finally
    close_writer(writer)
end