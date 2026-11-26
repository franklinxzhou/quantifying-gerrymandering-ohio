"""
Build independent MCT township and municipality edge multipliers.

Three modes, with separate w_T and w_MC:
  linear:      f_G(alpha) = 1 + (w_G-1)*alpha
  exponential: f_G(alpha) = exp(log(w_G)*alpha)
  full:        f_G = w_G if both nodes geometrically intersect one qualifying
               geography-G unit, else 1 (requires area-overlap incidence).

The score alpha for a geography type sums, over shared units, the products
P(unit ∩ node_i)/P(node_i) * P(unit ∩ node_j)/P(node_j).

Townships and geometrically non-nested incorporated municipalities have
separate channels. The joint edge multiplier is max(f_T, f_M), never larger than max(w_T,w_MC).
This avoids double weighting township/municipal geographic overlaps; it is an
intentional change from the additive convention of the previous MCT pilot.
The existing county multiplier is applied outside this function.

IMPORTANT: the geometric cross-place scope may include 30 municipalities
promoted to 'nested' in the POPULATION energy crosswalk. This is intentional:
municipal territory outside its parent township should not be silently ignored
by the edge proposal, even if its assigned block population is zero there.
"""
function build_mct_edge_multipliers(
    node_attributes,
    atoms_file::AbstractString;
    township_w::Real=1.0,
    cross_place_w::Real=1.0,
    edge_model::AbstractString="linear",
    geometric_cross_ids_file::Union{Nothing,AbstractString}=nothing,
    geometric_incidence_file::Union{Nothing,AbstractString}=nothing,
)
    for (label, w) in (("township_w", township_w), ("cross_place_w", cross_place_w))
        isfinite(w) && w >= 1.0 ||
            throw(ArgumentError("$label must be finite and >= 1"))
    end
    edge_model in ("linear", "exponential", "full") ||
        throw(ArgumentError("edge_model must be linear, exponential, or full"))
    if edge_model == "full" && geometric_incidence_file === nothing
        throw(ArgumentError("full mode requires a geometric incidence file; population atoms cannot detect zero-population overlaps"))
    end
    isfile(atoms_file) || throw(ArgumentError("MCT atom file is missing: $atoms_file"))

    geometric_cross = nothing
    if geometric_cross_ids_file !== nothing
        isfile(geometric_cross_ids_file) || throw(ArgumentError("Missing geometric cross-place IDs"))
        ids = Set{String}()
        for line in eachline(geometric_cross_ids_file)
            p = strip(line)
            isempty(p) && continue
            occursin(r"^39\d{5}$", p) || throw(ArgumentError("Invalid geometric place GEOID: $p"))
            p in ids && throw(ArgumentError("Duplicate geometric place GEOID: $p"))
            push!(ids, p)
        end
        isempty(ids) && throw(ArgumentError("Empty geometric cross-place ID file"))
        geometric_cross = ids
    end

    nn = length(node_attributes)
    node_index = Dict{String,Int}()
    pop = Vector{Int}(undef, nn)
    for (i, attrs) in enumerate(node_attributes)
        name = String(attrs["NAME"])
        haskey(node_index, name) && throw(ArgumentError("Duplicate graph NAME: $name"))
        node_index[name] = i
        pop[i] = Int(attrs["POP20"])
        pop[i] >= 0 || throw(ArgumentError("Negative population at $name"))
    end

    town = [Dict{String,Int}() for _ in 1:nn]
    cross = [Dict{String,Int}() for _ in 1:nn]
    observed = zeros(Int, nn)
    keys_seen = Set{Tuple{Int,String,String}}()
    place_classes = Dict{String,String}()
    nested_parents = Dict{String,String}()
    towns = Set{String}()
    cross_ids = Set{String}()
    nested_ids = Set{String}()
    n_atoms = 0

    open(atoms_file, "r") do io
        eof(io) && throw(ArgumentError("Empty MCT atom file"))
        expected = ["node_name", "township_geoid", "place_geoid", "place_class", "population"]
        split(chomp(readline(io)), '\t'; keepempty=true) == expected ||
            throw(ArgumentError("Unexpected MCT atom header"))
        for line in eachline(io)
            cols = split(chomp(line), '\t'; keepempty=true)
            length(cols) == 5 || throw(ArgumentError("MCT atom row $(n_atoms+2): expected five fields"))
            name, t, p, cls, pop_text = cols
            haskey(node_index, name) || throw(ArgumentError("Unknown node in atoms: $name"))
            i = node_index[name]
            cls in ("none", "cross", "nested") || throw(ArgumentError("Invalid place_class: $cls"))
            (isempty(p) == (cls == "none")) || throw(ArgumentError("Invalid place identifier/class combination"))
            isempty(t) || occursin(r"^39\d{8}$", t) || throw(ArgumentError("Invalid COUSUB GEOID: $t"))
            isempty(p) || occursin(r"^39\d{5}$", p) || throw(ArgumentError("Invalid PLACE GEOID: $p"))
            v = parse(Int, pop_text)
            v > 0 || throw(ArgumentError("Atom population must be positive"))
            key = (i,t,p)
            key in keys_seen && throw(ArgumentError("Repeated graph-node/unit intersection: $key"))
            push!(keys_seen, key)
            observed[i] += v
            if !isempty(t)
                town[i][t] = get(town[i], t, 0) + v
                push!(towns, t)
            end
            if !isempty(p)
                if haskey(place_classes, p) && place_classes[p] != cls
                    throw(ArgumentError("Inconsistent classification of PLACE $p"))
                end
                place_classes[p] = cls
                # Geometric scope deliberately reinstates a population-promoted
                # municipality when it is not literally contained in a township.
                eligible = geometric_cross === nothing ? cls == "cross" : p in geometric_cross
                if eligible
                    cross[i][p] = get(cross[i], p, 0) + v
                    push!(cross_ids, p)
                end
                if cls == "nested"
                    isempty(t) && throw(ArgumentError("Nested PLACE $p has population outside townships"))
                    if haskey(nested_parents, p) && nested_parents[p] != t
                        throw(ArgumentError("Nested PLACE $p has population in multiple townships"))
                    end
                    nested_parents[p] = t
                    push!(nested_ids, p)
                end
            end
            n_atoms += 1
        end
    end
    n_atoms > 0 || throw(ArgumentError("MCT atom file contains no population"))
    for i in 1:nn
        observed[i] == pop[i] || throw(ArgumentError(
            "MCT node population mismatch at $(node_attributes[i]["NAME"]): atoms=$(observed[i]), graph=$(pop[i])"))
    end
    if geometric_cross !== nothing
        # Every population-non-nested PLACE must also be geographically non-nested.
        for (p, cls) in place_classes
            if cls == "cross" && !(p in geometric_cross)
                throw(ArgumentError("Population-cross place $p missing from geometric cross scope"))
            end
        end
    end

    town_frac = Vector{Dict{String,Float64}}(undef, nn)
    cross_frac = Vector{Dict{String,Float64}}(undef, nn)
    for i in 1:nn
        denom = pop[i]
        town_frac[i] = denom == 0 ? Dict{String,Float64}() :
            Dict(k => v/denom for (k,v) in town[i])
        cross_frac[i] = denom == 0 ? Dict{String,Float64}() :
            Dict(k => v/denom for (k,v) in cross[i])
        sum(values(town_frac[i])) <= 1.0 + 1e-12 || throw(ArgumentError("Township fractions exceed 1 at node $i"))
        sum(values(cross_frac[i])) <= 1.0 + 1e-12 || throw(ArgumentError("Cross-place fractions exceed 1 at node $i"))
    end

    geo_town = [Set{String}() for _ in 1:nn]
    geo_cross = [Set{String}() for _ in 1:nn]
    incidence_count = 0
    if geometric_incidence_file !== nothing
        isfile(geometric_incidence_file) || throw(ArgumentError("Missing geometric incidence file"))
        open(geometric_incidence_file, "r") do io
            eof(io) && throw(ArgumentError("Empty geometric incidence file"))
            split(chomp(readline(io)), '\t') == ["node_name", "geography", "geoid"] ||
                throw(ArgumentError("Unexpected geometric incidence TSV header"))
            for line in eachline(io)
                fields = split(chomp(line), '\t'; keepempty=true)
                length(fields) == 3 || throw(ArgumentError("Invalid incidence row"))
                name, geography, geoid = fields
                i = get(node_index, name, 0)
                i != 0 || throw(ArgumentError("Unknown incidence graph node $name"))
                if geography == "township"
                    occursin(r"^39\d{8}$", geoid) || throw(ArgumentError("Invalid township incidence $geoid"))
                    geoid in geo_town[i] && throw(ArgumentError("Duplicate township incidence $name, $geoid"))
                    push!(geo_town[i], geoid)
                elseif geography == "place"
                    occursin(r"^39\d{5}$", geoid) || throw(ArgumentError("Invalid place incidence $geoid"))
                    geometric_cross === nothing || geoid in geometric_cross ||
                        throw(ArgumentError("Incidence for non-cross place $geoid"))
                    geoid in geo_cross[i] && throw(ArgumentError("Duplicate place incidence $name, $geoid"))
                    push!(geo_cross[i], geoid)
                else
                    throw(ArgumentError("Invalid incidence geography $geography"))
                end
                incidence_count += 1
            end
        end
        incidence_count > 0 || throw(ArgumentError("No incidence rows"))
        # A positive-population atom must not disappear from the polygon-derived
        # incidence. If it does, geometry/source vintages are incompatible.
        for i in 1:nn
            for t in keys(town[i])
                t in geo_town[i] || throw(ArgumentError("Missing geometric township overlap at node $(node_attributes[i]["NAME"]), $t"))
            end
            for p in keys(cross[i])
                p in geo_cross[i] || throw(ArgumentError("Missing geometric municipal overlap at node $(node_attributes[i]["NAME"]), $p"))
            end
        end
    end

    function shared_fraction(a::Dict{String,Float64}, b::Dict{String,Float64})
        if length(a) > length(b)
            a, b = b, a
        end
        score = 0.0
        for (id, x) in a
            score += x * get(b, id, 0.0)
        end
        -1e-12 <= score <= 1.0+1e-12 || throw(ArgumentError("Invalid MCT shared-population score $score"))
        return clamp(score, 0.0, 1.0)
    end
    function shared_incidence(a::Set{String}, b::Set{String})
        if length(a) > length(b)
            a, b = b, a
        end
        for key in a
            key in b && return true
        end
        return false
    end
    wt = Float64(township_w)
    wm = Float64(cross_place_w)
    function channel_gain(alpha::Float64, weight::Float64)
        if edge_model == "linear"
            return (weight-1)*alpha
        elseif edge_model == "exponential"
            return exp(log(weight)*alpha)-1
        end
        throw(ArgumentError("full mode uses geometric incidence"))
    end
    function multiplier(i::Int, j::Int)
        1 <= i <= nn && 1 <= j <= nn || throw(BoundsError(node_attributes, (i,j)))
        value = if edge_model == "full"
            t = shared_incidence(geo_town[i], geo_town[j]) ? wt : 1.0
            m = shared_incidence(geo_cross[i], geo_cross[j]) ? wm : 1.0
            max(t, m)
        else
            t = shared_fraction(town_frac[i], town_frac[j])
            m = shared_fraction(cross_frac[i], cross_frac[j])
            1.0 + max(channel_gain(t, wt), channel_gain(m, wm))
        end
        isfinite(value) && 1.0 <= value <= max(wt, wm) + 1e-10 ||
            throw(ArgumentError("Invalid MCT edge multiplier"))
        return value
    end

    return (
        multiplier=multiplier,
        n_atoms=n_atoms,
        n_townships=length(towns),
        n_cross_places=length(cross_ids),
        n_nested_places=length(nested_ids),
        n_geometric_cross_places=(geometric_cross === nothing ? length(cross_ids) : length(geometric_cross)),
        n_geometric_incidence=incidence_count,
        accounted_population=sum(pop),
        nesting_model=geometric_cross === nothing ? "population_effective" : "geometric_edge_population_energy",
        weighting_model="max_$(edge_model)_shared_membership",
    )
end
