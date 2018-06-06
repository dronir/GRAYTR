


const pattern_def = r"^(material|object|include|light_source|camera)\s+(\w+):"

const pattern_transform = r"^\s+(rotate|scale|translate)\s+(([xyzXYZ])\s+(-?[\d]+(\.[\d]*)?)|\[(\s*-?[\d]+(\.[\d]+)?,\s*-?[\d]+(\.[\d]+)?,\s*-?[\d]+(\.[\d]+)?\s*)\])"

const pattern_rotation = r"(rotate)\s+([xyzXYZ])\s+(-?[\d]+(\.[\d]+)?)"
const pattern_scale = r"(scale)\s+([xyzXYZ])\s+(-?[\d]+(\.[\d]+)?)"
const pattern_translate = r"(translate)\s+([xyzXYZ])\s+(-?[\d]+(\.[\d]+)?)"
const pattern_transform_vec = r"(translate|scale)\s+\[(\s*-?[\d]+(\.[\d]+)?,\s*-?[\d]+(\.[\d]+)?,\s*-?[\d]+(\.[\d]+)?\s*)\]"


vec_transformations = Dict(
    "translate" => translation,
    "scale" => scaling,
    "rotate" => rotation
)

VEC_AXES = Dict(
    "x" => Vector3(1.0, 0.0, 0.0),
    "y" => Vector3(0.0, 1.0, 0.0),
    "z" => Vector3(0.0, 0.0, 1.0)
)

str_to_vec(numbers) = Vector3(float.(split(numbers, ","))...)

function parse_rotation(line)
    M = match(pattern_rotation, line)
    ax = VEC_AXES[lowercase(M[2])]
    amount = float(M[3])
    return rotation(ax, amount)
end

function parse_scaling(line)
    M = match(pattern_scale, line)
    ax = VEC_AXES[lowercase(M[2])]
    amount = float(M[3])
    one = Vector3(1, 1, 1)
    return scaling(one - ax + amount * ax)
end

function parse_translation(line)
    M = match(pattern_translate, line)
    ax = VEC_AXES[lowercase(M[2])]
    amount = float(M[3])
    return translation(amount * ax)
end

function parse_vec_transform(line)
    M = match(pattern_transform_vec, line)
    f = vec_transformations[M[1]]
    v = str_to_vec(M[2])
    return f(v)
end

function get_transformation(line)
    if ismatch(pattern_rotation, line)
        return parse_rotation(line)
    elseif ismatch(pattern_scale, line)
        return parse_scaling(line)
    elseif ismatch(pattern_translate, line)
        return parse_translation(line)
    elseif ismatch(pattern_transform_vec, line)
        return parse_vec_transform(line)
    else
        error("Expected transformation, got: '$line'")
    end
end



"""
    skip(line)

Returns true/false whether line is skippable: only whitespace or starts with a # character.
"""
skip(line::AbstractString) = startswith(strip(line), "#") || length(strip(line)) == 0


"""
    noncomment(lines, state)

Advances in `lines` until the first non-skippable line. Returns that line and the iterator state.
"""
function noncomment(lines, state)
    line = ""
    while !done(lines, state) && skip(line)
        line, state = next(lines, state)
    end
    return line, state
end



"""
    parse_transformations(lines, state)

Gathers transformations until the block ends (i.e. a new block begins). Returns vector of strings.
"""
function parse_transformations(lines, state)
    T = String[]
    while true
        line, state = noncomment(lines, state)
        if ismatch(pattern_def, line) || done(lines, state)
            return T, state-1
        end
        if ismatch(r"^\s+(rotate|scale|translate)\s+(([xyzXYZ])\s+(-?[\d]+(\.[\d]*)?)|\[(\s*-?[\d]+(\.[\d]+)?,\s*-?[\d]+(\.[\d]+)?,\s*-?[\d]+(\.[\d]+)?\s*)\])", line)
            push!(T, String(strip(line)))
        else
            error("Expected a transformation, got: $line")
        end
    end
end

"""
    get_or_error(regex, pos, lines, state, expect, wrap=x->x)

Returns the group `pos` from math of `regex` from next non-comment non-empty line in lines.
If no match, throws an error.
The return value is mapped by function `wrap`.
"""
function get_or_error(regex, pos, lines, state, expect, wrap=x->x)
    line, state = noncomment(lines, state)
    m = match(regex, line)
    if m == nothing
        error("Expected $expect, got: $line")
    else
        result = wrap(m[pos])
    end
    return result, state
end

"""
    get_or_default(regex, pos, lines, state, default, wrap=x->x)
    
Returns the group `pos` from math of `regex` from next non-comment non-empty line in lines.
If no match, returns `default`.
The return value is mapped by function `wrap`.
"""
function get_or_default(regex, pos, lines, state, default, wrap=x->x)
    line, state = noncomment(lines, state)
    m = match(regex, line)
    if m == nothing
        return default, state-1
    else
        result = wrap(m[pos])
    end
    return result, state
end


"""
    parse_sphere(lines, state)

Read lines from an array of strings and parse a sphere definition.
"""
function parse_sphere(lines, state)
    radius, state = get_or_error(r"^\s+radius\s+(\d+(\.?\d*))", 1, lines, state, "sphere radius", float)
    material, state = get_or_default(r"^\s+material\s+(\w+)", 1, lines, state, "")
    T, state = parse_transformations(lines, state)
    obj = Dict{String,Any}(
        "type" => "sphere",
        "radius" => radius,
        "material" => material,
        "transformations" => T
    )
    return obj, state
end


"""
    parse_cylinder(lines, state)

Read lines from an array of strings and parse a cylinder definition.
"""
function parse_cylinder(lines, state)
    radius, state = get_or_error(r"^\s+radius\s+(\d+(\.?\d*))", 1, lines, state, "cylinder radius", float)
    height, state = get_or_error(r"^\s+height\s+(\d+(\.?\d*))", 1, lines, state, "cylinder height", float)
    material, state = get_or_default(r"^\s+material\s+(\w+)", 1, lines, state, "")
    T, state = parse_transformations(lines, state)
    obj = Dict{String,Any}(
        "type" => "cylinder",
        "radius" => radius,
        "height" => height,
        "material" => material,
        "transformations" => T
    )
    return obj, state
end

"""
    parse_object(lines, state)

Read lines from an array of strings and parse an object definition, calling sub-functions.
"""
function parse_object(lines, state)
    line, state = noncomment(lines, state)
    line = strip(line)
    if line == "sphere"
        obj, state = parse_sphere(lines, state)
    elseif line == "cylinder"
        obj, state = parse_cylinder(lines, state)
    else
        error("Not valid shape: $line")
    end
    return obj, state
end


function parse_material(lines, state)
    BRDF, state = get_or_error(r"^\s+BRDF\s+(.+)", 1, lines, state, "BRDF name")
    spec, state = get_or_error(r"^\s+spectrum\s+(from\s+(\w+\.?\w+))", 1, lines, state, "spectrum")
    mat = Dict{String,Any}(
        "BRDF" => BRDF,
        "spectrum" => spec
    )
    return mat, state
end


function parse_primitive(lines, state)
    material, state = get_or_default(r"^\s+material\s+(\w+)", 1, lines, state, "")
    T, state = parse_transformations(lines, state)
    prim = Dict{String,Any}(
        "material" => material,
        "transformations" => T
    )
    return prim, state
end

function parse_light(lines, state)
    light_type, state = get_or_error(r"^\s+(point_light|distant_light)", 1, lines, state, "light source type")
    spec, state = get_or_error(r"^\s+spectrum\s+(from\s+(\w+\.?\w+))", 1, lines, state, "spectrum")
    T, state = parse_transformations(lines, state)
    light = Dict{String,Any}(
        "type" => light_type,
        "spectrum" => spec,
        "transformations" => T
    )
    return light, state
end





parse(input::String) = parse(split(input, "\n"))

function parse(lines::Vector{T}) where T<:AbstractString
    contents = Dict{String,Any}()
    contents["materials"] = Dict{String,Any}()
    contents["objects"] = Dict{String,Any}()
    contents["light_sources"] = Dict{String,Any}()
    contents["primitives"] = []
    
    state = start(lines)
    
    while !done(lines, state)
        line, state = noncomment(lines, state)
        # Line starts a new definition block
        if ismatch(pattern_def, line)
            m = match(pattern_def, line)
            deftype = m[1]
            if deftype == "object"
               result, state = parse_object(lines, state)
               name = m[2]
               contents["objects"][name] = result
            elseif deftype == "material"
                result, state = parse_material(lines, state)
                name = m[2]
                contents["materials"][name] = result
            elseif deftype == "include"
                result, state = parse_primitive(lines, state)
                result["object"] = m[2]
                push!(contents["primitives"], result)
            elseif deftype == "light_source"
                name = m[2]
                result, state = parse_light(lines, state)
                contents["light_sources"][name] = result
            end
        end
    end
    return contents
end




function load_spectrum(line)
    if startswith(line, "from")
        fname = match(r"from ([\w\.\_\-]+)", line)[1]
        println(fname)
    end
    return SingleLine(532.0, 1.0)
end

function get_brdf(line, spec)
    M = match(r"((\w+)\s+(.+))", line)
    if M == nothing
        error("Unknown BRDF: $line")
    else
        brdf_type = M[2]
        params = M[3]
        if brdf_type == "Lambert"
            return Lambert(spec)
        elseif brdf_type == "LommelSeeliger"
            return LommelSeeliger(spec)
        elseif brdf_type == "AshkhminShirley"
            return AshkhminShirleySingle(spec, float(params))
        else:
            error("Unknown BRDF: $line")
        end
    end
    
end




function dict_to_objects(input::Dict)
    # Generate materials
    materials = Dict{String, MatteMaterial}()
    for (name, properties) in input["materials"]
        spectrum = load_spectrum(properties["spectrum"])
        BRDF = get_brdf(properties["BRDF"], spectrum)
        mat = MatteMaterial(BRDF)
        materials[name] = mat
    end
    
    # Generate objects
    primitives = GeometricPrimitive[]
    for (i, properties) in enumerate(input["primitives"])
        obj_ref = properties["object"]
        obj_properties = input["objects"][obj_ref]
        
        T = Transformation()
        for line in obj_properties["transformations"]
            T = get_transformation(line) * T
        end
        for line in properties["transformations"]
            T = get_transformation(line) * T
        end
        
        if obj_properties["type"] == "sphere"
            r = obj_properties["radius"]
            obj = Sphere(i, r, T)
        elseif obj_properties["type"] == "cylinder"
            r = obj_properties["radius"]
            h = obj_properties["height"]
            obj = Cylinder(i, r, 0.0, h, false, T, inv(T))
        end
        
        matname = properties["material"]
        obj_mat = obj_properties["material"]
        
        if matname != ""
            mat = materials[matname]
            if obj_mat != ""
                warn("Object $obj_ref material ($obj_mat) overridden in primitive by $matname.")
            end
        elseif obj_mat != ""
            mat = materials[obj_mat]
        else
            error("Material not given for $obj_ref nor in primitive.")
        end
        id = Int64(10^ceil(log(10, length(input["objects"]))) + i)
        push!(primitives, GeometricPrimitive(obj, mat, nothing, id))
    end
    
    accel = BVHAccelerator(primitives)
    
    # Generate lights
    lights = LightSource[]
    
    for (i, (name, properties)) in enumerate(input["light_sources"])
        spectrum = load_spectrum(properties["spectrum"])
        T = Transformation()
        for line in properties["transformations"]
            T = get_transformation(line) * T
        end
        ltype = properties["type"]
        if ltype == "point_light"
            L = PointLight(spectrum, T)
        elseif ltype == "distant_light"
            L = DistantLight(-Z_AXIS, spectrum, T)
        else
            error("Unknown light type: ltype")
        end
    end
    
    if length(lights) == 0
        warn("No light sources found.")
    end
    
    
    # Create camera
    return materials, primitives, lights
end

