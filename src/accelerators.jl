
using Statistics


"""
    PrimitiveInfo

A type holding information about a primitive (its index in an array, its centre point and
its bounding box). These are used when constructing the BVH accelerator tree structure.
    
"""
struct PrimitiveInfo
    number::Int64
    centroid::Point3
    BBox::BoundingBox
end


"""
    LinearBVHNode

A node of the flattened version of the BVH three. These nodes are used in the actual
ray-tracer acceleration.

# Fields

- `BBox::BoundingBox`: Bounding box of the primitive(s) contained in this node.
- `offset::Int64`: If this is a leaf node, contains the index of the corresponding primitive
  in the `primitives` array of the `BVHAccelerator`.
- `leaf::Bool`: Whether this is a leaf or branch node.
- `axis::Int64`: The axis along which the primitives are split if this is a branch node.
  Knowing this relative to the ray direction is relevant in the traversal of the tree.

"""
struct LinearBVHNode
    BBox::BoundingBox
    offset::Int64
    leaf::Bool
    axis::Int64
end




"""
    BVHBuildNode

A mutable type representing a node in the Bounding Volume Hierarchy tree. This is used when
creating the tree, which is then flattened into a more efficient structure using the
immutable `LinearBVHNode` types.

"""
mutable struct BVHBuildNode
    leaf::Bool
    split_axis::Int64
    first_offset::Int64
    BBox::BoundingBox
    childA::BVHBuildNode
    childB::BVHBuildNode
    BVHBuildNode() = new()
end


"""
    init_leaf!(node::BVHBuildNode, idx::Integer, BB::BoundingBox)

Initialize the `BVHBuildNode` as a leaf node of the BVH tree.

"""
function init_leaf!(node::BVHBuildNode, idx::Integer, BB::BoundingBox)
    node.leaf = true
    node.first_offset = idx
    node.BBox = BB
    node.split_axis = -1
    return node
end


"""
    init_interior!(node::BVHBuildNode, idx::Integer, BB::BoundingBox)

Initialize the `BVHBuildNode` as an interior node of the BVH tree.

"""

function init_interior!(node::BVHBuildNode, axis::Integer, childA::BVHBuildNode, childB::BVHBuildNode)
    node.leaf = false
    node.split_axis = axis
    node.BBox = BoundingBox(childA.BBox, childB.BBox)
    node.childA = childA
    node.childB = childB
    node.first_offset = -1
    return node
end



"""
    BVH_build_recursive(build_data, startN, endN, primitives, ordered_primitives,
                        total_nodes)

Recursively construct the BVH tree using `build_data` and `primitives`. This will populate
the array `ordered_primitives` with the primitives in the same order as to the tree leaf
nodes which refer to them.

"""
function BVH_build_recursive(build_data::Array{PrimitiveInfo,1}, startN::Integer, 
                             endN::Integer, primitives::Array{T,1}, ordered_primitives::Array{T,1}, 
                             total_nodes::Integer) where T <: Primitive
                            
    total_nodes += 1
    node = BVHBuildNode()
    
    if startN < 1 || endN < 1
        error()
    end
    
    if startN == endN
        # Only one primitive in the range, create leaf node
        n = build_data[startN].number
        push!(ordered_primitives, primitives[n])
        offset = length(ordered_primitives)
        init_leaf!(node, offset, build_data[startN].BBox)
    else
        # More than one primitive in the range, make partition and call recursively
        CBox = BoundingBox([build_data[i].centroid for i = startN:endN])
        dim = max_extent(CBox)
        pmid = mean(CBox.pMin, CBox.pMax)
        cmp(x,y) = x.centroid[dim] < y.centroid[dim]
        build_data[startN:endN] = sort(build_data[startN:endN] ; lt=cmp)
        mid = ceil(Int64, (startN+endN) / 2)
        
        # TODO: FIX THIS BETTER! THIS HAS SOME SUBTLE BUG
        #for i = startN:endN
        #    if (pmid[dim] <= build_data[i].centroid[dim])
        #        mid = i 
        #        break
        #    end
        #end
        
        childA, total_nodes = BVH_build_recursive(build_data, startN, mid-1, primitives, 
                                                  ordered_primitives, total_nodes)
        childB, total_nodes = BVH_build_recursive(build_data, mid, endN, primitives, 
                                                  ordered_primitives, total_nodes)
        init_interior!(node, dim, childA, childB)
    end
    return node, total_nodes
end



"""
    flatten_BVH!(node::BVHBuildNode, linear::Array{LinearBVHNode,1}, idx::Integer)

Recursively traverse the BVH tree starting from `node`, and produce the flattened version
of the tree in the array `linear`.

"""
function flatten_BVH!(node::BVHBuildNode, linear::Array{LinearBVHNode,1}, idx::Integer)
    idx0 = idx
    if node.leaf
        new_node = LinearBVHNode(node.BBox, node.first_offset, true, -1)
    else
        idx1 = flatten_BVH!(node.childA, linear, idx+1)
        idx2 = flatten_BVH!(node.childB, linear, idx1+1)
        new_node = LinearBVHNode(node.BBox, idx1+1, false, node.split_axis)
        idx = idx2
    end
    linear[idx0] = new_node
    return idx
end

# LEFT: your number is my number, plus 1. Also, give me the biggest number in your subtree.
# RIGHT: your number is the biggest number in the LEFT subtree, plus 1
# I got to give upwards the biggest number in my subtree. But that's just the highest number
# in my right subtree.





"""
    BVHAccelerator{P<:Primitive} <: Aggregate

The Bounding Volume Hierarchy accelerator keeps the list of primitives, and a list of
`LinearBVHNode` objects, which contain all the information about the tree structure.

"""
struct BVHAccelerator{P<:Primitive} <: Aggregate
    primitives::Array{P,1}
    nodes::Array{LinearBVHNode,1}
end


"""
    world_bounds(B::BVHAccelerator)

Return the bounding box containing all the primitives in the Accelerator.

"""
world_bounds(B::BVHAccelerator) = B.nodes[1].BBox



"""
    BVHAccelerator(prims::Array{T,1}) where T<:Primitive

Construct a `BVHAccelerator` from the list of primitives.

"""
function BVHAccelerator(prims::Array{T,1}) where T<:Primitive
    # Refine primitives to their fully intersectable components
    primtype = eltype(prims)
    primitives = primtype[]
    for p in prims
        fully_refine!(p, primitives)
    end
    
    # Extract relevant info from primitives and push into an array
    build_data = PrimitiveInfo[]
    for (i, obj) in enumerate(primitives)
        BB = world_bounds(obj)
        c = mean(BB.pMin, BB.pMax)
        info = PrimitiveInfo(i, c, BB)
        push!(build_data, info)
    end
    
    # Generate binary tree structure of bounding boxes
    ordered = primtype[]
    root_node, total_nodes = BVH_build_recursive(build_data, 1, length(primitives), primitives, 
                                                 ordered, 0)

    # Transform the tree structure into a more efficient form
    linear = Array{LinearBVHNode}(undef,total_nodes)
    flatten_BVH!(root_node, linear, 1)
    BVHAccelerator(ordered, linear)
end


# Two arrays that are reused by the `intersect` and `intersectP` methods to
# avoid reallocating them every time those functions are called (which is very
# often). A significant optimization.
const TODO_ARRAY = zeros(Int64, 64)
const dir_is_neg = zeros(Bool, 3)


"""
    update_isect(isect::Intersection, best_isect::Intersection)

Compare `isect` and `best_isect` and return the one with a smaller `tmin`.

"""
function update_isect(isect::Intersection, best_isect::Intersection)
    if isect.tmin < best_isect.tmin
        return isect
    else
        return best_isect
    end
end

"""
    update_isect(isect::Nothing, best_isect::Intersection)

Always return `best_isect` when `isect` is `nothing`.
"""
update_isect(isect::Nothing, best_isect::Intersection) = best_isect


"""
    update_isect(isect, best_isect::Nothing)

Always return `isect` when `best_isect` is `nothing`.
"""
update_isect(isect, best_isect::Nothing) = isect


"""
    intersect(ray::Ray, BVH::BVHAccelerator)

Find the nearest intersection of the ray and the primitives in the `BVHAccelerator`. If you
only need to know whether a ray intersects something or not, use [`intersectP`](@ref)
instead, since `intersect` must find _all possible_ intersections to determine the nearest
one and is therefore much slower.

"""
function intersect(ray::Ray, BVH::BVHAccelerator)
    if length(BVH.nodes) == 0
        return nothing
    end
    nodeN = 1
    todo_offset = 0
    todo = TODO_ARRAY
    for i = 1:64
        todo[i] = 0
    end
    for i = 1:3
        dir_is_neg[i] = ray.direction[i] < 0.0
    end
    
    tmin = Inf
    best_isect = nothing
    
    while true
        # Check the node given by the index nodeN
        node = BVH.nodes[nodeN]
        if intersectP(ray, node.BBox)
            # We hit the bounding box of this node.
            if node.leaf
                # This is a leaf node.
                # Check intersection with the primitive in the node and keep best intersection.
                isect = intersect(ray, BVH.primitives[node.offset])
                best_isect = update_isect(isect, best_isect)

                # If there's nothing in the todo queue, we exit the loop.
                if todo_offset == 0
                    break
                end
                # Othewise, set nodeN based on current todo value and shift todo queue left.
                nodeN = todo[todo_offset]
                todo_offset -= 1
            else
                # This is not a leaf node. 
                # We point nodeN to either the left or right child node, depending on the ray 
                # direction. We shift right in the todo queue and put the other node there.
                todo_offset += 1
                if dir_is_neg[node.axis]
                    # Put left child in todo queue and set nodeN to right child.
                    todo[todo_offset] = nodeN + 1
                    nodeN = node.offset
                else
                    # Put right child in todo queue and set nodeN to left child.
                    todo[todo_offset] = node.offset
                    nodeN = nodeN + 1
                end
            end
        else
            # We did not hit a bounding box. 
            # If there's nothing in the todo queue, we exit.
            if todo_offset == 0
                break
            end
            # Otherwise, set nodeN based on current todo value and shift left in the todo queue.
            nodeN = todo[todo_offset]
            todo_offset -= 1
            
        end
    end
    return best_isect
end


"""
    intersectP(ray::Ray, BVH::BVHAccelerator)

Returns `true` if the given ray intersects any primitive in the `BVHAccelerator`, else
returns false. This is much faster to do than [`intersect`](@ref), since it's enough to
find any intersection, not specifically the nearest one.

Possibly the most often called function in the entire ray-tracer.

"""
function intersectP(ray::Ray, BVH::BVHAccelerator)
    if length(BVH.nodes) == 0
        return false
    end
    nodeN = 1
    todo_offset = 0
    todo = TODO_ARRAY
    for i = 1:64
        todo[i] = 0
    end
    for i = 1:3
        dir_is_neg[i] = ray.direction[i] < 0.0
    end
    while true
        node = BVH.nodes[nodeN]
        if intersectP(ray, node.BBox)
            if node.leaf
                if intersectP(ray, BVH.primitives[node.offset])
                    return true
                end
                if todo_offset == 0
                    break
                end
                nodeN = todo[todo_offset]
                todo_offset -= 1
            else
                todo_offset += 1
                if dir_is_neg[node.axis]
                    todo[todo_offset] = nodeN + 1
                    nodeN = node.offset
                else
                    todo[todo_offset] = node.offset
                    nodeN = nodeN + 1
                end
            end
        else
            if todo_offset == 0
                break
            end
            nodeN = todo[todo_offset]
            todo_offset -= 1
            
        end
    end
    return false
end
