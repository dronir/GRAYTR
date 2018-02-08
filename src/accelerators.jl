

struct PrimitiveInfo
    number::Int32
    centroid::Point3
    BBox::BoundingBox
end

struct LinearBVHNode
    BBox::BoundingBox
    offset::Int32
    leaf::Bool
    axis::Int32
end


# BVH build process
#

mutable struct BVHBuildNode
    leaf::Bool
    split_axis::Int32
    first_offset::Int64
    BBox::BoundingBox
    childA::BVHBuildNode
    childB::BVHBuildNode
    BVHBuildNode() = new()
end

function init_leaf!(node::BVHBuildNode, idx::Integer, BB::BoundingBox)
    node.leaf = true
    node.first_offset = idx
    node.BBox = BB
    node.split_axis = -1
    return node
end


function init_interior!(node::BVHBuildNode, axis::Integer, childA::BVHBuildNode, childB::BVHBuildNode)
    node.leaf = false
    node.split_axis = axis
    node.BBox = union(childA.BBox, childB.BBox)
    node.childA = childA
    node.childB = childB
    node.first_offset = -1
    return node
end


function BVH_build_recursive(build_data::Array{PrimitiveInfo,1}, startN::Integer, endN::Integer, 
                             primitives::Array{Primitive,1}, ordered_primitives::Array{Primitive,1}, 
                             total_nodes::Integer)
                            
    total_nodes += 1
    node = BVHBuildNode()
    BB = BoundingBox([build_data[i].BBox for i = startN:endN])
    CBox = BoundingBox([build_data[i].centroid for i = startN:endN])
#    n_primitives = endN - startN
    dim = max_extent(CBox)
    if startN == endN #|| CBox.pMin[dim] == CBox.pMax[dim]
        # Only one primitive in the range, create leaf node
        n = build_data[startN].number
        push!(ordered_primitives, primitives[n])
        offset = length(ordered_primitives)
        init_leaf!(node, offset, BB)
    else
        # More than one primitive in the range, make partition and call recursively
        pmid = mean(CBox.pMin, CBox.pMax)
        cmp(x,y) = x.centroid[dim] < y.centroid[dim]
        build_data[startN:endN] = sort(build_data[startN:endN] ; lt=cmp)
        mid = div(startN+endN, 2)
        
        for i = startN:endN
            if (pmid[dim] < build_data[i].centroid[dim])
                mid = i 
                break
            end
        end
        childA, total_nodes = BVH_build_recursive(build_data, startN, mid-1, primitives, 
                                                  ordered_primitives, total_nodes)
        childB, total_nodes = BVH_build_recursive(build_data, mid, endN, primitives, 
                                                  ordered_primitives, total_nodes)
        init_interior!(node, dim, childA, childB)
    end
    return node, total_nodes
end

# This function takes the BVH binary tree structure and transforms it into a more
# efficient form as a linear array.
function flatten_BVH(node::BVHBuildNode, tree::Array{LinearBVHNode,1}, idx::Integer)
    idx0 = idx
    if node.leaf
        new_node = LinearBVHNode(node.BBox, node.first_offset, true, -1)
    else
        idx1 = flatten_BVH(node.childA, tree, idx+1)
        idx2 = flatten_BVH(node.childB, tree, idx1+1)
        new_node = LinearBVHNode(node.BBox, idx1+1, false, node.split_axis)
        idx = idx2
    end
    tree[idx0] = new_node
    return idx
end

# LEFT: your number is my number, plus 1. Also, give me the biggest number in your subtree.
# RIGHT: your number is the biggest number in the LEFT subtree, plus 1
# I got to give upwards the biggest number in my subtree. But that's just the highest number
# in my right subtree.



# BHVAccelerator constructor from a list of primitives
#

struct BVHAccelerator <: Aggregate
    primitives::Array{Primitive,1}
    nodes::Array{LinearBVHNode,1}
end

world_bounds(B::BVHAccelerator) = B.nodes[1].BBox

function BVHAccelerator(prims::Array{Primitive,1})
    # Refine primitives to their fully intersectable components
    primitives = Primitive[]
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
    ordered = Primitive[]
    root_node, total_nodes = BVH_build_recursive(build_data, 1, length(primitives), primitives, 
                                                 ordered, 0)

    # Transform the tree structure into a more efficient form
    linear = Array{LinearBVHNode}(total_nodes)
    flatten_BVH(root_node, linear, 1)
    BVHAccelerator(ordered, linear)
end

#function intersectP(ray::Ray, BB::BoundingBox, )


function intersect(ray::Ray, BVH::BVHAccelerator)
    if length(BVH.nodes) == 0
        return Nullable{Intersection}()
    end
    nodeN = 1
    todo_offset = 0
    todo = zeros(Int32, 64)
    dir_is_neg = [ray.direction[i] < 0.0 for i = 1:3]
    while true
        # Check the node given by the index nodeN
        node = BVH.nodes[nodeN]
        if intersectP(ray, node.BBox)
            # We hit the bounding box of this node.
            if node.leaf
                # This is a leaf node.
                # Check intersection with the primitive in the node and return if it hits.
                isect = intersect(ray, BVH.primitives[node.offset])
                if !isnull(isect)
                    return isect
                end
                # It didn't hit the primive.
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
    return Nullable{Intersection}()
end


function intersectP(ray::Ray, BVH::BVHAccelerator)
    if length(BVH.nodes) == 0
        return false
    end
    nodeN = 1
    todo_offset = 1
    todo = zeros(Int32, 64)
    dir_is_neg = [ray.direction[i] < 0.0 for i = 1:3]
    while true
        node = BVH.nodes[nodeN]
        if intersectP(ray, node.BBox)
            if node.leaf
                if intersectP(ray, BVH.primitives[node.offset])
                    return true
                end
                if todo_offset == 1
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
                    nodeN += 1
                end
            end
        else
            if todo_offset == 1
                break
            end
            nodeN = todo[todo_offset]
            todo_offset -= 1
            
        end
    end
    return false
end
