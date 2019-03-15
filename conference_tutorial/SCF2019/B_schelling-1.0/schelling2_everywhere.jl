# © Przemyslaw Szufel, 2018, run under Julia 1.0.0

struct ModelConfig
    # computational grid size
    GRID_SIZE :: Tuple{Int,Int}
    #data matrix size on a single process (all matrices are squares)
    dim_i::Int
    dim_j::Int
    #probabilities for empty, red and blue cells
    probs_w::Vector{Float64}
    #minimal required number of neighbours
    min_neighbours::Int
    slaves::Vector{Int64}
    master_id::Int64
end

mutable struct NodeData
    grid_i :: Int
    grid_j :: Int
    field :: Matrix{Int8}
    free_cells::Vector{Tuple{Int,Int}}
    rng::MersenneTwister
end

struct BorderData
    left :: Vector{Int8}
    right :: Vector{Int8}
    top :: Vector{Int8}
    bottom :: Vector{Int8}
end


struct CornerData
    topleft :: Int8
    topright :: Int8
    bottomleft :: Int8
    bottomright :: Int8
end

#The below const statements could be an alternative to barrier functions used
#const modelC  = ModelConfig((0,0),0,0,Vector{Float64}(),0)
#const node = NodeData(0,0,Matrix{Int8}(0,0),Vector{Tuple{Int,Int}}(),MersenneTwister(0),(0,0))
#const worker_grid = Matrix{Int}(0,0)
function get_field_interior():: Matrix{Int8}
    dim_i = modelC.dim_i
    dim_j = modelC.dim_j
    field = zeros(Int8,(dim_i,dim_j))
    for i = (1+1):(dim_i+1), j = (1+1):(dim_j+1)
        field[i-1,j-1]=node.field[i,j]
    end
    return field
end

function sample(rng::AbstractRNG, weights::Vector{Float64})::Int
    t = rand(rng) * sum(weights)
    n = length(weights)
    i = 1
    cw = weights[1]
    while cw < t && i < n
        i += 1
        @inbounds cw += weights[i]
    end
    return i
end

function init_node()::Bool
    dim_i = modelC.dim_i
    dim_j = modelC.dim_j
    values = Int8(0):Int8(length(modelC.probs_w)-1)
    node.field = zeros(Int8,(1+dim_i+1,1+dim_j+1))
    node.free_cells=Vector{Tuple{Int,Int}}()
    for i = (1+1):(dim_i+1), j = (1+1):(dim_j+1)
       val = values[sample(node.rng,modelC.probs_w)]
       if val == 0
           push!(node.free_cells,(i,j))
       end
       node.field[i,j] = val
    end
    return true
end


function get_border()::BorderData
    return get_border(modelC, node)
end
function get_border(modelC :: ModelConfig, node ::NodeData)::BorderData
    dim_i = modelC.dim_i
    dim_j = modelC.dim_j
    #Please note that 4 bytes are redundant, however this increases code readibility
    return BorderData(
        node.field[(1+1):(1+dim_i),1+1],
        node.field[(1+1):(1+dim_i),1+dim_j],
        node.field[1+1,(1+1):(1+dim_j)],
        node.field[1+dim_i,(1+1):(1+dim_j)] )
end

function get_border_top()::Vector{Int8}
    return get_border_top(modelC, node)
end
function get_border_bottom()::Vector{Int8}
    return get_border_bottom(modelC, node)
end
function get_border_left()::Vector{Int8}
    return get_border_left(modelC, node)
end
function get_border_right()::Vector{Int8}
    return get_border_right(modelC, node)
end

function get_border_top(modelC :: ModelConfig, node ::NodeData)::Vector{Int8}
    return node.field[(1+1),(1+1):(1+modelC.dim_j)]
end
function get_border_bottom(modelC :: ModelConfig, node ::NodeData)::Vector{Int8}
    return node.field[(1+modelC.dim_i),(1+1):(1+modelC.dim_j)]
end
function get_border_left(modelC :: ModelConfig, node ::NodeData)::Vector{Int8}
    return node.field[(1+1):(1+modelC.dim_i),(1+1)]
end
function get_border_right(modelC :: ModelConfig, node ::NodeData)::Vector{Int8}
    return node.field[(1+1):(1+modelC.dim_i),(1+modelC.dim_j)]
end

function get_corner_top_left()::Int8
    return get_corner_top_left(modelC, node)
end
function get_corner_top_right()::Int8
    return get_corner_top_right(modelC, node)
end
function get_corner_bottom_left()::Int8
    return get_corner_bottom_left(modelC, node)
end
function get_corner_bottom_right()::Int8
    return get_corner_bottom_right(modelC, node)
end


function get_corner_top_left(modelC :: ModelConfig, node ::NodeData)::Int8
    return node.field[1+1,1+1]
end
function get_corner_top_right(modelC :: ModelConfig, node ::NodeData)::Int8
    return node.field[1+1,1+modelC.dim_j]
end
function get_corner_bottom_left(modelC :: ModelConfig, node ::NodeData)::Int8
    return node.field[1+modelC.dim_i,1+1]
end
function get_corner_bottom_right(modelC :: ModelConfig, node ::NodeData)::Int8
    return node.field[1+modelC.dim_i,1+modelC.dim_j]
end


function step_setborder_nodes!()::Float64
    return step_setborder_nodes!(modelC, node)
end
function step_setborder_nodes!(modelC :: ModelConfig, node ::NodeData)::Bool
    start_t = time()
    dim_i = modelC.dim_i
    dim_j = modelC.dim_j
    #borders
    for j = 1:dim_j
        node.field[1,j+1]=border_data.top[j]
        node.field[1+dim_i+1,j+1]=border_data.bottom[j]
    end
    for i = 1:dim_i
        node.field[i+1,1]=border_data.left[i]
        node.field[i+1,1+dim_j+1]=border_data.right[i]
    end
    node.field[1,1] = corner_data.topleft
    node.field[1,1+dim_j+1] = corner_data.topright
    node.field[1+dim_i+1,1] = corner_data.bottomleft
    node.field[1+dim_i+1,1+dim_j+1] = corner_data.bottomright
    return time()-start_t
end

function step_setborder_nodes_p2p!()::Float64
    return step_setborder_nodes_p2p!(modelC, node, worker_grid)
end
function step_setborder_nodes_p2p!(modelC :: ModelConfig, node ::NodeData, worker_grid::Matrix{Int})::Float64
    start_t = time()
    j_left = node.grid_j-1
    if j_left < 1; j_left = modelC.GRID_SIZE[2]; end
    j_right = node.grid_j+1
    if j_right > modelC.GRID_SIZE[2]; j_right = 1; end
    i_top = node.grid_i-1
    if i_top < 1; i_top = modelC.GRID_SIZE[1]; end
    i_bottom = node.grid_i+1
    if i_bottom > modelC.GRID_SIZE[1]; i_bottom = 1; end
    dim_i = modelC.dim_i
    dim_j = modelC.dim_j
    #borders
    @sync begin
        #top
        @async node.field[1,(1+1):(dim_j+1)] = remotecall_fetch(get_border_bottom,worker_grid[i_top,node.grid_j])
        #bottom
        @async node.field[(1+dim_i+1),(1+1):(dim_j+1)] = remotecall_fetch(get_border_top,worker_grid[i_bottom,node.grid_j])
        #left
        @async node.field[(1+1):(dim_i+1),1] = remotecall_fetch(get_border_right,worker_grid[node.grid_i,j_left])
        #right
        @async node.field[(1+1):(dim_i+1),(1+dim_j+1)] = remotecall_fetch(get_border_left,worker_grid[node.grid_i,j_right])
        #top-left
        @async node.field[1,1] = remotecall_fetch(get_corner_bottom_right,worker_grid[i_top,j_left])
        #top-right
        @async node.field[1,(1+dim_j+1)] = remotecall_fetch(get_corner_bottom_left,worker_grid[i_top,j_right])
        #bottom-left
        @async node.field[(1+dim_i+1),1] = remotecall_fetch(get_corner_top_right,worker_grid[i_bottom,j_left])
        #bottom-right
        @async node.field[(1+dim_i+1),(1+dim_j+1)] = remotecall_fetch(get_corner_top_left,worker_grid[i_bottom,j_right])
    end
    return time()-start_t
end




struct LeaveLocData
    relocate_counts :: Vector{Int}
    free_cells_count :: Int
    elapsed_time :: Float64
end

function step_leave_location!()::LeaveLocData
    return step_leave_location!(modelC, node)::LeaveLocData
end
function step_leave_location!(modelC :: ModelConfig, node ::NodeData)::LeaveLocData
    start_t = time()
    field = node.field
    free_count = length(node.free_cells)
    dim_i = modelC.dim_i
    dim_j = modelC.dim_j
    ds = ((1,0),  (0,-1), (0,1), (-1,0), (1,1),  (1,-1), (-1,1), (-1,-1))
    for i = (1+1):(dim_i+1), j = (1+1):(dim_j+1)
        if field[i,j] > 0
            c = 0
            for (d1, d2) ∈ ds
                c += field[i,j] == field[i+d1,j+d2]
            end
            #print("c=",c," i=",i," j=",j," f=",field[i-1:i+1,j-1:j+1],"\r\n")
            if c < modelC.min_neighbours
                push!(node.free_cells,(i,j))
            end
        end
    end
    relocate_counts::Vector{Int} = zeros(Int,length(modelC.probs_w)-1)
    for (i, j) in node.free_cells[free_count+1:length(node.free_cells)]
        relocate_counts[field[i,j]] +=1
        field[i,j]=0
    end
    return LeaveLocData(relocate_counts,length(node.free_cells),time()-start_t)
end

function step_set_new_location!()::Float64
    return step_set_new_location!(modelC, node)
end
function step_set_new_location!(modelC :: ModelConfig, node ::NodeData)::Float64
    start_t = time()
    shuffle!(node.rng,node.free_cells)
    ix = length(node.free_cells)
    for z = 1:length(new_citizens)
        for a = 1:new_citizens[z]
            node.field[node.free_cells[ix]...]=z
            ix -=1
            pop!(node.free_cells)
        end
    end
    return time()-start_t
end
