# © Przemyslaw Szufel, 2018, run under Julia 1.0.0

function do_item(modelC :: ModelConfig,rng::MersenneTwister,worker_grid::Matrix{Int},worker_locations::Dict{Int,Tuple{Int,Int}},p2p::Bool) :: String
    GRID_SIZE=modelC.GRID_SIZE
    status=zeros(Float64,(GRID_SIZE))
    borders=Matrix{BorderData}(undef,GRID_SIZE...)
    leaveLocData=Matrix{LeaveLocData}(undef,GRID_SIZE...)
    dim_i=modelC.dim_i
    dim_j=modelC.dim_j

    avg_asyn_bord = 0.
    max_asyn_bord = 0.

    start = time()
    if p2p
        @sync for ix in modelC.slaves
            (i,j) = worker_locations[ix]
            @async status[i,j] = remotecall_fetch(step_setborder_nodes_p2p!,worker_grid[i,j])
        end
        avg_asyn_bord = mean(status)
        max_asyn_bord = maximum(status)
    else
        @sync for ix in modelC.slaves
            (i,j) = worker_locations[ix]
            @async borders[i,j] = remotecall_fetch(get_border,worker_grid[i,j])
        end


        @sync for ix in modelC.slaves
            (i,j) = worker_locations[ix]
            j_left = j-1
            if j_left < 1; j_left = GRID_SIZE[2]; end
            j_right = j+1
            if j_right > GRID_SIZE[2]; j_right = 1; end
            i_top = i-1
            if i_top < 1; i_top = GRID_SIZE[1]; end
            i_bottom = i+1
            if i_bottom > GRID_SIZE[1]; i_bottom = 1; end

            sendto(worker_grid[i,j],border_data=BorderData( #left right top bottom
                    borders[i,j_left].right,  borders[i,j_right].left,
                    borders[i_top,j].bottom,  borders[i_bottom,j].top),
                corner_data=CornerData( #topleft topright bottomleft bottomright
                    borders[i_top,j_left].bottom[modelC.dim_j],  borders[i_top,j_right].bottom[1],
                    borders[i_bottom,j_left].top[modelC.dim_j],  borders[i_bottom,j_right].top[1])
                    )
            @async status[i,j] = remotecall_fetch(step_setborder_nodes!,worker_grid[i,j])
        end
    end

    @sync for ix in modelC.slaves
        (i,j) = worker_locations[ix]
        @async begin
            leaveLocData[i,j] = remotecall_fetch(step_leave_location!,worker_grid[i,j])
            status[i,j] = leaveLocData[i,j].elapsed_time
        end
    end
    max_leave_loc_time = maximum(status)
    avg_leave_loc_time = mean(status)


    reloc_start_time = time()

    add_new_citizens=Matrix{Vector{Int}}(undef,GRID_SIZE...)

    # all units to be located on the board of each color type
    relocate_counts = zeros(Int,length(modelC.probs_w)-1)
    free_cells_count=zeros(Int,GRID_SIZE)

    for ix in modelC.slaves
        (i,j) = worker_locations[ix]
        relocate_counts += leaveLocData[i,j].relocate_counts
        add_new_citizens[i,j]=zeros(Int,length(modelC.probs_w)-1)
        free_cells_count[i,j]=leaveLocData[i,j].free_cells_count
    end


    while sum(relocate_counts) > 0
        for t in 1:length(relocate_counts)
            if relocate_counts[t] == 0
                continue
            end
            relocated = 0
            total_free_cells_count = sum(free_cells_count)
            for ix in shuffle(rng,modelC.slaves)
                (i,j) = worker_locations[ix]
                ρ = free_cells_count[i,j]/total_free_cells_count
                μ = ρ * relocate_counts[t]
                if μ < 2.0
                    μ = 2.0
                end
                σ = √(μ * (1-ρ))
                val = round(randn(rng)*σ+μ)
                t_rand = max(0, Int(min(val,relocate_counts[t]-relocated,free_cells_count[i,j]))  )
                free_cells_count[i,j] -= t_rand
                add_new_citizens[i,j][t] += t_rand
                relocated += t_rand
            end
            relocate_counts[t] -= relocated
        end
    end
    reloc_time = time() - reloc_start_time

    @sync for ix in modelC.slaves
        (i,j) = worker_locations[ix]

        sendto(worker_grid[i,j],new_citizens=add_new_citizens[i,j])
        @async status[i,j] = remotecall_fetch(step_set_new_location!,worker_grid[i,j])
    end
    max_new_loc_time = maximum(status)
    avg_new_loc_time = mean(status)
    elapsed_time = time()-start
    return string("$GRID_SIZE x $dim_i x $dim_j t: $(@sprintf("%.3f", elapsed_time  )) = loc $(@sprintf("%.3f", reloc_time  ))  avg/max_asyn_bord $(@sprintf("%.3f", avg_asyn_bord  ))/$(@sprintf("%.3f", max_asyn_bord  )) avg/max_leave $(@sprintf("%.3f", avg_leave_loc_time))/$(@sprintf("%.3f", max_leave_loc_time)) avg/max_new $(@sprintf("%.3f", avg_new_loc_time  ))/$(@sprintf("%.3f", max_new_loc_time  )) ")


end
function do_item() :: String
    return do_item(modelC,rng,worker_grid,worker_locations,p2p)
end
