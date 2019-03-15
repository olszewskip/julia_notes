# © Przemyslaw Szufel, 2018, run under Julia 1.0.0

#Usage
#julia -p 5 schelling2.jl 2x2 35X35 slurm=false plot_field=false plot_ascii=false max_step=100

using Random, Printf, Statistics
import Future

const cpu_grid = [2,2]
plot_field = false
slurm = false
plot_ascii = false
const field_size = [35,35]
max_step = 100
for arg in ARGS
    if startswith(arg,"max_step=")
        max_step = parse(Int,arg[length("max_step=")+1:end])
    elseif occursin("x", arg)
        ix = findfirst("x",arg)[1]
        cpu_grid[1]=parse(Int,arg[1:(ix-1)])
        cpu_grid[2]=parse(Int,arg[(ix+1):end])
    elseif occursin("X",arg)
        ix = findfirst("X",arg)[1]
        field_size[1]=parse(Int,arg[1:(ix-1)])
        field_size[2]=parse(Int,arg[(ix+1):end])
    elseif startswith(arg,"slurm=")
        slurm = parse(Bool,arg[length("slurm=")+1:end])
    elseif startswith(arg,"plot_field=")
        plot_field = parse(Bool,arg[length("plot_field=")+1:end])
    elseif startswith(arg,"plot_ascii=")
        plot_ascii = parse(Bool,arg[length("plot_ascii=")+1:end])
    else
        println("Sample usage\njulia  -p 5 schelling2.jl 2x2 35X35 slurm=false plot_field=false plot_ascii=false max_step=100")
        exit(1)
    end
end

const GRID_SIZE = (cpu_grid[1],cpu_grid[2])
plot_field=true
plot_ascii=true
#println(plot_field, plot_ascii, GRID_SIZE)


if slurm
    using ClusterManagers
    addprocs_slurm(GRID_SIZE[1]*GRID_SIZE[2]+1,job_name="schelling2", account="GC71-37", time="01:00:00", exename="/lustre/tetyda/home/pszufe/julia0.6.2_v002/usr/bin/julia")
    cd("/lustre/tetyda/home/pszufe")
else
    #aws version requires machinefile to be present, e.g.
    #julia --machinefile machinefile schelling2.jl .....
    #or use -p parameter
    #julia -p 17 schelling2.jl ............
    min_procs = GRID_SIZE[1]*GRID_SIZE[2]+1
    if nprocs() <= min_procs
        println("For this config at least ",min_procs, " workers are required")
        exit(1)
    end
end

if plot_field
    include("schelling2_plot.jl");
end

const p2p = true

const mywks = workers()[1:(GRID_SIZE[1]*GRID_SIZE[2]+1)]

@everywhere using Random, Printf, Statistics
@everywhere include("schelling2_everywhere.jl");
const modelC = ModelConfig(GRID_SIZE, field_size[1],field_size[2],[0.1,0.45,0.45],5,mywks[1:length(mywks)-1],mywks[length(mywks)])

@everywhere using ParallelDataTransfer
sendto(mywks, modelC = modelC)
sendto(modelC.master_id, p2p = p2p)

const n_rngs = length(modelC.slaves)+1
const rngs = accumulate(Future.randjump, [0; [big(10)^20 for i in 1:n_rngs-1]], init = MersenneTwister(0))

sendto(modelC.master_id, rng=rngs[length(rngs)])

const futures = Dict{Int,Distributed.Future}()

const worker_locations =  Dict{Int,Tuple{Int,Int}}()
const worker_grid = zeros(Int,GRID_SIZE)

function set_grid()
    ix=0
    for w in modelC.slaves
        i =  Int(floor(ix / GRID_SIZE[2])) +1
        j =  ix - ( (i-1)*GRID_SIZE[2] ) +1
        ix += 1
        #print("Creating worker ",w," ix=",ix," i=",i," j=",j,"\n")
        #creates a global variable at each worker
        sendto(w, node = NodeData(i, j, zeros(Int8,(0,0)),Vector{Tuple{Int,Int}}(),rngs[ix]))
        futures[w] = @spawnat w init_node()
        worker_locations[w] = (i,j)
        worker_grid[i,j] = w
    end
end
set_grid()

print("Waiting for ",length(modelC.slaves)," slaves\n")
flush(stdout)
sendto(modelC.master_id, worker_grid=worker_grid)
sendto(modelC.master_id, worker_locations=worker_locations)
for w in modelC.slaves
    init_ok = fetch(futures[w])
    sendto(w, worker_grid = worker_grid)
end
println("Completed creating workers max slave id ",maximum(modelC.slaves)," p2p=",p2p)
flush(stdout)

include("schelling2_master.jl")
fetch(@spawnat modelC.master_id include("schelling2_master.jl"))



function get_aggField(modelC::ModelConfig, GRID_SIZE::Tuple{Int,Int}, worker_grid::Matrix{Int},node_colors::Bool)::Matrix{Int8}
    dim_i=modelC.dim_i
    dim_j=modelC.dim_j
    aggField = zeros(Int8,modelC.dim_i*GRID_SIZE[1],modelC.dim_j*GRID_SIZE[2])
    for gi=1:GRID_SIZE[1], gj=1:GRID_SIZE[2]
        begin
            field = remotecall_fetch(get_field_interior,worker_grid[gi,gj])
            for i=1:modelC.dim_i, j=1:modelC.dim_j
                ii=dim_i*(gi-1)+i
                jj=dim_j*(gj-1)+j
                aggField[ii,jj]=field[i,j]
                if node_colors && aggField[ii,jj]==0
                    aggField[ii,jj]=length(modelC.probs_w)+GRID_SIZE[2]^(gi-1)+(gj-1)
                end
            end
        end
    end
    return aggField
end


function main_loop(modelC :: ModelConfig,worker_grid::Matrix{Int},worker_locations::Dict{Int,Tuple{Int,Int}},plotField::Bool,plotAscii::Bool,maxStep::Int)
    img = nothing
    for stepNo in 1:maxStep
        res = remotecall_fetch(do_item, modelC.master_id )
        println("[",stepNo,"] ",res)
        flush(stdout)

        if plotField
            aggField = get_aggField(modelC,modelC.GRID_SIZE,worker_grid,true)
            img = plot_schelling(aggField,stepNo,img)
        end
    end

    if plotAscii
        aggField = get_aggField(modelC,(min(GRID_SIZE[1],2),min(GRID_SIZE[2],2)),worker_grid,false)
        for i in 1:size(aggField)[1]
            for j=1:size(aggField)[2]
                if aggField[i,j]==0
                    print(".")
                else
                    print(aggField[i,j])
                end
            end
            print("\n")
        end
    end


end

main_loop(modelC,worker_grid,worker_locations,plot_field,plot_ascii,max_step)
if plot_field
    println("Press <ENTER>")
    readline()
end
