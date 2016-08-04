#!/usr/local/bin/julia

if !in(".", LOAD_PATH)
    push!(LOAD_PATH, ".")
end
using Disk
processes = []

for v in [:DSB] #, :LP, :RZ]
    #= push!(processes, @spawn Disk.simring(v)) =#
    #= push!(processes, @spawn Disk.simIllenseer(v)) =#
    #= push!(processes, @spawn Disk.simIllenseer(v, =#
    #=                   timestring=:(collect(linspace(1, 1e6, 1e6))), =#
    #=                   rstring=:(collect(linspace(-6, 25, 200))))) =#

    Disk.simIllenseer(v,
                      timestring=:(collect(linspace(1, 1e6, 1e6))),
                      rstring=:(collect(linspace(-6, 25, 200))))
    #= Disk.simring(v, =#
    #=                   timestring=:(collect(linspace(1, 1e6, 1e6))), =#
    #=                   rstring=:(collect(linspace(-6, 25, 200)))) =#
    #= push!(processes, @spawn Disk.simring(v, =#
    #=                   timestring=:(collect(linspace(1, 1e6, 1e6))), =#
    #=                   rstring=:(collect(linspace(-6, 25, 200))))) =#
end

#= for p in processes =#
#=     wait(p) =#
#= end =#
