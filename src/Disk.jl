"""
The purpose of this module is to simulate the time evolution of
self-gravitating accretion disks using a one-dimensional model.

# Usage examples

To run a simulation, call `include("Disk.jl")` and call, e.g.
`Disk.simIllenseer(:DSB)` to start the simulation. 

# Output files

The generated files will be put into the directory
`../sims/<name>-<viscositymodel>/`. The following files will be generated:

* `settings.txt`: Contains the parameters used for the simulation
* `Ωseries-<name>-<viscositymodel>.csv`: Table with four columns: ln(r),
ln(Ω(τ=1)), ln(Ω(τ=1000)), ln(Ω(τ=1000000))
* Various self-explanatory plotted images
"""
module Disk

using Sundials
#= using Documenter =#
using Loess
#= include("/Users/Nino/Google Drive/Knowledge/InlinePlot.jl") =#
#= using Plots =#
#= pyplot() =#

@deprecate setviscosityat! getviscosity
@deprecate Ωfixboundaries! x->nothing
@deprecate extendtoghostcellslinear! x->nothing
@deprecate setviscosityat! getviscosity


"Gravitational constant. Currently not in use."
const G = Float64(1.0)

"Mass of the central accreting object. Useful for deriving Ω from Σ."
const Mstar = Float64(1.0)

"""
Available options for viscosity models to pass as parameters to the simulation
code.
"""
const viscositymodels = [:DSB, :RZ, :LP]

"""
Simulation output path, relative to the main source file.
"""
const simsdirectory = "../sims"

"""
**Deprecated.**

Extrapolate values for `Ω` in the ghostcells and update the values in place. At
the origin, assume Keplerian rotation. In the large r limit, detect the slope
and extrapolate linearly.
"""
function Ωfixboundaries!(Ω::Array{Float64,1}, r::Array{Float64,1}, ghostcells::Int)
    for i in 1:ghostcells
        Ω[i] = -1.5 * (r[i] - r[ghostcells+1]) + Ω[ghostcells+1]
    end
    rightslope = (Ω[end-ghostcells-1] - Ω[end-ghostcells]) / (r[end-ghostcells-1] - r[end-ghostcells])
    for i in (length(r)-ghostcells+1):length(r)
        Ω[i] = rightslope * (r[i] - r[end-ghostcells]) + Ω[end-ghostcells]
    end
end

"""
**Deprecated**

Update ghostcells in place like `Ωfixboundaries()`, but extrapolate linearly at
the origin, too, instead of using the exponent -3/2.
"""
function extendtoghostcellslinear!(input::Array{Float64,1}, r::Array{Float64,1}, ghostcells::Int)
    leftslope  = (input[ghostcells+2] - input[ghostcells+1]) / (r[ghostcells+2] - r[ghostcells+1])
    rightslope = (input[end-ghostcells-1] - input[end-ghostcells]) / (r[end-ghostcells-1] - r[end-ghostcells])
    for i in 1:ghostcells
        input[i] = leftslope * (r[i] - r[ghostcells+1]) + input[ghostcells+1]
    end
    for i in (length(r)-ghostcells+1):length(r)
        input[i] = rightslope * (r[i] - r[end-ghostcells]) + input[end-ghostcells]
    end
end

"""
Calculate the right hand side of the main differential equation that defines
Ω. Write the resulting values for Ωdot into the passed variable Ωdot and return
zero.

Variables `x`, `z`, `f`, `dzdr` must be vectors of the same length as `r`. They
are used as persistent storage for their respective physical quantities in order
to avoid repeated memory allocations. This means `rhs()` will mutate all
arguments except `r`, `t`, and `viscositymodel`.

To monitor progress when run through Sundials, `rhs()` will occasionally output
the current time `t`.
"""
function rhs(Ω :: Array{Float64, 1}, Ωdot :: Array{Float64, 1}, t :: Float64,
             r :: Array{Float64, 1}, viscositymodel :: Symbol,
             x :: Array{Float64, 1}, z :: Array{Float64, 1},
             f :: Array{Float64, 1}, dxdr :: Array{Float64, 1})

    for i=1:(length(x)-1)
        x[i] = derivativeatinbounds(Ω, r, i)
        f[i] = getviscosity(viscositymodel, x[i])
        z[i] = getz(x[i], f[i])
    end
    x[end] = x[end-1] + (x[end-1]-x[end-2])
    f[end] = getviscosity(viscositymodel, x[end])
    z[end] = getz(x[end], f[end])
    dxdr[end] = dxdr[end-1]

    ghostx = -1.5 # default: for flow away from center
    if flowsleft(x[1]) # alternative: for flow into the center
        q = exp(1.5*(r[1]-r[2]))
        ghostx = q * x[1] + 1.5 * (q-1)
    end
    ghostf = getviscosity(viscositymodel, ghostx)
    ghostz = getz(ghostx, ghostf)
    ghostdxdr = (x[1]-ghostx)/(r[2]-r[1])
    up_x, up_z = upwindvalues_x_z(ghostx, x[1], ghostz, z[1], viscositymodel)
    ghostdzdr = ghostdxdr * (4*up_x + 3)

    Ωdot[1] = getΩdot(Ω[1], up_x, up_z, ghostdzdr)
    for i in 2:length(r)
        dxdr[i] = derivativeatinbounds(x, r, i-1)
        (up_x, up_z) = upwindvalues_x_z(x[i-1], x[i], z[i-1], z[i], viscositymodel)
        dzdx = (4*up_x + 3)
        Ωdot[i] = getΩdot(Ω[i], up_x, up_z, dzdx*dxdr[i])
    end
    return Int32(0)
end

@inline getz(x, f) = x*(2x+3)*f
@inline getΩdot(Ω, x, z, dzdr) = -exp(Ω) * (z * (5 + 4 * x) + dzdr)
@inline flowsright(x) = x< -1.25 #  -1.5 < x < -1.25
@inline flowsleft(x) = !flowsright(x)
@inline getζ(x, z) = (4x + 5)*z


@inline function upwindvalues_x_z(x1, x2, z1, z2, viscositymodel)
    if flowsright(x1) && flowsright(x2)
        return (x1, z1)
    elseif  flowsleft(x1) && flowsleft(x2)
        return (x2, z2)
    else
        return (mean([x1, x2]), mean([z1, z2]))
    end
end

"""
Calculate the derivative of `quantity` with respect to `var` and write the
result into `deriv`.
"""
function derivative!(quantity::Array{Float64,1}, var::Array{Float64,1}, deriv::Array{Float64,1})
    for i in 1:(length(var)-1)
        deriv[i] = (quantity[i+1]-quantity[i])/(var[i+1]-var[i])
    end
end

"""
Return the derivative of `quantity` with respect to `var` at the index `index`.
Return 0 if index is 1, `length(var)`, or out of bounds.
"""
function derivativeat(quantity::Array{Float64,1}, var::Array{Float64,1}, index::Int64)
    if index ≥ 1 && index < length(var)
        (quantity[index+1]-quantity[index])/(var[index+1]-var[index])
    else
        0.0
    end
end

"""
Like `derivativeat()` but does not check whether `index` is within bounds.
Very fast, but dangerous.
"""
@inline function derivativeatinbounds(quantity::Array{Float64,1}, var::Array{Float64,1}, index::Int64)
    return (quantity[index+1]-quantity[index])/(var[index+1]-var[index])
end

function deriv_secondorder(quantity, var, index)
    return (3quantity[index] - 4quantity[index-1] - quantity[index-2])/(var[i]-var[i-2])
end

"""
**Deprecated -- use `getviscosity()**

Calculate values for the viscosity function f(x) (as used in the paper by
Illenseer (2015)) for the `i`-th element and write the result into `f[i]`.
"""
function setviscosityat!(model::Symbol, x::Array{Float64,1}, i :: Int64, f::Array{Float64,1})
    if model == :DSB
        f[i] = 1.
    elseif model == :RZ
        f[i] = abs(x[i])
    elseif model == :LP
        f[i] = (2 * x[i] + 3)^2
    end
end

@inline function getviscosity(model::Symbol, x::Float64)
    1.0
    #= if model == :DSB =#
    #=     1. =#
    #= elseif model == :RZ =#
    #=     abs(x) =#
    #= else # if model == :LP =#
    #=     (2 * x + 3)^2 =#
    #= end =#
end

"""
Return the logarithmic power law exponent x (as used in the paper by
Illenseer (2015)) given Ω, for the entire solution space `r`.
"""
function getx(Ω, r)
    x = zeros(size(r))
    getx!(Ω, r, x)
    x
end

"""
Like `getx()`, but write the results into the variable `x` instead of returning.
"""
function getx!(Ω::Array{Float64,1}, r::Array{Float64,1}, x::Array{Float64,1})
    for i in 1:(length(x)-1)
        x[i] = derivativeatinbounds(Ω, r, i)
    end
end

"""
Return the logarithm of the surface density Σ(r) for a given angular velocity
distribution Ω.

Note: Values near the origin are prone to rounding errors since the formula
becomes 0⋅∞ in the limit. Applying `smoothcurve!()` to the result can mitigate
this effect.
"""
function getΣ(Ω, r)
    Σ = Array(Float64, length(r))
    x = getx(Ω, r)
    for i in 2:(length(r)-1)
        Σ[i] = exp(r[i] + 2Ω[i])*(2x[i]+3)/(2π)
        #= Σ[i] = (2 * exp(Ω[i])^2 * exp(r[i]) * derivativeatinbounds(Ω, r, i) + =#
        #=          3 * exp(r[i]) * exp(Ω[i])^2) =#
        Σ[i] = log(max(Σ[i], 1e-200))
    end
    Σ
end

function getΣlinear(Ω, r)
    Σ = Array(Float64, length(r))
    for i in 1:(length(r)-1)
        Σ[i] = (2 * exp(Ω[i])^2 * exp(r[i]) * derivativeatinbounds(Ω, r, i) +
                 3 * exp(r[i]) * exp(Ω[i])^2)
    end
    Σ[end] = Σ[end-1]
    Σ
end

"""
Return the accretion rate given an angular velocity distribution Ω(r) and a
choice of viscosity model.
"""
function getṀ(Ω, r, viscositymodel)
    Ω    = copy(Ω)
    Ωdot    = Array(Float64, size(r))
    x    = Array(Float64, size(r))
    z    = Array(Float64, size(r))
    f    = Array(Float64, size(r))
    dzdr = Array(Float64, size(r))
    #= global currenttime = Float64(0) =#
    rhs(Ω, Ωdot, Float64(1), r, viscositymodel, x, z, f, dzdr)
    β = 1000
    Ṁ = 2 .* exp(2Ω .+ 3r) .* Ωdot
end

"""
Return values of a normal distribution with mean `μ` and standard
deviation `σ` along the vector `x`. (This `x` has nothing to do with
the power law exponent `x`.)
"""
function gauss(x, μ, σ)
    1/(σ*sqrt(2π)) * exp(-(x.-μ).^2./(2σ^2))
end

"""
Apply sort of a gaussian blur filter to the vector `curve`. `σ` is the standard
deviation of the bell curve.
"""
function smoothcurve!(curve, σ)
    r = collect(eachindex(curve))
    len = length(curve)
    referencecurve = Array(Float64, len*3)
    longr = collect(eachindex(referencecurve))
    referencecurve[(len+1):(2*len)] = curve[:]
    leftslope  = (referencecurve[len+2] - referencecurve[len+1]) / (longr[len+2] - longr[len+1])
    rightslope = (referencecurve[end-len-1] - referencecurve[end-len]) / (longr[end-len-1] - longr[end-len])
    for i in 1:len
        referencecurve[i] = leftslope * (longr[i] - longr[len+1]) + referencecurve[len+1]
    end
    for i in (2len+1):3len
        referencecurve[i] = rightslope * (longr[i] - longr[end-len]) + referencecurve[end-len]
    end
    for i in 1:len
        curve[i] = sum(referencecurve .* gauss(eachindex(referencecurve), i+len, σ))
    end
end

"""
    conformtoarray(quantity, from, to)

Return an Array of values such that the mapping from `from` to `quantity` is
equivalent to the mapping from `to` to the returned value.

# Example

Say we have values of a quantity q(r), where `q` and `r` are vectors of the same
length, and we want the values for q(x) for another vector `x`.

```julia
q = [3, 4, 8, 6, 145, 3, 5]
r = [1, 3, 4, 5, 5.6, 6, 8]
x = [0, 1, 2.5, 4, 5, 6]
output = conformtoarray(q, r, x)
# => [2.5, 3.0, 3.75, 8.0, 6.0, 3.0]
```
"""
function conformtoarray(quantity, from, to)
    # create copies of the data
    from = copy(from)
    to   = copy(to)

    # Error handling:
    if from[1] > to[end] || to[1] > from[end]
        return zeros(size(to))
    end

    # Constrain `quantity` and `from` to the boundaries
    # of `out`.
    #= minindex = max(findlast(el->from[el] <= to[1], eachindex(from)), 1) =#
    #= maxindex = max(findfirst(el->from[el] >= to[end], eachindex(from)), length(from)) =#
    #= from     = from[minindex:maxindex] =#
    #= quantity = quantity[minindex:maxindex] =#

    # loess model
    model = loess(from, quantity, span=0.1)
    minindex = findfirst(e->e>=min(from...), to)
    maxindex = findlast(e-> e<=max(from...), to)

    # Populate resulting array
    result = Array(Float64, size(to))
    result[minindex:maxindex] = predict(model, to[minindex:maxindex])
    #= plot(to, result) =#
    #= plot!(from, quantity) =#
    #= png("tmp/cnfrm") =#
    #= for i in 1:(length(from)-1) =#
    #=     # Find all elements in `to` that lie between `from[i]` and `from[i+1]` =#
    #=     indices = filter(index->(to[index] >= from[i] && to[index] < from[i+1]), eachindex(to)) =#

    #=     # Populate `result` at all `indices` =#
    #=     for j in indices =#
    #=         fromdist = from[i+1] - from[i] =#
    #=         leftcomponent = quantity[i] * (from[i+1]-to[j])/fromdist =#
    #=         rightcomponent = quantity[i+1] * (to[j]-from[i])/fromdist =#
    #=         result[j] = leftcomponent + rightcomponent =#
    #=     end =#
    #= end =#

    # Extrapolate beyond the edges of `from`
    minindex_to = findlast(el->to[el] <= from[1], eachindex(to))
    maxindex_to = findfirst(el->to[el] >= from[end], eachindex(to))
    if minindex_to > 1 # TODO maybe fix this at some point
                       # (look at the if maxindex < length.. 
                       # thing to see what I mean)
        slope = (quantity[2]-quantity[1])/(from[2]-from[1])
        for j in (minindex_to-1):-1:1
            result[j] = quantity[1] - slope * (from[1] - to[j])
        end
    end
    if maxindex < length(to)
        slope = (result[maxindex]-result[maxindex-1])/(to[maxindex]-to[maxindex-1])
        for j in (maxindex+1):length(to)
            result[j] = result[maxindex] + slope * (to[j] - to[maxindex])
        end
    end
    return result
end

"""
    enclosedmass(Σ, r, Mstar)

Return an `Array` containing the values of the enclosed mass formula
M(r) for each cell in `r` by integrating over the provided surface
density `Σ`. `Mstar` is the mass of the central object and thus the
value of the first cell.
"""
function enclosedmass(Σ, r, Mstar)
    M = Array(Float64, length(r))
    ρ = exp(r)
    M[1] = Mstar
    for i in 2:length(r)
        M[i] = M[i-1] + 2π*(ρ[i]*Σ[i]*(ρ[i]-ρ[i-1]))
    end
    return M
end

"""
    ΩfromΣ(Σ, r, Mstar)

Return the angular velocity Ω(r) for all values of `r` from the
provided surface density `Σ`. This is useful for generating a starting
value for `Ω` that corresponds to a particular initial distribution of
mass.

* `Mstar`: The central object's mass
"""
function ΩfromΣ(Σ :: Array{Float64,1}, r :: Array{Float64,1}, Mstar :: Float64)
    Ω = Array(Float64, length(r))
    M = enclosedmass(Σ, r, Mstar)
    for i in 1:length(r)
        Ω[i] = log(sqrt(M[i]*exp(-3r[i])))
    end
    #= Ωfixboundaries!(Ω, r) =#
    return Ω
end

function cvodefun(t::Float64, Ω::Sundials.N_Vector, Ωdot::Sundials.N_Vector,
                  solverfun::Function)
    Ω = Sundials.asarray(Ω)
    Ωdot = Sundials.asarray(Ωdot)
    solverfun(t, Ω, Ωdot)
    return Int32(0)
end

"""
    runsimwithparams(name, comment, G, Mstar, viscositymodel,
                     r, rstring, time, timestring, Ω₀, Ω₀string)

Integrate `Ω` over time using provided initial conditions. Return `nothing`.

Write out 7 states of `Ω` into a `CSV` file and write plots of these states
as `PNG` images. The exported time steps will be spaced exponentially over
the given time range.

# Arguments

* `name [String]`: Name of the simulation. Will be used for directory and
  file names
* `comment [String]`: Additional comments to be put in the `settings.txt` file
* `G [Float64]`: Gravitational constant. Currently not being
  used -- TODO: Remove
* `Mstar [Float64]`: Mass of central object
* `viscositymodel [Symbol]`: Viscosity model to use. Can be `:DSB`, `:LP`,
  or `RZ`.
* `r [Array{Float64, 1}]`: Logarithmic radial coordinate vector
* `rstring [String|Symbol]`: String representation of `r`. Will be written
  into `settings.txt` to ensure replicability of simulations.
* `time [Array{Float64, 1}]`: Time steps for which the solver will generate
  solutions
* `timestring [String|Symbol]`: String representation of `time`, like `rstring`
* `Ω₀ [Array{Float64, 1}]`: Initial value for `Ω` at each `r`
* `Ω₀string [String|Symbol]`: String representation of `Ω₀`, like `rstring`

# Examples

For examples of how to use this function, see `simring()` and `simIllenseer()`.
"""
function runsimwithparams(name, comment, G, Mstar, viscositymodel,
                          r, rstring, time, timestring, Ω₀, Ω₀string)
    #= global currenttime = Float64(0) =#
    println(name)
    # Set up output files
    run(`mkdir -p $simsdirectory/$name`)
    settingsfile = open("$simsdirectory/$name/settings.txt", "w")
    write(settingsfile, """
    G = $G
    Mstar = $Mstar
    viscositymodel = $viscositymodel
    r = $rstring
    time = $timestring
    Ω₀ = $Ω₀string

    Comment:
    $comment
    """)
    close(settingsfile)
    #= p = plot(r, Ω₀) =#
    #= png(p, "$simsdirectory/$name/omega0-$name") =#

    x    = Array(Float64, size(r))
    z    = Array(Float64, size(r))
    f    = Array(Float64, size(r))
    dzdr = Array(Float64, size(r))

    ####### Solver setup #####
    solverfunction(t::Float64, om::Vector{Float64},
                   omdot::Vector{Float64}) = rhs(om, omdot, t, r,
                                                 viscositymodel, x, z, f, dzdr)

    rlength = length(r)

    ##### Create CVode objects #####
    mem = Sundials.CVodeCreate(Sundials.CV_BDF, Sundials.CV_NEWTON)
    if mem == C_NULL
        error("Failod to allocate CVODE solver object")
    end

    Ωseries = zeros(length(time), length(r))
    reltol = 1e-4
    abstol = 1e-6
    try
        flag = Sundials.CVodeInit(mem,
                                  cfunction(cvodefun, Int32,
                                            (Float64, Sundials.N_Vector,
                                             Sundials.N_Vector, Ref{Function})),
                                  time[1], Sundials.nvector(Ω₀))
        flag = Sundials.CVodeSetUserData(mem, solverfunction)
        flag = Sundials.CVodeSStolerances(mem, reltol, abstol)
        flag = Sundials.CVDense(mem, rlength)
        Ωseries[1,:] = Ω₀
        Ω = copy(Ω₀)
        tout = [0.0]
        #= omplot = plot(title="omega") =#
        #= gplot = plot(title="torque") =#
        #= xplot = plot(title="x") =#
        @time for k in 2:length(time)
            flag = Sundials.CVode(mem, time[k], Ω, tout, Sundials.CV_NORMAL)
            Ωseries[k, :] = Ω
            if mod(k, 10000) == 0
                println(k)
                indices = [round(Int, k^(i//6)) for i in 0:6]
                xseries = Array(Float64,size(Ωseries))
                Ωdotseries = Array(Float64, size(Ωseries))
                for i in indices
                    xseries[i, :]  = getx(collect(Ωseries[i,:]), r)
                    Ω2    = collect(Ωseries[i, :])
                    Ωdot2    = Array(Float64, size(r))
                    x    = Array(Float64, size(r))
                    z    = Array(Float64, size(r))
                    f    = Array(Float64, size(r))
                    dzdr = Array(Float64, size(r))
                    rhs(Ω2, Ωdot2, Float64(1), r, viscositymodel, x, z, f, dzdr)
                    Ωdotseries[i, :] = Ωdot2
                end
                Ωcsv = [r'; Ωseries[indices, :];
                         xseries[indices, :]; Ωdotseries[indices, :]]'
                writecsv("$simsdirectory/$name/omegaseries-$name-k$k.csv", Ωcsv)
                #= plot!(omplot, r, Ω, size=[2000, 1000]) =#
                #= plot!(xplot, r, getx(Ω, r)) =#
                #= plot!(gplot, r, gettorque(Ω, r)) =#
                #= plot(xplot, gplot, omplot, layout=@layout([a b; c])) =#
                #= png("tmp/tmp-$k") =#
            end
        end
    finally
        Sundials.CVodeFree([mem])
    end

    indices = [round(Int, (length(time))^(i//6)) for i in 0:6]
    Ωcsv = [r'; Ωseries[indices, :]]'
    writecsv("$simsdirectory/$name/omegaseries-$name.csv", Ωcsv)

    #= Ωplot = plot(xlabel="\$\\ln\\,r\$", ylabel="\$\\ln\\,\\Omega\$", title="Angular velocity") =#
    #= xplot = plot(xlabel="\$\\ln\\,r\$", ylabel="\$x\$", title="Local power law exponent", ylims=[-1.6, -0.9]) =#
    #= Σplot = plot(xlabel="\$\\ln\\,r\$", ylabel="\$\\ln\\,\\Sigma\$", title="Surface density", ylims=[-50, 20]) =#
    #= Ṁplot = plot(xlabel="\$\\ln\\,r\$", ylabel="\$\\dot M\$", title="Accretion rate") =#
    #= for i in indices =#
    #=     currentΩ = collect(Ωseries[i, :]) =#
    #=     plot!(Ωplot, r, currentΩ, label="\$\\tau = 10^{$(findfirst(indices,i)-1)}\$") =#
    #=     plot!(xplot, r, getx(currentΩ, r), label="\$\\tau = 10^{$(findfirst(indices,i)-1)}\$") =#
    #=     Σ = getΣ(currentΩ, r) =#
    #=     smoothcurve!(Σ, 1) =#
    #=     plot!(Σplot, r, Σ, label="\$\\tau = 10^{$(findfirst(indices,i)-1)}\$") =#
    #=     Ṁ = getṀ(currentΩ, r, viscositymodel) =#
    #=     smoothcurve!(Ṁ, 4) =#
    #=     plot!(Ṁplot, r, Ṁ, label="\$\\tau = 10^{$(findfirst(indices,i)-1)}\$") =#
    #= end =#

    #= plotsize = [16cm, 12cm] =#
    #= png(Ωplot, "$simsdirectory/$name/omega-$name") =#
    #= png(xplot, "$simsdirectory/$name/x-$name") =#
    #= png(Σplot, "$simsdirectory/$name/sigma-$name") =#
    #= png(Ṁplot, "$simsdirectory/$name/mdot-$name") =#
    nothing
end

function gettorque!(Ω, r, Σ, x, G)
    for i in 1:length(r)
        G[i] = exp(r[i]) * exp(Σ[i]) * exp(Ω[i]^2) * exp(r[i])
    end
end

function gettorque(Ω, r, Σ, x)
    G = Array(Float64, length(r))
    gettorque!(Ω, r, Σ, x, G)
    G
end

function gettorque(Ω, r)
    x = copy(r)
    derivative!(Ω, r, x)
    Ωl = exp(Ω)
    rl = exp(r)
    2π .* rl.^5 .* Ωl.^4 .* (2x .+ 3) .* x
end

# Sims ----------------------------------------

function simIllenseer(viscositymodel; rstring=:(collect(-6:0.15:25)),
                      timestring=:(collect(linspace(1, 1e6, 1e6))))
    r = eval(rstring)
    time = eval(timestring)
    ssdata = readdlm("../ssdisk/phys_t1.txt", '\t')
    ss_r = log(ssdata[:,1])
    ss_Ω = log(ssdata[:,2])
    Ω₀ = conformtoarray(ss_Ω, ss_r, r)
    #= Ω₀ = collect(r).*(-3/2) =#
    #= breakAt = round(Int16, findfirst(r, 0)) =#
    #= Ω₀[breakAt:end] = -1 .* (r[breakAt:end] .- r[breakAt]) .+ Ω₀[breakAt] =#
    #= smoothcurve!(Ω₀, 200) =#
    Ω₀string = "Startvalues taken from SSDisk"

    comment = """
    Simulation based on the simplest case in the Illenseer (2015) paper.
    This will be useful to verify the validity of the solver.
    Initial x is -3/2 until 0 on the r axis, then steps to -1.
    delta_x = 0.01, delta_t = 0.1, time from 1 till 1000.
    """

    runsimwithparams("illenseer-atzero-$viscositymodel", comment, G, Mstar,
                     viscositymodel, r, rstring, time,
                     timestring, Ω₀, Ω₀string)
end

# --------- sim2 - ring ----------

function simring(viscositymodel; rstring=:(collect(-6:0.15:25)),
                 timestring=:(collect(linspace(1, 1e6, 1e6))))
    r = eval(rstring)
    time = eval(timestring)
    amplitude = 1e0
    peakpos   = 2.0
    sdev      = 2.4
    Mstar     = 1.0
    Σ₀ = amplitude * gauss(r, peakpos, sdev)
    Ω₀ = ΩfromΣ(Σ₀, r, Mstar)
    Ω₀string = "Σ₀ = $amplitude * gauss(r, $peakpos, $sdev); Ω₀ = ΩfromΣ(Σ₀, r, $Mstar)"

    comment = """
    Σ₀ is a narrow bell curve, Ω₀ gets calculated based
    on that.
    """
    name = "ring-$viscositymodel"

    #= p = plot(r, Σ₀) =#
    #= png(p, "$simsdirectory/$name/sigma0") =#
    runsimwithparams("$name", comment, G, Mstar,
                     viscositymodel, r, rstring, time,
                     timestring, Ω₀, Ω₀string)
end

end

