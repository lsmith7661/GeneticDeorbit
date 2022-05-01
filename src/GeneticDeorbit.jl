module GeneticDeorbit

const μ = 3.986004418e14 # m^3/s^2
const Rₑ = 6378137 # m
const death_alt = 200000 # m
export
    μ,
    Rₑ,
    death_alt

using DifferentialEquations: ODEProblem, solve, Tsit5
include("propagator/gaussVOP.jl")
export
    COEs,
    dCOEs,
    RSW,
    gaussVOP,
    gaussVOPprop,
    initialstateCOES,
    force0,
    force1

using DifferentialEquations: ODEProblem, solve, Tsit5
using StaticArrays
using LinearAlgebra: norm, cross, dot, acos
include("propagator/cowellsprop.jl")
export
    rv_eci,
    drv_eci,
    cowells,
    cowellsprop,
    initialstate_rv_eci,
    thrust0_m_s_s,
    thrust1_m_s_s,
    exponential_atmosphere_density

using LinearAlgebra: norm, cross, dot, acos, cumsum
include("propagator/throttle.jl")
export
    throttle

using SatelliteDynamicsLCS
include("fitness/fitness.jl")
export
    spacecraft,
    trajectory,
    eval_fitness,
    eval_fitness_val

using Plots: plot3d, plot, plot!
include("plotting/orbitplots.jl")
export
    plotRV,
    plotCOES

end # module
