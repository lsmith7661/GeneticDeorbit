# Calculate the derivatives of the orbital paramaters for propagation
# Gauss Variation of Parameters: Vallado 4th ed, pg 628

struct COEs
    # Classical orbital elements
    a::Float64 # semi major axis, meters
    e::Float64 # eccentricity
    i::Float64 # inclination, rad
    Ω::Float64 # right ascenscion ascending node, rad
    ω::Float64 # argument of perigee, rad
    θ::Float64 # true anomaly, rad
end

struct dCOEs
    # time derivative of classical orbital elements
    dadt::Float64 # semi major axis, meters
    dedt::Float64 # eccentricity
    didt::Float64 # inclination, rad
    dΩdt::Float64 # right ascenscion ascending node, rad
    dωdt::Float64 # argument of perigee, rad
    dθdt::Float64 # true anomaly, rad
end

struct RSW
    # Perturbing force in gauss VOP (specific force, aka acceleration)
    # Force in the RSW frame as defined by vallado ch 3
    R::Float64 # m/s^2 acceleration in the R direction
    S::Float64 # m/s^2 acceleration in the S direction
    W::Float64 # m/s^2 acceleration in the F direction
end

function gaussVOP(state::COEs, force::RSW)
    # Derivative of the state as per Gauss VOP
    # Extract variables
    a = state.a
    e = state.e
    i = state.i
    Ω = state.Ω
    ω = state.ω
    θ = state.θ
    Fr = force.R
    Fs = force.S
    Fw = force.W

    # Helper variables
    n = sqrt(μ / a^3)
    p = a * (1 - e^2)
    r = p / (1 + e * cos(θ))
    u = ω + θ
    h = sqrt(μ * p)

    # Gauss VOP Equations
    dadt = 2 / (n * sqrt(1 - e^2)) * (e * sin(θ) * Fr + (p / r) * Fs)
    dedt = (sqrt(1 - e^2) / (n * a)) * (sin(θ) * Fr + ((cos(θ) + (e + cos(θ)) / (1 + e * cos(θ))) * Fs))
    didt = ((r * cos(u)) / (n * a^2 * sqrt(1 - e^2))) * Fw
    dΩdt = ((r * sin(u)) / (n * a^2 * sqrt(1 - e^2) * sin(i))) * Fw
    dωdt = ((sqrt(1 - e^2)) / (n * a * e)) * (-cos(θ) * Fr + sin(θ) * (1 + r / p) * Fs) - ((r * cot(i) * sin(u)) / h) * Fw
    dθdt = (h / r^2) + (1 / (e * h)) * (p * cos(θ)) * Fr - (p + r) * sin(θ) * Fs

    return dCOEs(dadt, dedt, didt, dΩdt, dωdt, dθdt)
end
function gaussVOP(state::AbstractArray{Float64}, force::AbstractArray{Float64}, t::Float64)
    s = COEs(state[1], state[2], state[3], state[4], state[5], state[6])
    f = RSW(force[1], force[2], force[3])
    dsdt = gaussVOP(s, f)
    return [dsdt.dadt, dsdt.dedt, dsdt.didt, dsdt.dΩdt, dsdt.dωdt, dsdt.dθdt]
end


function gaussVOPprop(state::COEs, force::RSW, tspan::Tuple{Float64,Float64})
    # propagate a state from tspan[0] to all instants of tspan
    # method gauss vop + ode45 

    # extract 
    a = state.a
    e = state.e
    i = state.i
    Ω = state.Ω
    ω = state.ω
    θ = state.θ
    Fr = force.R
    Fs = force.S
    Fw = force.W

    s0 = [a, e, i, Ω, ω, θ]
    f0 = [Fr, Fs, Fw]

    prob = ODEProblem(gaussVOP, s0, tspan, f0)
    sol = solve(prob, Tsit5())

    return sol
end

# Some initial variables to play with
initialstateCOES = COEs(6378000 + 600000, 0.003, 60 * 3.14 / 180, 0, 0, 0)
force0 = RSW(0, 0, 0)
force1 = RSW(0, 0.05 / 500, 0)
