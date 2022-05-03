# Calculate the derivatives of the orbital paramaters for propagation
# Gauss Variation of Parameters: Vallado 4th ed, pg 628

#=

struct rv_eci
    # cartesian state
    r_m::SVector{3,Float64}   # position vector, meters
    v_m_s::SVector{3,Float64} # velocity vector, meters per second

end

struct drv_eci
    # time derivative of state vector
    drdt_m_s::SVector{3,Float64}   # velocity vector, meters per second
    dvdt_m_s_s::SVector{3,Float64} # acceleration vector, meters per second per second
end

function exponential_atmosphere_density(h::Float64)
    # Calculate atmospheric density given an altitude, h, in meters
    # Table of values from Vallado pg 566
    #
    if h > 600e3
        h₀ = 600e3       # m
        ρ₀ = 1.454e-13   # kg/m^3
        H = 71.835e3     # m
    elseif h > 500e3
        h₀ = 500e3       # m
        ρ₀ = 6.967e-13   # kg/m^3
        H = 63.822e3     # m
    elseif h > 450e3
        h₀ = 450e3       # m
        ρ₀ = 1.585e-12   # kg/m^3
        H = 60.828e3     # m
    elseif h > 400e3
        h₀ = 400e3       # m
        ρ₀ = 3.725e-12   # kg/m^3
        H = 58.515e3     # m
    elseif h > 350e3
        h₀ = 350e3       # m
        ρ₀ = 9.518e-12   # kg/m^3
        H = 53.298e3     # m
    elseif h > 300e3
        h₀ = 300e3       # m
        ρ₀ = 2.418e-11   # kg/m^3
        H = 53.628e3     # m
    elseif h > 250e3
        h₀ = 250e3       # m
        ρ₀ = 7.248e-11   # kg/m^3
        H = 45.546e3     # m
    elseif h > 200e3
        h₀ = 200e3       # m
        ρ₀ = 2.789e-10   # kg/m^3
        H = 37.105e3     # m
    elseif h > 180e3
        h₀ = 180e3       # m
        ρ₀ = 5.464e-10   # kg/m^3
        H = 29.740e3     # m
    elseif h > 150e3
        h₀ = 150e3       # m
        ρ₀ = 2.070e-9    # kg/m^3
        H = 22.523e3     # m
    elseif h > 100e3
        h₀ = 100e3       # m
        ρ₀ = 5.297e-7    # kg/m^3
        H = 5.877e3      # m
    else
        h₀ = 0           # m
        ρ₀ = 1.225       # kg/m^3
        H = 7.249e3      # m
    end

    return ρ₀ * exp(-(h - h₀) / H)
    #

    #= Or curve fit so its a bit easier on the code
    a1 = 2.498e+33
    b1 = -1.119e+06
    c1 = 1.278e+05
    return a1 * exp(-((h - b1) / c1)^2) =#
    a = 4.628e+40
    b = -9.457
    return a * x^b

end

function accel_exp_drag(state::rv_eci, area2mass::Float64)
    # area to mass is m^2/kg
    r_m = state.r_m
    r_mag = norm(r_m)
    v_m_s = state.v_m_s
    v_mag = norm(v_m_s)
    v_hat = v_m_s / v_mag
    Cd = 2.2
    rho = exponential_atmosphere_density(r_mag - Rₑ)
    return -0.5 * Cd * area2mass * rho * v_mag^2 * v_hat
end
accel_exp_drag(state::rv_eci) = accel_exp_drag(state, 25.0 / 500.0)


function cowells(state::rv_eci, thrust_m_s_s::Float64)
    # Derivative of the state

    # Extract variables
    r_m = state.r_m
    v_m_s = state.v_m_s

    # Helper variables
    r_mag = norm(r_m)
    v_mag = norm(v_m_s)

    # perturbations
    thrust = thrust_m_s_s * -v_m_s / v_mag
    drag = accel_exp_drag(state)

    # derivative
    drdt_m_s = v_m_s
    dvdt_m_s_s = (-μ / r_mag^3) * r_m + thrust + drag
    return drv_eci(drdt_m_s, dvdt_m_s_s)
end
function cowells(state::SVector{6,Float64}, thrust_m_s_s::Float64, t::Float64)
    r = SVector{3,Float64}(state[1], state[2], state[3])
    v = SVector{3,Float64}(state[4], state[5], state[6])
    s = rv_eci(r, v)
    dsdt = cowells(s, thrust_m_s_s)
    return SVector{6,Float64}(dsdt.drdt_m_s[1], dsdt.drdt_m_s[2], dsdt.drdt_m_s[3], dsdt.dvdt_m_s_s[1], dsdt.dvdt_m_s_s[2], dsdt.dvdt_m_s_s[3])
end


function cowellsprop(state::rv_eci, thrust_m_s_s::Float64, tspan::Tuple{Float64,Float64})
    # propagate a state from tspan[0] to all instants of tspan
    # cowells method + ode45 

    # extract 
    r_m = state.r_m
    v_m_s = state.v_m_s

    # state
    s0 = SVector{6,Float64}(r_m[1], r_m[2], r_m[3], v_m_s[1], v_m_s[2], v_m_s[3])

    # solve ODE
    prob = ODEProblem(cowells, s0, tspan, thrust_m_s_s)
    sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6)

    return sol
end

# Some initial variables to play with
r = SVector{3,Float64}(1600000.0, 5310000.0, 3800000.0)
v = SVector{3,Float64}(-7350.0, 460.0, 2470.0)
initialstate_rv_eci = rv_eci(r, v)
thrust0_m_s_s = 0.0
thrust1_m_s_s = 0.05 / 500

=#