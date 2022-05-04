# Evaluate the fitness of a trajectory

struct spacecraft
    switching_alt::Real # altiitude for which the spacecraft switches from apogee centered to perigee centered burns
    burn_duration::Real # an angular width which controls the duration of the burn (uniform around burn center)
end

struct trajectory
    time::AbstractArray{<:Real,1}       # time points of the trajectory
    states::AbstractArray{<:Real,2}     # rv states of the trajectory
    throttle::AbstractArray{<:Real,1}   # throttle state of the trajectory
    fitness::Real                       # the fitness of the given trajectory
    time_fitness::Real                  # the fitness of the time objective only
    thrust_fitness::Real                # the fitness of the thrust objective only
end

function eval_fitness(sc::spacecraft, weights::AbstractArray{<:Real,1})

    # Declare simulation initial epoch and state
    epoch0 = Epoch(2015, 4, 25, 12, 0, 0, 0.0)
    epochf = epoch0 + 180 * 24 * 60 * 60
    r0 = SVector{3,Float64}(2712241.37, -358426.60, -6318958.88) # m
    v0 = SVector{3,Float64}(3891.70, 6403.93, 1265.18) # m/s
    s0_eci = [r0[1], r0[2], r0[3], v0[1], v0[2], v0[3]]

    #EarthInertialState orbit propagagator
    orb = EarthInertialState(epoch0, s0_eci, dt=60.0,
        mass=500.0, area_drag=25.0, n_grav=0, m_grav=0,
        drag=true, srp=false,
        moon=false, sun=false,
        relativity=false,
        thrust_m_s_s=0.100 / 500,
        throttle=x -> throttle(x, pi, sc.burn_duration, sc.switching_alt)
    )
    SatelliteDynamicsLCS.reinit!(orb)

    # Propagate
    t, epc, eci = sim!(orb, epochf)
    T = Array{Float64}(undef, 1)
    for i = 1:size(eci)[2]
        state = eci[:, i]
        push!(T, throttle(state,pi, sc.burn_duration, sc.switching_alt))
    end
    deleteat!(T, 1)
    traj = trajectory(t, copy(eci'), T, 0, 0, 0)

    # Calculate rmag for the states
    n = size(traj.states)[1]
    nskip = 500 # we dont need to check every state
    reentry = n # last entry if it s/c does not reenter
    for i = 1:nskip:n
        state = traj.states[i, :]
        r = state[1:3]
        v = state[4:6]
        rmag = norm(r)
        h = cross(r, v)
        hmag = norm(h)
        e = cross(v, h) / μ - r / rmag
        emag = norm(e)
        ra = hmag^2 / μ * (1 / 1 + emag)
        if ra < death_alt + Rₑ
            # state where apogee broke loop
            reentry = i
            break
        end
    end

    # Reentry time
    t_total = traj.time[reentry] - traj.time[1]

    # Total thrust
    T_total = sum(traj.throttle[1:reentry])

    # Weights
    w_t = weights[1] # weight on total time 
    w_T = weights[2] # weight on total thrust
    fitness1 = (w_t * t_total)
    fitness2 = (w_T * T_total)
    fitness = 1 * (fitness1 + fitness2)

    return trajectory(t, copy(eci'), T, fitness, fitness1, fitness2)
end
eval_fitness(sc::spacecraft) = eval_fitness(sc, [0.5, 0.5])

function eval_fitness_val(x::AbstractArray{<:Real,1}, weights::AbstractArray{<:Real,1})
    sc = spacecraft(x[1], x[2])
    # returns just the fitness parameter instead of the full trajectory
    fullout = eval_fitness(sc, weights)
    return fullout.fitness
end
eval_fitness_val(sc::spacecraft) = eval_fitness_val(sc, [0.5, 0.5])
