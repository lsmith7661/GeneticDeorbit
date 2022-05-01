# Some simple functions to plot orbits

nskip = 500

function plotRV(traj::trajectory)
    states = traj.states
    for i = 1:nskip:size(traj.states)[1]
        x = states[i, 1]
        y = states[i, 2]
        z = states[i, 3]
    end
    plot3d(x, y, z)
end

function plotCOES(traj::trajectory)
    t = []
    rp = []
    ra = []
    ecc = []
    θ = []
    nskip = 500
    for i = 1:nskip:size(traj.states)[1]
        state = traj.states[i, :]
        ti = traj.time[i]
        r = state[1:3]
        v = state[4:6]
        rmag = norm(r)
        h = cross(r, v)
        hmag = norm(h)
        e = cross(v, h) / μ - r / rmag
        emag = norm(e)
        rpi = hmag^2 / μ * (1 / 1 - emag)
        rai = hmag^2 / μ * (1 / 1 + emag)
        θi = acos(dot(e, r) / (emag * rmag))
        if dot(r, v) < 0
            θi = 2 * pi - θi
        end

        if rai < death_alt + Rₑ
            break
        end

        push!(t, ti)
        push!(θ, θi)
        push!(rp, rpi - Rₑ)
        push!(ra, rai - Rₑ)
        push!(ecc, emag)
    end
    p1 = plot(t, rp ./ 1e3, label="Perigee", xlabel="time", lw=3, title="Altitude")
    plot!(p1, t, ra ./ 1e3, label="Apogee", xlabel="time", lw=3, title="Altitude")
    p2 = plot(t, θ, label="True Anomaly", xlabel="time", lw=3, title="True Anomaly")
    p3 = plot(t, ecc, label="Eccentricity", xlabel="time", lw=3, title="Eccentricity")
    p4 = plot(traj.time, cumsum(traj.throttle), label="Cumulative Throttle", xlabel="time", lw=3, title="Cumulative Throttle")
    plot(p1, p2, p3, p4, layout=(2, 2), legend=false)

end