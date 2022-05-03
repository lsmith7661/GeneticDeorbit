# Some simple functions to plot orbits

nskip = 500

function plotRV(traj::trajectory)
    x = []
    y = []
    z = []
    states = traj.states
    for i = 1:1:size(traj.states)[1]
        push!(x , states[i, 1])
        push!(y , states[i, 2])
        push!(z , states[i, 3])
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
        cosθ = dot(e, r) / (emag * rmag)
        if cosθ > 1 # catch precision errors
            θi = 0
        elseif cosθ < -1
            θi = 0
        else
            θi = acos(cosθ)  
        end 
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
    p2 = plot(traj.time, cumsum(traj.throttle), label="Cumulative Throttle", xlabel="time", lw=3, title="Cumulative Throttle")
    p3 = scatter([traj.thrust_fitness], [traj.time_fitness], xlims = (0,1), ylims = (0,1), label="Optimal Point", xlabel="Thrust Objective", ylabel="Time Objective", lw=3, title="Objective Function")
    p4 = plot(t, ecc, label="Eccentricity", xlabel="time", lw=3, title="Eccentricity")
    plot(p1, p2, p3, p4, layout=(2, 2), legend=false)

end

function plotOptimal(traj::trajectory, genehist1::Vector{Any}, genehist2::Vector{Any})
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
        if rai < death_alt + Rₑ
            break
        end

        push!(t, ti)
        push!(rp, rpi - Rₑ)
        push!(ra, rai - Rₑ)
        push!(ecc, emag)
    end

    p1 = plot(t, rp ./ 1e3, label="Perigee", lw=3, title="Altitude", xlabel="Time", ylabel="Altitude (km)")
    plot!(p1, t, ra ./ 1e3, label="Apogee", lw=3, title="Altitude")
    p2 = scatter([traj.thrust_fitness], [traj.time_fitness], xlims = (0,1), ylims = (0,1), label="Optimal Point", xlabel="Thrust Objective", ylabel="Time Objective", lw=3, title="Objective Function")
    
    p3 = scatter(ones(size(genehist1[1])), genehist1[1]./1e3, label="Stopping Altitude", xlabel="Interation",  ylabel="Altitude (km)", lw=3, title="Perigee Stop Alt")
    for i = 2:length(genehist1)
        scatter!(p3, i.*ones(size(genehist1[i])), genehist1[i]./1e3)
    end
    
    p4 = scatter(ones(size(genehist2[1])), genehist2[1]./(pi), xlabel="Interation",  ylabel="Fraction of Orbit (%)", lw=3, title="Burn Fraction")
    for i = 2:length(genehist2)
        scatter!(p4, i.*ones(size(genehist2[i])), genehist2[i]./(pi))
    end    
    
    plot(p1, p2, p3, p4, layout=(2, 2), legend=false)

end


function plotCircularDeorbit()
    # Make Plots for a circular deorbit strategy
    stopping_alt = 200e3 # doesnt matter because always thrusting
    burn_dur = pi  # always be burnin
    sc = spacecraft(stopping_alt, burn_dur)
    traj = eval_fitness(sc)

    # Calculate elements to plot
    t = []
    rpvec = []
    ravec = []
    θvec = []
    n = size(traj.states)[1]
    nskip = 1  
    reentry = n # last entry if it s/c does not reenter
    for i = 1:nskip:n
        state = traj.states[i, :]
        ti = (traj.time[i] - traj.time[1])/60/60/24
        r = state[1:3]
        v = state[4:6]
        rmag = norm(r)
        h = cross(r, v)
        hmag = norm(h)
        e = cross(v, h) / μ - r / rmag
        emag = norm(e)
        ra = hmag^2 / μ * (1 / 1 + emag)
        rp = hmag^2 / μ * (1 / 1 - emag)
        cosθ = dot(e, r) / (emag * rmag)
        if cosθ > 1 # catch precision errors
            θ = 0
        elseif cosθ < -1
            θ = 0
        else
            θ = acos(cosθ)
        end
        if dot(r, v) < 0
            θ = 2 * pi - θ
        end
        
        if ra < death_alt + Rₑ
            # state where apogee broke loop
            reentry = i
            break
        end

        push!(t,ti)
        push!(rpvec, rp - Rₑ)
        push!(ravec, ra - Rₑ)
        push!(θvec,  θ)
    end

    p1 = plot(t, rpvec ./ 1e3, label="Perigee", xlims = (0,60), lw=3, title="Altitude", xlabel="Time (days)", ylabel="Altitude (km)")
    plot!(p1, t, ravec ./ 1e3, label="Apogee", xlims = (0,60), lw=3, title="Altitude")
    
    ii = 500
    perm = sortperm(θvec[1:ii])
    p2 = plot(θvec[perm].*180/pi, traj.throttle[perm], lw=3, ylim=(0,1),xlabel="True Anomaly (deg)", ylabel="Throttle", title="Throttle Profile", legend = false)

    plot(p1, p2, layout=(1, 2), legend=false)

end

function plotPerigeeDeorbit()
    # Make Plots for a perigee lowering aeroassist deorbit strategy
    stopping_alt = 75e3 # stop low for dramatic atmosphere affect
    burn_dur = pi/8  # burn around apogee
    sc = spacecraft(stopping_alt, burn_dur)
    traj = eval_fitness(sc)

    # Calculate elements to plot
    t = []
    rpvec = []
    ravec = []
    θvec = []
    n = size(traj.states)[1]
    nskip = 1  
    reentry = n # last entry if it s/c does not reenter
    for i = 1:nskip:n
        state = traj.states[i, :]
        ti = (traj.time[i] - traj.time[1])/60/60/24
        r = state[1:3]
        v = state[4:6]
        rmag = norm(r)
        h = cross(r, v)
        hmag = norm(h)
        e = cross(v, h) / μ - r / rmag
        emag = norm(e)
        ra = hmag^2 / μ * (1 / 1 + emag)
        rp = hmag^2 / μ * (1 / 1 - emag)
        cosθ = dot(e, r) / (emag * rmag)
        if cosθ > 1 # catch precision errors
            θ = 0
        elseif cosθ < -1
            θ = 0
        else
            θ = acos(cosθ)
        end
        if dot(r, v) < 0
            θ = 2 * pi - θ
        end
        
        if ra < death_alt + Rₑ
            # state where apogee broke loop
            reentry = i
            break
        end

        push!(t,ti)
        push!(rpvec, rp - Rₑ)
        push!(ravec, ra - Rₑ)
        push!(θvec,  θ)
    end

    p1 = plot(t, rpvec ./ 1e3, label="Perigee", xlims = (0,60), lw=3, title="Altitude", xlabel="Time (days)", ylabel="Altitude (km)")
    plot!(p1, t, ravec ./ 1e3, label="Apogee", xlims = (0,60),  lw=3, title="Altitude")
    
    ii = 500
    perm = sortperm(θvec[1:ii])
    p2 = plot(θvec[perm].*180/pi, traj.throttle[perm], lw=3, ylim=(0,1),xlabel="True Anomaly (deg)", ylabel="Throttle", title="Throttle Profile", legend = false)

    plot(p1, p2, layout=(1, 2), legend=false)

end

function plotCircularAndPerigee()
    # Make Plots for a circular deorbit strategy
    stopping_alt = 200e3 # doesnt matter because always thrusting
    burn_dur = pi  # always be burnin
    sc = spacecraft(stopping_alt, burn_dur)
    traj = eval_fitness(sc)

    # Calculate elements to plot
    t = []
    rpvec = []
    ravec = []
    θvec = []
    n = size(traj.states)[1]
    nskip = 1  
    reentry = n # last entry if it s/c does not reenter
    for i = 1:nskip:n
        state = traj.states[i, :]
        ti = (traj.time[i] - traj.time[1])/60/60/24
        r = state[1:3]
        v = state[4:6]
        rmag = norm(r)
        h = cross(r, v)
        hmag = norm(h)
        e = cross(v, h) / μ - r / rmag
        emag = norm(e)
        ra = hmag^2 / μ * (1 / 1 + emag)
        rp = hmag^2 / μ * (1 / 1 - emag)
        cosθ = dot(e, r) / (emag * rmag)
        if cosθ > 1 # catch precision errors
            θ = 0
        elseif cosθ < -1
            θ = 0
        else
            θ = acos(cosθ)
        end
        if dot(r, v) < 0
            θ = 2 * pi - θ
        end
        
        if ra < death_alt + Rₑ
            # state where apogee broke loop
            reentry = i
            break
        end

        push!(t,ti)
        push!(rpvec, rp - Rₑ)
        push!(ravec, ra - Rₑ)
        push!(θvec,  θ)
    end

    p1 = plot(t, rpvec ./ 1e3, label="Perigee", lw=3, title="Altitude", xlabel="Time (days)", ylabel="Altitude (km)", linecolor=:blue)
    plot!(p1, t, ravec ./ 1e3, label="Apogee", lw=3,linecolor=:blue)
    
    ii = 500
    perm = sortperm(θvec[1:ii])
    p2 = plot(θvec[perm].*180/pi, traj.throttle[perm], lw=3, ylim=(0,1),xlabel="True Anomaly (deg)", ylabel="Throttle", title="Throttle Profile", legend = false, linecolor=:blue)

    # Make Plots for a perigee lowering aeroassist deorbit strategy
    stopping_alt = 75e3 # stop low for dramatic atmosphere affect
    burn_dur = pi/8  # burn around apogee
    sc = spacecraft(stopping_alt, burn_dur)
    traj = eval_fitness(sc)

    # Calculate elements to plot
    t = []
    rpvec = []
    ravec = []
    θvec = []
    n = size(traj.states)[1]
    nskip = 1  
    reentry = n # last entry if it s/c does not reenter
    for i = 1:nskip:n
        state = traj.states[i, :]
        ti = (traj.time[i] - traj.time[1])/60/60/24
        r = state[1:3]
        v = state[4:6]
        rmag = norm(r)
        h = cross(r, v)
        hmag = norm(h)
        e = cross(v, h) / μ - r / rmag
        emag = norm(e)
        ra = hmag^2 / μ * (1 / 1 + emag)
        rp = hmag^2 / μ * (1 / 1 - emag)
        cosθ = dot(e, r) / (emag * rmag)
        if cosθ > 1 # catch precision errors
            θ = 0
        elseif cosθ < -1
            θ = 0
        else
            θ = acos(cosθ)
        end
        if dot(r, v) < 0
            θ = 2 * pi - θ
        end
        
        if ra < death_alt + Rₑ
            # state where apogee broke loop
            reentry = i
            break
        end

        push!(t,ti)
        push!(rpvec, rp - Rₑ)
        push!(ravec, ra - Rₑ)
        push!(θvec,  θ)
    end

    plot!(p1, t, rpvec ./ 1e3, label="Perigee", lw=3, linecolor=:red)
    plot!(p1, t, ravec ./ 1e3, label="Apogee", lw=3, linecolor=:red)
    
    ii = 500
    perm = sortperm(θvec[1:ii])
    plot!(p2, θvec[perm].*180/pi, traj.throttle[perm], lw=3, ylim=(0,1), linecolor=:red)

    plot(p1, p2, layout=(1, 2), legend=false)
end