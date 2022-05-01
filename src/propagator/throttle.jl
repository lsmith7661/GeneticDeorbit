# Throttle functions

function throttle(state::Array{<:Real,1}, center::Real, width::Real, switching_criteria::Real)
    # switching criteria is an altitude of perigee (meters)
    r = state[1:3]
    v = state[4:6]
    rmag = norm(r)
    h = cross(r, v)
    hmag = norm(h)
    e = cross(v, h) / μ - r / rmag
    emag = norm(e)
    θ = acos(dot(e, r) / (emag * rmag))
    if dot(r, v) < 0
        θ = 2 * pi - θ
    end
    rp = hmag^2 / μ * (1 / 1 - emag)
    ra = hmag^2 / μ * (1 / 1 + emag)

    # Switching Criteria 
    if (rp - Rₑ) < switching_criteria
        center = pi - center # swith burn 180 degrees
    end

    if rp > (350000 + Rₑ) # Only thrust if perigee is greater than value
        if abs(center - θ) <= width
            return 1.0
        else
            return 0.0
        end
    else
        return 0.0
    end
end
throttle(state::Array{<:Real,1}) = throttle(state, pi, 45 * pi / 180, 300000)