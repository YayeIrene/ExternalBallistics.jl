#module OneD
#using  ..Utils

#export timeOfFlight!


function trajectory!(p::AbstractPenetrator;dt=1e-6)
    ma = machNumber(p.velocity_m)
    cd0 = p.coeffAero(ma)
    ρ0 = density(0.0)
    drag = -0.125 * ρ0 * pi * p.calibre^2 * cd0 * p.velocity_m^2
    acceleration = drag / p.mass
    p.velocity_m += acceleration * dt
    p.position += p.velocity_m * dt
end

function timeOfFlight!(p::AbstractPenetrator, target::AbstractTarget;dt=1e-6)
#function flight(p::Projectile, ρ::Float64, dt::Float64)
    t = 0.0
    old_position = p.position
    old_velocity = p.velocity_m
    ρ = target.position

    while ρ > p.position
        old_position = p.position
        old_velocity = p.velocity_m
        t += dt
        #@process trajectory_projectile(sim, p, dt)
        trajectory!(p)
        #println(p.position)
    end
    flight_time = t-dt + dt * (ρ - old_position) / (p.position - old_position)
    p.position = old_position
    p.velocity_m = old_velocity + (p.velocity_m - old_velocity) * (flight_time - t + dt) / dt

    return flight_time


end

#end
