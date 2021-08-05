#module Fragmentation
#using ..scenario, ..Utils, ..Types, ..ExternalBallistics
#using LinearAlgebra
#export frag, fragFly!
function createFragment(z::Zone, s::Array{FragShapes,1},proj::Projectile1D,α::Float64)
    number = Int64(rand(Normal(z.number)))
    rads = rand(Uniform(0,360), number)
    axs = rand(Uniform(z.angle), number)
    f = Array{Fragment}(undef,number)
    c_bar=[cos(α), -sin(α)]
    y_bar = [0.0, 1.0]
    e_bar = cross(c_bar,y_bar)/norm(cross(c_bar,y_bar))
    d_bar = cross(e_bar,c_bar)
    for i=1:number
        velocity = rand(Normal(z.velocity))
        mass = abs(rand(Normal(z.mass)))
        shape = rand(s)
        cd = rand(Normal(shape.cd))
        γ = rand(LogNormal(s.cf))
        Aₚ = γ*(mass/z.density)^(2/3)
        #vx = velocity*cosd(axs[i])*cosd(rads[i])
        #vy = velocity*sind(axs[i])
        #vz = velocity*cosd(axs[i])*sind(rads[i])
        f_velocity = proj.velocity*[1.0,0.0] + velocity*[cosd(axs[i]), sind(axs[i])] #in the fragments axis (2D)
        f_position =proj.calibre*[1/tand(axs[i]), 1] #in the fragments axis (2D)
        f_position = proj.position+proj.calibre/2*[cosd(axs[i]), cosd(rads[i]), sind(rads[i])]
        f[i] = Fragment(rads[i], axs[i], mass, z.density,cd, Aₚ, f_position,f_velocity,0.0)
    end
    return f
end



function frag(proj::Projectile2D, zones::Array{Zone,1}, shapes::Array{FragShapes,1},cf::ShapeFactor, αₑ_bar::Array{Float64,1})
    θ  = asin(-proj.velocity[3]/norm(proj.velocity))
    Ψ = atan(proj.velocity[2]/proj.velocity[1])
    T = [cos(θ)*cos(Ψ) -sin(Ψ) sin(θ)*cos(Ψ); cos(θ)*sin(Ψ) cos(Ψ) sin(θ)*sin(Ψ); -sin(θ) 0.0 cos(θ)]
    frags = []
    for z in zones
        f = createFragment(z,shapes,cf, proj, αₑ_bar)

    #=    for i in f

            v = i.velocity#[i.vx, i.vy, i.vz]
            vrot = rotation(v, T)
            i.velocity = vrot
            #i.vy = vrot[2]
            #i.vz = vrot[3]
        end
        =#
        push!(frags,f)
    end
    return frags
end

function fragFly!(fragments, proj::Projectile2D, w_bar::Array{Float64,1}, target::Target2D, g₀::Float64, ω_bar::Array{Float64,1}, coeffAero;R=6.356766*1e6, tspan=(0.0,1000.0))
    for i=1:length(fragments)
        for j=1:length(fragments[i])

            u0=[proj.position[1], proj.position[2], proj.position[3],fragments[i][j].velocity[1], fragments[i][j].velocity[2], fragments[i][j].velocity[3]]

            p= [w_bar, fragments[i][j].Aₚ, R, g₀,ω_bar, fragments[i][j].mass,target.position,coeffAero,fragments[i][j]]
            ExternalBallistics.Pmm.trajectoryPMM!(u0, tspan, p, fragments[i][j])
        end

    end

end


#end
