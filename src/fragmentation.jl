#module Fragmentation
#using ..scenario, ..Utils, ..Types, ..ExternalBallistics
#using LinearAlgebra
#export frag, fragFly!
function createFragment(z::Zone, s::Array{FragShapes,1},proj::AbstractPenetrator)
    number = Int64(rand(Normal(z.number[1],z.number[2])))
    rads = rand(Uniform(0,360), number)
    axs = rand(Uniform(z.angle[1],z.angle[2]), number)
    f = Array{Fragment}(undef,number)
    c_bar = proj.velocity/norm(proj.velocity)*cos(norm(proj.αₑ))+proj.αₑ
    #c_bar=[cos(α), -sin(α)]
    y_bar = [0.0, 1.0,0.0]
    e_bar = cross(c_bar,y_bar)/norm(cross(c_bar,y_bar))
    d_bar = cross(e_bar,c_bar)
    for i=1:number
        velocity = rand(Normal(z.velocity[1],z.velocity[2]))
        mass = abs(rand(Normal(z.mass[1],z.mass[1])))
        shape = rand(s)
        cd = rand(Normal(shape.cd[1], shape.cd[2]))
        γ = rand(LogNormal(shape.cf[1], shape.cf[2]))
        Aₚ = γ*(mass/z.density)^(2/3)
        #vx = velocity*cosd(axs[i])*cosd(rads[i])
        #vy = velocity*sind(axs[i])
        #vz = velocity*cosd(axs[i])*sind(rads[i])
        #f_velocity = proj.velocity*[1.0,0.0] + velocity*[cosd(axs[i]), sind(axs[i])] #in the fragments axis (2D)
        rₛ_bar = proj.calibre/2*(e_bar*sind(rads[i])+d_bar*cosd(rads[i]))
        uₛₙ_bar =cross(proj.spin*c_bar,rₛ_bar)
        uₑⱼ_bar = velocity * (e_bar*sind(rads[i])+d_bar*cosd(rads[i]))
        f_velocity = proj.velocity + uₑⱼ_bar + uₛₙ_bar
        #f_position =proj.calibre*[1/tand(axs[i]), 1] #in the fragments axis (2D)
        f_position = proj.position + rₛ_bar
        #f_position = proj.position+proj.calibre/2*[cosd(axs[i]), cosd(rads[i]), sind(rads[i])]
        f[i] = Fragment(rads[i], axs[i], mass, z.density,cd, Aₚ, f_position,f_velocity,0.0)
    end
    return f
end

function createFragment(z::ZoneZdata, s::Array{FragShapes,1},proj::AbstractPenetrator)
    number = sum(z.number)
    rads = rand(Uniform(0,360), trunc(Int,number))
    axs = rand(Uniform(z.angle[1],z.angle[3]), trunc(Int,number))
    #f = Array{Fragment}(undef,number)
    f=[]
    c_bar = proj.velocity/norm(proj.velocity)*cos(norm(proj.αₑ))+proj.αₑ
    y_bar = [0.0, 1.0,0.0]
    e_bar = cross(c_bar,y_bar)/norm(cross(c_bar,y_bar))
    d_bar = cross(e_bar,c_bar)
    velEject=LinearInterpolation(z.angle,z.velocity)
    for i=1:length(z.number)
        for j=1:z.number[i]
            #fi =Cuas.createFragment(z.mass[i])
            fi = generateFragment(z.mass[i])
           # println(fi.mass)
            push!(f,fi)
        end
    end 
    for i=1:length(f)#trunc(Int,number)
        velocity = velEject(axs[i])
        #mass = abs(rand(Normal(z.mass[1],z.mass[1])))
        shape = rand(s)
        cd = rand(Normal(shape.cd[1], shape.cd[2]))
        γ = rand(LogNormal(shape.cf[1], shape.cf[2]))
        #println(f[i].mass)
        Aₚ = γ*(f[i].mass/z.density)^(2/3)
        #vx = velocity*cosd(axs[i])*cosd(rads[i])
        #vy = velocity*sind(axs[i])
        #vz = velocity*cosd(axs[i])*sind(rads[i])
        #f_velocity = proj.velocity*[1.0,0.0] + velocity*[cosd(axs[i]), sind(axs[i])] #in the fragments axis (2D)
        rₛ_bar = proj.calibre/2*(e_bar*sind(rads[i])+d_bar*cosd(rads[i]))
        uₛₙ_bar =cross(proj.spin*c_bar,rₛ_bar)
        uₑⱼ_bar = velocity * (e_bar*sind(rads[i])+d_bar*cosd(rads[i]))
        f_velocity = proj.velocity + uₑⱼ_bar + uₛₙ_bar
        #f_position =proj.calibre*[1/tand(axs[i]), 1] #in the fragments axis (2D)
        f_position = proj.position + rₛ_bar
        #f_position = proj.position+proj.calibre/2*[cosd(axs[i]), cosd(rads[i]), sind(rads[i])]
        #f[i] = Fragment(rads[i], axs[i], mass, z.density,cd, Aₚ, f_position,f_velocity,0.0)
        f[i].rad = rads[i]
        f[i].ax = axs[i]
        f[i].density = z.density 
        f[i].cd_sub = cd
        f[i].Aₚ = Aₚ 
        f[i].position = f_position
        f[i].velocity = f_velocity
        f[i].tof = 0.0
    end
    return f
end


"""
    frag(proj, zones, shapes)
Returns the a vector of fragments:
* proj : the fragmenting projectile with projectile velocity and position
* zones : contains information on the static fragmentation
* shapes: contains information on fragments shapes
"""
function frag(proj::AbstractPenetrator, zones::ZoneList, shapes::Array{FragShapes,1})
    #θ  = asin(-proj.velocity[3]/norm(proj.velocity))
    #Ψ = atan(proj.velocity[2]/proj.velocity[1])
    #T = [cos(θ)*cos(Ψ) -sin(Ψ) sin(θ)*cos(Ψ); cos(θ)*sin(Ψ) cos(Ψ) sin(θ)*sin(Ψ); -sin(θ) 0.0 cos(θ)]
    frags = []
    for z in zones.list
        f = createFragment(z,shapes, proj)

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


function fragFly!(fragments, proj::AbstractPenetrator, w_bar::Array{Float64,1}, target::AbstractTarget, g₀::Float64, ω_bar::Array{Float64,1}, coeffAero;R=6.356766*1e6, tspan=(0.0,1000.0))
    for i=1:length(fragments)
        for j=1:length(fragments[i])

            u0=[proj.position[1], proj.position[2], proj.position[3],fragments[i][j].velocity[1], fragments[i][j].velocity[2], fragments[i][j].velocity[3]]

            p= [w_bar, fragments[i][j].Aₚ, R, g₀,ω_bar, fragments[i][j].mass,target.position,coeffAero,fragments[i][j]]
            ExternalBallistics.Pmm.trajectoryPMM!(u0, tspan, p, fragments[i][j])
        end

    end

end


#end
