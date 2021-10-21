

#using aeroCoeff
#include("aerodynamic_coeff.jl")
#export trajectoryMPMM!, QEfinderMPMM!, iniCond

function mpmm(du,u,p,t)
    proj,w_bar, R, g₀,ω_bar,target,dist,tp,aero = p

    Iₓ=proj.Ix
    d =proj.calibre
    m=proj.mass

    vx= u[4]
    vy=u[5]
    vz=u[6]
    X₁ = u[1]
    X₂ = u[2]
    X₃ = u[3]
    pₛ = u[7]
   u_bar =[u[4], u[5], u[6]]
   v_bar= u_bar - w_bar
    v = norm(v_bar)
    du_bar = [du[4], du[5],du[6]]
    #println(X₂)
    ρ =  density(X₂)
    #println("ρ", " ", ρ)
    #α = atan(vz/vx)
    #β = asin(vy/v)
    α = atan(vy/vx)
    β = asin(vz/v)

    #mach = Utils.machNumber(v,X₂)
    if tp==nothing
        tp = temperature(X₂)
    end

    mach = machNumber(v,tp)
    CMα₀ = inter(aero.mach, aero.cma)#cma
    #CMα₂ = CM_α2_inter(mach)
    CMα = CMα₀(mach) #+ CMα₂*(sin(αₜ))^2
    #println(du_bar)

   αₑ_bar = yawOfRepose(v_bar, du_bar,pₛ, t, Iₓ, ρ,d,CMα,v)
   #println("αₑ", " ", αₑ_bar)
    αₑ = norm(αₑ_bar)
    #αₜ = asin(αₑ)
    αₜ = αₑ
    #println("1", " ", αₜ)
    #αₜ = asin(sqrt((sin(α)*cos(β))^2+(sin(β)^2)))
    #println("2", " ", αₜ)
    #println("αₜ", " ", αₜ )
    #println("αₜ2", " ", sin(αₑ))

    CX0 = inter(aero.mach,aero.cx)
    CX2 = inter(aero.mach,aero.cx2)
    CDα⁴ = 0.0
    CNa = inter(aero.mach,aero.cna)
    CLα³ = 0.0
    CLα⁵ = 0.0
    cnpα_matrix = matrix(aero.cnpα_1,aero.cnpα_2, aero.cnpα_5,aero.cnpα_10)
    cnpα = inter_matrix(aero.mach, cnpα_matrix, [1.0,2.0,5.0,10.0])

    #Cmag = CN_pα_inter(mach,(rad2deg(αₜ))^2)#cnpa
    Cmag = cnpα(mach,rad2deg(αₜ))
    clp = inter(aero.mach,aero.clp)
    Cₛₚᵢₙ = clp(mach)#clp

    CD = drag(CX0,CX2,αₜ,mach,CNa)
    #CD = drag(CX0,CX2,αₜ,α,β,mach,CNa)
    #println("α", " ",αₜ)
    #println("mach", " ", mach)
    #println("CD", " ", CD)
    CL = lift(CX0,CX2,αₜ,mach,CNa)
    #println("CL", " ", CL)
    FDx = - pi*ρ*d^2/(8)*(CD)*v*v_bar[1] #drag
    FDy = - pi*ρ*d^2/(8)*(CD)*v*v_bar[2] #drag
    FDz = - pi*ρ*d^2/(8)*(CD)*v*v_bar[3] #drag
    FLx = pi*ρ*d^2/(8)*(CL)*v^2*αₑ_bar[1] #lift
    FLy = pi*ρ*d^2/(8)*(CL)*v^2*αₑ_bar[2] #lift
    FLz = pi*ρ*d^2/(8)*(CL)*v^2*αₑ_bar[3] #lift
    FMx = -pi*ρ*d^3*pₛ*Cmag/8*cross(αₑ_bar,v_bar)[1]
    FMy = -pi*ρ*d^3*pₛ*Cmag/8*cross(αₑ_bar,v_bar)[2]
    FMz = -pi*ρ*d^3*pₛ*Cmag/8*cross(αₑ_bar,v_bar)[3]
    #println("cmag", " ", Cmag)
    FGx = - g₀*X₁/R*m
    FGy = -g₀*(1-2*X₂/R)*m
    FGz = - g₀*X₃/R*m
    #FGz = 0.0
    FCx = -2*cross(ω_bar,u_bar)[1]*m
    FCy = -2*cross(ω_bar,u_bar)[2]*m
    FCz = -2*cross(ω_bar,u_bar)[3]*m
    #println("FDx", " ", FDx)
    #println("FLx", " ", FLx)
    #println("FMx", " ", FMx)
    #println("FGx", " ", FGx)
    #println("FCx", " ", FCx)
    Fx =  FDx + FLx +FGx +FMx+FCx
    Fy =  FDy + FLy +FGy +FMy+FCy
    Fz =  FDz + FLz +FGz +FMz+FCz
    #println("Fx", " ", Fx)
    #println("FDx", " ", FDx)
    #println("FLx", " ", FLx)
    #println("FMx", " ", FMx)
    #println("FMy", " ", FMy)
    #println("FMz", " ", FMz)
    #println("FGx", " ", FGx)
    #println("FCx", " ", FCx)
    du[1] = u[4] #x
    du[2] = u[5] #y
    du[3] = u[6] #z
    du[4] =Fx/m #ax
    du[5] = Fy/m#ay
    du[6] = Fz/m#az
    du[7] = pi*ρ*d^4*pₛ*v*Cₛₚᵢₙ/(8*Iₓ)#pdot
    global yaw = αₑ_bar

end

function yawOfRepose(v_bar::Array{Float64,1}, du_bar::Array{Float64,1},pₛ::Float64, t::Float64, Iₓ::Float64, ρ::Float64,d::Float64,CMα::Float64,v::Float64)
    if t==0
        return [0.0,0.0,0.0]
    else
       αₑ_bar = -8*Iₓ*pₛ*cross(v_bar,du_bar)/(pi*ρ*d^3*(CMα)*v^4)

        return αₑ_bar
    end
end



#u₀ = norm()




#canon = gun(0.0,0.0,0.0)

#sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6)
#sol = solve(prob, AutoTsit5(Rosenbrock23()))
#sol = solve(prob)

function condition(u,t,integrator) # Event when event_f(u,t) == 0

  norm(integrator.p[6])-norm([u[1],u[2],u[3]])
  #integrator.p[8][1]-u[1]
  #integrator.p[8][1]>u[1]
end

"""
    conditionHeight(u, t, integrator)

    ContinousCallback function is used to stop algorithm when the Height reach 0.00
"""
function conditionHeight(u,t,integrator)

  if u[1] > 0.1 && u[2] < 0.00
    return norm(integrator.p[10])
  else
    return 1
  end

end


function timeFuze(u,t,integrator)
    #euclidean([u[1],u[2],u[3]],integrator.p[8])-integrator.p[9]

    #close=(euclidean([u[1],u[2],u[3]],integrator.p[8])-integrator.p[9])<0.0
    close=(euclidean([u[1],u[2],u[3]],integrator.p[8])-integrator.p[9])-1e-6
    #println((euclidean([u[1],u[2],u[3]],integrator.p[8])-integrator.p[9]))
    #println(close)
    return close
end

affect!(integrator) = terminate!(integrator)
affect_tf!(integrator) = terminate!(integrator)

"""
    trajectoryMPMM(proj, target, gun, aero; args ... )

Returns the impact point (vector), the velocity at impact (vector) and the time of flight of the projectile (proj).
The stop condition for the trajectory computation is the range of the target. Launcher conditions are specified by the gun
    parameters. The projectile aerodynamic coefficient are specified in the aero variable.
Optional arguments are:

* tspan = the computational time interval. Default is (0.0,1000.0)
* R = the earth radius. Default is 6.356766*1e6
* Ω = the earth spin rate. Default is 7.292115*1e-5
* w_bar = the wind velocity. Default is [0.0,0.0,0.0]
* atm =the atmosphere characteristics. Default is nothing

The projectile trajectory is computed using MPMM

"""
function trajectoryMPMM(proj::AbstractPenetrator, target::AbstractTarget, gun::Gun, aero::DataFrame; tspan = (0.0,1000.0), R = 6.356766*1e6, Ω = 7.292115*1e-5, w_bar=[0.0,0.0,0.0],atm=nothing )
    #global dist = target.position[1]
    #append!(p,dist)
    ω_bar = [Ω*cosd(gun.lat)*cosd(gun.AZ), Ω*sind(gun.lat), -Ω*cosd(gun.lat)*sind(gun.AZ)]
    g₀=grav0(gun.lat)
    if atm==nothing
        tp =nothing

    else
        tp = temperature(atm.t, atm.rh, atm.p)
    end
    p = [proj,w_bar, R, g₀, ω_bar,target.position,5.0,tp,aero]

    p₀ = 2*pi*gun.u₀/(gun.tc*proj.calibre)
    u0=[proj.position[1],proj.position[2],proj.position[3],proj.velocity[1],proj.velocity[2],proj.velocity[3],p₀]
    prob = ODEProblem(mpmm,u0,tspan, p)
    #cb = DiscreteCallback(condition,affect!)
#println(target.x)
    cb = ContinuousCallback(condition,affect!)
    #cb_tf = ContinuousCallback(timeFuze,affect!)
    #cb_tf=DiscreteCallback(timeFuze,affect_tf!)
    #cbs = CallbackSet(cb,cb_tf)
    sol = solve(prob,Tsit5(),callback=cb, reltol=1e-8, abstol=1e-8)
#sol = solve(prob,Tsit5(),callback=cb)
#sol = solve(prob,BS5(),callback=cb)
#sol = solve(prob,AutoTsit5(Rosenbrock23()),callback=cb)
 #proj.position[1] = sol.u[end][1]
 #proj.position[2] = sol.u[end][2]
 #proj.position[3] = sol.u[end][3]
 #proj.velocity[1] = sol.u[end][4]
 #proj.velocity[2] = sol.u[end][5]
 #proj.velocity[3] = sol.u[end][6]
 #proj.tof = sol.t[end]
 #proj.spin = sol.u[end][7]
 return [sol.u[end][1], sol.u[end][2], sol.u[end][3]], [sol.u[end][4], sol.u[end][5], sol.u[end][6]],sol.t[end]
end

function iniCond(gun::Gun, calibre::Float64)
    p₀ = 2*pi*gun.u₀/(gun.tc*calibre)
    u₀_bar = [gun.u₀*cosd(gun.QE)*cosd(gun.AZ), gun.u₀*sind(gun.QE), gun.u₀*cosd(gun.QE)*sind(gun.AZ)]
    X₀_bar = [gun.lw*cosd(gun.QE)*cosd(gun.AZ), gun.X2w + gun.lw *sind(gun.QE), gun.lw*cosd(gun.QE)*sind(gun.AZ)]
    u0 = [X₀_bar[1], X₀_bar[2], X₀_bar[3], u₀_bar[1], u₀_bar[2], u₀_bar[3],p₀]
    return u0
end

"""
    QEfinderMPMM!(drone, proj, gun,aero;args ...)

Returns the elevation and azimuth departure angles of the projectile (proj) with target (drone).
Launcher conditions are specified by the gun parameters, projectile aerodynamic coefficient are specified in the aero variable.
Optional arguments are:

* tspan = the computational time interval. Default is (0.0,1000.0)
* w_bar = the wind velocity. Default is [0.0,0.0,0.0]
* atm =the atmosphere characteristics. Default is nothing

The trajectory is computed using MPMM

"""
function QEfinderMPMM!(drone::AbstractTarget, proj::AbstractPenetrator, gun::Gun,aero::DataFrame;w_bar=[0.0,0.0,0.0],atmosphere=nothing, maxiters=1e6)

    epsilonAz = 1e6
    epsilonQE = 1e6
    precisie = 0.001
    g₀=grav0(gun.lat)
    ddoel = euclidean(proj.position, drone.position)
    tdoel = sqrt(drone.position[1]^2+drone.position[3]^2)/gun.u₀
    #gun.QE = (drone.position[2] - proj.position[2] + g₀ /2 *tdoel^2)*tdoel/gun.u₀
    #gun.AZ = 0.0
    gun.QE=deg2rad(gun.QE)
    gun.AZ=deg2rad(gun.AZ)
    tspan = (0.0,1000.0)
    #R = 6.356766*1e6 #m
    #Ω = 7.292115*1e-5 #rad/s
    #ω_bar = [Ω*cosd(gun.lat)*cosd(gun.AZ), Ω*sind(gun.lat), -Ω*cosd(gun.lat)*sind(gun.AZ)]
    #p = [proj.inertia[1],w_bar, proj.calibre, R, g₀, ω_bar, proj.mass,drone]
    #αₑ_bar = [0.0,0.0,0.0]
    #p = [proj.inertia[1],w_bar, proj.calibre, R, g₀, ω_bar, proj.mass,drone.position,0.0]
    global n=0

    while abs(epsilonAz)>precisie || abs(epsilonQE)>precisie

        #proj.velocity = [gun.u₀*cos(gun.QE)*cos(gun.AZ), gun.u₀*sin(gun.QE), gun.u₀*cos(gun.QE)*sin(gun.AZ)]
        #proj.position = [gun.lw*cos(gun.QE)*cos(gun.AZ), gun.X2w + gun.lw *sin(gun.QE), gun.lw*cos(gun.QE)*sin(gun.AZ)]


        #proj = createProjectile(projectile.mass,projectile.calibre, velocity=muzzle_velocity,  position=muzzle_position, Ix=projectile.Ix )
        #target = ExternalBallistics.createTarget(drone_position)
        #gun = createGun(weapon.u₀,latitude,weapon.θₜ+tourelle.θ,weapon.ξₜ+tourelle.ξ,weapon.tc)

        #u0 = iniCond(gun, proj.calibre)
        #global proj = projectile(u0[1],u0[2],u0[3],u0[4],u0[5],u0[6],0.0)
        #proj = projectile(u0[1],u0[2],u0[3],u0[4],u0[5],u0[6],0.0)
        #proj.position = [u0[1],u0[2],u0[3]]
        #proj.velocity = [u0[4],u0[5],u0[6]]
        #proj.tof = 0.0
        impactP =trajectoryMPMM(proj, drone, gun,aero,atm=atmosphere)[1]
        #trajectoryMPMM!(u0, tspan, p, proj, drone)

        #global epsilonQE = proj.y - drone.y
        epsilonQE = impactP[2] - drone.position[2]
        #println("epsilonQE", " ", epsilonQE)
        #QE = QE0 + (accuracy)/(range_/QE0)
        #global epsilonAz = (sqrt((proj.z-drone.z)^2+(proj.x-drone.x)^2)*sign(atan(proj.z)/proj.x)-atan(drone.z/drone.x))
        epsilonAz = (sqrt((impactP[3]-drone.position[3])^2+(impactP[1]-drone.position[1])^2)*sign(atan(impactP[3])/impactP[1])-atan(drone.position[3]/drone.position[1]))
        #println("epsilonAz", " ", epsilonAz)
        #global AZ = AZ - epsilonAz/ddoel
        gun.AZ = gun.AZ - epsilonAz/ddoel
        #global QE = QE - epsilonQE/ddoel
        gun.QE = gun.QE - epsilonQE/ddoel

        proj.velocity = [gun.u₀*cos(gun.QE)*cos(gun.AZ), gun.u₀*sin(gun.QE), gun.u₀*cos(gun.QE)*sin(gun.AZ)]
        proj.position = [gun.lw*cos(gun.QE)*cos(gun.AZ), gun.X2w + gun.lw *sin(gun.QE), gun.lw*cos(gun.QE)*sin(gun.AZ)]
        if n>=maxiters
            break
        end
        global n+=1


        #trajectory!(u0, tspan, p, proj, drone)
        #calcRange = euclidean([0.0,0.0,0.0], [proj.x,proj.y, proj.z])
        #global accuracy = range_  - calcRange
        #println("AZ", " ", gun.AZ)
        #println("QE", " ", gun.QE)


    end

return gun.QE,gun.AZ
end


function DirectFireTableMPMM(proj::AbstractPenetrator, target::AbstractTarget, gun::Gun, aero::DataFrame, stop_condition::String; tspan = (0.0,1000.0), R = 6.356766*1e6, Ω = 7.292115*1e-5, w_bar=[0.0,0.0,0.0],atm=nothing )

    ω_bar = [Ω*cosd(gun.lat)*cosd(gun.AZ), Ω*sind(gun.lat), -Ω*cosd(gun.lat)*sind(gun.AZ)]
    g₀=grav0(gun.lat)
    if atm==nothing
        tp = nothing
    else
        tp = temperature(atm.t, atm.rh, atm.p)
    end

    p = [proj,w_bar, R, g₀, ω_bar,target.position,5.0,tp,aero, 0.00]
    p₀ = 2*pi*gun.u₀/(gun.tc*proj.calibre)

    u0=[proj.position[1],proj.position[2],proj.position[3],proj.velocity[1],proj.velocity[2],proj.velocity[3],p₀]
    
    prob = ODEProblem(mpmm,u0,tspan, p)
    
    if stop_condition == "Range"
        cb  = ContinuousCallback(condition,affect!)    
    elseif stop_condition == "Height"
        cb  = ContinuousCallback(conditionHeight, affect!)
    end

    sol = solve(prob, Tsit5(), callback=cb, reltol=1e-8, abstol=1e-8)
    return prob, cb, sol

end

"""
```julia
DirecFireTableQEfinderMPMM!(drone, proj, gun,aero, stop_condition;
                            w_bar=[0.0,0.0,0.0],
                            atmosphere=nothing,
                            maxiters=1e6)
```
Returns the elevation and azimuth departure angles of the projectile (proj) with target (drone).
Launcher conditions are specified by the gun parameters, projectile aerodynamic coefficient are specified in the aero variable.

# Arguments
- ´drone´         : AbstractTarget
- `proj`          : AbstractPenetrator
- `gun`           : Gun
- `aero`          : DataFrame
- `stop_condition`: String
# Optional Arguments
- `w_bar`     : the wind velocity. Default is [0.0,0.0,0.0]
- `atmosphere`: the atmosphere characteristics. Defalut is nothing
- `maxiters`  :
"""

function DirecFireTableQEfinderMPMM!(drone::AbstractTarget, proj::AbstractPenetrator, gun::Gun,aero::DataFrame, stop_condition::String;w_bar=[0.0,0.0,0.0],atmosphere=nothing, maxiters=1e6)
    # error precision
    epsilonAz = 1e6
    epsilonQE = 1e6
    precision = 0.001
    
    # out condition initialization
    global n = 0

    # Calibration and Conversion
    g₀     = grav0(gun.lat)
    gun.QE = deg2rad(gun.QE)
    gun.AZ = deg2rad(gun.AZ)

    ddoel = euclidean(proj.position, drone.position)
    tdoel = sqrt(drone.position[1]^2+drone.position[3]^2) / gun.u₀

    tspan = (0.0,1000.0)

    # data array to export value
    arrProb = []
    arrCb   = []
    arrSol  = []

    while abs(epsilonAz)>precision || abs(epsilonQE)>precision

        # calculate the MPMM Trajectory
        prob, cb, sol = DirectFireTableMPMM(proj, drone, gun,aero,stop_condition, atm=atmosphere)

        # Exporting 
        epsilonQE = sol.u[end][2] - drone.position[2]
        epsilonAz = (sqrt((sol.u[end][3]-drone.position[3])^2+(sol.u[end][1]-drone.position[1])^2)*sign(atan(sol.u[end][3])/sol.u[end][1])-atan(drone.position[3]/drone.position[1]))

        gun.AZ = gun.AZ - epsilonAz/ddoel
        gun.QE = gun.QE - epsilonQE/ddoel

        proj.velocity = [gun.u₀*cos(gun.QE)*cos(gun.AZ), gun.u₀*sin(gun.QE), gun.u₀*cos(gun.QE)*sin(gun.AZ)]
        proj.position = [gun.lw*cos(gun.QE)*cos(gun.AZ), gun.X2w + gun.lw *sin(gun.QE), gun.lw*cos(gun.QE)*sin(gun.AZ)]
        
        # out condition in case of infinity loop
        if n>=maxiters
            break
        end # end if maxiters

        push!(arrProb, prob)
        push!(arrCb, cb)
        push!(arrSol, sol)

        global n+=1

    end # end while precision
    return gun.QE,gun.AZ, arrProb, arrCb, arrSol
end # end DirecFireTableQEfinderMPMM!