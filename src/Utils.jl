#module Utils
#using ..Types
#using LinearAlgebra

#export machNumber, density


function machNumber(v::Float64)
    ma = v / 340.0
    return ma
end


function machNumber(v::Float64,t::Float64)
    #tv = temperature(t, RH, P)
    a = speedOfSound(t)
    ma = v / a
    return ma
end


function temperature(y::Float64;y0 = 0.0, t0 = 288.15, β = -6.5)
    t = (t0+β * (y-y0)*1e-3)-273.15 #deg Celsius
    return t
end
function pressvap(t::Float64;A₀=6.1117675,A₁=0.443986062,A₂=1.43053301*1e-2,A₃=2.65027242*1e-4,A₄=3.02246994*1e-6,A₅=2.03886313*1e-8,A₆=6.38780966*1e-11)
    p = (A₀+A₁*t+A₂*t^2+A₃*t^3+A₄*t^4+A₅*t^5+A₆*t^6)*100
end
function temperature(t::Float64, RH::Float64, P::Float64;MadMav=0.622004944)
    ta = t + 273.15 #kelvin
    Ps = pressvap(t)
    t=ta/(1-(Ps*RH)/(100*P)*(1-MadMav))
end

function speedOfSound(t::Float64; γ = 1.4, R = 287.05287)
    #a = sqrt(γ*R*t)
    a=20.046796*sqrt(t+273.15)
    return a
end



function gravity(y::Float64; g0 = 9.80665, r = 6356766.0)
    g = g0*(r/(r+y))^2
    return g
end
function pressure(y::Float64;p0 = 101325.0, y0 = 0.0, t0 = 288.15, β = -6.5, R = 287.05287)
    p = p0*(1+β/t0*(y-y0)*1e-3)^(-gravity(y)/(β*R))
    return p
end
function density(y::Float64; R = 287.05287)
    t = temperature(y)
    p = pressure(y)
    #ρ = p/(R*t)
    ρ = 0.003483678761*p/(t+273.15)
    return ρ
end



function grav0(lat::Float64)
    g₀ =  9.80665*(1-0.0026*cosd(2*lat))
    return  g₀
end


#end
