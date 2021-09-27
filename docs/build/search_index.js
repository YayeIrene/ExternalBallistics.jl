var documenterSearchIndex = {"docs":
[{"location":"#ExternalBallistics.jl-Documentation","page":"Home","title":"ExternalBallistics.jl Documentation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Trajectory computation.","category":"page"},{"location":"#Introduction","page":"Home","title":"Introduction","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"using ExternalBallistics","category":"page"},{"location":"#Table-of-contents","page":"Home","title":"Table of contents","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Pages = [\"index.md\"]","category":"page"},{"location":"#Package-Features","page":"Home","title":"Package Features","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"MPMM direct and indirect","category":"page"},{"location":"#Functions-Documentation","page":"Home","title":"Functions Documentation","text":"","category":"section"},{"location":"#Objects-generation","page":"Home","title":"Objects generation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"createTargetRect","category":"page"},{"location":"#ExternalBallistics.createTargetRect","page":"Home","title":"ExternalBallistics.createTargetRect","text":"createTargetRect(a,b, position)     Return a rectangular target by specifying the length 'a', the height 'b' and the 'position' (range)\n\n\n\n\n\n","category":"function"},{"location":"","page":"Home","title":"Home","text":"createTargetCirc","category":"page"},{"location":"#ExternalBallistics.createTargetCirc","page":"Home","title":"ExternalBallistics.createTargetCirc","text":"createTargetCirc(size, position)\n\nReturns a circular target by specifying the radius (size) and the position (range)\n\n\n\n\n\n","category":"function"},{"location":"","page":"Home","title":"Home","text":"createTarget","category":"page"},{"location":"#ExternalBallistics.createTarget","page":"Home","title":"ExternalBallistics.createTarget","text":"createTarget()\n\nRetursns an object of type 'Target'. Optional arguments are:\n\nρ the size of the target\nposition, the range to the target\n\n\n\n\n\n","category":"function"},{"location":"","page":"Home","title":"Home","text":"createGun","category":"page"},{"location":"#ExternalBallistics.createGun","page":"Home","title":"ExternalBallistics.createGun","text":"createGun(u₀,lat,QE,AZ,tc; args ...)\n\nReturns an object of type 'Gun' by specifying the muzzle velocity (u₀) its latitude (lat), the elevation angle (QE), azimuth angle (AZ) and the twist (tc). Optional arguments are\n\nlw, the gun length\nX2w, the height of the gun\n\n\n\n\n\n","category":"function"},{"location":"#Trajectory-computation","page":"Home","title":"Trajectory computation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"trajectoryMPMM","category":"page"},{"location":"#ExternalBallistics.trajectoryMPMM","page":"Home","title":"ExternalBallistics.trajectoryMPMM","text":"trajectoryMPMM(proj, target, gun, aero; args ... )\n\nReturns the impact point (vector), the velocity at impact (vector) and the time of flight of the projectile (proj). The stop condition for the trajectory computation is the range of the target. Launcher conditions are specified by the gun     parameters. The projectile aerodynamic coefficient are specified in the aero variable. Optional arguments are:\n\ntspan = the computational time interval. Default is (0.0,1000.0)\nR = the earth radius. Default is 6.356766*1e6\nΩ = the earth spin rate. Default is 7.292115*1e-5\nw_bar = the wind velocity. Default is [0.0,0.0,0.0]\natm =the atmosphere characteristics. Default is nothing\n\nThe projectile trajectory is computed using MPMM\n\n\n\n\n\n","category":"function"},{"location":"","page":"Home","title":"Home","text":"QEfinderMPMM!","category":"page"},{"location":"#ExternalBallistics.QEfinderMPMM!","page":"Home","title":"ExternalBallistics.QEfinderMPMM!","text":"QEfinderMPMM!(drone, proj, gun,aero;args ...)\n\nReturns the elevation and azimuth departure angles of the projectile (proj) with target (drone). Launcher conditions are specified by the gun parameters, projectile aerodynamic coefficient are specified in the aero variable. Optional arguments are:\n\ntspan = the computational time interval. Default is (0.0,1000.0)\nw_bar = the wind velocity. Default is [0.0,0.0,0.0]\natm =the atmosphere characteristics. Default is nothing\n\nThe trajectory is computed using MPMM\n\n\n\n\n\n","category":"function"},{"location":"#Types-Documentation","page":"Home","title":"Types Documentation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"AbstractTarget","category":"page"},{"location":"#ExternalBallistics.AbstractTarget","page":"Home","title":"ExternalBallistics.AbstractTarget","text":"abstract type under which target types are defined.\n\n\n\n\n\n","category":"type"},{"location":"","page":"Home","title":"Home","text":"AbstractPenetrator","category":"page"},{"location":"#ExternalBallistics.AbstractPenetrator","page":"Home","title":"ExternalBallistics.AbstractPenetrator","text":"abstract type under which projectiles types are defined.\n\n\n\n\n\n","category":"type"},{"location":"","page":"Home","title":"Home","text":"Wind","category":"page"},{"location":"#ExternalBallistics.Wind","page":"Home","title":"ExternalBallistics.Wind","text":"Wind(cross, range)\n\nWind is defined by the:\n\ncross wind magnitude\nrange wind magnitude\n\n\n\n\n\n","category":"type"},{"location":"","page":"Home","title":"Home","text":"Projectile","category":"page"},{"location":"#ExternalBallistics.Projectile","page":"Home","title":"ExternalBallistics.Projectile","text":"Projectile(mass,calibre,position,tof,Ix,Iy,Xcg)\n\nProjectile is defined by the following parameters:\n\nmass (kg)\ncalibre (m)\npositionsion: coordinate of the 'aimpoint'\nvelocity: velocity vector\ntof: time of flight\nIx: axial moment of inertia\nIy: transversal moment of inertia\nXcg: distance to the center of gravity\n\n\n\n\n\n","category":"type"},{"location":"","page":"Home","title":"Home","text":"TargetRect","category":"page"},{"location":"#ExternalBallistics.TargetRect","page":"Home","title":"ExternalBallistics.TargetRect","text":"TargetRect(a,b,position)\n\ndefines a rectangular target by specifying the lenght a, the height b and the position (range)\n\n\n\n\n\n","category":"type"},{"location":"","page":"Home","title":"Home","text":"TargetCirc","category":"page"},{"location":"#ExternalBallistics.TargetCirc","page":"Home","title":"ExternalBallistics.TargetCirc","text":"TargetCirc(ρ,position)\n\ndefines a circular target by specifying the radius ρ and the position (range)\n\n\n\n\n\n","category":"type"},{"location":"","page":"Home","title":"Home","text":"Air","category":"page"},{"location":"#ExternalBallistics.Air","page":"Home","title":"ExternalBallistics.Air","text":"Air(p,t,rh)\n\ndefines the atmosphere by specifying:\n\np: the pressure\nt: the temperature\nrh: the relative humidity\n\n\n\n\n\n","category":"type"},{"location":"","page":"Home","title":"Home","text":"Gun","category":"page"},{"location":"#ExternalBallistics.Gun","page":"Home","title":"ExternalBallistics.Gun","text":"Gun(u₀,lat,tc,lw,X2w,QE,AZ)\n\nGun is defined by the following parameters:\n\nu₀: the muzzle velocity\nlat: the latitude\ntc: the twist\nlw: the gun length\nX2w: the height\nQE: the the elevation angle\nAZ: the azimuth\n\n\n\n\n\n","category":"type"}]
}