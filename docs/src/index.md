# ExternalBallistics.jl Documentation

*Trajectory computation.*
## Introduction

```@repl
using ExternalBallistics
```
## Table of contents


```@contents
Pages = ["index.md"]
```

## Package Features
- MPMM direct and indirect
## Functions Documentation

### Objects generation
```@docs
createTargetRect
```

```@docs
createTargetCirc
```

```@docs
createTarget
```

```@docs
createGun
```
### Trajectory computation
```@docs
trajectoryMPMM
```

```@docs
QEfinderMPMM!
```

## Types Documentation

```@docs
AbstractTarget
```

```@docs
AbstractPenetrator
```
```@docs
Wind
```
```@docs
Projectile
```
```@docs
TargetRect
```

```@docs
TargetCirc
```

```@docs
Air
```

```@docs
Gun
```
