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
Air
```

```@docs
Gun
```
