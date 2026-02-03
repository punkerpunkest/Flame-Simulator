# Flame Simulator


![Flame_Simulator_C__SFML](https://github.com/user-attachments/assets/f4ad4f36-0168-411c-815d-fa11101fcda8)


We once heard that once you learn to simulate a liquid, you can also simulate gases and fire since they share a lot of the same underlying principles. We decided to test this claim by building on the work of Stam and Nguyen et al., who did wonders for liquid simulation in games. To make the simulation as realistic as possible, we have chosen to include the following:


## Physics Features

### Navier-Stokes Fluid Dynamics
- **Diffusion**: Smoothing of fluid motion with viscocity
- **Advection**: Advection of velocity through a flow field
- **Mass Conservation**: Incompressibility of a fluid

### Thermal Dynamics
- **Temperature Field**: Separate temperature field diffused and advected
- **Buoyancy**: Hot gases rise, cold gases sink
- **Cooling**: Cooling towards ambient temperature

### Additional Effects
- **Vorticity Confinement**: Captures the small-scale swirls you see in flames
- **Smoke**: Smoke is advected and diffused around the flame

This project also gave us the chance to experiment with a lot of optimizations using multi-threading and SIMD that we benchmarked and found to be quite effective. SIMD was found to increase the performance of diffusion and projection by 30%! To try out this project, the following are necessary:

## Prerequisites

- C++20 compatible compiler
- SFML
- OpenMP (libomp on macOS)

## Installation

### macOS
```bash
brew install sfml libomp
```

### Linux (Ubuntu/Debian)
```bash
sudo apt-get install libsfml-dev libomp-dev
```

## To run

```bash
make
make run
```


