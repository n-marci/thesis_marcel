# Master thesis repository of Marcel Neugebauer

## Technical Goals of this work

- [x] estimate energy consumption of LEDs (biggest energy consumer in traditional vertical farms with no natural lighting)
  - [x] get detailed data on irradiance by the sun according to the physical location and orientation of the farm
  - [x] gather information on the plant to be cultivated (lettuce for now)
- [ ] DOING estimate energy consumption of HVAC system (second biggest energy consumer)
  - [x] erect thermodynamic energy balances
  - [ ] DOING simulate in modelica - supposes homogeneous environment inside farm (1d simulation)
  - [ ] 3d gas simulation using openfoam (finite volume method) or dualsphysics (smoothed particle hydrodynamics) in unity/blender/freecad
  - [ ] combine simulations with fmu and co-simulation
- [ ] estimate energy consumption of remaining components (irrigation, sensors, mechanical actors)
- [ ] DOING develop system architecture
  - [x] use case analysis
  - [ ] design functional architecture
  - [ ] assign modules to functions
- [ ] DOING build 3d model to convey idea better
- [ ] calculate appropriate solar installation
- [ ] calculate appropriate energy storage
- [ ] calculate insulation value of "naked farm" (in winter it will most likely be to energy intensive to operate)

## Overview

- irradiance model is a python script which focuses on the lighting for the plants
- full model is a openmodelica model where I eventually want to put everything. At the moment it is focused on the thermodynamic examination
- visualizations contains some useful visualizations which I make to communicate my results
