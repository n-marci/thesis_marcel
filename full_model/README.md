### Quick Start

- These steps are subject to change as I may change from a model to a package as the simulation becomes more complicated
- I am using [OpenModelica](https://openmodelica.org/#) for my simulation
- Install OpenModelica
- Install the Buildings library
  - On the welcome page click on install library
  - Click on name and choose Buildings library
  - Click on OK
- Open up the vertical_farm.mo file from this repository
- Open up the file airvolume from the left sidebar
- Click on the weaDat Block all the way on the left (the little sun and cloud)
- Change the weather data file to the appropriate path for you
  - The file is located in irradiance_model/tmy/DEU_BY_Nurnberg.AP.107630_TMYx.mos
- Click on Simulate (top bar)
- Hopefully everything works

### Status

- seperating the led calculation into its own model
