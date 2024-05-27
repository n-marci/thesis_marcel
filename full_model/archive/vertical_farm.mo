package vertical_farm
  model airvolume
  Buildings.HeatTransfer.Sources.FixedTemperature outside_air(T = 268.15)  annotation(
      Placement(transformation(origin = {-190, -10}, extent = {{-10, -10}, {10, 10}})));
  Buildings.HeatTransfer.Sources.FixedTemperature building_air(T = 293.15)  annotation(
      Placement(transformation(origin = {170, -10}, extent = {{10, -10}, {-10, 10}})));
  Buildings.Fluid.MixingVolumes.MixingVolume farm_air(m_flow_nominal = 0.1, V = 594, redeclare package Medium = Buildings.Media.Air "Moist air")  annotation(
      Placement(transformation(origin = {-10, -30}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  Buildings.HeatTransfer.Conduction.SingleLayer building_glass_conduction(A = 581.04, material = glass)  annotation(
      Placement(transformation(origin = {70, 70}, extent = {{-10, -10}, {10, 10}})));
  Buildings.HeatTransfer.Conduction.SingleLayer envelope_glass_conduction(A = 1236.9, material = glass)  annotation(
      Placement(transformation(origin = {-90, -10}, extent = {{-10, -10}, {10, 10}})));
  Buildings.HeatTransfer.Convection.Interior building_glass_building_convection(A = 581.04, conMod = Buildings.HeatTransfer.Types.InteriorConvection.Temperature, til = 1.5707963267948966)  annotation(
      Placement(transformation(origin = {100, 70}, extent = {{-10, -10}, {10, 10}})));
  Buildings.HeatTransfer.Convection.Interior building_facade_building_convection(A = 523.2, conMod = Buildings.HeatTransfer.Types.InteriorConvection.Temperature, til = 1.5707963267948966)  annotation(
      Placement(transformation(origin = {100, -50}, extent = {{-10, -10}, {10, 10}})));
  Buildings.HeatTransfer.Convection.Exterior building_glass_farm_convection(A = 581.04, conMod = Buildings.HeatTransfer.Types.ExteriorConvection.TemperatureWind, azi = 0, til = 1.5707963267948966)  annotation(
      Placement(transformation(origin = {40, 70}, extent = {{10, -10}, {-10, 10}})));
  Buildings.HeatTransfer.Convection.Exterior building_facade_farm_convection(A = 523.2, conMod = Buildings.HeatTransfer.Types.ExteriorConvection.TemperatureWind, roughness = Buildings.HeatTransfer.Types.SurfaceRoughness.Medium, azi = 0, til = 1.5707963267948966)  annotation(
      Placement(transformation(origin = {40, -50}, extent = {{10, -10}, {-10, 10}})));
  Buildings.HeatTransfer.Convection.Exterior envelope_glass_farm_convection(A = 1236.9, conMod = Buildings.HeatTransfer.Types.ExteriorConvection.TemperatureWind, roughness = Buildings.HeatTransfer.Types.SurfaceRoughness.VerySmooth, azi = 0, til = 1.5707963267948966)  annotation(
      Placement(transformation(origin = {-60, -10}, extent = {{-10, -10}, {10, 10}})));
  Buildings.HeatTransfer.Convection.Exterior envelope_glass_outside_convection(conMod = Buildings.HeatTransfer.Types.ExteriorConvection.TemperatureWind, roughness = Buildings.HeatTransfer.Types.SurfaceRoughness.VerySmooth, azi = 0, til = 1.5707963267948966, A = 1236.9)  annotation(
      Placement(transformation(origin = {-120, -10}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
  Modelica.Blocks.Sources.Constant windspeed_and_direction_envelope(k = 0)  annotation(
      Placement(transformation(origin = {-90, 50}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  Modelica.Blocks.Sources.Constant windspeed_and_direction_building(k = 0)  annotation(
      Placement(transformation(origin = {40, -2}, extent = {{-10, -10}, {10, 10}})));
  parameter Buildings.HeatTransfer.Data.Solids.Concrete concrete(x = 0.18, k = 2, c = 550, d = 3500)  annotation(
      Placement(transformation(origin = {50, -104}, extent = {{-10, -10}, {10, 10}})));
  parameter Buildings.HeatTransfer.Data.Solids.InsulationBoard eps(x = 0.05, k = 0.035, c = 1200, d = 21)  annotation(
      Placement(transformation(origin = {90, -104}, extent = {{-10, -10}, {10, 10}})));
  parameter Buildings.HeatTransfer.Data.Solids.Glass glass(x(displayUnit = "mm") = 0.025)  annotation(
      Placement(transformation(origin = {70, 42}, extent = {{-10, -10}, {10, 10}})));
  parameter Buildings.HeatTransfer.Data.OpaqueConstructions.Generic facade_material(nLay = 2, material = {concrete, eps})  annotation(
      Placement(transformation(origin = {70, -76}, extent = {{-10, -10}, {10, 10}})));
  Buildings.HeatTransfer.Conduction.MultiLayer building_facade_conduction(A = 523.2, layers = facade_material) annotation(
      Placement(transformation(origin = {70, -50}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
  Buildings.HeatTransfer.Windows.Window win(fFra = 0, A = 1236.9, glaSys = datGlaSys1, til = 1.5707963267948966)  annotation(
      Placement(transformation(origin = {-40, -80}, extent = {{-20, -20}, {20, 20}})));
  parameter Buildings.HeatTransfer.Data.GlazingSystems.SingleClear3 datGlaSys annotation(
      Placement(transformation(origin = {-170, -70}, extent = {{-10, -10}, {10, 10}})));
  parameter Buildings.HeatTransfer.Data.GlazingSystems.DoubleClearAir13Clear datGlaSys1 annotation(
      Placement(transformation(origin = {-170, -90}, extent = {{-10, -10}, {10, 10}})));
  parameter Buildings.HeatTransfer.Data.GlazingSystems.Generic datGlaSys2 annotation(
      Placement(transformation(origin = {-170, -48}, extent = {{-10, -10}, {10, 10}})));
  Buildings.HeatTransfer.Windows.BaseClasses.CenterOfGlass centerOfGlass annotation(
      Placement(transformation(origin = {-110, -130}, extent = {{-10, -10}, {10, 10}})));
  Buildings.HeatTransfer.Windows.BaseClasses.WindowRadiation windowRadiation(haveExteriorShade = false, haveInteriorShade = false, N = 1, xGla = 0.025, AWin = 1236.9, tauGlaSol = 0.89)  annotation(
      Placement(transformation(origin = {-70, -130}, extent = {{-10, -10}, {10, 10}})));
  equation
    connect(outside_air.port, envelope_glass_outside_convection.fluid) annotation(
      Line(points = {{-180, -10}, {-130, -10}}, color = {191, 0, 0}));
    connect(envelope_glass_outside_convection.solid, envelope_glass_conduction.port_a) annotation(
      Line(points = {{-110, -10}, {-100, -10}}, color = {191, 0, 0}));
    connect(envelope_glass_conduction.port_b, envelope_glass_farm_convection.solid) annotation(
      Line(points = {{-80, -10}, {-70, -10}}, color = {191, 0, 0}));
    connect(windspeed_and_direction_envelope.y, envelope_glass_outside_convection.v) annotation(
      Line(points = {{-90, 40}, {-90, 12}, {-102, 12}, {-102, 0}, {-108, 0}}, color = {0, 0, 127}));
    connect(windspeed_and_direction_envelope.y, envelope_glass_outside_convection.dir) annotation(
      Line(points = {{-90, 40}, {-90, 12}, {-102, 12}, {-102, -4}, {-108, -4}}, color = {0, 0, 127}));
    connect(windspeed_and_direction_envelope.y, envelope_glass_farm_convection.v) annotation(
      Line(points = {{-90, 40}, {-90, 12}, {-78, 12}, {-78, 0}, {-72, 0}}, color = {0, 0, 127}));
    connect(windspeed_and_direction_envelope.y, envelope_glass_farm_convection.dir) annotation(
      Line(points = {{-90, 40}, {-90, 12}, {-78, 12}, {-78, -4}, {-72, -4}}, color = {0, 0, 127}));
    connect(envelope_glass_farm_convection.fluid, farm_air.heatPort) annotation(
      Line(points = {{-50, -10}, {-10, -10}, {-10, -20}}, color = {191, 0, 0}));
    connect(envelope_glass_farm_convection.fluid, building_facade_farm_convection.fluid) annotation(
      Line(points = {{-50, -10}, {20, -10}, {20, -50}, {30, -50}}, color = {191, 0, 0}));
    connect(envelope_glass_farm_convection.fluid, building_glass_farm_convection.fluid) annotation(
      Line(points = {{-50, -10}, {20, -10}, {20, 70}, {30, 70}}, color = {191, 0, 0}));
    connect(building_facade_building_convection.fluid, building_glass_building_convection.fluid) annotation(
      Line(points = {{110, -50}, {120, -50}, {120, 70}, {110, 70}}, color = {191, 0, 0}));
    connect(building_facade_building_convection.fluid, building_air.port) annotation(
      Line(points = {{110, -50}, {120, -50}, {120, -10}, {160, -10}}, color = {191, 0, 0}));
    connect(building_glass_building_convection.solid, building_glass_conduction.port_b) annotation(
      Line(points = {{90, 70}, {80, 70}}, color = {191, 0, 0}));
    connect(building_glass_conduction.port_a, building_glass_farm_convection.solid) annotation(
      Line(points = {{60, 70}, {50, 70}}, color = {191, 0, 0}));
    connect(windspeed_and_direction_building.y, building_glass_farm_convection.dir) annotation(
      Line(points = {{52, -2}, {58, -2}, {58, 76}, {52, 76}}, color = {0, 0, 127}));
    connect(windspeed_and_direction_building.y, building_glass_farm_convection.v) annotation(
      Line(points = {{52, -2}, {58, -2}, {58, 80}, {52, 80}}, color = {0, 0, 127}));
    connect(windspeed_and_direction_building.y, building_facade_farm_convection.v) annotation(
      Line(points = {{52, -2}, {58, -2}, {58, -40}, {52, -40}}, color = {0, 0, 127}));
    connect(windspeed_and_direction_building.y, building_facade_farm_convection.dir) annotation(
      Line(points = {{52, -2}, {58, -2}, {58, -44}, {52, -44}}, color = {0, 0, 127}));
  connect(building_facade_farm_convection.solid, building_facade_conduction.port_b) annotation(
      Line(points = {{50, -50}, {60, -50}}, color = {191, 0, 0}));
  connect(building_facade_conduction.port_a, building_facade_building_convection.solid) annotation(
      Line(points = {{80, -50}, {90, -50}}, color = {191, 0, 0}));
    annotation(
      Diagram(coordinateSystem(extent = {{-200, 100}, {180, -100}}), graphics = {Rectangle(origin = {70, -75}, extent = {{-44, 45}, {44, -45}}), Rectangle(origin = {71, 58}, extent = {{-45, 32}, {45, -32}}), Rectangle(origin = {-90, 21}, extent = {{-46, 47}, {46, -47}}), Text(origin = {50, -118}, extent = {{-18, 2}, {18, -2}}, textString = "https://doi.org/10.4028/www.scientific.net/SSP.321.113"), Text(origin = {92, -118}, extent = {{-18, 2}, {18, -2}}, textString = "https://doi.org/10.1088/1742-6596%2F2069%2F1%2F012090")}));
end airvolume;
  annotation(
    uses(Buildings(version = "11.0.0"), Modelica(version = "4.0.0")));
end vertical_farm;