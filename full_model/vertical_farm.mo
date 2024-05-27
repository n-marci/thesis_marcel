package vertical_farm
model airvolume
  Buildings.HeatTransfer.Sources.FixedTemperature building_air(T = 293.15) annotation(
    Placement(transformation(origin = {170, -10}, extent = {{10, -10}, {-10, 10}})));
  Buildings.Fluid.MixingVolumes.MixingVolume farm_air(m_flow_nominal = 0.1, V = 594, redeclare package Medium = Buildings.Media.Air "Moist air") annotation(
    Placement(transformation(origin = {-10, -30}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  Buildings.HeatTransfer.Conduction.SingleLayer building_glass_conduction(A = 581.04, material = glass) annotation(
    Placement(transformation(origin = {70, 70}, extent = {{-10, -10}, {10, 10}})));
  Buildings.HeatTransfer.Conduction.SingleLayer envelope_glass_conduction(A = 1236.9, material = glass) annotation(
    Placement(transformation(origin = {-90, -10}, extent = {{-10, -10}, {10, 10}})));
  Buildings.HeatTransfer.Convection.Interior building_glass_building_convection(A = 581.04, conMod = Buildings.HeatTransfer.Types.InteriorConvection.Temperature, til = 1.5707963267948966) annotation(
    Placement(transformation(origin = {100, 70}, extent = {{-10, -10}, {10, 10}})));
  Buildings.HeatTransfer.Convection.Interior building_facade_building_convection(A = 523.2, conMod = Buildings.HeatTransfer.Types.InteriorConvection.Temperature, til = 1.5707963267948966) annotation(
    Placement(transformation(origin = {100, -50}, extent = {{-10, -10}, {10, 10}})));
  Buildings.HeatTransfer.Convection.Exterior building_glass_farm_convection(A = 581.04, conMod = Buildings.HeatTransfer.Types.ExteriorConvection.TemperatureWind, azi = 0, til = 1.5707963267948966) annotation(
    Placement(transformation(origin = {40, 70}, extent = {{10, -10}, {-10, 10}})));
  Buildings.HeatTransfer.Convection.Exterior building_facade_farm_convection(A = 523.2, conMod = Buildings.HeatTransfer.Types.ExteriorConvection.TemperatureWind, roughness = Buildings.HeatTransfer.Types.SurfaceRoughness.Medium, azi = 0, til = 1.5707963267948966) annotation(
    Placement(transformation(origin = {40, -50}, extent = {{10, -10}, {-10, 10}})));
  Buildings.HeatTransfer.Convection.Exterior envelope_glass_farm_convection(A = 1236.9, conMod = Buildings.HeatTransfer.Types.ExteriorConvection.TemperatureWind, roughness = Buildings.HeatTransfer.Types.SurfaceRoughness.VerySmooth, azi = 0, til = 1.5707963267948966) annotation(
    Placement(transformation(origin = {-60, -10}, extent = {{-10, -10}, {10, 10}})));
  Buildings.HeatTransfer.Convection.Exterior envelope_glass_outside_convection(conMod = Buildings.HeatTransfer.Types.ExteriorConvection.TemperatureWind, roughness = Buildings.HeatTransfer.Types.SurfaceRoughness.VerySmooth, azi = 0, til = 1.5707963267948966, A = 1236.9) annotation(
    Placement(transformation(origin = {-120, -10}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
  Modelica.Blocks.Sources.Constant windspeed_and_direction_envelope(k = 0) annotation(
    Placement(transformation(origin = {-78, 50}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  Modelica.Blocks.Sources.Constant windspeed_and_direction_building(k = 0) annotation(
    Placement(transformation(origin = {40, -2}, extent = {{-10, -10}, {10, 10}})));
  parameter Buildings.HeatTransfer.Data.Solids.Concrete concrete(x = 0.18, k = 2, c = 550, d = 3500) annotation(
    Placement(transformation(origin = {50, -104}, extent = {{-10, -10}, {10, 10}})));
  parameter Buildings.HeatTransfer.Data.Solids.InsulationBoard eps(x = 0.05, k = 0.035, c = 1200, d = 21) annotation(
    Placement(transformation(origin = {90, -104}, extent = {{-10, -10}, {10, 10}})));
  parameter Buildings.HeatTransfer.Data.Solids.Glass glass(x(displayUnit = "mm") = 0.025) annotation(
    Placement(transformation(origin = {70, 42}, extent = {{-10, -10}, {10, 10}})));
  parameter Buildings.HeatTransfer.Data.OpaqueConstructions.Generic facade_material(nLay = 2, material = {concrete, eps}) annotation(
    Placement(transformation(origin = {70, -76}, extent = {{-10, -10}, {10, 10}})));
  Buildings.HeatTransfer.Conduction.MultiLayer building_facade_conduction(A = 523.2, layers = facade_material) annotation(
    Placement(transformation(origin = {70, -50}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
  Modelica.Blocks.Math.Gain gaiWin(k = 1236.9*0.03) annotation(
    Placement(transformation(origin = {-110, 80}, extent = {{-10, -10}, {10, 10}})));
  Buildings.HeatTransfer.Sources.PrescribedHeatFlow preHeaFlo annotation(
    Placement(transformation(origin = {-70, 80}, extent = {{-10, -10}, {10, 10}})));
  Buildings.HeatTransfer.Sources.PrescribedTemperature preTem1 annotation(
    Placement(transformation(origin = {-196, -34}, extent = {{-10, -10}, {10, 10}})));
  Buildings.BoundaryConditions.WeatherData.ReaderTMY3 weaDat(filNam = "/home/marci/dev/ma/irradiance_model/tmy/DEU_BY_Nurnberg.AP.107630_TMYx.mos") annotation(
    Placement(transformation(origin = {-320, -30}, extent = {{-10, -10}, {10, 10}})));
  Buildings.BoundaryConditions.WeatherData.Bus weaBus annotation(
    Placement(transformation(origin = {-268, -30}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {-274, -32}, extent = {{-10, -10}, {10, 10}})));

  model testing
    Buildings.BoundaryConditions.WeatherData.ReaderTMY3 weaDat annotation(
      Placement(transformation(origin = {-170, 8}, extent = {{-10, -10}, {10, 10}})));
    Buildings.BoundaryConditions.WeatherData.Bus weaBus annotation(
      Placement(transformation(origin = {-122, 10}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {-122, 10}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Math.Gain gain annotation(
      Placement(transformation(origin = {-66, 12}, extent = {{-10, -10}, {10, 10}})));
  equation
    connect(weaBus, weaDat.weaBus) annotation(
      Line(points = {{-122, 10}, {-152, 10}, {-152, 8}, {-160, 8}}, thickness = 0.5));
  annotation(
      Diagram(coordinateSystem(extent = {{-180, 40}, {-40, -20}})));
end testing;

    Buildings.BoundaryConditions.WeatherData.ReaderTMY3 weaDat1(filNam = "/home/marci/dev/ma/irradiance_model/tmy/DEU_BY_Nurnberg.AP.107630_TMYx.mos") annotation(
      Placement(transformation(origin = {-260, -118}, extent = {{-60, 20}, {-40, 40}})));
  Buildings.BoundaryConditions.SolarIrradiation.DirectTiltedSurface HDirRoo(azi = Buildings.Types.Azimuth.S, til (displayUnit = "rad")= 90) annotation(
      Placement(transformation(origin = {-260, -118}, extent = {{20, 20}, {40, 40}})));
  Buildings.BoundaryConditions.SolarIrradiation.DiffuseIsotropic HDifTilIso(til = 1.5707963267948966)  annotation(
      Placement(transformation(origin = {-230, -108}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Math.Add add annotation(
      Placement(transformation(origin = {-190, -98}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Math.Gain solar_to_par(k = 2.02/1000000)  annotation(
      Placement(transformation(origin = {-152, -98}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Continuous.Integrator integrator(use_reset = true)  annotation(
      Placement(transformation(origin = {-72, -98}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Interaction.Show.RealValue realValue annotation(
      Placement(transformation(origin = {-38, -70}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Discrete.Sampler sampler(samplePeriod = 3600)  annotation(
      Placement(transformation(origin = {-118, -98}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.SampleTrigger sampleTrigger(period = 86400)  annotation(
      Placement(transformation(origin = {-90, -130}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Math.Mean mean(f(displayUnit = "nHz") = 3.73e-7, yGreaterOrEqualZero = true)  annotation(
      Placement(transformation(origin = {20, -98}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.Constant optimal_dli(k = 14.4)  annotation(
      Placement(transformation(origin = {-16, -156}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Math.Add missing_dli(k1 = -1)  annotation(
      Placement(transformation(origin = {20, -130}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Discrete.Sampler sampler1(samplePeriod = 86400, startTime = 72000)  annotation(
      Placement(transformation(origin = {-14, -98}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.Constant hourly_light_integral(k = 218.28*60*60/1000000)  annotation(
      Placement(transformation(origin = {20, -166}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Math.Division hours_of_suppl_lighting annotation(
      Placement(transformation(origin = {70, -150}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.Constant led_watts(k = 77.76)  annotation(
      Placement(transformation(origin = {70, -186}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Math.Product led_power_draw annotation(
      Placement(transformation(origin = {120, -170}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.Constant led_efficiency(k = 0.622)  annotation(
      Placement(transformation(origin = {120, -204}, extent = {{-10, -10}, {10, 10}})));
  equation
    connect(envelope_glass_outside_convection.solid, envelope_glass_conduction.port_a) annotation(
      Line(points = {{-110, -10}, {-100, -10}}, color = {191, 0, 0}));
    connect(envelope_glass_conduction.port_b, envelope_glass_farm_convection.solid) annotation(
      Line(points = {{-80, -10}, {-70, -10}}, color = {191, 0, 0}));
    connect(windspeed_and_direction_envelope.y, envelope_glass_farm_convection.v) annotation(
      Line(points = {{-78, 39}, {-78, 0}, {-72, 0}}, color = {0, 0, 127}));
    connect(windspeed_and_direction_envelope.y, envelope_glass_farm_convection.dir) annotation(
      Line(points = {{-78, 39}, {-78, -4}, {-72, -4}}, color = {0, 0, 127}));
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
    connect(gaiWin.y, preHeaFlo.Q_flow) annotation(
      Line(points = {{-99, 80}, {-81, 80}}, color = {0, 0, 127}));
    connect(preHeaFlo.port, farm_air.heatPort) annotation(
      Line(points = {{-60, 80}, {-10, 80}, {-10, -20}}, color = {191, 0, 0}));
    connect(preTem1.port, envelope_glass_outside_convection.fluid) annotation(
      Line(points = {{-186, -34}, {-160, -34}, {-160, -10}, {-130, -10}}, color = {191, 0, 0}));
    connect(weaDat.weaBus, weaBus) annotation(
      Line(points = {{-310, -30}, {-268, -30}}, color = {255, 204, 51}, thickness = 0.5));
    connect(weaBus.TDryBul, preTem1.T) annotation(
      Line(points = {{-268, -30}, {-232, -30}, {-232, -34}, {-208, -34}}, color = {0, 0, 127}));
    connect(gaiWin.u, weaBus.HDirNor);
//  connect(HGloTil.HDifTil.incAng.weaBus, weaBus) annotation(
//    Line(points = {{-210, 80}, {-268, 80}, {-268, -30}}, color = {255, 204, 51}, thickness = 0.5));
    connect(envelope_glass_outside_convection.v, weaBus.winSpe) annotation(
      Line(points = {{-108, 0}, {-102, 0}, {-102, 14}, {-268, 14}, {-268, -30}}, color = {0, 0, 127}));
    connect(envelope_glass_outside_convection.dir, weaBus.winDir) annotation(
      Line(points = {{-108, -4}, {-100, -4}, {-100, 16}, {-268, 16}, {-268, -30}}, color = {0, 0, 127}));
    connect(weaDat1.weaBus, HDirRoo.weaBus) annotation(
      Line(points = {{-300, -88}, {-240, -88}}, color = {255, 204, 51}, thickness = 0.5));
    connect(HDifTilIso.weaBus, weaDat1.weaBus) annotation(
      Line(points = {{-240, -108}, {-270, -108}, {-270, -88}, {-300, -88}}, color = {255, 204, 51}, thickness = 0.5));
    connect(HDirRoo.H, add.u1) annotation(
      Line(points = {{-218, -88}, {-211, -88}, {-211, -92}, {-202, -92}}, color = {0, 0, 127}));
    connect(HDifTilIso.H, add.u2) annotation(
      Line(points = {{-218, -108}, {-214, -108}, {-214, -104}, {-202, -104}}, color = {0, 0, 127}));
    connect(add.y, solar_to_par.u) annotation(
      Line(points = {{-178, -98}, {-164, -98}}, color = {0, 0, 127}));
    connect(solar_to_par.y, sampler.u) annotation(
      Line(points = {{-140, -98}, {-130, -98}}, color = {0, 0, 127}));
    connect(sampler.y, integrator.u) annotation(
      Line(points = {{-106, -98}, {-84, -98}}, color = {0, 0, 127}));
    connect(integrator.y, realValue.numberPort) annotation(
      Line(points = {{-60, -98}, {-58, -98}, {-58, -70}, {-50, -70}}, color = {0, 0, 127}));
  connect(sampleTrigger.y, integrator.reset) annotation(
      Line(points = {{-78, -130}, {-66, -130}, {-66, -110}}, color = {255, 0, 255}));
  connect(optimal_dli.y, missing_dli.u2) annotation(
      Line(points = {{-5, -156}, {-1, -156}, {-1, -136}, {7, -136}}, color = {0, 0, 127}));
  connect(sampler1.u, integrator.y) annotation(
      Line(points = {{-26, -98}, {-60, -98}}, color = {0, 0, 127}));
  connect(sampler1.y, mean.u) annotation(
      Line(points = {{-2, -98}, {8, -98}}, color = {0, 0, 127}));
  connect(sampler1.y, missing_dli.u1) annotation(
      Line(points = {{-2, -98}, {0, -98}, {0, -124}, {8, -124}}, color = {0, 0, 127}));
  connect(missing_dli.y, hours_of_suppl_lighting.u1) annotation(
      Line(points = {{32, -130}, {44, -130}, {44, -144}, {58, -144}}, color = {0, 0, 127}));
  connect(hourly_light_integral.y, hours_of_suppl_lighting.u2) annotation(
      Line(points = {{31, -166}, {52.5, -166}, {52.5, -156}, {58, -156}}, color = {0, 0, 127}));
  connect(hours_of_suppl_lighting.y, led_power_draw.u1) annotation(
      Line(points = {{82, -150}, {92, -150}, {92, -164}, {108, -164}}, color = {0, 0, 127}));
  connect(led_watts.y, led_power_draw.u2) annotation(
      Line(points = {{82, -186}, {88, -186}, {88, -176}, {108, -176}}, color = {0, 0, 127}));
    annotation(
    Diagram(coordinateSystem(extent = {{-340, 100}, {180, -120}}), graphics = {Rectangle(origin = {70, -77}, extent = {{-44, 47}, {44, -47}}), Rectangle(origin = {71, 58}, extent = {{-45, 32}, {45, -32}}), Rectangle(origin = {-90, 37}, extent = {{-46, 63}, {46, -63}}), Text(origin = {50, -118}, extent = {{-18, 2}, {18, -2}}, textString = "https://doi.org/10.4028/www.scientific.net/SSP.321.113"), Text(origin = {92, -118}, extent = {{-18, 2}, {18, -2}}, textString = "https://doi.org/10.1088/1742-6596%2F2069%2F1%2F012090")}));
  annotation(
    Diagram(coordinateSystem(extent = {{-340, 100}, {180, -120}})));
end airvolume;

  model testing3
    parameter Modelica.Units.SI.Area A=1 "Window surface area";
    parameter Real fFra=0.2
      "Fraction of frame, = frame area divided by total area";
    final parameter Modelica.Units.SI.Area AFra=fFra*A "Frame area";
    final parameter Modelica.Units.SI.Area AGla=A - AFra "Glass area";
    parameter Boolean linearize = false "Set to true to linearize emissive power";
    parameter Modelica.Units.SI.Angle azi=0 "Surface azimuth";
    parameter Modelica.Units.SI.Angle til=1.5707963267949 "Surface tilt";
  Buildings.BoundaryConditions.WeatherData.ReaderTMY3 weaDat(filNam = Modelica.Utilities.Files.loadResource("modelica://Buildings/Resources/weatherdata/USA_IL_Chicago-OHare.Intl.AP.725300_TMY3.mos"))  annotation(
      Placement(transformation(origin = {-70, 96}, extent = {{-10, -10}, {10, 10}})));
  Buildings.HeatTransfer.Windows.SideFins fin(h = 0.2, dep = 0.5, gap = 0.1, hWin = 1.0, wWin = 1.0)  annotation(
      Placement(transformation(origin = {10, 86}, extent = {{-10, -10}, {10, 10}})));
  Buildings.BoundaryConditions.SolarIrradiation.DirectTiltedSurface HDirRoo(azi = 0.78539816339745, til = Buildings.Types.Tilt.Ceiling) annotation(
      Placement(transformation(origin = {-60, 86}, extent = {{20, 20}, {40, 40}})));
  equation
    connect(weaDat.weaBus, HDirTil.incAng.weaBus) annotation(
      Line(points = {{-60, 96}, {-52, 96}, {-52, 112}, {-40, 112}}, color = {255, 204, 51}, thickness = 0.5));
    connect(weaDat.weaBus, HDirTil.incAng.weaBus) annotation(
      Line(points = {{-60, 96}, {-48, 96}, {-48, 112}, {-34, 112}}, color = {255, 204, 51}, thickness = 0.5));
  connect(weaDat.weaBus, HDirRoo.incAng.weaBus) annotation(
      Line(points = {{-60, 96}, {-54, 96}, {-54, 116}, {-40, 116}}, color = {255, 204, 51}, thickness = 0.5));
  connect(HDirRoo.H, fin.HDirTilUns) annotation(
      Line(points = {{-18, 116}, {-12, 116}, {-12, 92}, {-2, 92}}, color = {0, 0, 127}));
  connect(HDirRoo.inc, fin.incAng) annotation(
      Line(points = {{-18, 112}, {-14, 112}, {-14, 80}, {-2, 80}}, color = {0, 0, 127}));
  connect(fin.incAng, HDirRoo.inc) annotation(
      Line(points = {{-2, 80}, {-16, 80}, {-16, 112}, {-18, 112}}, color = {0, 0, 127}));
  connect(fin.weaBus, weaDat.weaBus) annotation(
      Line(points = {{0, 86}, {-40, 86}, {-40, 96}, {-60, 96}}, color = {255, 204, 51}, thickness = 0.5));
    annotation(
      Diagram(coordinateSystem(extent = {{-100, 180}, {320, -80}})));
end testing3;

  model testing4
  parameter Real rho=0.2 "Ground reflectance";
  Buildings.BoundaryConditions.WeatherData.ReaderTMY3 weaDat(filNam = Modelica.Utilities.Files.loadResource("modelica://Buildings/Resources/weatherdata/USA_IL_Chicago-OHare.Intl.AP.725300_TMY3.mos")) annotation(
      Placement(transformation(extent = {{-60, 20}, {-40, 40}})));
  Buildings.BoundaryConditions.SolarIrradiation.DirectTiltedSurface HDirRoo(azi = 0.78539816339745, til = Buildings.Types.Tilt.Ceiling) annotation(
      Placement(transformation(extent = {{20, 20}, {40, 40}})));
  equation
    connect(weaDat.weaBus, HDirRoo.weaBus) annotation(
      Line(points = {{-40, 30}, {20, 30}}, color = {255, 204, 51}, thickness = 0.5));
  end testing4;
  annotation(
    uses(Buildings(version = "11.0.0"), Modelica(version = "4.0.0")));
end vertical_farm;