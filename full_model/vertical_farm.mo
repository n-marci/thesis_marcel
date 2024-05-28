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
    Placement(transformation(origin = {156, -76}, extent = {{-10, -10}, {10, 10}})));
  parameter Buildings.HeatTransfer.Data.Solids.InsulationBoard eps(x = 0.05, k = 0.035, c = 1200, d = 21) annotation(
    Placement(transformation(origin = {196, -76}, extent = {{-10, -10}, {10, 10}})));
  parameter Buildings.HeatTransfer.Data.Solids.Glass glass(x(displayUnit = "mm") = 0.025) annotation(
    Placement(transformation(origin = {70, 42}, extent = {{-10, -10}, {10, 10}})));
  parameter Buildings.HeatTransfer.Data.OpaqueConstructions.Generic facade_material(nLay = 2, material = {concrete, eps}) annotation(
    Placement(transformation(origin = {176, -48}, extent = {{-10, -10}, {10, 10}})));
  Buildings.HeatTransfer.Conduction.MultiLayer building_facade_conduction(A = 523.2, layers = facade_material) annotation(
    Placement(transformation(origin = {70, -50}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
  Modelica.Blocks.Math.Gain gaiWin(k = 1236.9*0.03) annotation(
    Placement(transformation(origin = {-110, 80}, extent = {{-10, -10}, {10, 10}})));
  Buildings.HeatTransfer.Sources.PrescribedHeatFlow preHeaFlo annotation(
    Placement(transformation(origin = {-70, 80}, extent = {{-10, -10}, {10, 10}})));
  Buildings.HeatTransfer.Sources.PrescribedTemperature preTem1 annotation(
    Placement(transformation(origin = {-196, -34}, extent = {{-10, -10}, {10, 10}})));
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
      Placement(transformation(origin = {-520, -104}, extent = {{-60, 20}, {-40, 40}})));
  Buildings.BoundaryConditions.SolarIrradiation.DirectTiltedSurface HDirRoo(azi = Buildings.Types.Azimuth.S, til (displayUnit = "rad")= 90) annotation(
      Placement(transformation(origin = {-520, -104}, extent = {{20, 20}, {40, 40}})));
  Buildings.BoundaryConditions.SolarIrradiation.DiffuseIsotropic HDifTilIso(til = 1.5707963267948966)  annotation(
      Placement(transformation(origin = {-490, -94}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Math.Add add annotation(
      Placement(transformation(origin = {-450, -84}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Math.Gain solar_to_par(k = 2.02/1000000)  annotation(
      Placement(transformation(origin = {-412, -84}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Continuous.Integrator integrator(use_reset = true)  annotation(
      Placement(transformation(origin = {-332, -84}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Interaction.Show.RealValue realValue annotation(
      Placement(transformation(origin = {-294, -56}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.SampleTrigger daily_trigger(period = 86400)  annotation(
      Placement(transformation(origin = {-350, -116}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Math.Mean mean(f(displayUnit = "nHz") = 3.73e-7, yGreaterOrEqualZero = true)  annotation(
      Placement(transformation(origin = {-240, -84}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.Constant optimal_dli(k = 14.4)  annotation(
      Placement(transformation(origin = {-276, -142}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Math.Add missing_dli(k1 = -1)  annotation(
      Placement(transformation(origin = {-240, -116}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Discrete.Sampler sampler1(samplePeriod = 86400, startTime = 72000)  annotation(
      Placement(transformation(origin = {-274, -84}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.Constant hourly_light_integral(k = 218.28*60*60/1000000)  annotation(
      Placement(transformation(origin = {-240, -152}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Math.Division hours_of_suppl_lighting annotation(
      Placement(transformation(origin = {-190, -136}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.Constant led_watts(k = 77.76)  annotation(
      Placement(transformation(origin = {-158, -180}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Math.Product led_power_draw annotation(
      Placement(transformation(origin = {-108, -156}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.Constant led_efficiency(k = 0.622)  annotation(
      Placement(transformation(origin = {-138, -228}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.Constant full_efficiency(k = 1)  annotation(
      Placement(transformation(origin = {-138, -264}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Math.Add led_heat_loss(k1 = -1)  annotation(
      Placement(transformation(origin = {-94, -242}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Math.Product led_heat_generation annotation(
      Placement(transformation(origin = {-54, -162}, extent = {{-10, -10}, {10, 10}})));
  Buildings.HeatTransfer.Sources.PrescribedHeatFlow led_prescriped_heat_flow annotation(
      Placement(transformation(origin = {-22, -162}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Nonlinear.Limiter limiter(uMin = 0, uMax = Modelica.Constants.inf)  annotation(
      Placement(transformation(origin = {-152, -140}, extent = {{-10, -10}, {10, 10}})));
  Buildings.BoundaryConditions.SolarGeometry.IncidenceAngle incAng(azi = 0, til = 1.5707963267948966)  annotation(
      Placement(transformation(origin = {-490, -122}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Continuous.Integrator led_cumulative_power_draw annotation(
      Placement(transformation(origin = {-60, -110}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Math.Sum sum1 annotation(
      Placement(transformation(origin = {-448, -126}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Math.MultiSum multiSum annotation(
      Placement(transformation(origin = {-572, -134}, extent = {{-10, -10}, {10, 10}})));
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
      Line(points = {{-560, -74}, {-500, -74}}, color = {255, 204, 51}, thickness = 0.5));
    connect(HDifTilIso.weaBus, weaDat1.weaBus) annotation(
      Line(points = {{-500, -94}, {-530, -94}, {-530, -74}, {-560, -74}}, color = {255, 204, 51}, thickness = 0.5));
    connect(HDirRoo.H, add.u1) annotation(
      Line(points = {{-479, -74}, {-472, -74}, {-472, -78}, {-463, -78}}, color = {0, 0, 127}));
    connect(HDifTilIso.H, add.u2) annotation(
      Line(points = {{-479, -94}, {-475, -94}, {-475, -90}, {-463, -90}}, color = {0, 0, 127}));
    connect(add.y, solar_to_par.u) annotation(
      Line(points = {{-439, -84}, {-425, -84}}, color = {0, 0, 127}));
    connect(integrator.y, realValue.numberPort) annotation(
      Line(points = {{-321, -84}, {-314.25, -84}, {-314.25, -56}, {-305.5, -56}}, color = {0, 0, 127}));
    connect(daily_trigger.y, integrator.reset) annotation(
      Line(points = {{-339, -116}, {-327, -116}, {-327, -96}}, color = {255, 0, 255}));
    connect(optimal_dli.y, missing_dli.u2) annotation(
      Line(points = {{-265, -142}, {-261, -142}, {-261, -122}, {-253, -122}}, color = {0, 0, 127}));
    connect(sampler1.u, integrator.y) annotation(
      Line(points = {{-286, -84}, {-320, -84}}, color = {0, 0, 127}));
    connect(sampler1.y, mean.u) annotation(
      Line(points = {{-263, -84}, {-253, -84}}, color = {0, 0, 127}));
    connect(sampler1.y, missing_dli.u1) annotation(
      Line(points = {{-263, -84}, {-261, -84}, {-261, -110}, {-253, -110}}, color = {0, 0, 127}));
    connect(missing_dli.y, hours_of_suppl_lighting.u1) annotation(
      Line(points = {{-229, -116}, {-217, -116}, {-217, -130}, {-203, -130}}, color = {0, 0, 127}));
    connect(hourly_light_integral.y, hours_of_suppl_lighting.u2) annotation(
      Line(points = {{-229, -152}, {-207.5, -152}, {-207.5, -142}, {-202, -142}}, color = {0, 0, 127}));
    connect(led_watts.y, led_power_draw.u2) annotation(
      Line(points = {{-147, -180}, {-141, -180}, {-141, -162}, {-121, -162}}, color = {0, 0, 127}));
    connect(full_efficiency.y, led_heat_loss.u2) annotation(
      Line(points = {{-127, -264}, {-121, -264}, {-121, -248}, {-107, -248}}, color = {0, 0, 127}));
    connect(led_efficiency.y, led_heat_loss.u1) annotation(
      Line(points = {{-127, -228}, {-121, -228}, {-121, -236}, {-107, -236}}, color = {0, 0, 127}));
    connect(led_heat_loss.y, led_heat_generation.u2) annotation(
      Line(points = {{-83, -242}, {-76, -242}, {-76, -168}, {-67, -168}}, color = {0, 0, 127}));
    connect(led_power_draw.y, led_heat_generation.u1) annotation(
      Line(points = {{-97, -156}, {-67, -156}}, color = {0, 0, 127}));
    connect(led_heat_generation.y, led_prescriped_heat_flow.Q_flow) annotation(
      Line(points = {{-42, -162}, {-32, -162}}, color = {0, 0, 127}));
    connect(led_prescriped_heat_flow.port, farm_air.heatPort) annotation(
      Line(points = {{-12, -162}, {-2, -162}, {-2, -66}, {-34, -66}, {-34, -10}, {-10, -10}, {-10, -20}}, color = {191, 0, 0}));
    connect(hours_of_suppl_lighting.y, limiter.u) annotation(
      Line(points = {{-178, -136}, {-174, -136}, {-174, -140}, {-164, -140}}, color = {0, 0, 127}));
    connect(limiter.y, led_power_draw.u1) annotation(
      Line(points = {{-140, -140}, {-134, -140}, {-134, -150}, {-120, -150}}, color = {0, 0, 127}));
    connect(weaDat1.weaBus, weaBus) annotation(
      Line(points = {{-560, -74}, {-530, -74}, {-530, -30}, {-268, -30}}, color = {255, 204, 51}, thickness = 0.5));
    connect(weaDat1.weaBus, incAng.weaBus) annotation(
      Line(points = {{-560, -74}, {-530, -74}, {-530, -122}, {-500, -122}}, color = {255, 204, 51}, thickness = 0.5));
    connect(led_power_draw.y, led_cumulative_power_draw.u) annotation(
      Line(points = {{-96, -156}, {-84, -156}, {-84, -110}, {-72, -110}}, color = {0, 0, 127}));
  connect(solar_to_par.y, integrator.u) annotation(
      Line(points = {{-400, -84}, {-344, -84}}, color = {0, 0, 127}));
    annotation(
      Diagram(coordinateSystem(extent = {{-580, 100}, {220, -280}})));
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

  model led "model for the power consumption and heat generation of an led module"
  //  parameter Modelica.Units.SI.Angle azi "Surface azimuth";
  //  parameter Modelica.Units.SI.Angle til "Surface tilt";
  //  Real HDirNor(final unit="W/m2") "Direct normal solar irradiation";
    parameter Modelica.Units.SI.Power power "led power draw";
    parameter Real efficiency=0.622 "efficiency of the led";
    
    Real heat_loss "the heat loss of the led";
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b port annotation (Placement(transformation(extent={{90,
          -10},{110,10}})));
    
  equation
    heat_loss = (1 - efficiency) * power;
    port.Q_flow = heat_loss;
  end led;
  annotation(
    uses(Buildings(version = "11.0.0"), Modelica(version = "4.0.0")));
end vertical_farm;