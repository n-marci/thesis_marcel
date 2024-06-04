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
      Placement(transformation(origin = {50, -106}, extent = {{-10, -10}, {10, 10}})));
    parameter Buildings.HeatTransfer.Data.Solids.InsulationBoard eps(x = 0.05, k = 0.035, c = 1200, d = 21) annotation(
      Placement(transformation(origin = {90, -106}, extent = {{-10, -10}, {10, 10}})));
    parameter Buildings.HeatTransfer.Data.Solids.Glass glass(x(displayUnit = "mm") = 0.025) annotation(
      Placement(transformation(origin = {70, 42}, extent = {{-10, -10}, {10, 10}})));
    parameter Buildings.HeatTransfer.Data.OpaqueConstructions.Generic facade_material(nLay = 2, material = {concrete, eps}) annotation(
      Placement(transformation(origin = {70, -78}, extent = {{-10, -10}, {10, 10}})));
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
      Placement(transformation(origin = {-310, -102}, extent = {{-60, 20}, {-40, 40}})));
    Buildings.BoundaryConditions.SolarIrradiation.DirectTiltedSurface HDirRoo(azi = Buildings.Types.Azimuth.S, til(displayUnit = "rad") = 90) annotation(
      Placement(transformation(origin = {-322, -102}, extent = {{20, 20}, {40, 40}})));
    Buildings.BoundaryConditions.SolarIrradiation.DiffuseIsotropic HDifTilIso(til = 1.5707963267948966) annotation(
      Placement(transformation(origin = {-292, -92}, extent = {{-10, -10}, {10, 10}})));
    Buildings.BoundaryConditions.SolarGeometry.IncidenceAngle incAng(azi = 0, til = 1.5707963267948966) annotation(
      Placement(transformation(origin = {-292, -120}, extent = {{-10, -10}, {10, 10}})));
  led led_module annotation(
      Placement(transformation(origin = {-150, -90}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Continuous.Integrator led_cumulative_power_draw(k = 1/(60*60))  annotation(
      Placement(transformation(origin = {-110, -66}, extent = {{-10, -10}, {10, 10}})));
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
      Line(points = {{-350, -72}, {-302, -72}}, color = {255, 204, 51}, thickness = 0.5));
    connect(HDifTilIso.weaBus, weaDat1.weaBus) annotation(
      Line(points = {{-302, -92}, {-332, -92}, {-332, -72}, {-350, -72}}, color = {255, 204, 51}, thickness = 0.5));
    connect(weaDat1.weaBus, weaBus) annotation(
      Line(points = {{-350, -72}, {-321, -72}, {-321, -30}, {-268, -30}}, color = {255, 204, 51}, thickness = 0.5));
    connect(weaDat1.weaBus, incAng.weaBus) annotation(
      Line(points = {{-350, -72}, {-332, -72}, {-332, -120}, {-302, -120}}, color = {255, 204, 51}, thickness = 0.5));
  connect(led_module.y, led_cumulative_power_draw.u) annotation(
      Line(points = {{-138, -86}, {-132, -86}, {-132, -66}, {-122, -66}}, color = {0, 0, 127}));
  connect(led_module.port, farm_air.heatPort) annotation(
      Line(points = {{-140, -94}, {-30, -94}, {-30, -10}, {-10, -10}, {-10, -20}}, color = {191, 0, 0}));
  connect(HDifTilIso.H, led_module.h_dif) annotation(
      Line(points = {{-280, -92}, {-194, -92}, {-194, -94}, {-162, -94}}, color = {0, 0, 127}));
  connect(HDirRoo.H, led_module.h_dir) annotation(
      Line(points = {{-280, -72}, {-216, -72}, {-216, -86}, {-162, -86}}, color = {0, 0, 127}));
    annotation(
      Diagram(coordinateSystem(extent = {{-380, 100}, {200, -140}}), graphics = {Rectangle(origin = {-90, 35}, extent = {{-50, 63}, {50, -63}}), Rectangle(origin = {70, 59}, extent = {{-48, 35}, {48, -35}}), Rectangle(origin = {70, -79}, extent = {{48, -47}, {-48, 47}})}));
  end airvolume;

  model testing3
    parameter Modelica.Units.SI.Area A = 1 "Window surface area";
    parameter Real fFra = 0.2 "Fraction of frame, = frame area divided by total area";
    final parameter Modelica.Units.SI.Area AFra = fFra*A "Frame area";
    final parameter Modelica.Units.SI.Area AGla = A - AFra "Glass area";
    parameter Boolean linearize = false "Set to true to linearize emissive power";
    parameter Modelica.Units.SI.Angle azi = 0 "Surface azimuth";
    parameter Modelica.Units.SI.Angle til = 1.5707963267949 "Surface tilt";
    Buildings.BoundaryConditions.WeatherData.ReaderTMY3 weaDat(filNam = Modelica.Utilities.Files.loadResource("modelica://Buildings/Resources/weatherdata/USA_IL_Chicago-OHare.Intl.AP.725300_TMY3.mos")) annotation(
      Placement(transformation(origin = {-70, 96}, extent = {{-10, -10}, {10, 10}})));
    Buildings.HeatTransfer.Windows.SideFins fin(h = 0.2, dep = 0.5, gap = 0.1, hWin = 1.0, wWin = 1.0) annotation(
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
    parameter Real rho = 0.2 "Ground reflectance";
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
    parameter Modelica.Units.SI.Power power=77.76 "led power draw";
    parameter Modelica.Units.SI.Area agricultural_area = 523.2 "planting area";
    parameter Real efficiency = 0.622 "efficiency of the led";
    parameter Real opti_dli(final unit="mol/(m2*day)") = 14.4 "optimal dli for the chosen plant";
    parameter Real ppfd_led(final unit="µmol/(m2*s)") = 218.28 "ppfd of the led module";
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b port annotation(
      Placement(transformation(origin = {0, -30}, extent = {{90, -10}, {110, 10}}), iconTransformation(origin = {0, -40}, extent = {{90, -10}, {110, 10}})));
  Modelica.Blocks.Math.Add total_solar_irradiance annotation(
      Placement(transformation(origin = {-310, 2}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Math.Gain solar_to_par(k = 2.02)  annotation(
      Placement(transformation(origin = {-280, 2}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Continuous.Integrator dli(use_reset = true, k = 1/1000000)  annotation(
      Placement(transformation(origin = {-248, 2}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.SampleTrigger daily_trigger(period = 86400)  annotation(
      Placement(transformation(origin = {-270, -30}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Discrete.Sampler dli_end_of_day(samplePeriod = 86400, startTime = 86400 - 1)  annotation(
      Placement(transformation(origin = {-216, 2}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Math.Add missing_dli(k1 = -1)  annotation(
      Placement(transformation(origin = {-180, 2}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.Constant optimal_dli(k = opti_dli)  annotation(
      Placement(transformation(origin = {-216, -32}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Logical.Timer timer annotation(
      Placement(transformation(origin = {-104, -90}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Logical.OnOffController onOffController(bandwidth = 0)  annotation(
      Placement(transformation(origin = {-70, 2}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.Constant ppfd(k = ppfd_led)  annotation(
      Placement(transformation(origin = {-180, -50}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Math.Gain gain(k = 1/1000000)  annotation(
      Placement(transformation(origin = {-150, -50}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Math.Division time_to_achieve_opt_dli annotation(
      Placement(transformation(origin = {-116, 2}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Math.BooleanToReal booleanToReal(realTrue = power*agricultural_area)  annotation(
      Placement(transformation(origin = {-34, 2}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Math.Gain heat_loss_efficiency(k = (1 - efficiency))  annotation(
      Placement(transformation(origin = {10, -50}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Interfaces.RealOutput y annotation(
      Placement(transformation(origin = {110, 30}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {110, 40}, extent = {{-10, -10}, {10, 10}})));
  Buildings.HeatTransfer.Sources.PrescribedHeatFlow preHeaFlo annotation(
      Placement(transformation(origin = {50, -50}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Interfaces.RealInput h_dir annotation(
      Placement(transformation(origin = {-400, 30}, extent = {{-20, -20}, {20, 20}}), iconTransformation(origin = {-120, 40}, extent = {{-20, -20}, {20, 20}})));
  Modelica.Blocks.Interfaces.RealInput h_dif annotation(
      Placement(transformation(origin = {-400, -28}, extent = {{-20, -20}, {20, 20}}), iconTransformation(origin = {-120, -40}, extent = {{-20, -20}, {20, 20}})));
  Modelica.Blocks.Nonlinear.Limiter limiter(uMax = Modelica.Constants.inf, uMin = 0)  annotation(
      Placement(transformation(origin = {-150, 2}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Math.RealToBoolean realToBoolean(threshold = 0.0001)  annotation(
      Placement(transformation(origin = {-228, -88}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Logical.Xor xor annotation(
      Placement(transformation(origin = {-150, -90}, extent = {{-10, -10}, {10, 10}})));
  equation
    //heat_loss = (1 - efficiency)*power;
    //port.Q_flow = heat_loss;
    connect(total_solar_irradiance.y, solar_to_par.u) annotation(
      Line(points = {{-298, 2}, {-292, 2}}, color = {0, 0, 127}));
    connect(solar_to_par.y, dli.u) annotation(
      Line(points = {{-268, 2}, {-260, 2}}, color = {0, 0, 127}));
    connect(daily_trigger.y, dli.reset) annotation(
      Line(points = {{-258, -30}, {-242, -30}, {-242, -10}}, color = {255, 0, 255}));
    connect(dli.y, dli_end_of_day.u) annotation(
      Line(points = {{-236, 2}, {-228, 2}}, color = {0, 0, 127}));
    connect(optimal_dli.y, missing_dli.u2) annotation(
      Line(points = {{-204, -32}, {-200, -32}, {-200, -4}, {-192, -4}}, color = {0, 0, 127}));
    connect(dli_end_of_day.y, missing_dli.u1) annotation(
      Line(points = {{-204, 2}, {-200, 2}, {-200, 8}, {-192, 8}}, color = {0, 0, 127}));
    connect(timer.y, onOffController.u) annotation(
      Line(points = {{-93, -90}, {-88.5, -90}, {-88.5, -4}, {-82, -4}}, color = {0, 0, 127}));
    connect(ppfd.y, gain.u) annotation(
      Line(points = {{-169, -50}, {-162, -50}}, color = {0, 0, 127}));
    connect(gain.y, time_to_achieve_opt_dli.u2) annotation(
      Line(points = {{-139, -50}, {-132, -50}, {-132, -4}, {-128, -4}}, color = {0, 0, 127}));
    connect(time_to_achieve_opt_dli.y, onOffController.reference) annotation(
      Line(points = {{-104, 2}, {-92, 2}, {-92, 8}, {-82, 8}}, color = {0, 0, 127}));
    connect(onOffController.y, booleanToReal.u) annotation(
      Line(points = {{-58, 2}, {-46, 2}}, color = {255, 0, 255}));
    connect(booleanToReal.y, y) annotation(
      Line(points = {{-22, 2}, {52, 2}, {52, 30}, {110, 30}}, color = {0, 0, 127}));
    connect(booleanToReal.y, heat_loss_efficiency.u) annotation(
      Line(points = {{-22, 2}, {-16, 2}, {-16, -50}, {-2, -50}}, color = {0, 0, 127}));
    connect(heat_loss_efficiency.y, preHeaFlo.Q_flow) annotation(
      Line(points = {{22, -50}, {40, -50}}, color = {0, 0, 127}));
    connect(preHeaFlo.port, port) annotation(
      Line(points = {{60, -50}, {78, -50}, {78, -30}, {100, -30}}, color = {191, 0, 0}));
  connect(h_dir, total_solar_irradiance.u1) annotation(
      Line(points = {{-400, 30}, {-360, 30}, {-360, 8}, {-322, 8}}, color = {0, 0, 127}));
  connect(h_dif, total_solar_irradiance.u2) annotation(
      Line(points = {{-400, -28}, {-360, -28}, {-360, -4}, {-322, -4}}, color = {0, 0, 127}));
  connect(missing_dli.y, limiter.u) annotation(
      Line(points = {{-168, 2}, {-162, 2}}, color = {0, 0, 127}));
  connect(limiter.y, time_to_achieve_opt_dli.u1) annotation(
      Line(points = {{-138, 2}, {-136, 2}, {-136, 8}, {-128, 8}}, color = {0, 0, 127}));
  connect(limiter.y, realToBoolean.u) annotation(
      Line(points = {{-138, 2}, {-136, 2}, {-136, -68}, {-252, -68}, {-252, -88}, {-240, -88}}, color = {0, 0, 127}));
  connect(daily_trigger.y, xor.u2) annotation(
      Line(points = {{-258, -30}, {-256, -30}, {-256, -106}, {-209, -106}, {-209, -98}, {-162, -98}}, color = {255, 0, 255}));
  connect(realToBoolean.y, xor.u1) annotation(
      Line(points = {{-216, -88}, {-176, -88}, {-176, -90}, {-162, -90}}, color = {255, 0, 255}));
  connect(xor.y, timer.u) annotation(
      Line(points = {{-138, -90}, {-116, -90}}, color = {255, 0, 255}));
    annotation(
      Diagram(coordinateSystem(extent = {{-420, 60}, {120, -120}})),
  Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}})));
end led;
  annotation(
    uses(Buildings(version = "11.0.0"), Modelica(version = "4.0.0")));
end vertical_farm;