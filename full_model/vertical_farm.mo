package vertical_farm
  model airvolume
    Buildings.HeatTransfer.Sources.FixedTemperature building_air(T = 293.15) annotation(
      Placement(transformation(origin = {170, -10}, extent = {{10, -10}, {-10, 10}})));
  Buildings.Fluid.MixingVolumes.MixingVolume farm_air(m_flow_nominal = 0.001, V = 594, redeclare package Medium = Buildings.Media.Air "Moist air", nPorts = 2) annotation(
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
      Placement(transformation(origin = {-150, -88}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Continuous.Integrator led_cumulative_power_draw(k = 1/(60*60)) annotation(
      Placement(transformation(origin = {-110, -66}, extent = {{-10, -10}, {10, 10}})));
    Buildings.Fluid.Sources.Outside_CpLowRise outside_air_pressure(nPorts = 4, redeclare package Medium = Buildings.Media.Air "Moist air", azi = 0, s = 5) annotation(
      Placement(transformation(origin = {-190, -150}, extent = {{-10, -10}, {10, 10}})));
    Buildings.Airflow.Multizone.MediumColumn farm_air_column_bottom(redeclare package Medium = Buildings.Media.Air "Moist air", h = 22.5/2, densitySelection = Buildings.Airflow.Multizone.Types.densitySelection.actual) annotation(
      Placement(transformation(origin = {-20, -90}, extent = {{-10, -10}, {10, 10}})));
    Buildings.Airflow.Multizone.MediumColumn farm_air_column_top(redeclare package Medium = Buildings.Media.Air "Moist air", h = 22.5/2, densitySelection = Buildings.Airflow.Multizone.Types.densitySelection.actual) annotation(
      Placement(transformation(origin = {-20, 34}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Sources.Constant window_bottom_opening_percent(k = 0) annotation(
      Placement(transformation(origin = {-204, -220}, extent = {{-10, -10}, {10, 10}})));
  Buildings.Airflow.Multizone.DoorOperable window_bottom(redeclare package Medium = Buildings.Media.Air "Moist air", wOpe = 52.8, hOpe = 0.2, LClo = 0.001)  annotation(
      Placement(transformation(origin = {-70, -150}, extent = {{-10, -10}, {10, 10}})));
  Buildings.Airflow.Multizone.DoorOperable window_top(LClo = 0.001, redeclare package Medium = Buildings.Media.Air "Moist air", hOpe = 0.2, wOpe = 52.8) annotation(
      Placement(transformation(origin = {-50, 130}, extent = {{-10, -10}, {10, 10}})));
  plant plant1 annotation(
      Placement(transformation(origin = {-250, 50}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.Constant windspeed_and_direction_envelope1(k = 0) annotation(
      Placement(transformation(origin = {-276, 82}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
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
      Line(points = {{-139, -84}, {-132, -84}, {-132, -66}, {-122, -66}}, color = {0, 0, 127}));
    connect(led_module.port, farm_air.heatPort) annotation(
      Line(points = {{-140, -92}, {-140, -93}, {-118, -93}, {-118, -94}, {-30, -94}, {-30, -10}, {-10, -10}, {-10, -20}}, color = {191, 0, 0}));
    connect(HDifTilIso.H, led_module.h_dif) annotation(
      Line(points = {{-280, -92}, {-162, -92}}, color = {0, 0, 127}));
    connect(HDirRoo.H, led_module.h_dir) annotation(
      Line(points = {{-280, -72}, {-216, -72}, {-216, -84}, {-162, -84}}, color = {0, 0, 127}));
    connect(farm_air_column_bottom.port_a, farm_air.ports[1]) annotation(
      Line(points = {{-20, -80}, {-20, -30}}, color = {0, 127, 255}));
    connect(farm_air.ports[2], farm_air_column_top.port_b) annotation(
      Line(points = {{-20, -30}, {-20, 24}}, color = {0, 127, 255}));
    connect(outside_air_pressure.weaBus, weaDat1.weaBus) annotation(
      Line(points = {{-200, -150}, {-332, -150}, {-332, -72}, {-350, -72}}, color = {255, 204, 51}, thickness = 0.5));
    connect(window_bottom.port_b1, farm_air_column_bottom.port_b) annotation(
      Line(points = {{-60, -144}, {-20, -144}, {-20, -100}}, color = {0, 127, 255}));
    connect(window_bottom.port_a2, farm_air_column_bottom.port_b) annotation(
      Line(points = {{-60, -156}, {-20, -156}, {-20, -100}}, color = {0, 127, 255}));
    connect(window_bottom.y, window_bottom_opening_percent.y) annotation(
      Line(points = {{-80, -150}, {-119.5, -150}, {-119.5, -220}, {-193, -220}}, color = {0, 0, 127}));
    connect(window_bottom.port_a1, outside_air_pressure.ports[1]) annotation(
      Line(points = {{-80, -144}, {-160, -144}, {-160, -150}, {-180, -150}}, color = {0, 127, 255}));
    connect(window_bottom.port_b2, outside_air_pressure.ports[2]) annotation(
      Line(points = {{-80, -156}, {-160, -156}, {-160, -150}, {-180, -150}}, color = {0, 127, 255}));
    connect(window_top.port_b1, farm_air_column_top.port_a) annotation(
      Line(points = {{-40, 136}, {-20, 136}, {-20, 44}}, color = {0, 127, 255}));
    connect(window_top.port_a2, farm_air_column_top.port_a) annotation(
      Line(points = {{-40, 124}, {-20, 124}, {-20, 44}}, color = {0, 127, 255}));
    connect(window_top.y, window_bottom_opening_percent.y) annotation(
      Line(points = {{-60, 130}, {-406, 130}, {-406, -250}, {-178, -250}, {-178, -220}, {-192, -220}}, color = {0, 0, 127}));
    connect(window_top.port_b2, outside_air_pressure.ports[3]) annotation(
      Line(points = {{-60, 124}, {-396, 124}, {-396, -180}, {-160, -180}, {-160, -150}, {-180, -150}}, color = {0, 127, 255}));
    connect(window_top.port_a1, outside_air_pressure.ports[4]) annotation(
      Line(points = {{-60, 136}, {-396, 136}, {-396, -180}, {-160, -180}, {-160, -150}, {-180, -150}}, color = {0, 127, 255}));
    connect(plant1.weaBus, weaDat1.weaBus) annotation(
      Line(points = {{-260, 50}, {-340, 50}, {-340, -72}, {-350, -72}}, color = {255, 204, 51}, thickness = 0.5));
    connect(windspeed_and_direction_envelope1.y, plant1.v) annotation(
      Line(points = {{-276, 71}, {-276, 60}, {-262, 60}}, color = {0, 0, 127}));
    annotation(
      Diagram(coordinateSystem(extent = {{-380, 100}, {200, -140}}), graphics = {Rectangle(origin = {-90, 35}, extent = {{-50, 63}, {50, -63}}), Rectangle(origin = {70, 59}, extent = {{-48, 35}, {48, -35}}), Rectangle(origin = {70, -79}, extent = {{48, -47}, {-48, 47}})}),
      experiment(StartTime = 0, StopTime = 3.07584e+07, Tolerance = 1e-06, Interval = 307.584));
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
    parameter Modelica.Units.SI.Power power = 77.76 "led power draw";
    parameter Modelica.Units.SI.Area agricultural_area = 523.2 "planting area";
    parameter Real efficiency = 0.622 "efficiency of the led";
    parameter Real opti_dli(final unit = "mol/(m2*day)") = 14.4 "optimal dli for the chosen plant";
    parameter Real ppfd_led(final unit = "µmol/(m2*s)") = 218.28 "ppfd of the led module";
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b port annotation(
      Placement(transformation(origin = {0, -30}, extent = {{90, -10}, {110, 10}}), iconTransformation(origin = {0, -40}, extent = {{90, -10}, {110, 10}})));
    Modelica.Blocks.Math.Add total_solar_irradiance annotation(
      Placement(transformation(origin = {-310, 2}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Math.Gain solar_to_par(k = 2.02) annotation(
      Placement(transformation(origin = {-280, 2}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Continuous.Integrator dli(use_reset = true, k = 1/1000000) annotation(
      Placement(transformation(origin = {-248, 2}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Sources.SampleTrigger daily_trigger(period = 86400) annotation(
      Placement(transformation(origin = {-270, -30}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Discrete.Sampler dli_end_of_day(samplePeriod = 86400, startTime = 86400 - 1) annotation(
      Placement(transformation(origin = {-216, 2}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Math.Add missing_dli(k1 = -1) annotation(
      Placement(transformation(origin = {-180, 2}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Sources.Constant optimal_dli(k = opti_dli) annotation(
      Placement(transformation(origin = {-216, -32}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Logical.Timer timer annotation(
      Placement(transformation(origin = {-104, -90}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Logical.OnOffController onOffController(bandwidth = 0) annotation(
      Placement(transformation(origin = {-70, 2}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Sources.Constant ppfd(k = ppfd_led) annotation(
      Placement(transformation(origin = {-180, -50}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Math.Gain gain(k = 1/1000000) annotation(
      Placement(transformation(origin = {-150, -50}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Math.Division time_to_achieve_opt_dli annotation(
      Placement(transformation(origin = {-116, 2}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Math.BooleanToReal booleanToReal(realTrue = power*agricultural_area) annotation(
      Placement(transformation(origin = {-34, 2}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Math.Gain heat_loss_efficiency(k = (1 - efficiency)) annotation(
      Placement(transformation(origin = {10, -50}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Interfaces.RealOutput y annotation(
      Placement(transformation(origin = {110, 30}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {110, 40}, extent = {{-10, -10}, {10, 10}})));
    Buildings.HeatTransfer.Sources.PrescribedHeatFlow preHeaFlo annotation(
      Placement(transformation(origin = {50, -50}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Interfaces.RealInput h_dir annotation(
      Placement(transformation(origin = {-400, 30}, extent = {{-20, -20}, {20, 20}}), iconTransformation(origin = {-120, 40}, extent = {{-20, -20}, {20, 20}})));
    Modelica.Blocks.Interfaces.RealInput h_dif annotation(
      Placement(transformation(origin = {-400, -28}, extent = {{-20, -20}, {20, 20}}), iconTransformation(origin = {-120, -40}, extent = {{-20, -20}, {20, 20}})));
    Modelica.Blocks.Nonlinear.Limiter limiter(uMax = Modelica.Constants.inf, uMin = 0) annotation(
      Placement(transformation(origin = {-150, 2}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Math.RealToBoolean realToBoolean(threshold = 0.0001) annotation(
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

  model orifice_operable
    //equation
    extends Buildings.Airflow.Multizone.Coefficient_V_flow(m = 0.5, final C = CD*A*sqrt(2.0/rho_default));
    parameter Modelica.Units.SI.Area A "Area of orifice" annotation(
      Dialog(group = "Orifice characteristics"));
    parameter Real CD = 0.65 "Discharge coefficient" annotation(
      Dialog(group = "Orifice characteristics"));
    Modelica.Blocks.Interfaces.RealInput y(min = 0.0001, max = 1, unit = "1") "Opening signal, 0.0001=closed, 1=open" annotation(
      Placement(transformation(origin = {0, 120}, extent = {{-20, -20}, {20, 20}}, rotation = -90), iconTransformation(origin = {0, 120}, extent = {{-20, -20}, {20, 20}}, rotation = -90)));
    Modelica.Units.SI.Velocity v(nominal = 1) = V_flow/(A*y) "Average velocity";
    annotation(
      Icon(graphics = {Rectangle(extent = {{-100, 8}, {100, -8}}, lineColor = {0, 0, 255}, pattern = LinePattern.None, fillColor = {0, 127, 0}, fillPattern = FillPattern.Solid), Rectangle(extent = {{-20, 100}, {20, 20}}, lineColor = {0, 0, 255}, pattern = LinePattern.None, fillColor = {0, 0, 0}, fillPattern = FillPattern.Solid), Rectangle(extent = {{-20, -20}, {20, -100}}, lineColor = {0, 0, 255}, pattern = LinePattern.None, fillColor = {0, 0, 0}, fillPattern = FillPattern.Solid), Text(extent = {{24, -24}, {96, -100}}, textColor = {0, 0, 255}, fillColor = {0, 127, 0}, fillPattern = FillPattern.Solid, textString = "A=%A")}),
      defaultComponentName = "ori",
      Documentation(info = "<html>
  <p>
  This model describes the mass flow rate and pressure difference relation
  of an orifice in the form
  </p>
  <p align=\"center\" style=\"font-style:italic;\">
  V&#775; = C &Delta;p<sup>m</sup>,
  </p>
  <p>
  where
  <i>V&#775;</i> is the volume flow rate,
  <i>C</i> is a flow coefficient and
  <i>m</i> is the flow exponent.
  The flow coefficient is
  </p>
  <p align=\"center\" style=\"font-style:italic;\">
  C = C<sub>D</sub> A (2/&rho;<sub>0</sub>)<sup>0.5</sup>,
  </p>
  <p>
  where
  <i>C<sub>D</sub></i> is the discharge coefficient,
  <i>A</i> is the cross section area and
  <i>&rho;<sub>0</sub></i> is the mass density at the medium default pressure, temperature and humidity.
  </p>
  <p>
  For turbulent flow, set <i>m=1/2</i> and
  for laminar flow, set <i>m=1</i>.
  Large openings are characterized by values close to <i>0.5</i>,
  while values near <i>0.65</i> have been found for small
  crack-like openings (Dols and Walton, 2002).
  </p>
  <h4>References</h4>
  <ul>
  <li>
  W. Stuart Dols and George N. Walton, <i>CONTAMW 2.0 User Manual,
  Multizone Airflow and Contaminant Transport Analysis Software</i>,
  Building and Fire Research Laboratory,
  National Institute of Standards and Technology,
  Tech. Report NISTIR 6921,
  November, 2002.
  </li>
  <li>Michael Wetter.
  <a href=\"modelica://Buildings/Resources/Images/Airflow/Multizone/Wetter-airflow-2006.pdf\">Multizone Airflow Model in Modelica.</a>
  Proc. of the 5th International Modelica Conference, p. 431-440. Vienna, Austria, September 2006.
  </li>
  </ul>
  </html>", revisions = "<html>
  <ul>
  <li>
  February 2, 2022, by Michael Wetter:<br/>
  Revised implementation.<br/>
  This is for
  <a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/1436\">IBPSA, #1436</a>.
  </li>
  <li>
  Apr 6, 2021, by Klaas De Jonge:<br/>
  Changes due to changes in the baseclass, velocity is now a top-level variable.
  </li>
  <li>
  June 27, 2018, by Michael Wetter:<br/>
  Corrected old parameter annotation.
  </li>
  <li>
  June 24, 2018, by Michael Wetter:<br/>
  Removed parameter <code>lWet</code> as it is only used to compute
  the Reynolds number, and the Reynolds number is not used by this model.
  Also removed the variable <code>Re</code> for the Reynolds number.<br/>
  This change is non-backward compatible.<br/>
  This is for
  <a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/932\">IBPSA, #932</a>.
  </li>
  <li>
  May 30, 2018, by Michael Wetter:<br/>
  Improved documentation for
  <a href=\"https://github.com/ibpsa/modelica-ibpsa/issues/546\">IBPSA, #546</a>.
  </li>
  <li>
  October 8, 2013 by Michael Wetter:<br/>
  Changed the parameter <code>useConstantDensity</code> to
  <code>useDefaultProperties</code> to use consistent names within this package.
  A conversion script can be used to update this parameter.
  </li>
  <li>
  December 6, 2011 by Michael Wetter:<br/>
  Replaced <code>rho</code> with <code>rho_nominal</code> because
  <code>rho</code> is computed in an <code>equation</code> section and not
  in the <code>initial equation</code> section.
  </li>
  <li>
  July 20, 2010 by Michael Wetter:<br/>
  Migrated model to Modelica 3.1 and integrated it into the Buildings library.
  </li>
  <li>
  February 4, 2005 by Michael Wetter:<br/>
  Released first version.
  </li>
  </ul>
  </html>"));
  end orifice_operable;

  model plant
  //Modelica.Blocks.Math.Add total_solar_irradiance annotation(
  //    Placement(transformation(origin = {-310, 2}, extent = {{-10, -10}, {10, 10}})));
  //Modelica.Blocks.Interfaces.RealInput h_dir annotation(
  //    Placement(transformation(origin = {-400, 30}, extent = {{-20, -20}, {20, 20}}), iconTransformation(origin = {-120, 40}, extent = {{-20, -20}, {20, 20}})));
  //Modelica.Blocks.Interfaces.RealInput h_dif annotation(
  //    Placement(transformation(origin = {-400, -28}, extent = {{-20, -20}, {20, 20}}), iconTransformation(origin = {-120, -40}, extent = {{-20, -20}, {20, 20}})));
  //Modelica.Blocks.Interfaces.RealInput pAtm annotation(
  //    Placement(transformation(origin = {-400, -90}, extent = {{-20, -20}, {20, 20}}), iconTransformation(origin = {-120, -82}, extent = {{-20, -20}, {20, 20}})));
  Buildings.BoundaryConditions.WeatherData.Bus weaBus annotation(
      Placement(transformation(extent = {{-110, -10}, {-90, 10}}), iconTransformation(extent = {{-120, -18}, {-80, 22}})));
  Modelica.Blocks.Interfaces.RealInput v(unit="m/s") "Wind speed" annotation (
      Placement(transformation(extent={{-140,80},{-100,120}})));
  Buildings.BoundaryConditions.SolarIrradiation.DiffuseIsotropic HDifTilIso(til = 1.5707963267948966) annotation(
      Placement(transformation(origin = {-50, 50}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Interfaces.RealOutput y annotation(
      Placement(transformation(origin = {110, 0}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {110, 40}, extent = {{-10, -10}, {10, 10}})));
  Real et_r(final unit = "mm/h") "reference evapotranspiration";
  Real delta(final unit = "kPa/°C") "hourly slope of the vapor pressure curve";
  Real r_n(final unit = "MJ/(m2*h)") "hourly net radiation flux";
  Real gamma(final unit = "kPa/°C") "psychrometric constant";
  Real c_n(final unit = "1") = 37 "numerator constant for the reference crop type and time step";//@lopez_mora2024
  Real e_s(final unit = "kPa") "saturation vapor pressure";
  Real e_a(final unit = "kPa") "actual vapor pressure";
  Real vpd(final unit = "kPa") "vapor pressure deficit";
  Real c_d(final unit = "1") "denominator constant for the reference crop type and time step";
  
  equation
//connect(h_dir, total_solar_irradiance.u1) annotation(
//    Line(points = {{-400, 30}, {-340, 30}, {-340, 8}, {-322, 8}}, color = {0, 0, 127}));
//connect(h_dif, total_solar_irradiance.u2) annotation(
//    Line(points = {{-400, -28}, {-340, -28}, {-340, -4}, {-322, -4}}, color = {0, 0, 127}));
    connect(HDifTilIso.weaBus, weaBus) annotation(
      Line(points = {{-60, 50}, {-80, 50}, {-80, 0}, {-100, 0}}, color = {255, 204, 51}, thickness = 0.5));
  delta = 4098 * (0.6108 * Modelica.Math.exp((17.27 * (weaBus.TDryBul-273.15))/((weaBus.TDryBul-273.15) + 237.3))) / (((weaBus.TDryBul-273.15) + 237.3)^2);
  r_n = 200;
  gamma = 0.000665 * weaBus.pAtm/1000;
  e_s = 0.6108 * Modelica.Math.exp((17.27 * (weaBus.TDryBul-273.15))/((weaBus.TDryBul-273.15) + 237.3));
  e_a = 0.6108 * Modelica.Math.exp((17.27 * (weaBus.TDewPoi-273.15))/((weaBus.TDewPoi-273.15) + 237.3));
  vpd = e_s - e_a;
  //c_d = 0.24;
  if (weaBus.HGloHor == 0) then
    c_d = 0.96;
  else
    c_d = 0.24;
  end if;
  et_r = (0.408 * delta * r_n + gamma * c_n * v * (e_s - e_a) / weaBus.TDryBul) / (delta + gamma * (1 + c_d * v));
  y = et_r;
  end plant;
  annotation(
    uses(Buildings(version = "11.0.0"), Modelica(version = "4.0.0")));
end vertical_farm;