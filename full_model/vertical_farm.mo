package vertical_farm
  model airvolume
    Buildings.HeatTransfer.Sources.FixedTemperature building_air(T = 293.15) annotation(
      Placement(transformation(origin = {170, -10}, extent = {{10, -10}, {-10, 10}})));
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
    Buildings.HeatTransfer.Sources.PrescribedTemperature preTem1 annotation(
      Placement(transformation(origin = {-196, -34}, extent = {{-10, -10}, {10, 10}})));
    Buildings.BoundaryConditions.WeatherData.Bus weaBus annotation(
      Placement(transformation(origin = {-268, -30}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {-274, -32}, extent = {{-10, -10}, {10, 10}})));

    Buildings.BoundaryConditions.WeatherData.ReaderTMY3 weaDat1(filNam = "/home/marci/Downloads/DEU_BY_Nurnberg.AP.107630_TMYx.2007-2021/DEU_BY_Nurnberg.AP.107630_TMYx.2007-2021.mos") annotation(
      Placement(transformation(origin = {-320, -102}, extent = {{-60, 20}, {-40, 40}})));
    Buildings.BoundaryConditions.SolarIrradiation.DirectTiltedSurface HDirRoo(azi (displayUnit = "deg")= 0.7330382858376184, til(displayUnit = "deg") = 1.5707963267948966) annotation(
      Placement(transformation(origin = {-322, -102}, extent = {{20, 20}, {40, 40}})));
    Buildings.BoundaryConditions.SolarIrradiation.DiffuseIsotropic HDifTilIso(til = 1.5707963267948966) annotation(
      Placement(transformation(origin = {-292, -92}, extent = {{-10, -10}, {10, 10}})));
    Buildings.BoundaryConditions.SolarGeometry.IncidenceAngle incAng(azi = 0, til = 1.5707963267948966) annotation(
      Placement(transformation(origin = {-292, -120}, extent = {{-10, -10}, {10, 10}})));
    led led_module annotation(
      Placement(transformation(origin = {-150, -88}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Continuous.Integrator led_daily_power_draw(k = 1/(60*60), use_reset = true) annotation(
      Placement(transformation(origin = {-110, -66}, extent = {{-10, -10}, {10, 10}})));
    Buildings.Airflow.Multizone.MediumColumn farm_air_column_bottom(redeclare package Medium = Buildings.Media.Air "Moist air", h = 22.5/2, densitySelection = Buildings.Airflow.Multizone.Types.densitySelection.actual) annotation(
      Placement(transformation(origin = {-20, -90}, extent = {{-10, -10}, {10, 10}})));
    Buildings.Airflow.Multizone.MediumColumn farm_air_column_top(redeclare package Medium = Buildings.Media.Air "Moist air", h = 22.5/2, densitySelection = Buildings.Airflow.Multizone.Types.densitySelection.actual) annotation(
      Placement(transformation(origin = {-20, 34}, extent = {{-10, -10}, {10, 10}})));
  Buildings.Airflow.Multizone.DoorOperable window_bottom(redeclare package Medium = Buildings.Media.Air "Moist air", wOpe = 52.8, hOpe = 0.2, LClo = 0.001)  annotation(
      Placement(transformation(origin = {-70, -150}, extent = {{-10, -10}, {10, 10}})));
  Buildings.Airflow.Multizone.DoorOperable window_top(LClo = 0.001, redeclare package Medium = Buildings.Media.Air "Moist air", hOpe = 0.2, wOpe = 52.8) annotation(
      Placement(transformation(origin = {-60, 260}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.Constant windspeed_and_direction_envelope1(k = 2) annotation(
      Placement(transformation(origin = {-276, 82}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  Buildings.Fluid.Sources.Outside outside_air_pressure(nPorts = 4, redeclare package Medium = Buildings.Media.Air "Moist air")  annotation(
      Placement(transformation(origin = {-190, -150}, extent = {{-10, -10}, {10, 10}})));
  Buildings.Fluid.MixingVolumes.MixingVolume farm_air(m_flow_nominal = 0.001, V = 594, redeclare package Medium = Buildings.Media.Air "Moist air", nPorts = 3) annotation(
      Placement(transformation(origin = {-10, -30}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  Modelica.Blocks.Sources.Constant basal_crop_coefficient(k = 0.1) annotation(
      Placement(transformation(origin = {-310, 82}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor farm_temperature_sensor annotation(
      Placement(transformation(origin = {-190, -230}, extent = {{-10, -10}, {10, 10}})));
  Buildings.Controls.Continuous.LimPID window_control annotation(
      Placement(transformation(origin = {-150, -226}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.Constant farm_temperature_setpoint(k = 273.15 + 20) annotation(
      Placement(transformation(origin = {-190, -264}, extent = {{-10, -10}, {10, 10}})));
  pump irrigation_pump annotation(
      Placement(transformation(origin = {-150, -110}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Continuous.Integrator pump_daily_power_draw(k = 1/(60*60), use_reset = true) annotation(
      Placement(transformation(origin = {-110, -110}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.SampleTrigger daily_trigger(period = 31536000)  annotation(
      Placement(transformation(origin = {-150, -132}, extent = {{-10, -10}, {10, 10}})));
  Buildings.Electrical.DC.Sources.PVSimple pVSimple(A = 1212, eta = 0.2, V_nominal = 24)  annotation(
      Placement(transformation(origin = {-290, -250}, extent = {{-10, -10}, {10, 10}})));
  Buildings.BoundaryConditions.SolarIrradiation.DirectTiltedSurface HDirRoo1(azi = Buildings.Types.Azimuth.S, til(displayUnit = "deg") = 0.5235987755982988) annotation(
      Placement(transformation(origin = {-380, -240}, extent = {{20, 20}, {40, 40}})));
  Buildings.BoundaryConditions.SolarIrradiation.DiffuseIsotropic HDifTilIso1(til = 0.5235987755982988) annotation(
      Placement(transformation(origin = {-350, -230}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Math.Add add annotation(
      Placement(transformation(origin = {-310, -220}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Electrical.Analog.Basic.Ground ground annotation(
      Placement(transformation(origin = {-340, -310}, extent = {{-92, -40}, {-72, -20}})));
  Buildings.Electrical.DC.Loads.Resistor res(R = 0.5, V_nominal = 12) annotation(
      Placement(transformation(origin = {-340, -310}, extent = {{-2, -10}, {18, 10}})));
  Buildings.Electrical.DC.Sources.ConstantVoltage sou(V = 24) annotation(
      Placement(transformation(origin = {-340, -310}, extent = {{-82, -10}, {-62, 10}})));
  Buildings.Electrical.DC.Lines.TwoPortResistance lin(R = 0.05) annotation(
      Placement(transformation(origin = {-340, -310}, extent = {{-38, 30}, {-18, 50}})));
  Buildings.Electrical.DC.Sensors.GeneralizedSensor sen annotation(
      Placement(transformation(origin = {-340, -310}, extent = {{0, 30}, {20, 50}})));
  Modelica.Blocks.Continuous.Integrator solar_daily_power_output(k = 1/(60*60), use_reset = true) annotation(
      Placement(transformation(origin = {-256, -250}, extent = {{-10, -10}, {10, 10}})));
  plant_evapotranspiration plant_evapotranspiration1(a_plant = 523.2)  annotation(
      Placement(transformation(origin = {-250, 50}, extent = {{-10, -10}, {10, 10}})));
  plant_yield plant_yield1(a_plant = 523.2)  annotation(
      Placement(transformation(origin = {-250, 130}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.Constant co2_concentration(k = 365) annotation(
      Placement(transformation(origin = {-272, 178}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  Buildings.BoundaryConditions.SolarIrradiation.DirectTiltedSurface HDirRoo11(azi = 0.7330382858376184, til(displayUnit = "deg") = 1.5707963267948966) annotation(
      Placement(transformation(origin = {-400, 102}, extent = {{20, 20}, {40, 40}})));
  Buildings.BoundaryConditions.SolarIrradiation.DiffuseIsotropic HDifTilIso11(til = 1.5707963267948966) annotation(
      Placement(transformation(origin = {-370, 112}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Math.Add add1 annotation(
      Placement(transformation(origin = {-330, 122}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Thermal.HeatTransfer.Sensors.TemperatureSensor farm_temperature_sensor1 annotation(
      Placement(transformation(origin = {-290, 138}, extent = {{-10, -10}, {10, 10}})));
  Buildings.Controls.OBC.CDL.Reals.MovingAverage farm_temp_average(delta = 86400)  annotation(
      Placement(transformation(origin = {-170, 170}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Thermal.HeatTransfer.Celsius.FromKelvin fromKelvin annotation(
      Placement(transformation(origin = {-210, 170}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Math.Gain gaiWin(k = 1236.9*0.2*0.8) annotation(
      Placement(transformation(origin = {-110, 80}, extent = {{-10, -10}, {10, 10}})));
  Buildings.HeatTransfer.Sources.PrescribedHeatFlow preHeaFlo annotation(
      Placement(transformation(origin = {-70, 80}, extent = {{-10, -10}, {10, 10}})));
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
    connect(preTem1.port, envelope_glass_outside_convection.fluid) annotation(
      Line(points = {{-186, -34}, {-160, -34}, {-160, -10}, {-130, -10}}, color = {191, 0, 0}));
    connect(weaBus.TDryBul, preTem1.T) annotation(
      Line(points = {{-268, -30}, {-232, -30}, {-232, -34}, {-208, -34}}, color = {0, 0, 127}));
//  connect(gaiWin.u, weaBus.HDirNor);
//  connect(HGloTil.HDifTil.incAng.weaBus, weaBus) annotation(
//    Line(points = {{-210, 80}, {-268, 80}, {-268, -30}}, color = {255, 204, 51}, thickness = 0.5));
    connect(envelope_glass_outside_convection.v, weaBus.winSpe) annotation(
      Line(points = {{-108, 0}, {-102, 0}, {-102, 14}, {-268, 14}, {-268, -30}}, color = {0, 0, 127}));
    connect(envelope_glass_outside_convection.dir, weaBus.winDir) annotation(
      Line(points = {{-108, -4}, {-100, -4}, {-100, 16}, {-268, 16}, {-268, -30}}, color = {0, 0, 127}));
    connect(weaDat1.weaBus, HDirRoo.weaBus) annotation(
      Line(points = {{-360, -72}, {-302, -72}}, color = {255, 204, 51}, thickness = 0.5));
    connect(HDifTilIso.weaBus, weaDat1.weaBus) annotation(
      Line(points = {{-302, -92}, {-332, -92}, {-332, -72}, {-360, -72}}, color = {255, 204, 51}, thickness = 0.5));
    connect(weaDat1.weaBus, weaBus) annotation(
      Line(points = {{-360, -72}, {-321, -72}, {-321, -30}, {-268, -30}}, color = {255, 204, 51}, thickness = 0.5));
    connect(weaDat1.weaBus, incAng.weaBus) annotation(
      Line(points = {{-360, -72}, {-332, -72}, {-332, -120}, {-302, -120}}, color = {255, 204, 51}, thickness = 0.5));
    connect(led_module.y, led_daily_power_draw.u) annotation(
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
    connect(window_bottom.port_b1, farm_air_column_bottom.port_b) annotation(
      Line(points = {{-60, -144}, {-20, -144}, {-20, -100}}, color = {0, 127, 255}));
    connect(window_bottom.port_a2, farm_air_column_bottom.port_b) annotation(
      Line(points = {{-60, -156}, {-20, -156}, {-20, -100}}, color = {0, 127, 255}));
    connect(window_top.port_b1, farm_air_column_top.port_a) annotation(
      Line(points = {{-50, 266}, {-20, 266}, {-20, 44}}, color = {0, 127, 255}));
    connect(window_top.port_a2, farm_air_column_top.port_a) annotation(
      Line(points = {{-50, 254}, {-20, 254}, {-20, 44}}, color = {0, 127, 255}));
    connect(outside_air_pressure.weaBus, weaDat1.weaBus) annotation(
      Line(points = {{-200, -150}, {-332, -150}, {-332, -72}, {-360, -72}}, color = {255, 204, 51}, thickness = 0.5));
    connect(window_bottom.port_a1, outside_air_pressure.ports[1]) annotation(
      Line(points = {{-80, -144}, {-140, -144}, {-140, -150}, {-180, -150}}, color = {0, 127, 255}));
    connect(window_bottom.port_b2, outside_air_pressure.ports[2]) annotation(
      Line(points = {{-80, -156}, {-140, -156}, {-140, -150}, {-180, -150}}, color = {0, 127, 255}));
    connect(window_top.port_a1, outside_air_pressure.ports[3]) annotation(
      Line(points = {{-70, 266}, {-402, 266}, {-402, -174}, {-180, -174}, {-180, -150}}, color = {0, 127, 255}));
    connect(window_top.port_b2, outside_air_pressure.ports[4]) annotation(
      Line(points = {{-70, 254}, {-398, 254}, {-398, -174}, {-180, -174}, {-180, -150}}, color = {0, 127, 255}));
    connect(farm_temperature_sensor.T, window_control.u_s) annotation(
      Line(points = {{-179, -230}, {-170.5, -230}, {-170.5, -226}, {-162, -226}}, color = {0, 0, 127}));
    connect(farm_temperature_setpoint.y, window_control.u_m) annotation(
      Line(points = {{-179, -264}, {-150, -264}, {-150, -238}}, color = {0, 0, 127}));
    connect(window_control.y, window_bottom.y) annotation(
      Line(points = {{-139, -226}, {-120, -226}, {-120, -150}, {-80, -150}}, color = {0, 0, 127}));
    connect(window_control.y, window_top.y) annotation(
      Line(points = {{-139, -226}, {-120, -226}, {-120, -190}, {-406, -190}, {-406, 260}, {-71, 260}}, color = {0, 0, 127}));
    connect(farm_temperature_sensor.port, farm_air.heatPort) annotation(
      Line(points = {{-200, -230}, {-220, -230}, {-220, -202}, {-30, -202}, {-30, -10}, {-10, -10}, {-10, -20}}, color = {191, 0, 0}));
    connect(irrigation_pump.weaBus, weaDat1.weaBus) annotation(
      Line(points = {{-160, -110}, {-220, -110}, {-220, -150}, {-332, -150}, {-332, -72}, {-360, -72}}, color = {255, 204, 51}, thickness = 0.5));
    connect(irrigation_pump.y, pump_daily_power_draw.u) annotation(
      Line(points = {{-138, -110}, {-122, -110}}, color = {0, 0, 127}));
    connect(daily_trigger.y, pump_daily_power_draw.reset) annotation(
      Line(points = {{-139, -132}, {-104, -132}, {-104, -122}}, color = {255, 0, 255}));
    connect(daily_trigger.y, led_daily_power_draw.reset) annotation(
      Line(points = {{-139, -132}, {-128, -132}, {-128, -88}, {-104, -88}, {-104, -78}}, color = {255, 0, 255}));
    connect(HDirRoo1.H, add.u1) annotation(
      Line(points = {{-338, -210}, {-332, -210}, {-332, -214}, {-322, -214}}, color = {0, 0, 127}));
    connect(HDifTilIso1.H, add.u2) annotation(
      Line(points = {{-338, -230}, {-332, -230}, {-332, -226}, {-322, -226}}, color = {0, 0, 127}));
    connect(HDirRoo1.weaBus, weaDat1.weaBus) annotation(
      Line(points = {{-360, -210}, {-370, -210}, {-370, -150}, {-332, -150}, {-332, -72}, {-360, -72}}, color = {255, 204, 51}, thickness = 0.5));
    connect(HDifTilIso1.weaBus, weaDat1.weaBus) annotation(
      Line(points = {{-360, -230}, {-370, -230}, {-370, -150}, {-332, -150}, {-332, -72}, {-360, -72}}, color = {255, 204, 51}, thickness = 0.5));
    connect(add.y, pVSimple.G) annotation(
      Line(points = {{-298, -220}, {-290, -220}, {-290, -238}}, color = {0, 0, 127}));
    connect(sou.terminal, res.terminal) annotation(
      Line(points = {{-402, -310}, {-342, -310}}, color = {0, 0, 255}));
    connect(lin.terminal_n, res.terminal) annotation(
      Line(points = {{-378, -270}, {-390, -270}, {-390, -310}, {-342, -310}, {-342, -310}}, color = {0, 0, 255}));
    connect(lin.terminal_p, sen.terminal_n) annotation(
      Line(points = {{-358, -270}, {-340, -270}}, color = {0, 0, 255}));
    connect(sou.n, ground.p) annotation(
      Line(points = {{-422, -310}, {-422, -330}}, color = {0, 0, 255}));
    connect(sen.terminal_p, pVSimple.terminal) annotation(
      Line(points = {{-320, -270}, {-312, -270}, {-312, -250}, {-300, -250}}, color = {0, 0, 255}));
    connect(pVSimple.P, solar_daily_power_output.u) annotation(
      Line(points = {{-278, -242}, {-276, -242}, {-276, -250}, {-268, -250}}, color = {0, 0, 127}));
    connect(daily_trigger.y, solar_daily_power_output.reset) annotation(
      Line(points = {{-138, -132}, {-128, -132}, {-128, -188}, {-232, -188}, {-232, -278}, {-250, -278}, {-250, -262}}, color = {255, 0, 255}));
    connect(plant_evapotranspiration1.port, farm_air.heatPort) annotation(
      Line(points = {{-240, 46}, {-178, 46}, {-178, 114}, {-10, 114}, {-10, -20}}, color = {191, 0, 0}));
    connect(windspeed_and_direction_envelope1.y, plant_evapotranspiration1.v) annotation(
      Line(points = {{-276, 72}, {-276, 60}, {-262, 60}}, color = {0, 0, 127}));
    connect(basal_crop_coefficient.y, plant_evapotranspiration1.k_cb) annotation(
      Line(points = {{-310, 72}, {-310, 56}, {-262, 56}}, color = {0, 0, 127}));
    connect(plant_evapotranspiration1.weaBus, weaDat1.weaBus) annotation(
      Line(points = {{-260, 50}, {-332, 50}, {-332, -72}, {-360, -72}}, color = {255, 204, 51}, thickness = 0.5));
    connect(co2_concentration.y, plant_yield1.co2_concentration) annotation(
      Line(points = {{-272, 167}, {-272, 130}, {-262, 130}}, color = {0, 0, 127}));
    connect(HDirRoo11.H, add1.u1) annotation(
      Line(points = {{-359, 132}, {-353, 132}, {-353, 128}, {-343, 128}}, color = {0, 0, 127}));
    connect(HDifTilIso11.H, add1.u2) annotation(
      Line(points = {{-359, 112}, {-353, 112}, {-353, 116}, {-343, 116}}, color = {0, 0, 127}));
    connect(add1.y, plant_yield1.par) annotation(
      Line(points = {{-318, 122}, {-262, 122}}, color = {0, 0, 127}));
    connect(HDifTilIso11.weaBus, weaDat1.weaBus) annotation(
      Line(points = {{-380, 112}, {-386, 112}, {-386, -30}, {-332, -30}, {-332, -72}, {-360, -72}}, color = {255, 204, 51}, thickness = 0.5));
    connect(HDirRoo11.weaBus, weaDat1.weaBus) annotation(
      Line(points = {{-380, 132}, {-386, 132}, {-386, -30}, {-332, -30}, {-332, -72}, {-360, -72}}, color = {255, 204, 51}, thickness = 0.5));
    connect(farm_temperature_sensor1.T, plant_yield1.air_temperature) annotation(
      Line(points = {{-278, 138}, {-262, 138}}, color = {0, 0, 127}));
    connect(farm_temperature_sensor1.port, farm_air.heatPort) annotation(
      Line(points = {{-300, 138}, {-308, 138}, {-308, 202}, {-10, 202}, {-10, -20}}, color = {191, 0, 0}));
    connect(fromKelvin.Celsius, farm_temp_average.u) annotation(
      Line(points = {{-198, 170}, {-182, 170}}, color = {0, 0, 127}));
    connect(farm_temperature_sensor1.T, fromKelvin.Kelvin) annotation(
      Line(points = {{-278, 138}, {-276, 138}, {-276, 146}, {-222, 146}, {-222, 170}}, color = {0, 0, 127}));
    connect(gaiWin.y, preHeaFlo.Q_flow) annotation(
      Line(points = {{-99, 80}, {-81, 80}}, color = {0, 0, 127}));
    connect(preHeaFlo.port, farm_air.heatPort) annotation(
      Line(points = {{-60, 80}, {-10, 80}, {-10, -20}}, color = {191, 0, 0}));
    connect(HDirRoo11.H, gaiWin.u) annotation(
      Line(points = {{-358, 132}, {-352, 132}, {-352, 142}, {-312, 142}, {-312, 102}, {-154, 102}, {-154, 80}, {-122, 80}}, color = {0, 0, 127}));
    annotation(
      Diagram(coordinateSystem(extent = {{-440, 280}, {180, -360}}), graphics = {Rectangle(origin = {-90, 35}, extent = {{-50, 63}, {50, -63}}), Rectangle(origin = {70, 59}, extent = {{-48, 35}, {48, -35}}), Rectangle(origin = {70, -79}, fillColor = {153, 193, 241}, extent = {{48, -47}, {-48, 47}})}),
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
  Buildings.HeatTransfer.Conduction.MultiLayer heaCon(A = 1, layers(nLay = 2, material = {glass, concrete}))  annotation(
      Placement(transformation(origin = {8, 0}, extent = {{-10, -10}, {10, 10}})));
  parameter Buildings.HeatTransfer.Data.Solids.Glass glass(x = 0.1)  annotation(
      Placement(transformation(origin = {-10, -50}, extent = {{-10, -10}, {10, 10}})));
  parameter Buildings.HeatTransfer.Data.Solids.Concrete concrete(x = 0.1)  annotation(
      Placement(transformation(origin = {20, -50}, extent = {{-10, -10}, {10, 10}})));
  Buildings.HeatTransfer.Sources.FixedTemperature inside(T = 293.15)  annotation(
      Placement(transformation(origin = {58, 0}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
  Buildings.HeatTransfer.Sources.FixedTemperature outside(T = 283.15)  annotation(
      Placement(transformation(origin = {-62, 0}, extent = {{-10, -10}, {10, 10}})));
  Buildings.HeatTransfer.Convection.Interior con(A = 1, conMod = Buildings.HeatTransfer.Types.InteriorConvection.Temperature, til = 1.5707963267948966)  annotation(
      Placement(transformation(origin = {-26, 0}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
  greenshell.building.facade facade1 annotation(
      Placement(transformation(origin = {0, 54}, extent = {{-10, -10}, {10, 10}})));
  equation
  connect(inside.port, heaCon.port_b) annotation(
      Line(points = {{48, 0}, {18, 0}}, color = {191, 0, 0}));
  connect(heaCon.port_a, con.solid) annotation(
      Line(points = {{-2, 0}, {-16, 0}}, color = {191, 0, 0}));
  connect(con.fluid, outside.port) annotation(
      Line(points = {{-36, 0}, {-52, 0}}, color = {191, 0, 0}));
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
  Modelica.Blocks.Sources.BooleanPulse booleanPulse(width = 0.1, period = 86400)  annotation(
      Placement(transformation(origin = {-270, -110}, extent = {{-10, -10}, {10, 10}})));
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
    connect(realToBoolean.y, xor.u1) annotation(
      Line(points = {{-216, -88}, {-176, -88}, {-176, -90}, {-162, -90}}, color = {255, 0, 255}));
  connect(xor.y, timer.u) annotation(
      Line(points = {{-138, -90}, {-116, -90}}, color = {255, 0, 255}));
  connect(booleanPulse.y, xor.u2) annotation(
      Line(points = {{-258, -110}, {-176, -110}, {-176, -98}, {-162, -98}}, color = {255, 0, 255}));
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

  model plant_evapotranspiration
    Buildings.BoundaryConditions.WeatherData.Bus weaBus annotation(
        Placement(transformation(extent = {{-110, -10}, {-90, 10}}), iconTransformation(extent = {{-120, -18}, {-80, 22}})));
    Modelica.Blocks.Interfaces.RealInput v(unit="m/s") "Wind speed" annotation (
        Placement(transformation(extent={{-140,80},{-100,120}})));
    Modelica.Blocks.Interfaces.RealInput k_cb(unit="1") "basal crop coefficient - 0.01 for the initial stage and 0.21 for the midseason and late stages" annotation (
        Placement(transformation(origin = {0, -32}, extent = {{-140, 80}, {-100, 120}}), iconTransformation(origin = {0, -50}, extent = {{-140, 80}, {-100, 120}})));
    Buildings.BoundaryConditions.SolarIrradiation.DiffuseIsotropic HDifTilIso(til = 1.5707963267948966) annotation(
        Placement(transformation(origin = {-50, 50}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Interfaces.RealOutput y annotation(
        Placement(transformation(origin = {110, 0}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {110, 40}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b port annotation(
      Placement(transformation(origin = {0, -30}, extent = {{90, -10}, {110, 10}}), iconTransformation(origin = {0, -40}, extent = {{90, -10}, {110, 10}})));
    
    parameter Modelica.Units.SI.Area a_plant "planting area";
    Real et_r(final unit = "mm") "reference evapotranspiration";
    Real et_a(final unit = "L") "actual evapotranspiration";
    Real delta(final unit = "kPa/°C") "hourly slope of the vapor pressure curve";
    Real r_n(final unit = "MJ/m2") "hourly net radiation flux";
    Real gamma(final unit = "kPa/°C") "psychrometric constant";
    Real c_n(final unit = "1") = 37 "numerator constant for the reference crop type and time step";//@lopez_mora2024
    Real e_s(final unit = "kPa") "saturation vapor pressure";
    Real e_a(final unit = "kPa") "actual vapor pressure";
    Real vpd(final unit = "kPa") "vapor pressure deficit";
    Real c_d(final unit = "1") "denominator constant for the reference crop type and time step";
    
    Real l_v(final unit = "kJ/kg") = 2260 "latent heat of vaporization for water";
    Modelica.Units.SI.HeatFlowRate q_flow "heat absorbed during evapotranspiration";
    Modelica.Units.SI.Density density_water = 1000 "density of water = 1000 kg/m3";
  
  equation
    connect(HDifTilIso.weaBus, weaBus) annotation(
      Line(points = {{-60, 50}, {-80, 50}, {-80, 0}, {-100, 0}}, color = {255, 204, 51}, thickness = 0.5));
    delta = 4098 * (0.6108 * Modelica.Math.exp((17.27 * (weaBus.TDryBul-273.15))/((weaBus.TDryBul-273.15) + 237.3))) / (((weaBus.TDryBul-273.15) + 237.3)^2);
    r_n = 500;
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
    et_r = (1/3600) * (0.408 * delta * r_n + gamma * c_n * v * (e_s - e_a) / weaBus.TDryBul) / (delta + gamma * (1 + c_d * v));
    et_a = a_plant * k_cb * et_r;
    q_flow = et_a * density_water/1000 * l_v;
    port.Q_flow = q_flow;
    y = et_a;
  end plant_evapotranspiration;

  model analytics
  equation

  end analytics;

  model pump
  
    Buildings.BoundaryConditions.WeatherData.Bus weaBus annotation(
        Placement(transformation(extent = {{-110, -10}, {-90, 10}}), iconTransformation(extent = {{-120, -18}, {-80, 22}})));
    Modelica.Blocks.Interfaces.RealOutput y annotation(
        Placement(transformation(origin = {110, 0}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {110, 0}, extent = {{-10, -10}, {10, 10}})));
  
  equation
  if (weaBus.HGloHor > 0) then
      y = 4000;
    else
      y = 0;
    end if;
  end pump;

  model plant_yield
  //  Buildings.BoundaryConditions.WeatherData.Bus weaBus annotation(
  //      Placement(transformation(extent = {{-110, -10}, {-90, 10}}), iconTransformation(extent = {{-120, -18}, {-80, 22}})));
  
  // model taken from doi.org/10.1016/S0308-521X(94)90280-1
  
    Modelica.Blocks.Interfaces.RealInput air_temperature(unit="K") "Air temperature" annotation (
        Placement(transformation(origin = {0, -80}, extent = {{-140, 80}, {-100, 120}}), iconTransformation(origin = {0, -20}, extent = {{-140, 80}, {-100, 120}})));
    Modelica.Blocks.Interfaces.RealInput co2_concentration(unit="ppm") "CO2 concentration" annotation (
        Placement(transformation(origin = {0, -120}, extent = {{-140, 80}, {-100, 120}}), iconTransformation(origin = {0, -100}, extent = {{-140, 80}, {-100, 120}})));
    Modelica.Blocks.Interfaces.RealInput par(unit="W/m2") "Photosynthetically active radiation" annotation (
        Placement(transformation(origin = {0, -160}, extent = {{-140, 80}, {-100, 120}}), iconTransformation(origin = {0, -180}, extent = {{-140, 80}, {-100, 120}})));
    Modelica.Blocks.Interfaces.RealOutput y annotation(
        Placement(transformation(origin = {110, 0}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {110, 0}, extent = {{-10, -10}, {10, 10}})));
    
    parameter Modelica.Units.SI.Area a_plant "planting area";
    
    Real x_nsdw(start=0.001, final unit = "g/m2") "state variable non-structural dry weight";
    Real x_sdw(start=0.001, final unit = "g/m2") "state variable structural dry weight";
    Real dw(final unit = "kg") "total dry weight";
    
    Real c_alpha(final unit = "1") = 0.68 "Conversion rate of CO2 to CH2O";
    Real c_beta(final unit = "1") = 0.8 "Yield factor";
    
    Real f_phot(final unit = "g/(m2s") "gross canopy photosynthesis";
    Real f_resp(final unit = "g/(m2s") "maintenance respiration";
    Real c_resp_sht(final unit = "1/s") = 3.47e-7 "Shoot maintenance respiration coefﬁcient at 25°C";
    Real c_tau(final unit = "1") = 0.14 "Ratio of root dry mass to total plant dry mass (hydroponic) @abedi2003";
    Real c_resp_rt(final unit = "1/s") = 1.16e-7 "Root maintenance respiration coefﬁcient at 25°C";
    Real c_q_10_resp(final unit = "1") = 2.0 "Sensitivity of maintenance respiration to the canopy temperature";
    
    Real c_k(final unit = "1") = 0.9 "Extinction coefﬁcient";
    Real c_lar(final unit = "m2/g") = 75e-3 "Structural leaf area ratio";
    Real f_phot_max(final unit = "g/(m2s") "gross CO2 assimilation rate for a canopy with 1 square meter of effective surface area";
    
    Real epsilon(final unit = "g/J") "light use efficiency";
    Real g_co2(final unit = "m/s") "canopy conductance to CO2 diffusion";
    Real c_omega(final unit = "g/m3") = 1.83e-3 "Density of CO2";
    Real capital_gamma(final unit = "ppm") "CO2 compensation point, accounting for the impact of the temperature on photosynthesis rate";
    
    Real c_cap_gamma(final unit = "ppm") = 40 "CO2 compensation point at 20°C";
    Real c_q_10_cap_gamma(final unit = "1") = 2.0 "Sensitivity of CO2 compensation with canopy temperature";
    
    Real c_epsilon(final unit = "g/J") = 17e-6 "Quantum use efﬁciency as energy required for a reduction of one molecule of CO2";
    
    Real g_bnd(final unit = "m/s") = 0.007 "Boundary layer conductance";
    Real g_stm(final unit = "m/s") = 0.005 "Stomatal conductance";
    Real g_car(final unit = "m/s") "Carboxylation conductance";
    
    Real r_gr(final unit = "1/s") "specific growth rate";
    Real c_gr_max(final unit = "1/s") = 5e-6 "saturation growth rate at 20°C";
    Real c_gamma(final unit = "1") = 1.0 "growth rate coefficient";
    Real c_q_10_gr(final unit = "1") = 1.6 "growth rate sensitivity to the canopy temperature";
    
  equation
    
    der(x_nsdw) = c_alpha * f_phot - r_gr * x_sdw - f_resp - ((1 - c_beta) / c_beta) * r_gr * x_sdw;
    f_resp = (c_resp_sht * (1 - c_tau) * x_sdw + c_resp_rt * c_tau * x_sdw) * c_q_10_resp ^ (((air_temperature-273.15) - 25) / 10);
    f_phot = (1 - exp(-c_k * c_lar * (1 - c_tau) * x_sdw)) * f_phot_max;
    
    f_phot_max = (epsilon * par * g_co2 * c_omega * (co2_concentration - capital_gamma)) / (epsilon * par + g_co2 * c_omega * (co2_concentration - capital_gamma));
    capital_gamma = c_cap_gamma * c_q_10_cap_gamma ^ (((air_temperature-273.15) - 20) / 10);
    epsilon = c_epsilon * (co2_concentration - capital_gamma) / (co2_concentration + 2 * capital_gamma);
    g_co2 = 1 / ((1/g_bnd) + (1/g_stm) + (1/g_car));
    g_car = -1.32e-5 * air_temperature^2 + 5.94e-4 * air_temperature - 2.64e-3;
  
    der(x_sdw) = r_gr * x_sdw;
    r_gr = c_gr_max * (x_nsdw / (c_gamma * x_sdw + x_nsdw)) * c_q_10_gr ^ (((air_temperature-273.15) - 20) / 10 );
    
    dw = a_plant * (x_nsdw + x_sdw) / 1000;
    y = dw;
  
  end plant_yield;

  model system_analysis
  plant_yield plant_yield1(a_plant = 1, x_nsdw(start = 20*0.015*0.25), x_sdw(start = 20*0.015*0.75))  annotation(
      Placement(transformation(origin = {10, 30}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.Constant const(k = 30+ 273.15)  annotation(
      Placement(transformation(origin = {-50, 50}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.Constant const1(k = 1200)  annotation(
      Placement(transformation(origin = {-50, 30}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.Constant const2(k = 300/2.02)  annotation(
      Placement(transformation(origin = {-50, 10}, extent = {{-10, -10}, {10, 10}})));
  equation
  connect(const.y, plant_yield1.air_temperature) annotation(
      Line(points = {{-38, 50}, {-20, 50}, {-20, 38}, {-2, 38}}, color = {0, 0, 127}));
  connect(const1.y, plant_yield1.co2_concentration) annotation(
      Line(points = {{-38, 30}, {-2, 30}}, color = {0, 0, 127}));
  connect(const2.y, plant_yield1.par) annotation(
      Line(points = {{-38, 10}, {-20, 10}, {-20, 22}, {-2, 22}}, color = {0, 0, 127}));
  end system_analysis;
  annotation(
    uses(Buildings(version = "11.0.0"), Modelica(version = "4.0.0")));
end vertical_farm;