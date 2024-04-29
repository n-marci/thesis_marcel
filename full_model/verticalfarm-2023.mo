model verticalfarm
  Greenhouses.Components.Greenhouse.Floor plant_panels(A = surface.k, V = 100, c_p = 1300, rho = 1.35e6) annotation(
    Placement(transformation(origin = {-188, -80}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  Greenhouses.Components.Greenhouse.Canopy canopy(A = surface.k, LAI = 0.5) annotation(
    Placement(transformation(origin = {-146, -74}, extent = {{10, -10}, {-10, 10}})));
  Greenhouses.Components.Greenhouse.Illumination illumination(A = surface.k, P_el = 7000, p_el = 70, power_input = false) annotation(
    Placement(transformation(origin = {-90, -90}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  Greenhouses.Components.Greenhouse.Air air(A = surface.k) annotation(
    Placement(transformation(origin = {-96, -34}, extent = {{10, -10}, {-10, 10}})));
  Greenhouses.Components.Greenhouse.Cover cover(A = surface.k, phi = 1.5707963267948966) annotation(
    Placement(transformation(origin = {-80, 70}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  Greenhouses.Flows.HeatTransfer.Radiation_T4 radiation_T4(A = surface.k, epsilon_a = 0.84, epsilon_b = 1) annotation(
    Placement(transformation(origin = {-30, 90}, extent = {{-10, -10}, {10, 10}})));
  Greenhouses.Flows.HeatTransfer.OutsideAirConvection outsideAirConvection(A = surface.k, phi = 1.5707963267948966) annotation(
    Placement(transformation(origin = {-30, 70}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Thermal.HeatTransfer.Celsius.PrescribedTemperature prescribedTemperature annotation(
    Placement(transformation(origin = {28, 70}, extent = {{8, -8}, {-8, 8}})));
  Modelica.Thermal.HeatTransfer.Celsius.PrescribedTemperature prescribedTemperature1 annotation(
    Placement(transformation(origin = {28, 90}, extent = {{8, -8}, {-8, 8}})));
  Modelica.Blocks.Sources.CombiTimeTable combiTimeTable1(columns = 1:2, fileName = "lut.txt", tableName = "tab", tableOnFile = true) annotation(
    Placement(transformation(origin = {90, 70}, extent = {{10, -10}, {-10, 10}})));
  Modelica.Blocks.Sources.CombiTimeTable combiTimeTable2(columns = 1:7, fileName = "lut.txt", tableName = "tab", tableOnFile = true) annotation(
    Placement(transformation(origin = {90, 90}, extent = {{10, -10}, {-10, 10}})));
  Modelica.Blocks.Sources.CombiTimeTable combiTimeTable3(columns = 1:10, fileName = "lut.txt", tableName = "tab", tableOnFile = true) annotation(
    Placement(transformation(origin = {-170, -170}, extent = {{-10, -10}, {10, 10}})));
  Greenhouses.Flows.HeatAndVapourTransfer.Ventilation ventilation(A = surface.k, U_vents = 1) annotation(
    Placement(transformation(origin = {-80, 30}, extent = {{-10, -10}, {10, 10}})));
  Greenhouses.Flows.Sources.Vapour.PrescribedPressure prescribedPressure annotation(
    Placement(transformation(origin = {8, 28}, extent = {{8, -8}, {-8, 8}})));
  Modelica.Blocks.Sources.CombiTimeTable combiTimeTable4(columns = 1:3, fileName = "lut.txt", tableName = "tab", tableOnFile = true) annotation(
    Placement(transformation(origin = {90, 30}, extent = {{10, -10}, {-10, 10}})));
  Modelica.Blocks.Sources.RealExpression realExpression(y = Greenhouses.Functions.WaterVapourPressure(combiTimeTable4.y[2], combiTimeTable4.y[3])) annotation(
    Placement(transformation(origin = {50, 28}, extent = {{10, -10}, {-10, 10}})));
  Greenhouses.Flows.HeatTransfer.Radiation_T4 radiation_T41(A = surface.k, FFb = canopy.FF, epsilon_a = 0.89, epsilon_b = 1) annotation(
    Placement(transformation(origin = {-172, -80}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.Constant surface(k = 100) annotation(
    Placement(transformation(origin = {-200, 70}, extent = {{-10, -10}, {10, 10}})));
  Greenhouses.Flows.HeatTransfer.Radiation_T4 radiation_T42(A = surface.k, FFa = canopy.FF, FFb = 1, epsilon_a = 1, epsilon_b = 0.85) annotation(
    Placement(transformation(origin = {-110, 70}, extent = {{-10, -10}, {10, 10}})));
  Greenhouses.Flows.HeatTransfer.CanopyFreeConvection canopyFreeConvection(A = surface.k, LAI = canopy.LAI) annotation(
    Placement(transformation(origin = {-96, -64}, extent = {{-10, -10}, {10, 10}})));
  Greenhouses.Flows.HeatTransfer.FreeConvection freeConvection(A = surface.k, floor = true, phi = 1.5707963267948966) annotation(
    Placement(transformation(origin = {-170, -38}, extent = {{-10, -10}, {10, 10}})));
  Greenhouses.Flows.VapourMassTransfer.MV_CanopyTranspiration mV_CanopyTranspiration(A = surface.k, LAI = canopy.LAI, T_can = canopy.T) annotation(
    Placement(transformation(origin = {-116, -10}, extent = {{-10, -10}, {10, 10}})));
equation
  connect(illumination.R_IluFlr_Glob, plant_panels.R_Flr_Glob[1]) annotation(
    Line(points = {{-97, -84}, {-114.5, -84}, {-114.5, -56}, {-162.25, -56}, {-162.25, -74}, {-184, -74}}, color = {255, 207, 14}));
  connect(illumination.R_IluCan_Glob, canopy.R_Can_Glob[1]) annotation(
    Line(points = {{-97, -90}, {-97, -89}, {-120.5, -89}, {-120.5, -66}, {-143, -66}}, color = {255, 207, 14}));
  connect(illumination.R_IluAir_Glob, air.R_Air_Glob[1]) annotation(
    Line(points = {{-97, -96}, {-122, -96}, {-122, -28}, {-91, -28}}, color = {255, 207, 14}));
  connect(outsideAirConvection.port_a, cover.heatPort) annotation(
    Line(points = {{-40, 70}, {-80, 70}}, color = {191, 0, 0}));
  connect(radiation_T4.port_a, cover.heatPort) annotation(
    Line(points = {{-40, 90}, {-40, 91}, {-66, 91}, {-66, 69.5}, {-80, 69.5}, {-80, 70}}, color = {191, 0, 0}));
  connect(prescribedTemperature.port, outsideAirConvection.port_b) annotation(
    Line(points = {{20, 70}, {-20, 70}}, color = {191, 0, 0}));
  connect(prescribedTemperature1.port, radiation_T4.port_b) annotation(
    Line(points = {{20, 90}, {-20, 90}}, color = {191, 0, 0}));
  connect(combiTimeTable1.y[2], prescribedTemperature.T) annotation(
    Line(points = {{80, 70}, {38, 70}}, color = {0, 0, 127}));
  connect(combiTimeTable2.y[7], prescribedTemperature1.T) annotation(
    Line(points = {{80, 90}, {38, 90}}, color = {0, 0, 127}));
  connect(combiTimeTable3.y[10], illumination.switch) annotation(
    Line(points = {{-159, -170}, {-159, -170.062}, {-147, -170.062}, {-147, -124}, {-67.5, -124}, {-67.5, -102}, {-67.75, -102}, {-67.75, -88}, {-84, -88}}, color = {0, 0, 127}));
  connect(ventilation.HeatPort_a, air.heatPort) annotation(
    Line(points = {{-90, 32}, {-90, 31}, {-94, 31}, {-94, -34}}, color = {191, 0, 0}));
  connect(ventilation.MassPort_a, air.massPort) annotation(
    Line(points = {{-90, 28}, {-92, 28}, {-92, -8}, {-98, -8}, {-98, -34}}, color = {0, 0, 255}));
  connect(prescribedPressure.port, ventilation.MassPort_b) annotation(
    Line(points = {{0, 28}, {-70, 28}}, color = {0, 0, 255}));
  connect(prescribedTemperature.port, ventilation.HeatPort_b) annotation(
    Line(points = {{20, 70}, {11, 70}, {11, 50}, {-70, 50}, {-70, 32}}, color = {191, 0, 0}));
  connect(realExpression.y, prescribedPressure.VP) annotation(
    Line(points = {{39, 28}, {18, 28}}, color = {0, 0, 127}));
  connect(plant_panels.heatPort, radiation_T41.port_a) annotation(
    Line(points = {{-188, -80}, {-182, -80}}, color = {191, 0, 0}));
  connect(radiation_T41.port_b, canopy.heatPort) annotation(
    Line(points = {{-162, -80}, {-159.5, -80}, {-159.5, -81}, {-146, -81}}, color = {191, 0, 0}));
  connect(radiation_T42.port_b, cover.heatPort) annotation(
    Line(points = {{-100, 70}, {-80, 70}}, color = {191, 0, 0}));
  connect(radiation_T42.port_a, canopy.heatPort) annotation(
    Line(points = {{-120, 70}, {-130, 70}, {-130, -82}, {-146, -82}}, color = {191, 0, 0}));
  connect(canopyFreeConvection.port_a, canopy.heatPort) annotation(
    Line(points = {{-106, -64}, {-106, -82}, {-146, -82}}, color = {191, 0, 0}));
  connect(canopyFreeConvection.port_b, air.heatPort) annotation(
    Line(points = {{-86, -64}, {-86, -34}, {-94, -34}}, color = {191, 0, 0}));
  connect(radiation_T41.port_a, freeConvection.port_a) annotation(
    Line(points = {{-182, -80}, {-180, -80}, {-180, -38}}, color = {191, 0, 0}));
  connect(freeConvection.port_b, air.heatPort) annotation(
    Line(points = {{-160, -38}, {-94, -38}, {-94, -34}}, color = {191, 0, 0}));
  connect(canopy.massPort, mV_CanopyTranspiration.port_a) annotation(
    Line(points = {{-146, -84}, {-132, -84}, {-132, -10}, {-126, -10}}, color = {0, 0, 255}));
  annotation(
    Diagram(graphics = {Text(origin = {-50, -85}, extent = {{-30, 3}, {30, -3}}, textString = "Model for high pressure sodium (HPS) lamps.
See documentation"), Rectangle(origin = {-140, -20}, lineColor = {119, 118, 123}, lineThickness = 1, extent = {{-60, 120}, {60, -120}}, radius = 1), Rectangle(origin = {-240, -20}, fillColor = {152, 106, 68}, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{-40, 120}, {40, -120}}, radius = 1), Polygon(origin = {-240, 120}, fillColor = {165, 29, 45}, fillPattern = FillPattern.Solid, lineThickness = 0.5, points = {{-40, -20}, {0, 20}, {40, -20}, {40, -20}, {-40, -20}}), Text(origin = {-350, 20}, extent = {{-50, 80}, {50, -80}}, textString = "DOING
- define system
\t- list state variables
\t- list flow variables
\t- list input variables
\t- order the variables into functional subsystems

TODO
- make model of eei-tower
\t- get geodata of the location
\t- get building data

DONE", horizontalAlignment = TextAlignment.Left)}, coordinateSystem(extent = {{-400, 140}, {100, -180}})),
    version = "",
    uses,
    experiment(StartTime = 0, StopTime = 29804400, Tolerance = 1e-6, Interval = 3599.57));end verticalfarm;