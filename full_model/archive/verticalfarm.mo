within ;
package verticalfarm
  extends Modelica.Icons.Package;

  model canuse
  Buildings.Airflow.Multizone.MediumColumnDynamic MediumColumnDynamic1 annotation(
      Placement(transformation(origin = {-70, 10}, extent = {{-10, -10}, {10, 10}})));
  Buildings.Airflow.Multizone.Orifice orifice1 annotation(
      Placement(transformation(origin = {-32, -30}, extent = {{-10, -10}, {10, 10}})));
  Buildings.Airflow.Multizone.Examples.NaturalVentilation naturalVentilation annotation(
      Placement(transformation(origin = {50, 50}, extent = {{-10, -10}, {10, 10}})));
  Buildings.Airflow.Multizone.Examples.CO2TransportStep cO2TransportStep annotation(
      Placement(transformation(origin = {70, 30}, extent = {{-10, -10}, {10, 10}})));
  Buildings.BoundaryConditions.GroundTemperature.UndisturbedSoilTemperature undisturbedSoilTemperature annotation(
      Placement(transformation(origin = {-50, -90}, extent = {{-10, -10}, {10, 10}})));
  equation

  annotation(
      Diagram(graphics = {Text(origin = {-11, -33}, extent = {{-13, 3}, {13, -3}}, textString = "vent"), Text(origin = {-47, 10}, extent = {{-11, 4}, {11, -4}}, textString = "inside air")}));
end canuse;

  model environment
  Buildings.Fluid.Sources.Boundary_pT outside_air(nPorts = 2, redeclare package Medium = Buildings.Media.Air "Moist air", T = 273.15)  annotation(
      Placement(transformation(origin = {102, 10}, extent = {{10, -10}, {-10, 10}})));
  Buildings.Airflow.Multizone.MediumColumn col(redeclare package Medium = Buildings.Media.Air "Moist air", densitySelection = Buildings.Airflow.Multizone.Types.densitySelection.fromBottom)  annotation(
      Placement(transformation(origin = {70, 30}, extent = {{-10, -10}, {10, 10}})));
  Buildings.Airflow.Multizone.MediumColumn col1(redeclare package Medium = Buildings.Media.Air "Moist air", densitySelection = Buildings.Airflow.Multizone.Types.densitySelection.fromTop)  annotation(
      Placement(transformation(origin = {70, -10}, extent = {{-10, -10}, {10, 10}})));
  Buildings.Airflow.Multizone.Orifice ori(redeclare package Medium = Buildings.Media.Air "Moist air", A = 1)  annotation(
      Placement(transformation(origin = {0, -30}, extent = {{-10, -10}, {10, 10}})));
  Buildings.Airflow.Multizone.Orifice ori1(redeclare package Medium = Buildings.Media.Air "Moist air", A = 1)  annotation(
      Placement(transformation(origin = {0, 50}, extent = {{-10, -10}, {10, 10}})));
  Buildings.Airflow.Multizone.MediumColumn col2(redeclare package Medium = Buildings.Media.Air "Moist air", densitySelection = Buildings.Airflow.Multizone.Types.densitySelection.fromTop)  annotation(
      Placement(transformation(origin = {-30, -10}, extent = {{-10, -10}, {10, 10}})));
  Buildings.Airflow.Multizone.MediumColumn col3(redeclare package Medium = Buildings.Media.Air "Moist air", densitySelection = Buildings.Airflow.Multizone.Types.densitySelection.fromBottom)  annotation(
      Placement(transformation(origin = {-30, 30}, extent = {{-10, -10}, {10, 10}})));
  Buildings.Fluid.MixingVolumes.MixingVolume inside_air1(nPorts = 2, redeclare package Medium = Buildings.Media.Air "Moist air", m_flow_nominal = 0.001, V = 1)  annotation(
      Placement(transformation(origin = {-70, 34}, extent = {{-10, -10}, {10, 10}})));
  equation
    connect(ori1.port_b, col.port_a) annotation(
      Line(points = {{10, 50}, {70, 50}, {70, 40}}, color = {0, 127, 255}));
    connect(ori.port_b, col1.port_b) annotation(
      Line(points = {{10, -30}, {70, -30}, {70, -20}}, color = {0, 127, 255}));
    connect(outside_air.ports[1], col.port_b) annotation(
      Line(points = {{92, 10}, {70, 10}, {70, 20}}, color = {0, 127, 255}));
    connect(col1.port_a, outside_air.ports[2]) annotation(
      Line(points = {{70, 0}, {70, 10}, {92, 10}}, color = {0, 127, 255}));
    connect(ori1.port_a, col3.port_a) annotation(
      Line(points = {{-10, 50}, {-30, 50}, {-30, 40}}, color = {0, 127, 255}));
    connect(col2.port_b, ori.port_a) annotation(
      Line(points = {{-30, -20}, {-30, -30}, {-10, -30}}, color = {0, 127, 255}));
  connect(col3.port_b, inside_air1.ports[1]) annotation(
      Line(points = {{-30, 20}, {-70, 20}, {-70, 24}}, color = {0, 127, 255}));
  connect(col2.port_a, inside_air1.ports[2]) annotation(
      Line(points = {{-30, 0}, {-70, 0}, {-70, 24}}, color = {0, 127, 255}));
    annotation(
      Diagram(graphics = {Rectangle(origin = {-50, 10}, lineColor = {192, 191, 188}, lineThickness = 1, extent = {{-50, 50}, {50, -50}})}, coordinateSystem(extent = {{-100, 60}, {120, -60}})));
end environment;
  annotation(
    uses(Buildings(version = "10.0.0")));
end verticalfarm;