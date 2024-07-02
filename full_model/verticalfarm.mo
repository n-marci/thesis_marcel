package greenshell
    model plant
    model plant_transpiration
      Buildings.BoundaryConditions.WeatherData.Bus weaBus annotation(
          Placement(transformation(extent = {{-110, -10}, {-90, 10}}), iconTransformation(extent = {{-120, -18}, {-80, 22}})));
      Modelica.Blocks.Interfaces.RealInput v(unit="m/s") "Wind speed" annotation (
          Placement(transformation(extent={{-140,80},{-100,120}})));
      Modelica.Blocks.Interfaces.RealInput k_cb(unit="1") "basal crop coefficient - 0.01 for the initial stage and 0.21 for the midseason and late stages" annotation (
          Placement(transformation(origin = {0, -32}, extent = {{-140, 80}, {-100, 120}}), iconTransformation(origin = {0, -38}, extent = {{-140, 80}, {-100, 120}})));
      Modelica.Blocks.Interfaces.RealInput t(unit="K") "Air temperature" annotation (
          Placement(transformation(origin = {0, -130}, extent = {{-120, 100}, {-100, 120}}), iconTransformation(origin = {0, -150}, extent = {{-120, 100}, {-100, 120}})));
      Modelica.Blocks.Interfaces.RealInput p(unit="Pa") "Air pressure" annotation (
          Placement(transformation(origin = {0, -146}, extent = {{-130, 90}, {-100, 120}}), iconTransformation(origin = {0, -170}, extent = {{-120, 100}, {-100, 120}})));
      Modelica.Blocks.Interfaces.RealInput par(unit="umol/m2*s") "total radiation hitting the plant" annotation (
          Placement(transformation(origin = {0, -90}, extent = {{-120, 100}, {-100, 120}}), iconTransformation(origin = {0, -90}, extent = {{-120, 100}, {-100, 120}})));
      Modelica.Blocks.Interfaces.RealInput rh(unit="%") "relative humidity" annotation (
          Placement(transformation(origin = {0, -70}, extent = {{-120, 100}, {-100, 120}}), iconTransformation(origin = {0, -70}, extent = {{-120, 100}, {-100, 120}})));
      Modelica.Blocks.Interfaces.RealOutput y annotation(
          Placement(transformation(origin = {110, 0}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {110, 40}, extent = {{-10, -10}, {10, 10}})));
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b port annotation(
        Placement(transformation(origin = {0, -30}, extent = {{90, -10}, {110, 10}}), iconTransformation(origin = {0, -40}, extent = {{90, -10}, {110, 10}})));
    //  Modelica.Blocks.Interfaces.
      
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
      delta = 4098*(0.6108*Modelica.Math.exp((17.27*(weaBus.TDryBul - 273.15))/((weaBus.TDryBul - 273.15) + 237.3)))/(((weaBus.TDryBul - 273.15) + 237.3)^2);
      r_n = 500;
      gamma = 0.000665*weaBus.pAtm/1000;
      e_s = 0.6108*Modelica.Math.exp((17.27*(weaBus.TDryBul - 273.15))/((weaBus.TDryBul - 273.15) + 237.3));
      e_a = 0.6108*Modelica.Math.exp((17.27*(weaBus.TDewPoi - 273.15))/((weaBus.TDewPoi - 273.15) + 237.3));
      vpd = e_s - e_a;
//c_d = 0.24;
      if (weaBus.HGloHor == 0) then
        c_d = 0.96;
      else
        c_d = 0.24;
      end if;
      et_r = (1/3600)*(0.408*delta*r_n + gamma*c_n*v*(e_s - e_a)/weaBus.TDryBul)/(delta + gamma*(1 + c_d*v));
      et_a = a_plant*k_cb*et_r;
      q_flow = et_a*density_water/1000*l_v;
      port.Q_flow = q_flow;
      y = et_a;
    end plant_transpiration;

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
  equation

  end plant;

  model led
  equation

  end led;

  model pump
  equation

  end pump;

  model farm
    Buildings.Fluid.MixingVolumes.MixingVolumeMoistAir vol(redeclare package Medium = Modelica.Media.Air.MoistAir, m_flow_nominal = 0.001, V = 100, nPorts = 3)  annotation(
      Placement(transformation(origin = {-10, 30}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.Constant const(k = 0.1)  annotation(
      Placement(transformation(origin = {-54, 38}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor(G = 10000)  annotation(
      Placement(transformation(origin = {-54, 8}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Thermal.HeatTransfer.Sources.FixedTemperature fixedTemperature(T = 274.15)  annotation(
      Placement(transformation(origin = {-94, 8}, extent = {{-10, -10}, {10, 10}})));
  Buildings.Fluid.Sources.Boundary_pT bou(nPorts = 3, redeclare package Medium = Buildings.Media.Air "Moist air")  annotation(
      Placement(transformation(origin = {-10, -76}, extent = {{-10, -10}, {10, 10}})));
  Buildings.Airflow.Multizone.DoorOperable doo(LClo = 0.1, redeclare package Medium = Buildings.Media.Air "Moist air")  annotation(
      Placement(transformation(origin = {58, -4}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.Constant const11(k = 1) annotation(
      Placement(transformation(origin = {28, 26}, extent = {{-10, -10}, {10, 10}})));
  equation
    connect(const.y, vol.mWat_flow) annotation(
      Line(points = {{-42, 38}, {-22, 38}}, color = {0, 0, 127}));
    connect(vol.heatPort, thermalConductor.port_b) annotation(
      Line(points = {{-20, 30}, {-32, 30}, {-32, 8}, {-44, 8}}, color = {191, 0, 0}));
    connect(thermalConductor.port_a, fixedTemperature.port) annotation(
      Line(points = {{-64, 8}, {-84, 8}}, color = {191, 0, 0}));
    connect(doo.port_a1, vol.ports[2]) annotation(
      Line(points = {{48, 2}, {-10, 2}, {-10, 20}}, color = {0, 127, 255}));
    connect(doo.port_b2, vol.ports[3]) annotation(
      Line(points = {{48, -10}, {-10, -10}, {-10, 20}}, color = {0, 127, 255}));
    connect(doo.port_a2, bou.ports[2]) annotation(
      Line(points = {{68, -10}, {90, -10}, {90, -76}, {0, -76}}, color = {0, 127, 255}));
    connect(doo.port_b1, bou.ports[3]) annotation(
      Line(points = {{68, 2}, {90, 2}, {90, -76}, {0, -76}}, color = {0, 127, 255}));
    connect(const11.y, doo.y) annotation(
      Line(points = {{40, 26}, {42, 26}, {42, -4}, {48, -4}}, color = {0, 0, 127}));
    annotation(
      Diagram(coordinateSystem(extent = {{-120, 60}, {40, -100}})));
end farm;
  annotation(
    uses(Buildings(version = "11.0.0"), Modelica(version = "4.0.0")));
end greenshell;