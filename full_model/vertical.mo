package vertical
  model plant
    replaceable package Medium = Modelica.Media.Interfaces.PartialMedium "Medium in the component" annotation(
      choices(choice(redeclare package Medium = Buildings.Media.Air "Moist air")));
    Modelica.Fluid.Interfaces.FluidPort_a farm_air(redeclare package Medium = Medium, p(start = Medium.p_default)) "Fluid connector a (positive design flow direction is from port_a to port_b)" annotation(
      Placement(transformation(extent = {{-10, 90}, {10, 110}}, rotation = 90), iconTransformation(extent = {{-10, 90}, {10, 110}}, rotation = 90)));
    Modelica.Blocks.Interfaces.RealInput par(unit = "W/m2") "Photosynthetically active radiation" annotation(
      Placement(transformation(origin = {0, -20}, extent = {{-140, 80}, {-100, 120}}), iconTransformation(origin = {0, -20}, extent = {{-140, 80}, {-100, 120}})));
    Modelica.Blocks.Interfaces.RealInput v_air(unit = "m/s") "Air Speed" annotation(
      Placement(transformation(origin = {0, -160}, extent = {{-140, 80}, {-100, 120}}), iconTransformation(origin = {0, -180}, extent = {{-140, 80}, {-100, 120}})));
    Modelica.Blocks.Interfaces.RealOutput y annotation(
      Placement(transformation(origin = {110, 0}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {110, 0}, extent = {{-10, -10}, {10, 10}})));
    parameter Modelica.Units.SI.Area a_plant "planting area";

    model plant_transpiration
      //  Buildings.BoundaryConditions.WeatherData.Bus weaBus annotation(
      //      Placement(transformation(extent = {{-110, -10}, {-90, 10}}), iconTransformation(extent = {{-120, -18}, {-80, 22}})));
      Modelica.Blocks.Interfaces.RealInput v_air(unit = "m/s") "Air speed" annotation(
        Placement(transformation(extent = {{-140, 80}, {-100, 120}})));
      Modelica.Blocks.Interfaces.RealInput k_cb(unit = "1") "basal crop coefficient - 0.01 for the initial stage and 0.21 for the midseason and late stages" annotation(
        Placement(transformation(origin = {0, -40}, extent = {{-140, 80}, {-100, 120}}), iconTransformation(origin = {0, -40}, extent = {{-140, 80}, {-100, 120}})));
      Modelica.Blocks.Interfaces.RealInput air_temperature(unit = "K") "Air temperature" annotation(
        Placement(transformation(origin = {0, -170}, extent = {{-120, 100}, {-100, 120}}), iconTransformation(origin = {0, -170}, extent = {{-120, 100}, {-100, 120}})));
      Modelica.Blocks.Interfaces.RealInput air_pressure(unit = "Pa") "Air pressure" annotation(
        Placement(transformation(origin = {0, -186}, extent = {{-130, 90}, {-100, 120}}), iconTransformation(origin = {0, -190}, extent = {{-120, 100}, {-100, 120}})));
      Modelica.Blocks.Interfaces.RealInput par(unit = "umol/m2*s") "total radiation hitting the plant" annotation(
        Placement(transformation(origin = {0, -110}, extent = {{-120, 100}, {-100, 120}}), iconTransformation(origin = {0, -110}, extent = {{-120, 100}, {-100, 120}})));
      Modelica.Blocks.Interfaces.RealInput relative_humidity(unit = "%") "relative humidity" annotation(
        Placement(transformation(origin = {0, -150}, extent = {{-120, 100}, {-100, 120}}), iconTransformation(origin = {0, -150}, extent = {{-120, 100}, {-100, 120}})));
      Modelica.Blocks.Interfaces.RealOutput y annotation(
        Placement(transformation(origin = {110, 0}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {110, 40}, extent = {{-10, -10}, {10, 10}})));
      //  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b port annotation(
      //    Placement(transformation(origin = {0, -30}, extent = {{90, -10}, {110, 10}}), iconTransformation(origin = {0, -40}, extent = {{90, -10}, {110, 10}})));
      //  Modelica.Blocks.Interfaces.
      parameter Modelica.Units.SI.Area a_plant "planting area";
      Real et_r(final unit = "mm") "reference evapotranspiration";
      Real et_a(final unit = "L") "actual evapotranspiration";
      Real delta(final unit = "kPa/°C") "hourly slope of the vapor pressure curve";
      Real r_n(final unit = "MJ/m2") "hourly net radiation flux";
      Real gamma(final unit = "kPa/°C") "psychrometric constant";
      Real c_n(final unit = "1") = 37 "numerator constant for the reference crop type and time step";
      //@lopez_mora2024
      Real e_s(final unit = "kPa") "saturation vapor pressure";
      Real e_a(final unit = "kPa") "actual vapor pressure";
      Real vpd(final unit = "kPa") "vapor pressure deficit";
      Real c_d(final unit = "1") "denominator constant for the reference crop type and time step";
      Real l_v(final unit = "kJ/kg") = 2260 "latent heat of vaporization for water";
      //  Modelica.Units.SI.HeatFlowRate q_flow "heat absorbed during evapotranspiration";
      Modelica.Units.SI.Density density_water = 1000 "density of water = 1000 kg/m3";
      Real t_dew(final unit = "°C") "dew point temperature";
    equation
      delta = 4098*(0.6108*Modelica.Math.exp((17.27*(air_temperature - 273.15))/((air_temperature - 273.15) + 237.3)))/(((air_temperature - 273.15) + 237.3)^2);
      r_n = 500;
      gamma = 0.000665*air_pressure/1000;
      e_s = 0.6108*Modelica.Math.exp((17.27*(air_temperature - 273.15))/((air_temperature - 273.15) + 237.3));
      e_a = 0.6108*Modelica.Math.exp((17.27*t_dew)/(t_dew + 237.3));
      t_dew = (air_temperature - 273.15) - ((100 - relative_humidity)/5);
      vpd = e_s - e_a;
//c_d = 0.24;
      if (par > 0.1) then
        c_d = 0.96;
      else
        c_d = 0.24;
      end if;
      et_r = (1/3600)*(0.408*delta*r_n + gamma*c_n*v_air*(e_s - e_a)/air_temperature)/(delta + gamma*(1 + c_d*v_air));
      et_a = a_plant*k_cb*et_r;
//  q_flow = et_a*density_water/1000*l_v;
//  port.Q_flow = q_flow;
      y = et_a;
    end plant_transpiration;

    model plant_yield
      //  Buildings.BoundaryConditions.WeatherData.Bus weaBus annotation(
      //      Placement(transformation(extent = {{-110, -10}, {-90, 10}}), iconTransformation(extent = {{-120, -18}, {-80, 22}})));
      // model taken from doi.org/10.1016/S0308-521X(94)90280-1
      Modelica.Blocks.Interfaces.RealInput air_temperature(unit = "K") "Air temperature" annotation(
        Placement(transformation(origin = {0, -80}, extent = {{-140, 80}, {-100, 120}}), iconTransformation(origin = {0, -20}, extent = {{-140, 80}, {-100, 120}})));
      Modelica.Blocks.Interfaces.RealInput co2_concentration(unit = "ppm") "CO2 concentration" annotation(
        Placement(transformation(origin = {0, -120}, extent = {{-140, 80}, {-100, 120}}), iconTransformation(origin = {0, -100}, extent = {{-140, 80}, {-100, 120}})));
      Modelica.Blocks.Interfaces.RealInput par(unit = "W/m2") "Photosynthetically active radiation" annotation(
        Placement(transformation(origin = {0, -160}, extent = {{-140, 80}, {-100, 120}}), iconTransformation(origin = {0, -180}, extent = {{-140, 80}, {-100, 120}})));
      Modelica.Blocks.Interfaces.BooleanInput reset "Connector of reset signal" annotation(
        Placement(transformation(extent = {{-20, -20}, {20, 20}}, rotation = 90, origin = {60, -120})));
      Modelica.Blocks.Interfaces.RealOutput dw(unit = "kg") annotation(
        Placement(transformation(origin = {110, 0}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {110, 0}, extent = {{-10, -10}, {10, 10}})));
      parameter Modelica.Units.SI.Area a_plant "planting area";
      parameter Modelica.Units.SI.Area x_nsdw_start = 0.25*0.72 "start value for non-structural dry weight";
      parameter Modelica.Units.SI.Area x_sdw_start = 0.75*0.72 "start value for structural dry weight";
      Real x_nsdw(final unit = "g/m2") "state variable non-structural dry weight";
      Real x_sdw(final unit = "g/m2") "state variable structural dry weight";
      //"start = 0.001," "start = 0.001,"
      //  Real dw(final unit = "g") "total dry weight";
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
      Real g_bnd = /*(final unit = "m/s")*/0.007 "Boundary layer conductance";
      Real g_stm = /*(final unit = "m/s")*/0.005 "Stomatal conductance";
      Real g_car "Carboxylation conductance";
      /*(final unit = "m/s")*/
      Real r_gr(final unit = "1/s") "specific growth rate";
      Real c_gr_max(final unit = "1/s") = 5e-6 "saturation growth rate at 20°C";
      Real c_gamma(final unit = "1") = 1.0 "growth rate coefficient";
      Real c_q_10_gr(final unit = "1") = 1.6 "growth rate sensitivity to the canopy temperature";
    initial equation
      x_nsdw = x_nsdw_start;
      x_sdw = x_sdw_start;
    equation
      when reset then
        reinit(x_nsdw, x_nsdw_start);
        reinit(x_sdw, x_sdw_start);
      end when;
      der(x_nsdw) = c_alpha*f_phot - r_gr*x_sdw - f_resp - ((1 - c_beta)/c_beta)*r_gr*x_sdw;
      der(x_sdw) = r_gr*x_sdw;
      r_gr = c_gr_max*(x_nsdw/(c_gamma*x_sdw + x_nsdw))*c_q_10_gr^(((air_temperature - 273.15) - 20)/10);
      f_resp = (c_resp_sht*(1 - c_tau)*x_sdw + c_resp_rt*c_tau*x_sdw)*c_q_10_resp^(((air_temperature - 273.15) - 25)/10);
      f_phot = (1 - exp(-c_k*c_lar*(1 - c_tau)*x_sdw))*f_phot_max;
      f_phot_max = (epsilon*par*g_co2*c_omega*(co2_concentration - capital_gamma))/(epsilon*par + g_co2*c_omega*(co2_concentration - capital_gamma));
      capital_gamma = c_cap_gamma*c_q_10_cap_gamma^(((air_temperature - 273.15) - 20)/10);
      epsilon = c_epsilon*(co2_concentration - capital_gamma)/(co2_concentration + 2*capital_gamma);
      g_co2 = 1/((1/g_bnd) + (1/g_stm) + (1/g_car));
      g_car = max(-1.32e-5*(air_temperature - 273.15)^2 + 5.94e-4*(air_temperature - 273.15) - 2.64e-3, Modelica.Constants.small);
      dw = a_plant*(x_nsdw + x_sdw)/1000;
    end plant_yield;

    plant_transpiration transpiration(a_plant = a_plant) annotation(
      Placement(transformation(origin = {10, -30}, extent = {{-10, -10}, {10, 10}})));
    plant_yield yield(a_plant = a_plant, x_nsdw_start = 0.25*2.7, x_sdw_start = 0.75*2.7) annotation(
      Placement(transformation(origin = {10, 50}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Fluid.Sensors.Pressure pressure(redeclare package Medium = Buildings.Media.Air "Moist air") annotation(
      Placement(transformation(origin = {-70, 10}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Fluid.Sensors.Temperature temperature(redeclare package Medium = Buildings.Media.Air "Moist air") annotation(
      Placement(transformation(origin = {-90, 30}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Sources.Constant co2_concentration(k = 416) annotation(
      Placement(transformation(origin = {-20, 80}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Blocks.Sources.Constant basal_crop_coefficient(k = 0.1) annotation(
      Placement(transformation(origin = {-20, 2}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Blocks.Sources.Constant relative_humidity(k = 0.8) annotation(
      Placement(transformation(origin = {-20, -82}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
    Modelica.Blocks.Math.Gain ppfd_to_par(k = 1/4.56) annotation(
      Placement(transformation(origin = {-24, 42}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Interfaces.BooleanInput reset annotation(
      Placement(transformation(origin = {56, -120}, extent = {{-20, -20}, {20, 20}}, rotation = 90), iconTransformation(origin = {60, -120}, extent = {{-20, -20}, {20, 20}}, rotation = 90)));
  equation
    connect(pressure.port, farm_air) annotation(
      Line(points = {{-70, 0}, {-100, 0}}, color = {0, 127, 255}));
    connect(temperature.port, farm_air) annotation(
      Line(points = {{-90, 20}, {-90, 0}, {-100, 0}}, color = {0, 127, 255}));
    connect(transpiration.air_temperature, temperature.T) annotation(
      Line(points = {{0, -36}, {-48, -36}, {-48, 30}, {-82, 30}}, color = {0, 0, 127}));
    connect(transpiration.air_pressure, pressure.p) annotation(
      Line(points = {{0, -38}, {-50, -38}, {-50, 10}, {-58, 10}}, color = {0, 0, 127}));
    connect(co2_concentration.y, yield.co2_concentration) annotation(
      Line(points = {{-20, 70}, {-20, 50}, {-2, 50}}, color = {0, 0, 127}));
    connect(transpiration.par, par) annotation(
      Line(points = {{0, -30}, {-40, -30}, {-40, 80}, {-120, 80}}, color = {0, 0, 127}));
    connect(basal_crop_coefficient.y, transpiration.k_cb) annotation(
      Line(points = {{-20, -9}, {-20, -24}, {-2, -24}}, color = {0, 0, 127}));
    connect(transpiration.v_air, v_air) annotation(
      Line(points = {{-2, -20}, {-10, -20}, {-10, -60}, {-120, -60}}, color = {0, 0, 127}));
    connect(yield.dw, y) annotation(
      Line(points = {{22, 50}, {60, 50}, {60, 0}, {110, 0}}, color = {0, 0, 127}));
    connect(relative_humidity.y, transpiration.relative_humidity) annotation(
      Line(points = {{-20, -70}, {-20, -34}, {0, -34}}, color = {0, 0, 127}));
    connect(ppfd_to_par.y, yield.par) annotation(
      Line(points = {{-12, 42}, {-2, 42}}, color = {0, 0, 127}));
    connect(ppfd_to_par.u, par) annotation(
      Line(points = {{-36, 42}, {-40, 42}, {-40, 80}, {-120, 80}}, color = {0, 0, 127}));
    connect(temperature.T, yield.air_temperature) annotation(
      Line(points = {{-82, 30}, {-48, 30}, {-48, 58}, {-2, 58}}, color = {0, 0, 127}));
    connect(reset, yield.reset) annotation(
      Line(points = {{56, -120}, {56, 32}, {16, 32}, {16, 38}}, color = {255, 0, 255}));
    annotation(
      Icon(graphics = {Bitmap(extent = {{-100, -100}, {100, 100}}, imageSource = "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")}));
  end plant;

  model led
  equation

  end led;

  model farm
    model envelope
      Buildings.BoundaryConditions.WeatherData.Bus weaBus annotation(
        Placement(transformation(origin = {-40, 0}, extent = {{-110, -10}, {-90, 10}}), iconTransformation(extent = {{-120, -18}, {-80, 22}})));
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b interior annotation(
        Placement(transformation(extent = {{90, -10}, {110, 10}}), iconTransformation(extent = {{90, -10}, {110, 10}})));
      Modelica.Blocks.Interfaces.RealInput v_air(unit = "m/s") annotation(
        Placement(transformation(origin = {0, -40}, extent = {{140, 80}, {100, 120}}), iconTransformation(origin = {-10, -40}, extent = {{140, 80}, {100, 120}}, rotation = -0)));
      parameter Modelica.Units.SI.Area area "envelope area";
      parameter Modelica.Units.SI.Angle azimuth "Surface azimuth";
      Buildings.HeatTransfer.Convection.Exterior convection_out(A = area, azi = azimuth, roughness = Buildings.HeatTransfer.Types.SurfaceRoughness.VerySmooth, til = 1.5707963267948966) annotation(
        Placement(transformation(origin = {-40, 0}, extent = {{10, -10}, {-10, 10}})));
      Buildings.HeatTransfer.Conduction.SingleLayer conduction_glass(A = area, material = glass) annotation(
        Placement(transformation(extent = {{-10, -10}, {10, 10}})));
      Buildings.HeatTransfer.Data.Solids.Glass glass(x = 25) annotation(
        Placement(transformation(origin = {2, -40}, extent = {{-10, -10}, {10, 10}})));
      Modelica.Blocks.Sources.Constant wind_direction(k = 3.14159 - azimuth) annotation(
        Placement(transformation(origin = {0, 50}, extent = {{-10, -10}, {10, 10}}, rotation = -0)));
      Modelica.Blocks.Math.Gain radiation_heat(k = area*0.2) annotation(
        Placement(transformation(origin = {10, -100}, extent = {{-10, -10}, {10, 10}})));
      Buildings.HeatTransfer.Sources.PrescribedHeatFlow radiation_preHeaFlo annotation(
        Placement(transformation(origin = {50, -100}, extent = {{-10, -10}, {10, 10}})));
      Buildings.HeatTransfer.Convection.Exterior convection_in(A = area, azi = azimuth, roughness = Buildings.HeatTransfer.Types.SurfaceRoughness.VerySmooth, til = 1.5707963267948966) annotation(
        Placement(transformation(origin = {40, 0}, extent = {{-10, -10}, {10, 10}}, rotation = -0)));
      Buildings.HeatTransfer.Sources.PrescribedTemperature city_temperature annotation(
        Placement(transformation(origin = {-80, 0}, extent = {{-10, -10}, {10, 10}})));
      Buildings.BoundaryConditions.SolarIrradiation.DiffuseIsotropic HDifTilIso(til = 1.5707963267948966) annotation(
        Placement(transformation(origin = {-70, -90}, extent = {{-10, -10}, {10, 10}})));
      Buildings.BoundaryConditions.SolarIrradiation.DirectTiltedSurface HDirTil(til = 1.5707963267948966, azi = azimuth) annotation(
        Placement(transformation(origin = {-70, -110}, extent = {{-10, -10}, {10, 10}})));
      Modelica.Blocks.Math.Add add annotation(
        Placement(transformation(origin = {-30, -100}, extent = {{-10, -10}, {10, 10}})));
    equation
      connect(HDifTilIso.H, add.u1) annotation(
        Line(points = {{-59, -90}, {-51, -90}, {-51, -94}, {-43, -94}}, color = {0, 0, 127}));
      connect(HDirTil.H, add.u2) annotation(
        Line(points = {{-59, -110}, {-51, -110}, {-51, -106}, {-43, -106}}, color = {0, 0, 127}));
      connect(add.y, radiation_heat.u) annotation(
        Line(points = {{-19, -100}, {-3, -100}}, color = {0, 0, 127}));
      connect(radiation_heat.y, radiation_preHeaFlo.Q_flow) annotation(
        Line(points = {{21, -100}, {39, -100}}, color = {0, 0, 127}));
      connect(wind_direction.y, convection_in.dir) annotation(
        Line(points = {{12, 50}, {20, 50}, {20, 6}, {28, 6}}, color = {0, 0, 127}));
      connect(convection_in.v, v_air) annotation(
        Line(points = {{28, 10}, {26, 10}, {26, 60}, {120, 60}}, color = {0, 0, 127}));
      connect(weaBus.TDryBul, city_temperature.T) annotation(
        Line(points = {{-140, 0}, {-92, 0}}, color = {0, 0, 127}));
      connect(convection_out.solid, conduction_glass.port_a) annotation(
        Line(points = {{-30, 0}, {-10, 0}}, color = {191, 0, 0}));
      connect(conduction_glass.port_b, convection_in.solid) annotation(
        Line(points = {{10, 0}, {30, 0}}, color = {191, 0, 0}));
      connect(convection_in.fluid, interior) annotation(
        Line(points = {{50, 0}, {100, 0}}, color = {191, 0, 0}));
      connect(convection_out.fluid, city_temperature.port) annotation(
        Line(points = {{-50, 0}, {-70, 0}}, color = {191, 0, 0}));
      connect(convection_out.v, weaBus.winSpe) annotation(
        Line(points = {{-28, 10}, {-24, 10}, {-24, 20}, {-140, 20}, {-140, 0}}, color = {0, 0, 127}));
      connect(convection_out.dir, weaBus.winDir) annotation(
        Line(points = {{-28, 6}, {-20, 6}, {-20, 24}, {-140, 24}, {-140, 0}}, color = {0, 0, 127}));
      connect(HDifTilIso.weaBus, weaBus) annotation(
        Line(points = {{-80, -90}, {-140, -90}, {-140, 0}}, color = {255, 204, 51}, thickness = 0.5));
      connect(HDirTil.weaBus, weaBus) annotation(
        Line(points = {{-80, -110}, {-140, -110}, {-140, 0}}, color = {255, 204, 51}, thickness = 0.5));
      connect(radiation_preHeaFlo.port, interior) annotation(
        Line(points = {{60, -100}, {100, -100}, {100, 0}}, color = {191, 0, 0}));
      annotation(
        Diagram(coordinateSystem(extent = {{-160, 80}, {140, -120}})),
        Icon(graphics = {Rectangle(lineColor = {119, 118, 123}, fillColor = {225, 246, 255}, fillPattern = FillPattern.Solid, lineThickness = 8, extent = {{-90, 100}, {90, -100}})}));
    end envelope;

    model illumination_control
      Modelica.Blocks.Interfaces.RealOutput y annotation(
        Placement(transformation(origin = {110, 0}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {110, 0}, extent = {{-10, -10}, {10, 10}})));

      model shading
        Buildings.BoundaryConditions.WeatherData.Bus weaBus annotation(
          Placement(transformation(extent = {{-110, -10}, {-90, 10}}), iconTransformation(extent = {{-120, -18}, {-80, 22}})));
        Modelica.Blocks.Interfaces.RealOutput y annotation(
          Placement(transformation(origin = {110, 0}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {110, 0}, extent = {{-10, -10}, {10, 10}})));
        parameter Real light_setpoint "setpoint for natural lighting in umol/(m2*s)";
        parameter Modelica.Units.SI.Angle azimuth "Surface azimuth";
        Buildings.BoundaryConditions.SolarIrradiation.DiffuseIsotropic HDifTilIso(til = 1.5707963267948966) annotation(
          Placement(transformation(origin = {-70, 10}, extent = {{-10, -10}, {10, 10}})));
        Buildings.BoundaryConditions.SolarIrradiation.DirectTiltedSurface HDirTil(azi = azimuth, til = 1.5707963267948966) annotation(
          Placement(transformation(origin = {-70, -10}, extent = {{-10, -10}, {10, 10}})));
        Modelica.Blocks.Math.Add add annotation(
          Placement(transformation(origin = {-30, 0}, extent = {{-10, -10}, {10, 10}})));
        Buildings.Controls.Continuous.LimPID shading_control(controllerType = Modelica.Blocks.Types.SimpleController.P, yMax = 1, yMin = 0.6, k = 1) annotation(
          Placement(transformation(origin = {10, -44}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
        Modelica.Blocks.Sources.Constant light_set(k = light_setpoint) annotation(
          Placement(transformation(origin = {50, -44}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
        Modelica.Blocks.Math.Product shade annotation(
          Placement(transformation(origin = {10, 0}, extent = {{-10, -10}, {10, 10}})));
        Modelica.Blocks.Math.Gain irradiance_to_ppfd(k = 2.02) annotation(
          Placement(transformation(origin = {50, 0}, extent = {{-10, -10}, {10, 10}})));
      equation
        connect(shading_control.y, shade.u2) annotation(
          Line(points = {{0, -44}, {-10, -44}, {-10, -6}, {-2, -6}}, color = {0, 0, 127}));
        connect(add.y, shade.u1) annotation(
          Line(points = {{-18, 0}, {-14, 0}, {-14, 6}, {-2, 6}}, color = {0, 0, 127}));
        connect(HDifTilIso.H, add.u1) annotation(
          Line(points = {{-58, 10}, {-50, 10}, {-50, 6}, {-42, 6}}, color = {0, 0, 127}));
        connect(HDirTil.H, add.u2) annotation(
          Line(points = {{-58, -10}, {-50, -10}, {-50, -6}, {-42, -6}}, color = {0, 0, 127}));
        connect(weaBus, HDifTilIso.weaBus) annotation(
          Line(points = {{-100, 0}, {-88, 0}, {-88, 10}, {-80, 10}}, thickness = 0.5));
        connect(weaBus, HDirTil.weaBus) annotation(
          Line(points = {{-100, 0}, {-88, 0}, {-88, -10}, {-80, -10}}, thickness = 0.5));
        connect(shading_control.u_s, light_set.y) annotation(
          Line(points = {{22, -44}, {39, -44}}, color = {0, 0, 127}));
        connect(shade.y, irradiance_to_ppfd.u) annotation(
          Line(points = {{22, 0}, {38, 0}}, color = {0, 0, 127}));
        connect(irradiance_to_ppfd.y, y) annotation(
          Line(points = {{61, 0}, {110, 0}}, color = {0, 0, 127}));
        connect(irradiance_to_ppfd.y, shading_control.u_m) annotation(
          Line(points = {{62, 0}, {80, 0}, {80, -68}, {10, -68}, {10, -56}}, color = {0, 0, 127}));
        annotation(
          Icon(graphics = {Bitmap(extent = {{-100, 100}, {100, -100}}, imageSource = "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"), Text(origin = {0, -90}, extent = {{-100, 10}, {100, -10}}, textString = "Shading")}));
      end shading;

      model lighting
        parameter Modelica.Units.SI.Power power = 77.76 "led power draw";
        parameter Modelica.Units.SI.Area agricultural_area = 523.2 "planting area";
        parameter Real efficiency = 0.622 "efficiency of the led";
        parameter Real opti_dli(final unit = "mol/(m2*day)") = 14.4 "optimal dli for the chosen plant";
        parameter Real ppfd_led(final unit = "µmol/(m2*s)") = 218.28 "ppfd of the led module";
        Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b port annotation(
          Placement(transformation(origin = {60, -32}, extent = {{90, -10}, {110, 10}}), iconTransformation(origin = {0, -60}, extent = {{90, -10}, {110, 10}})));
        Modelica.Blocks.Continuous.Integrator dli(k = 1/1000000, use_reset = true) annotation(
          Placement(transformation(origin = {-90, 10}, extent = {{-10, -10}, {10, 10}})));
        Modelica.Blocks.Sources.SampleTrigger daily_trigger(period = 86400) annotation(
          Placement(transformation(origin = {-110, -30}, extent = {{-10, -10}, {10, 10}})));
        Modelica.Blocks.Discrete.Sampler dli_end_of_day(samplePeriod = 86400, startTime = 86400 - 1) annotation(
          Placement(transformation(origin = {-56, 10}, extent = {{-10, -10}, {10, 10}})));
        Modelica.Blocks.Math.Add missing_dli(k1 = -1) annotation(
          Placement(transformation(origin = {-20, 10}, extent = {{-10, -10}, {10, 10}})));
        Modelica.Blocks.Sources.Constant optimal_dli(k = opti_dli) annotation(
          Placement(transformation(origin = {-58, -30}, extent = {{-10, -10}, {10, 10}})));
        Modelica.Blocks.Logical.Timer timer annotation(
          Placement(transformation(origin = {40, -82}, extent = {{-10, -10}, {10, 10}})));
        Modelica.Blocks.Logical.OnOffController onOffController(bandwidth = 0) annotation(
          Placement(transformation(origin = {74, 10}, extent = {{-10, -10}, {10, 10}})));
        Modelica.Blocks.Sources.Constant ppfd(k = ppfd_led) annotation(
          Placement(transformation(origin = {-20, -30}, extent = {{-10, -10}, {10, 10}})));
        Modelica.Blocks.Math.Gain gain(k = 1/1000000) annotation(
          Placement(transformation(origin = {10, -30}, extent = {{-10, -10}, {10, 10}})));
        Modelica.Blocks.Math.Division time_to_achieve_opt_dli annotation(
          Placement(transformation(origin = {42, 10}, extent = {{-10, -10}, {10, 10}})));
        Modelica.Blocks.Math.BooleanToReal booleanToReal(realTrue = power*agricultural_area) annotation(
          Placement(transformation(origin = {110, 10}, extent = {{-10, -10}, {10, 10}})));
        Modelica.Blocks.Math.Gain heat_loss_efficiency(k = (1 - efficiency)) annotation(
          Placement(transformation(origin = {90, -32}, extent = {{-10, -10}, {10, 10}})));
        Modelica.Blocks.Interfaces.RealOutput p annotation(
          Placement(transformation(origin = {170, 10}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {110, 0}, extent = {{-10, -10}, {10, 10}})));
        Modelica.Blocks.Interfaces.RealOutput ppfd_out annotation(
          Placement(transformation(origin = {170, 50}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {110, 60}, extent = {{-10, -10}, {10, 10}})));
        Buildings.HeatTransfer.Sources.PrescribedHeatFlow preHeaFlo annotation(
          Placement(transformation(origin = {124, -32}, extent = {{-10, -10}, {10, 10}})));
        Modelica.Blocks.Interfaces.RealInput shaded_radiation annotation(
          Placement(transformation(origin = {-140, 10}, extent = {{-20, -20}, {20, 20}}), iconTransformation(origin = {-120, 0}, extent = {{-20, -20}, {20, 20}})));
        //Modelica.Blocks.Interfaces.RealInput h_dif annotation(
        //    Placement(transformation(origin = {-218, -20}, extent = {{-20, -20}, {20, 20}}), iconTransformation(origin = {-120, -40}, extent = {{-20, -20}, {20, 20}})));
        Modelica.Blocks.Nonlinear.Limiter limiter(uMax = Modelica.Constants.inf, uMin = 0) annotation(
          Placement(transformation(origin = {10, 10}, extent = {{-10, -10}, {10, 10}})));
        Modelica.Blocks.Math.RealToBoolean realToBoolean(threshold = 0.0001) annotation(
          Placement(transformation(origin = {-38, -76}, extent = {{-10, -10}, {10, 10}})));
        Modelica.Blocks.Logical.Xor xor annotation(
          Placement(transformation(origin = {4, -82}, extent = {{-10, -10}, {10, 10}})));
        Modelica.Blocks.Sources.BooleanPulse booleanPulse(period = 86400, width = 0.1) annotation(
          Placement(transformation(origin = {-38, -104}, extent = {{-10, -10}, {10, 10}})));
        Modelica.Blocks.Math.BooleanToReal booleanToReal1(realTrue = ppfd_led) annotation(
          Placement(transformation(origin = {110, 50}, extent = {{-10, -10}, {10, 10}})));
      equation
        connect(daily_trigger.y, dli.reset) annotation(
          Line(points = {{-99, -30}, {-85, -30}, {-85, -2}}, color = {255, 0, 255}));
        connect(dli.y, dli_end_of_day.u) annotation(
          Line(points = {{-79, 10}, {-68, 10}}, color = {0, 0, 127}));
        connect(optimal_dli.y, missing_dli.u2) annotation(
          Line(points = {{-47, -30}, {-41, -30}, {-41, 4}, {-33, 4}}, color = {0, 0, 127}));
        connect(dli_end_of_day.y, missing_dli.u1) annotation(
          Line(points = {{-45, 10}, {-41, 10}, {-41, 16}, {-33, 16}}, color = {0, 0, 127}));
        connect(timer.y, onOffController.u) annotation(
          Line(points = {{51, -82}, {57.5, -82}, {57.5, 4}, {62, 4}}, color = {0, 0, 127}));
        connect(ppfd.y, gain.u) annotation(
          Line(points = {{-9, -30}, {-2, -30}}, color = {0, 0, 127}));
        connect(gain.y, time_to_achieve_opt_dli.u2) annotation(
          Line(points = {{21, -30}, {26, -30}, {26, 3}, {30, 3}, {30, 4}}, color = {0, 0, 127}));
        connect(time_to_achieve_opt_dli.y, onOffController.reference) annotation(
          Line(points = {{53, 10}, {56.5, 10}, {56.5, 16}, {62, 16}}, color = {0, 0, 127}));
        connect(onOffController.y, booleanToReal.u) annotation(
          Line(points = {{85, 10}, {97, 10}}, color = {255, 0, 255}));
        connect(booleanToReal.y, p) annotation(
          Line(points = {{121, 10}, {170, 10}}, color = {0, 0, 127}));
        connect(booleanToReal.y, heat_loss_efficiency.u) annotation(
          Line(points = {{121, 10}, {143, 10}, {143, -10}, {67, -10}, {67, -32}, {78, -32}}, color = {0, 0, 127}));
        connect(heat_loss_efficiency.y, preHeaFlo.Q_flow) annotation(
          Line(points = {{101, -32}, {114, -32}}, color = {0, 0, 127}));
        connect(preHeaFlo.port, port) annotation(
          Line(points = {{134, -32}, {160, -32}}, color = {191, 0, 0}));
        connect(missing_dli.y, limiter.u) annotation(
          Line(points = {{-9, 10}, {-2, 10}}, color = {0, 0, 127}));
        connect(limiter.y, time_to_achieve_opt_dli.u1) annotation(
          Line(points = {{21, 10}, {25.5, 10}, {25.5, 16}, {30, 16}}, color = {0, 0, 127}));
        connect(limiter.y, realToBoolean.u) annotation(
          Line(points = {{21, 10}, {23, 10}, {23, -12}, {-35, -12}, {-35, -54}, {-60.5, -54}, {-60.5, -76}, {-50, -76}}, color = {0, 0, 127}));
        connect(realToBoolean.y, xor.u1) annotation(
          Line(points = {{-27, -76}, {-17, -76}, {-17, -82}, {-9, -82}}, color = {255, 0, 255}));
        connect(xor.y, timer.u) annotation(
          Line(points = {{15, -82}, {28, -82}}, color = {255, 0, 255}));
        connect(booleanPulse.y, xor.u2) annotation(
          Line(points = {{-27, -104}, {-18, -104}, {-18, -90}, {-9, -90}}, color = {255, 0, 255}));
        connect(dli.u, shaded_radiation) annotation(
          Line(points = {{-102, 10}, {-140, 10}}, color = {0, 0, 127}));
        connect(onOffController.y, booleanToReal1.u) annotation(
          Line(points = {{85, 10}, {91, 10}, {91, 50}, {97, 50}}, color = {255, 0, 255}));
        connect(booleanToReal1.y, ppfd_out) annotation(
          Line(points = {{121, 50}, {170, 50}}, color = {0, 0, 127}));
        annotation(
          Icon(graphics = {Bitmap(origin = {0, 100}, extent = {{-100, -200}, {100, 0}}, imageSource = "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")}),
          Diagram(coordinateSystem(extent = {{-160, 60}, {300, -120}})));
      end lighting;
    equation

    end illumination_control;

    model atmosphere_control
      model window
        extends Buildings.Fluid.Interfaces.PartialTwoPortTransport;
        replaceable package Medium = Modelica.Media.Interfaces.PartialMedium "Medium in the component" annotation(
          choices(choice(redeclare package Medium = Buildings.Media.Air "Moist air")));
        //  parameter Integer nPorts=4 "Number of ports"
        //    annotation(Evaluate=true, Dialog(connectorSizing=true, tab="General",group="Ports"));
        Modelica.Blocks.Interfaces.RealInput T annotation(
          Placement(transformation(origin = {0, 60}, extent = {{140, -20}, {100, 20}}, rotation = -0), iconTransformation(origin = {0, 80}, extent = {{140, -20}, {100, 20}}, rotation = -0)));
        parameter Modelica.Units.SI.Length window_width "width of window area";
        Buildings.Airflow.Multizone.MediumColumn air_column_top(redeclare package Medium = Buildings.Media.Air "Moist air", densitySelection = Buildings.Airflow.Multizone.Types.densitySelection.actual, h = 22.5/2) annotation(
          Placement(transformation(origin = {50, 50}, extent = {{-10, -10}, {10, 10}})));
        Buildings.Airflow.Multizone.DoorOperable window_top(redeclare package Medium = Buildings.Media.Air "Moist air", wOpe = window_width, hOpe = 0.2, LClo = 0.01) annotation(
          Placement(transformation(origin = {-10, 70}, extent = {{-10, -10}, {10, 10}})));
        Buildings.Airflow.Multizone.MediumColumn air_column_bottom(redeclare package Medium = Buildings.Media.Air "Moist air", densitySelection = Buildings.Airflow.Multizone.Types.densitySelection.actual, h = 22.5/2) annotation(
          Placement(transformation(origin = {50, -50}, extent = {{-10, -10}, {10, 10}})));
        Buildings.Airflow.Multizone.DoorOperable window_bottom(redeclare package Medium = Buildings.Media.Air "Moist air", wOpe = window_width, hOpe = 0.2, LClo = 0.01) annotation(
          Placement(transformation(origin = {-10, -70}, extent = {{-10, -10}, {10, 10}})));
        Buildings.Controls.Continuous.LimPID conPID(controllerType = Modelica.Blocks.Types.SimpleController.P) annotation(
          Placement(transformation(origin = {-10, 0}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
        Modelica.Blocks.Sources.Constant temperature_setpoint(k = 273.15 + 24) annotation(
          Placement(transformation(origin = {-10, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
        //protected
        //  Medium.ThermodynamicState sta_a=Medium.setState_phX(
        //      farm_air.p,
        //      actualStream(farm_air.h_outflow),
        //      actualStream(farm_air.Xi_outflow)) "Medium properties in port_a";
        //  Medium.MassFraction Xi[Medium.nXi] "Mass fraction used to compute density";
        //initial equation
        //  // The next assert tests for all allowed values of the enumeration.
        //  // Testing against densitySelection > 0 gives an error in OpenModelica as enumerations start with 1.
        //  assert(densitySelection == Buildings.Airflow.Multizone.Types.densitySelection.fromTop
        //     or densitySelection == Buildings.Airflow.Multizone.Types.densitySelection.fromBottom
        //     or densitySelection == Buildings.Airflow.Multizone.Types.densitySelection.actual,
        //    "You need to set the parameter \"densitySelection\" for the model MediumColumn.");
      equation
        connect(temperature_setpoint.y, conPID.u_m) annotation(
          Line(points = {{-10, -19}, {-10, -12}}, color = {0, 0, 127}));
        connect(conPID.y, window_top.y) annotation(
          Line(points = {{-20, 0}, {-40, 0}, {-40, 70}, {-20, 70}}, color = {0, 0, 127}));
        connect(conPID.y, window_bottom.y) annotation(
          Line(points = {{-20, 0}, {-40, 0}, {-40, -70}, {-20, -70}}, color = {0, 0, 127}));
        connect(window_top.port_b1, air_column_top.port_a) annotation(
          Line(points = {{0, 76}, {50, 76}, {50, 60}}, color = {0, 127, 255}));
        connect(window_top.port_a2, air_column_top.port_a) annotation(
          Line(points = {{0, 64}, {50, 64}, {50, 60}}, color = {0, 127, 255}));
        connect(window_bottom.port_a2, air_column_bottom.port_b) annotation(
          Line(points = {{0, -76}, {50, -76}, {50, -60}}, color = {0, 127, 255}));
        connect(window_bottom.port_b1, air_column_bottom.port_b) annotation(
          Line(points = {{0, -64}, {50, -64}, {50, -60}}, color = {0, 127, 255}));
        connect(conPID.u_s, T) annotation(
          Line(points = {{2, 0}, {30, 0}, {30, 20}, {80, 20}, {80, 60}, {120, 60}}, color = {0, 0, 127}));
        connect(port_b, air_column_top.port_b) annotation(
          Line(points = {{100, 0}, {50, 0}, {50, 40}}));
        connect(port_b, air_column_bottom.port_a) annotation(
          Line(points = {{100, 0}, {50, 0}, {50, -40}}));
        connect(port_a, window_top.port_a1) annotation(
          Line(points = {{-100, 0}, {-60, 0}, {-60, 76}, {-20, 76}}));
        connect(window_top.port_b2, port_a) annotation(
          Line(points = {{-20, 64}, {-60, 64}, {-60, 0}, {-100, 0}}, color = {0, 127, 255}));
        connect(window_bottom.port_a1, port_a) annotation(
          Line(points = {{-20, -64}, {-60, -64}, {-60, 0}, {-100, 0}}, color = {0, 127, 255}));
        connect(window_bottom.port_b2, port_a) annotation(
          Line(points = {{-20, -76}, {-60, -76}, {-60, 0}, {-100, 0}}, color = {0, 127, 255}));
        annotation(
          Icon(graphics = {Bitmap(extent = {{100, 100}, {-100, -100}}, imageSource = "iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAAAXNSR0IArs4c6QAABEBJREFUeF7tXdtx2zAQPCiNOJUkbsR6/bCL2FWEH5FMupE4lViNxEyoSJFsEST3eAQIeTXjGX8AuLvdu8VDEOmEn6gIuKjWaVxIQOQkIAEkIDICkc2zAkhAZAQim2cFkIDICEQ2zwogAZERiGyeFZA6AYtFdvM6k4Vz8kVEvo4TT50n1ThD46PuRGT/537LQ1Hk9f/qz6AKuFtl987JN7V1tGPNwSCPUYOd7XdVJeXTNr/vbOlpoA5nvs5+XmT8AaDp4aSFp3e/53KT3/ZufdZQRcBilS0qJ49TAHoKPtR4VpU8aCoBJqDW/GomLxOTAk3ywX06yN45J8viR/6MDAwTEFz3kWiM2vqA7qo2TRXABMzX2ePfFcCiMdYuD00AslsRjeAuPBdoCHgRkZu3WNqBcj7uCACZpMD/Qc4d/Pf/rtzknxEjGgIOC/Jm0MtNDo+JODxfZ60bgtTsw2BdGwAI+XVb6/hJAMgACaAEXZcGgwVACbKWABIAIkACJqDBbfuDMMtQ/74Htc9VUOQKJAEkAEOAc8D7OeCdIKMaiMFvvxOF7a+yqu0oHo3/qiSozoWnsc+iSMB1bQQHVUDTchAtQVgCJrAM9vvspNx8hzCFGo9xGnhdBIigCUgCwAywXoWRAFMCKEGwBID4xz8NvVtnVVvZoBoYG4DY9odL0EfbiBmvwtIlwHMkmloFDifgXQ2nBoBWgnxH4mj8JABkYNrL0MnfpALRbmreEWP0CjAIMekhwhBglelW4/goizB+GAIuAh7nbqi+FAB/jEkKSAAQpB7J5HoGJCABbIyzuzXig60wBIQMLAGez10MQ0BioOzdbUiaMfIoOgGoAyiX1huh2Pa5EwYZ8CZA0DmgxWlWAPYLIVaAqgJ4N9QLW2oVyApQVYC/E5oAJIAEYAhwGWr8nSgG/wQu5xrHTwkCM8C6As0IOG7r0UkIjN/8Xk5s+2YEHM9byi22EYkNQGz7tgQIfjk1NgA6+9yIcSPmQ4BzACbBlCBQg6axCmr5JoMVEKwC+MCmpuJBE5ASlKQEHZzmj/RqIN4qAStg7N8JT+IsiJOw2T6Ec0DKc8DJ95MOohoIxs/DuPk6a3hwKwpjR/uOG1NWF6qsxjmPBk1AjQRdPrbei2eKF3gH+VyUm3yJpCNMgOrh3WOkGhJl37YD/Qzy8O7ux9cfMiiQjPTFNkQ75+R29MfX14EcX+DQO6ggZPSTjoFJ7hfbSpbFNi96Y3JoCEvQ0UDrY+xRL1Jtf2ITfmz9MWQ1AfUAp/mgX/Zd4LwPQNnXlDS9DxrdP3d9EAF7OapfafJp/yal+p0C+tdY+bRhLM1oIrBfQuxfYVVV8mv2KkXU11iZJuEHHWxwBXxQ3MzCJgFmUOoGIgE63Mx6kQAzKHUDkQAdbma9SIAZlLqBSIAON7NeJMAMSt1AJECHm1kvEmAGpW4gEqDDzazXH04uF45wgVW9AAAAAElFTkSuQmCC")}),
          experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-06, Interval = 0.002));
      end window;
    equation

    end atmosphere_control;

    model irrigation
    equation

    end irrigation;
  equation

    annotation(
      experiment(StartTime = 0, StopTime = 864000, Tolerance = 1e-06, Interval = 86.4));
  end farm;

  model system_analysis
    plant.plant_yield plant_yield1(a_plant = 1, x_nsdw(start = 0.25*0.72), x_sdw(start = 0.75*0.72)) annotation(
      Placement(transformation(origin = {10, 10}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Sources.Constant const(k = 20 + 273.15) annotation(
      Placement(transformation(origin = {-150, 130}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Sources.Constant const1(k = 600) annotation(
      Placement(transformation(origin = {-150, 100}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Sources.Constant const2(k = 200/2.02) annotation(
      Placement(transformation(origin = {-150, 70}, extent = {{-10, -10}, {10, 10}})));
    Buildings.BoundaryConditions.WeatherData.ReaderTMY3 weaDat1(filNam = "/home/marci/Downloads/DEU_BY_Nurnberg.AP.107630_TMYx.2007-2021/DEU_BY_Nurnberg.AP.107630_TMYx.2007-2021.mos") annotation(
      Placement(transformation(origin = {-80, -20}, extent = {{-60, 20}, {-40, 40}})));
    Buildings.BoundaryConditions.SolarIrradiation.DirectTiltedSurface HDirRoo11(azi = 0, til(displayUnit = "deg") = 1.5707963267948966) annotation(
      Placement(transformation(origin = {-120, -70}, extent = {{20, 20}, {40, 40}})));
    Buildings.BoundaryConditions.SolarIrradiation.DiffuseIsotropic HDifTilIso11(til = 1.5707963267948966) annotation(
      Placement(transformation(origin = {-90, -60}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Math.Add add1 annotation(
      Placement(transformation(origin = {-50, -50}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Math.Gain gain(k = 1/18) annotation(
      Placement(transformation(origin = {50, 10}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Math.Gain gain1(k = 1/7) annotation(
      Placement(transformation(origin = {-20, -50}, extent = {{-10, -10}, {10, 10}})));
  equation
    connect(HDirRoo11.H, add1.u1) annotation(
      Line(points = {{-79, -40}, {-73, -40}, {-73, -44}, {-63, -44}}, color = {0, 0, 127}));
    connect(HDifTilIso11.H, add1.u2) annotation(
      Line(points = {{-79, -60}, {-73, -60}, {-73, -56}, {-63, -56}}, color = {0, 0, 127}));
    connect(weaDat1.weaBus, HDirRoo11.weaBus) annotation(
      Line(points = {{-120, 10}, {-110, 10}, {-110, -40}, {-100, -40}}, color = {255, 204, 51}, thickness = 0.5));
    connect(weaDat1.weaBus, HDifTilIso11.weaBus) annotation(
      Line(points = {{-120, 10}, {-110, 10}, {-110, -60}, {-100, -60}}, color = {255, 204, 51}, thickness = 0.5));
    connect(const1.y, plant_yield1.co2_concentration) annotation(
      Line(points = {{-138, 100}, {-44, 100}, {-44, 10}, {-2, 10}}, color = {0, 0, 127}));
    connect(weaDat1.weaBus.TDryBul, plant_yield1.air_temperature) annotation(
      Line(points = {{-120, 10}, {-110, 10}, {-110, 18}, {-2, 18}}, color = {0, 0, 127}));
    connect(plant_yield1.y, gain.u) annotation(
      Line(points = {{22, 10}, {38, 10}}, color = {0, 0, 127}));
    connect(add1.y, gain1.u) annotation(
      Line(points = {{-38, -50}, {-32, -50}}, color = {0, 0, 127}));
    connect(gain1.y, plant_yield1.par) annotation(
      Line(points = {{-8, -50}, {-6, -50}, {-6, 2}, {-2, 2}}, color = {0, 0, 127}));
    annotation(
      Diagram(coordinateSystem(extent = {{-160, 140}, {60, -80}})));
  end system_analysis;

  model transpiration_validation
    plant.plant_transpiration plant_transpiration1(a_plant = 1) annotation(
      Placement(transformation(origin = {30, 10}, extent = {{-10, -10}, {10, 10}})));
  equation

  end transpiration_validation;

  model building
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a side_sw annotation(
      Placement(transformation(origin = {-20, -80}, extent = {{-110, -10}, {-90, 10}}), iconTransformation(origin = {0, -80}, extent = {{-110, -10}, {-90, 10}})));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a side_nw annotation(
      Placement(transformation(origin = {-20, 40}, extent = {{-110, -10}, {-90, 10}}), iconTransformation(origin = {0, 40}, extent = {{-110, -10}, {-90, 10}})));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a side_ne annotation(
      Placement(transformation(origin = {20, 40}, extent = {{110, -10}, {90, 10}}, rotation = -0), iconTransformation(origin = {0, 40}, extent = {{110, -10}, {90, 10}})));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a side_se annotation(
      Placement(transformation(origin = {20, -80}, extent = {{110, -10}, {90, 10}}, rotation = -0), iconTransformation(origin = {0, -82}, extent = {{110, -10}, {90, 10}})));
    Modelica.Blocks.Interfaces.RealInput v_air_sw(unit = "m/s") "Air Speed for the south-west side" annotation(
      Placement(transformation(origin = {-20, -140}, extent = {{-140, 80}, {-100, 120}}), iconTransformation(origin = {0, -140}, extent = {{-140, 80}, {-100, 120}})));
    Modelica.Blocks.Interfaces.RealInput v_air_nw(unit = "m/s") "Air Speed for the north-west side" annotation(
      Placement(transformation(origin = {-20, -20}, extent = {{-140, 80}, {-100, 120}}), iconTransformation(origin = {0, -20}, extent = {{-140, 80}, {-100, 120}})));
    Modelica.Blocks.Interfaces.RealInput v_air_ne(unit = "m/s") "Air Speed for the north-east side" annotation(
      Placement(transformation(origin = {20, -20}, extent = {{140, 80}, {100, 120}}, rotation = -0), iconTransformation(origin = {0, -20}, extent = {{140, 80}, {100, 120}})));
    Modelica.Blocks.Interfaces.RealInput v_air_se(unit = "m/s") "Air Speed for the south-east side" annotation(
      Placement(transformation(origin = {20, -140}, extent = {{140, 80}, {100, 120}}, rotation = -0), iconTransformation(origin = {0, -140}, extent = {{140, 80}, {100, 120}})));
    Buildings.HeatTransfer.Sources.FixedTemperature air_temperature(T = 293.15) annotation(
      Placement(transformation(extent = {{-10, -10}, {10, 10}})));

    model facade
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a exterior annotation(
        Placement(transformation(extent = {{-110, -10}, {-90, 10}})));
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b interior annotation(
        Placement(transformation(extent = {{90, -10}, {110, 10}})));
      Modelica.Blocks.Interfaces.RealInput v_air(unit = "m/s") "Air Speed" annotation(
        Placement(transformation(origin = {0, -40}, extent = {{-140, 80}, {-100, 120}}), iconTransformation(origin = {10, -40}, extent = {{-140, 80}, {-100, 120}})));
      parameter Modelica.Units.SI.Area area "facade area";
      parameter Modelica.Units.SI.Angle azimuth "Surface azimuth";
      Buildings.HeatTransfer.Convection.Exterior convection_out(A = area, azi = azimuth, roughness = Buildings.HeatTransfer.Types.SurfaceRoughness.Medium, til = 1.5707963267948966) annotation(
        Placement(transformation(origin = {-40, 0}, extent = {{10, -10}, {-10, 10}})));
      Buildings.HeatTransfer.Convection.Interior convection_in(A = area, conMod = Buildings.HeatTransfer.Types.InteriorConvection.Temperature, til = 1.5707963267948966) annotation(
        Placement(transformation(origin = {40, 0}, extent = {{-10, -10}, {10, 10}})));
      Buildings.HeatTransfer.Conduction.MultiLayer conduction(A = area, layers(nLay = 2, material = {concrete, eps})) annotation(
        Placement(transformation(extent = {{10, -10}, {-10, 10}})));
      Modelica.Blocks.Sources.Constant wind_direction(k = 3.14159 - azimuth) annotation(
        Placement(transformation(origin = {0, 50}, extent = {{10, -10}, {-10, 10}})));
      Buildings.HeatTransfer.Data.Solids.Concrete concrete(c = 550, d = 3500, k = 2, x = 0.18) annotation(
        Placement(transformation(origin = {0, -40}, extent = {{-10, -10}, {10, 10}})));
      Buildings.HeatTransfer.Data.Solids.InsulationBoard eps(c = 1200, d = 21, k = 0.035, x = 0.05) annotation(
        Placement(transformation(origin = {40, -40}, extent = {{-10, -10}, {10, 10}})));
    equation
      connect(convection_out.solid, conduction.port_b) annotation(
        Line(points = {{-30, 0}, {-10, 0}}, color = {191, 0, 0}));
      connect(conduction.port_a, convection_in.solid) annotation(
        Line(points = {{10, 0}, {30, 0}}, color = {191, 0, 0}));
      connect(wind_direction.y, convection_out.dir) annotation(
        Line(points = {{-11, 50}, {-20, 50}, {-20, 6}, {-28, 6}}, color = {0, 0, 127}));
      connect(convection_in.fluid, interior) annotation(
        Line(points = {{50, 0}, {100, 0}}, color = {191, 0, 0}));
      connect(convection_out.fluid, exterior) annotation(
        Line(points = {{-50, 0}, {-100, 0}}, color = {191, 0, 0}));
      connect(convection_out.v, v_air) annotation(
        Line(points = {{-28, 10}, {-24, 10}, {-24, 60}, {-120, 60}}, color = {0, 0, 127}));
      annotation(
        Diagram,
        Icon(graphics = {Rectangle(lineColor = {154, 153, 150}, fillColor = {246, 245, 244}, fillPattern = FillPattern.Solid, lineThickness = 1.5, extent = {{-90, 100}, {90, -100}})}));
    end facade;

    model window
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a exterior annotation(
        Placement(transformation(extent = {{-110, -10}, {-90, 10}})));
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b interior annotation(
        Placement(transformation(extent = {{90, -10}, {110, 10}})));
      Modelica.Blocks.Interfaces.RealInput v_air(unit = "m/s") "Air Speed" annotation(
        Placement(transformation(origin = {0, -40}, extent = {{-140, 80}, {-100, 120}}), iconTransformation(origin = {10, -40}, extent = {{-140, 80}, {-100, 120}})));
      parameter Modelica.Units.SI.Area area "window area";
      parameter Modelica.Units.SI.Angle azimuth "Surface azimuth";
      Buildings.HeatTransfer.Convection.Exterior convection_out_glass_sw(A = area, azi = azimuth, roughness = Buildings.HeatTransfer.Types.SurfaceRoughness.VerySmooth, til = 1.5707963267948966) annotation(
        Placement(transformation(origin = {-40, 0}, extent = {{10, -10}, {-10, 10}})));
      Buildings.HeatTransfer.Convection.Interior convection_in_glass_sw(A = area, conMod = Buildings.HeatTransfer.Types.InteriorConvection.Temperature, til = 1.5707963267948966) annotation(
        Placement(transformation(origin = {40, 0}, extent = {{-10, -10}, {10, 10}})));
      Buildings.HeatTransfer.Conduction.SingleLayer conduction_glass_sw(A = area, material = glass) annotation(
        Placement(transformation(extent = {{-10, -10}, {10, 10}})));
      Buildings.HeatTransfer.Data.Solids.Glass glass(x = 25) annotation(
        Placement(transformation(origin = {2, -40}, extent = {{-10, -10}, {10, 10}})));
      Modelica.Blocks.Sources.Constant wind_direction(k = 3.14159 - azimuth) annotation(
        Placement(transformation(origin = {0, 50}, extent = {{10, -10}, {-10, 10}})));
    equation
      connect(wind_direction.y, convection_out_glass_sw.dir) annotation(
        Line(points = {{-10, 50}, {-20, 50}, {-20, 6}, {-28, 6}}, color = {0, 0, 127}));
      connect(convection_out_glass_sw.solid, conduction_glass_sw.port_a) annotation(
        Line(points = {{-30, 0}, {-10, 0}}, color = {191, 0, 0}));
      connect(conduction_glass_sw.port_b, convection_in_glass_sw.solid) annotation(
        Line(points = {{10, 0}, {30, 0}}, color = {191, 0, 0}));
      connect(convection_in_glass_sw.fluid, interior) annotation(
        Line(points = {{50, 0}, {100, 0}}, color = {191, 0, 0}));
      connect(convection_out_glass_sw.fluid, exterior) annotation(
        Line(points = {{-50, 0}, {-100, 0}}, color = {191, 0, 0}));
      connect(convection_out_glass_sw.v, v_air) annotation(
        Line(points = {{-28, 10}, {-24, 10}, {-24, 60}, {-120, 60}}, color = {0, 0, 127}));
      annotation(
        Icon(graphics = {Rectangle(lineColor = {245, 194, 17}, fillColor = {225, 246, 255}, fillPattern = FillPattern.Solid, lineThickness = 8, extent = {{-90, 100}, {90, -100}})}));
    end window;

    facade facade_sw(area = 523.2, azimuth = 0.7330382858376184) annotation(
      Placement(transformation(origin = {-50, -50}, extent = {{-10, -10}, {10, 10}})));
    window window_sw(area = 581.04, azimuth = 0.7330382858376184) annotation(
      Placement(transformation(origin = {-50, -30}, extent = {{-10, -10}, {10, 10}})));
    greenshell.building.facade facade_sw1(area = 369.72, azimuth = 2.303834612632515, redeclare Buildings.HeatTransfer.Data.OpaqueConstructions.Insulation100Concrete200 layers) annotation(
      Placement(transformation(origin = {-70, 50}, extent = {{-10, -10}, {10, 10}})));
    greenshell.building.window window_sw1(area = 118.8, azimuth = 2.303834612632515) annotation(
      Placement(transformation(origin = {-70, 70}, extent = {{-10, -10}, {10, 10}})));
    greenshell.building.facade facade_sw2(area = 523.2, azimuth = 3.8746309394274117) annotation(
      Placement(transformation(origin = {50, 30}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
    greenshell.building.window window_sw2(area = 581.04, azimuth = 3.8746309394274117) annotation(
      Placement(transformation(origin = {50, 50}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
    greenshell.building.facade facade_sw3(area = 369.71, azimuth = 5.445427266222308) annotation(
      Placement(transformation(origin = {70, -70}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
    greenshell.building.window window_sw3(area = 118.8, azimuth = 5.445427266222308) annotation(
      Placement(transformation(origin = {70, -50}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));

    model facade_test
      Buildings.HeatTransfer.Convection.Interior building_facade_building_convection(A = 523.2, conMod = Buildings.HeatTransfer.Types.InteriorConvection.Temperature, til = 1.5707963267948966) annotation(
        Placement(transformation(origin = {40, 0}, extent = {{-10, -10}, {10, 10}})));
      Buildings.HeatTransfer.Convection.Exterior building_facade_farm_convection(A = 523.2, azi = 0, conMod = Buildings.HeatTransfer.Types.ExteriorConvection.TemperatureWind, roughness = Buildings.HeatTransfer.Types.SurfaceRoughness.Medium, til = 1.5707963267948966) annotation(
        Placement(transformation(origin = {-40, 0}, extent = {{10, -10}, {-10, 10}})));
      Buildings.HeatTransfer.Conduction.MultiLayer building_facade_conduction(A = 523.2, layers = datOpaCon) annotation(
        Placement(transformation(extent = {{10, -10}, {-10, 10}})));
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a exterior annotation(
        Placement(transformation(extent = {{-110, -10}, {-90, 10}}), iconTransformation(extent = {{-110, -10}, {-90, 10}})));
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b interior annotation(
        Placement(transformation(extent = {{90, -10}, {110, 10}}), iconTransformation(extent = {{90, -10}, {110, 10}})));
      Modelica.Blocks.Interfaces.RealInput v_air(unit = "m/s") annotation(
        Placement(transformation(origin = {0, -40}, extent = {{-140, 80}, {-100, 120}}), iconTransformation(origin = {10, -40}, extent = {{-140, 80}, {-100, 120}})));
      Modelica.Blocks.Sources.Constant wind_direction(k = 3.14159 - azimuth) annotation(
        Placement(transformation(origin = {0, 50}, extent = {{10, -10}, {-10, 10}})));
      parameter Buildings.HeatTransfer.Data.OpaqueConstructions.Insulation100Concrete200 datOpaCon annotation(
        Placement(transformation(origin = {0, -42}, extent = {{-10, -10}, {10, 10}})));
      Buildings.HeatTransfer.Data.Solids.Concrete concrete(c = 550, d = 3500, k = 2, x = 0.18) annotation(
        Placement(transformation(origin = {0, -68}, extent = {{-10, -10}, {10, 10}})));
      Buildings.HeatTransfer.Data.Solids.InsulationBoard eps(c = 1200, d = 21, k = 0.035, x = 0.05) annotation(
        Placement(transformation(origin = {40, -68}, extent = {{-10, -10}, {10, 10}})));
      Buildings.HeatTransfer.Data.OpaqueConstructions.Generic facade(material = {concrete, eps}, nLay = 2, roughness_a = Buildings.HeatTransfer.Types.SurfaceRoughness.Medium) annotation(
        Placement(transformation(origin = {-40, -68}, extent = {{-10, -10}, {10, 10}})));
    equation
      connect(building_facade_farm_convection.solid, building_facade_conduction.port_b) annotation(
        Line(points = {{-30, 0}, {-10, 0}}, color = {191, 0, 0}));
      connect(building_facade_conduction.port_a, building_facade_building_convection.solid) annotation(
        Line(points = {{10, 0}, {30, 0}}, color = {191, 0, 0}));
      connect(building_facade_building_convection.fluid, interior) annotation(
        Line(points = {{50, 0}, {100, 0}}, color = {191, 0, 0}));
      connect(building_facade_farm_convection.fluid, exterior) annotation(
        Line(points = {{-50, 0}, {-100, 0}}, color = {191, 0, 0}));
      connect(wind_direction.y, building_facade_farm_convection.dir) annotation(
        Line(points = {{-10, 50}, {-20, 50}, {-20, 6}, {-28, 6}}, color = {0, 0, 127}));
      connect(building_facade_farm_convection.v, v_air) annotation(
        Line(points = {{-28, 10}, {-22, 10}, {-22, 60}, {-120, 60}}, color = {0, 0, 127}));
      annotation(
        Diagram,
        Icon(graphics = {Rectangle(lineColor = {154, 153, 150}, fillColor = {246, 245, 244}, fillPattern = FillPattern.Solid, lineThickness = 1.5, extent = {{-90, 100}, {90, -100}})}));
    end facade_test;

    model facade2
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a exterior annotation(
        Placement(transformation(extent = {{-110, -10}, {-90, 10}}), iconTransformation(extent = {{-110, -10}, {-90, 10}})));
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b interior annotation(
        Placement(transformation(extent = {{90, -10}, {110, 10}}), iconTransformation(extent = {{90, -10}, {110, 10}})));
      Modelica.Blocks.Interfaces.RealInput v_air(unit = "m/s") annotation(
        Placement(transformation(origin = {0, -40}, extent = {{-140, 80}, {-100, 120}}), iconTransformation(origin = {10, -40}, extent = {{-140, 80}, {-100, 120}})));
      parameter Modelica.Units.SI.Area area "facade area";
      parameter Modelica.Units.SI.Angle azimuth "Surface azimuth";
      Buildings.HeatTransfer.Convection.Exterior convection_out(A = area, azi = azimuth, roughness = Buildings.HeatTransfer.Types.SurfaceRoughness.Medium, til = 1.5707963267948966) annotation(
        Placement(transformation(origin = {-40, 0}, extent = {{10, -10}, {-10, 10}})));
      Buildings.HeatTransfer.Convection.Interior convection_in(A = area, conMod = Buildings.HeatTransfer.Types.InteriorConvection.Temperature, til = 1.5707963267948966) annotation(
        Placement(transformation(origin = {40, 0}, extent = {{-10, -10}, {10, 10}})));
      Modelica.Blocks.Sources.Constant wind_direction(k = 3.14159 - azimuth) annotation(
        Placement(transformation(origin = {0, 50}, extent = {{10, -10}, {-10, 10}})));
      Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor(G = 346.4900662251656) annotation(
        Placement(transformation(origin = {2, 0}, extent = {{-10, -10}, {10, 10}})));
    equation
      connect(wind_direction.y, convection_out.dir) annotation(
        Line(points = {{-11, 50}, {-20, 50}, {-20, 6}, {-28, 6}}, color = {0, 0, 127}));
      connect(convection_in.fluid, interior) annotation(
        Line(points = {{50, 0}, {100, 0}}, color = {191, 0, 0}));
      connect(convection_out.fluid, exterior) annotation(
        Line(points = {{-50, 0}, {-100, 0}}, color = {191, 0, 0}));
      connect(convection_out.v, v_air) annotation(
        Line(points = {{-28, 10}, {-24, 10}, {-24, 60}, {-120, 60}}, color = {0, 0, 127}));
      connect(convection_out.solid, thermalConductor.port_a) annotation(
        Line(points = {{-30, 0}, {-8, 0}}, color = {191, 0, 0}));
      connect(thermalConductor.port_b, convection_in.solid) annotation(
        Line(points = {{12, 0}, {30, 0}}, color = {191, 0, 0}));
    end facade2;

    model window2
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a exterior annotation(
        Placement(transformation(extent = {{-110, -10}, {-90, 10}}), iconTransformation(extent = {{-110, -10}, {-90, 10}})));
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b interior annotation(
        Placement(transformation(extent = {{90, -10}, {110, 10}}), iconTransformation(extent = {{90, -10}, {110, 10}})));
      Modelica.Blocks.Interfaces.RealInput v_air(unit = "m/s") annotation(
        Placement(transformation(origin = {0, -40}, extent = {{-140, 80}, {-100, 120}}), iconTransformation(origin = {10, -40}, extent = {{-140, 80}, {-100, 120}})));
      parameter Modelica.Units.SI.Area area "window area";
      parameter Modelica.Units.SI.Angle azimuth "Surface azimuth";
      Buildings.HeatTransfer.Convection.Exterior convection_out_glass_sw(A = area, azi = azimuth, roughness = Buildings.HeatTransfer.Types.SurfaceRoughness.VerySmooth, til = 1.5707963267948966) annotation(
        Placement(transformation(origin = {-40, 0}, extent = {{10, -10}, {-10, 10}})));
      Buildings.HeatTransfer.Convection.Interior convection_in_glass_sw(A = area, conMod = Buildings.HeatTransfer.Types.InteriorConvection.Temperature, til = 1.5707963267948966) annotation(
        Placement(transformation(origin = {40, 0}, extent = {{-10, -10}, {10, 10}})));
      Modelica.Blocks.Sources.Constant wind_direction(k = 3.14159 - azimuth) annotation(
        Placement(transformation(origin = {0, 50}, extent = {{10, -10}, {-10, 10}})));
      Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor(G = 1154.3841059602646) annotation(
        Placement(transformation(extent = {{-10, -10}, {10, 10}})));
    equation
      connect(wind_direction.y, convection_out_glass_sw.dir) annotation(
        Line(points = {{-10, 50}, {-20, 50}, {-20, 6}, {-28, 6}}, color = {0, 0, 127}));
      connect(convection_in_glass_sw.fluid, interior) annotation(
        Line(points = {{50, 0}, {100, 0}}, color = {191, 0, 0}));
      connect(convection_out_glass_sw.fluid, exterior) annotation(
        Line(points = {{-50, 0}, {-100, 0}}, color = {191, 0, 0}));
      connect(convection_out_glass_sw.v, v_air) annotation(
        Line(points = {{-28, 10}, {-24, 10}, {-24, 60}, {-120, 60}}, color = {0, 0, 127}));
      connect(convection_out_glass_sw.solid, thermalConductor.port_a) annotation(
        Line(points = {{-30, 0}, {-10, 0}}, color = {191, 0, 0}));
      connect(thermalConductor.port_b, convection_in_glass_sw.solid) annotation(
        Line(points = {{10, 0}, {30, 0}}, color = {191, 0, 0}));
    end window2;
  equation
    connect(window_sw.interior, facade_sw.interior) annotation(
      Line(points = {{-40, -30}, {-20, -30}, {-20, -50}, {-40, -50}}, color = {191, 0, 0}));
    connect(window_sw3.interior, facade_sw3.interior) annotation(
      Line(points = {{60, -50}, {40, -50}, {40, -70}, {60, -70}}, color = {191, 0, 0}));
    connect(window_sw2.interior, facade_sw2.interior) annotation(
      Line(points = {{40, 50}, {20, 50}, {20, 30}, {40, 30}}, color = {191, 0, 0}));
    connect(window_sw1.interior, facade_sw1.interior) annotation(
      Line(points = {{-60, 70}, {-40, 70}, {-40, 50}, {-60, 50}}, color = {191, 0, 0}));
    connect(facade_sw1.interior, air_temperature.port) annotation(
      Line(points = {{-60, 50}, {0, 50}, {0, 20}, {20, 20}, {20, 0}, {10, 0}}, color = {191, 0, 0}));
    connect(facade_sw2.interior, air_temperature.port) annotation(
      Line(points = {{40, 30}, {20, 30}, {20, 0}, {10, 0}}, color = {191, 0, 0}));
    connect(window_sw3.interior, air_temperature.port) annotation(
      Line(points = {{60, -50}, {20, -50}, {20, 0}, {10, 0}}, color = {191, 0, 0}));
    connect(window_sw.interior, air_temperature.port) annotation(
      Line(points = {{-40, -30}, {20, -30}, {20, 0}, {10, 0}}, color = {191, 0, 0}));
    connect(facade_sw1.exterior, side_nw) annotation(
      Line(points = {{-80, 50}, {-100, 50}, {-100, 40}, {-120, 40}}, color = {191, 0, 0}));
    connect(window_sw1.exterior, side_nw) annotation(
      Line(points = {{-80, 70}, {-100, 70}, {-100, 40}, {-120, 40}}, color = {191, 0, 0}));
    connect(window_sw1.v_air, v_air_nw) annotation(
      Line(points = {{-80, 76}, {-90, 76}, {-90, 80}, {-140, 80}}, color = {0, 0, 127}));
    connect(facade_sw1.v_air, v_air_nw) annotation(
      Line(points = {{-80, 56}, {-90, 56}, {-90, 80}, {-140, 80}}, color = {0, 0, 127}));
    connect(facade_sw2.exterior, side_ne) annotation(
      Line(points = {{60, 30}, {90, 30}, {90, 40}, {120, 40}}, color = {191, 0, 0}));
    connect(window_sw2.exterior, side_ne) annotation(
      Line(points = {{60, 50}, {90, 50}, {90, 40}, {120, 40}}, color = {191, 0, 0}));
    connect(facade_sw2.v_air, v_air_ne) annotation(
      Line(points = {{62, 36}, {80, 36}, {80, 80}, {140, 80}}, color = {0, 0, 127}));
    connect(window_sw2.v_air, v_air_ne) annotation(
      Line(points = {{62, 56}, {80, 56}, {80, 80}, {140, 80}}, color = {0, 0, 127}));
    connect(facade_sw3.exterior, side_se) annotation(
      Line(points = {{80, -70}, {100, -70}, {100, -80}, {120, -80}}, color = {191, 0, 0}));
    connect(window_sw3.exterior, side_se) annotation(
      Line(points = {{80, -50}, {100, -50}, {100, -80}, {120, -80}}, color = {191, 0, 0}));
    connect(facade_sw3.v_air, v_air_se) annotation(
      Line(points = {{82, -64}, {92, -64}, {92, -40}, {140, -40}}, color = {0, 0, 127}));
    connect(window_sw3.v_air, v_air_se) annotation(
      Line(points = {{82, -44}, {92, -44}, {92, -40}, {140, -40}}, color = {0, 0, 127}));
    connect(facade_sw.exterior, side_sw) annotation(
      Line(points = {{-60, -50}, {-90, -50}, {-90, -80}, {-120, -80}}, color = {191, 0, 0}));
    connect(window_sw.exterior, side_sw) annotation(
      Line(points = {{-60, -30}, {-90, -30}, {-90, -80}, {-120, -80}}, color = {191, 0, 0}));
    connect(facade_sw.v_air, v_air_sw) annotation(
      Line(points = {{-60, -44}, {-80, -44}, {-80, -40}, {-140, -40}}, color = {0, 0, 127}));
    connect(window_sw.v_air, v_air_sw) annotation(
      Line(points = {{-60, -24}, {-80, -24}, {-80, -40}, {-140, -40}}, color = {0, 0, 127}));
    annotation(
      Diagram(graphics = {Polygon(origin = {9, 3}, points = {{-115, 11}, {11, -103}, {99, -15}, {-29, 97}, {-115, 11}})}, coordinateSystem(extent = {{-120, 100}, {120, -100}})),
      Icon(graphics = {Bitmap(extent = {{-100, -100}, {100, 100}}, imageSource = "iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAAAXNSR0IArs4c6QAABKpJREFUeF7tnHFS2zwQxdfhIhyFXqQ4MJ3JLQqXaDPTAgkXaY+SHqSodUhwbGQ/rXZlyWTz1zefrJX0ftqntfC0IvtlVaDKOroNTgYg8yYwAAYgswKZh7cMMACZFcg8vGWAAcisQObhLQMMQGYFMg9vGWAA4hW4vl09OUd/nh/Xd/FR8vacZQbU9erSXdATEV0d5NtVf+nTZrPe5ZWTP/rsANQ3X2pXXTTi938752g7t2yYFYDGcoioHttnztH9nCDMAoDHclCuz8aSigdQ36xqV+39nvubhSUVDSDEchCV0i2pSAARloM4FGtJxQE4tRxH5P2T3f7/DzUOoyjSkooCoGE5KBVKs6QiAERZjttnx/KF6LKq6CsSvtdejCVlBxBZ5eyqipabn+vfjbAHgL+I6JIBoghLygog0nI224f1si90A+FlQTU3G3JbUhYAMZbTnLkLR8vN43oztss/36zuuBD+J1E2S5ocgIblIJuZkyVNCqC1nGbYZk+P/Q7PLGiz/fHeckIgxFnS4v758dtk19uTAIixnEbgKsByEIjSLSk5gBDL8bxTdaocJDJq51lSRY5cU+JOUiUlBdCtckJsZy+lt8pBIqP2UqukJABay6musNe30mlYDgLBs6S3TZOsSlIHEGI5rzc8nUNY1XIQBJ4lvUVLYkmqADRfrJCI0vZSLEkFAL/Kec2AKSwHgeJZUpsNWh8BiAF0LCf8inhSy0EQYizJEe1I4SMAMYDr29XezBl39EmqHCQyao+1pO3DWqShqHOzqCMAtECtF6uQcSTPcC1pLgCKshwEiGNJcwBQpOWEQPDeJfXOuSIBnM5ROkEkVOp2ZLHS9SU/A6QTTC0wim8AkEKJ2/cARspr6QazDAAALQMS73AU3gAghRK3+wG0l4lmQVkAtIOWBcBzWEknmFhfGP40A06Xd/xv6frsEI48hA0A3Ls6D9ghrKNjdBQDEC2dTsezA4AWPCZryN+DuIcmmg83Xn/+xR3C6NVfuq+5gp0nAKnKI/0NQHDZF/whFwuXAQgGwNI1+GEDkBVARduH76xzb4Z3QafWob/g4K0+8OCZZYAWAKns8ZdnVgX1tEc7WFsw7XjzeA8QlJHagmnHmwSA5KsI6YKl/fsCDcX7sLehUgGl/UMBHJ9DlohOL1ZJ5gtW2oIH53PYslzBtNc3iQWdDjL1grUF045nAJhfMxsAK0PHjxHtHSKNJ+0/20N46I8hdgaMb+CCqyD/dTQCmjwDyvw8ffjuHgk2lvLtWsO/REsOgHkm2XuAUDBtoFaGWhmatqrS3rHa8c4wA7rnk+RM8m0dbrzJAaBDiNuOFox2LHe8oeeLvg0N+UDKvzD8JYQMwCF+/ATfTRvNBwFP/h6AJjDW7tMJLbiTAZ4A8drHvZeg9RcNIMZzfRbE+GcUkF7lZoB/Z2FL4a44NAPidzpvRmg+KNpZZAASQdJeLoBEWxAteKoq6AgNzQfBnTADJHakdxeEBOG2ZwfAnbA931VAnAEmqEwBAyDTT9zbAIgllAUwADL9xL0NgFhCWQADINNP3NsAiCWUBTAAMv3EvQ2AWEJZAAMg00/c2wCIJZQF+AcXnnmOvw3zCwAAAABJRU5ErkJggg==")}));
  end building;

  model building2
    building.facade2 facade_sw(area = 523.2, azimuth = 0.7330382858376184) annotation(
      Placement(transformation(origin = {-50, -50}, extent = {{-10, -10}, {10, 10}})));
    Buildings.HeatTransfer.Sources.FixedTemperature air_temperature(T = 293.15) annotation(
      Placement(transformation(extent = {{-10, -10}, {10, 10}})));
    building.facade2 facade_ne(area = 523.2, azimuth = 3.8746309394274117) annotation(
      Placement(transformation(origin = {50, 30}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
    building.facade2 facade_nw(area = 369.7, azimuth = 2.303834612632515) annotation(
      Placement(transformation(origin = {-50, 50}, extent = {{-10, -10}, {10, 10}})));
    building.facade2 facade_se(area = 369.7, azimuth = 5.445427266222308) annotation(
      Placement(transformation(origin = {50, -50}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
    building.window2 window_sw(area = 581.04, azimuth = 0.7330382858376184) annotation(
      Placement(transformation(origin = {-50, -28}, extent = {{-10, -10}, {10, 10}})));
    building.window2 window_nw(area = 118.8, azimuth = 2.303834612632515) annotation(
      Placement(transformation(origin = {-50, 72}, extent = {{-10, -10}, {10, 10}})));
    building.window2 window_ne(area = 581.04, azimuth = 3.8746309394274117) annotation(
      Placement(transformation(origin = {52, 54}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
    building.window2 window_se(area = 118.8, azimuth = 5.445427266222308) annotation(
      Placement(transformation(origin = {50, -26}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a side_sw annotation(
      Placement(transformation(origin = {-20, -80}, extent = {{-110, -10}, {-90, 10}}), iconTransformation(origin = {0, -80}, extent = {{-110, -10}, {-90, 10}})));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a side_nw annotation(
      Placement(transformation(origin = {-20, 40}, extent = {{-110, -10}, {-90, 10}}), iconTransformation(origin = {0, 40}, extent = {{-110, -10}, {-90, 10}})));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a side_ne annotation(
      Placement(transformation(origin = {20, 40}, extent = {{110, -10}, {90, 10}}), iconTransformation(origin = {0, 40}, extent = {{110, -10}, {90, 10}})));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a side_se annotation(
      Placement(transformation(origin = {20, -80}, extent = {{110, -10}, {90, 10}}), iconTransformation(origin = {0, -82}, extent = {{110, -10}, {90, 10}})));
    Modelica.Blocks.Interfaces.RealInput v_air_sw(unit = "m/s") annotation(
      Placement(transformation(origin = {-20, -140}, extent = {{-140, 80}, {-100, 120}}), iconTransformation(origin = {0, -140}, extent = {{-140, 80}, {-100, 120}})));
    Modelica.Blocks.Interfaces.RealInput v_air_nw(unit = "m/s") annotation(
      Placement(transformation(origin = {-20, -20}, extent = {{-140, 80}, {-100, 120}}), iconTransformation(origin = {0, -20}, extent = {{-140, 80}, {-100, 120}})));
    Modelica.Blocks.Interfaces.RealInput v_air_ne(unit = "m/s") annotation(
      Placement(transformation(origin = {20, -20}, extent = {{140, 80}, {100, 120}}), iconTransformation(origin = {0, -20}, extent = {{140, 80}, {100, 120}})));
    Modelica.Blocks.Interfaces.RealInput v_air_se(unit = "m/s") annotation(
      Placement(transformation(origin = {20, -140}, extent = {{140, 80}, {100, 120}}), iconTransformation(origin = {0, -140}, extent = {{140, 80}, {100, 120}})));
  equation
    connect(air_temperature.port, facade_ne.interior) annotation(
      Line(points = {{10, 0}, {20, 0}, {20, 30}, {40, 30}}, color = {191, 0, 0}));
    connect(air_temperature.port, window_ne.interior) annotation(
      Line(points = {{10, 0}, {20, 0}, {20, 54}, {42, 54}}, color = {191, 0, 0}));
    connect(air_temperature.port, window_se.interior) annotation(
      Line(points = {{10, 0}, {20, 0}, {20, -26}, {40, -26}}, color = {191, 0, 0}));
    connect(air_temperature.port, facade_se.interior) annotation(
      Line(points = {{10, 0}, {20, 0}, {20, -50}, {40, -50}}, color = {191, 0, 0}));
    connect(air_temperature.port, facade_nw.interior) annotation(
      Line(points = {{10, 0}, {20, 0}, {20, 50}, {-40, 50}}, color = {191, 0, 0}));
    connect(air_temperature.port, window_nw.interior) annotation(
      Line(points = {{10, 0}, {20, 0}, {20, 72}, {-40, 72}}, color = {191, 0, 0}));
    connect(air_temperature.port, window_sw.interior) annotation(
      Line(points = {{10, 0}, {20, 0}, {20, -28}, {-40, -28}}, color = {191, 0, 0}));
    connect(air_temperature.port, facade_sw.interior) annotation(
      Line(points = {{10, 0}, {20, 0}, {20, -50}, {-40, -50}}, color = {191, 0, 0}));
    connect(window_nw.exterior, side_nw) annotation(
      Line(points = {{-60, 72}, {-100, 72}, {-100, 40}, {-120, 40}}, color = {191, 0, 0}));
    connect(facade_nw.exterior, side_nw) annotation(
      Line(points = {{-60, 50}, {-100, 50}, {-100, 40}, {-120, 40}}, color = {191, 0, 0}));
    connect(window_sw.exterior, side_sw) annotation(
      Line(points = {{-60, -28}, {-100, -28}, {-100, -80}, {-120, -80}}, color = {191, 0, 0}));
    connect(facade_sw.exterior, side_sw) annotation(
      Line(points = {{-60, -50}, {-100, -50}, {-100, -80}, {-120, -80}}, color = {191, 0, 0}));
    connect(window_ne.exterior, side_ne) annotation(
      Line(points = {{62, 54}, {100, 54}, {100, 40}, {120, 40}}, color = {191, 0, 0}));
    connect(facade_ne.exterior, side_ne) annotation(
      Line(points = {{60, 30}, {100, 30}, {100, 40}, {120, 40}}, color = {191, 0, 0}));
    connect(window_se.exterior, side_se) annotation(
      Line(points = {{60, -26}, {100, -26}, {100, -80}, {120, -80}}, color = {191, 0, 0}));
    connect(facade_se.exterior, side_se) annotation(
      Line(points = {{60, -50}, {100, -50}, {100, -80}, {120, -80}}, color = {191, 0, 0}));
    connect(facade_nw.v_air, v_air_nw) annotation(
      Line(points = {{-60, 56}, {-80, 56}, {-80, 80}, {-140, 80}}, color = {0, 0, 127}));
    connect(window_nw.v_air, v_air_nw) annotation(
      Line(points = {{-60, 78}, {-80, 78}, {-80, 80}, {-140, 80}}, color = {0, 0, 127}));
    connect(window_sw.v_air, v_air_sw) annotation(
      Line(points = {{-60, -22}, {-80, -22}, {-80, -40}, {-140, -40}}, color = {0, 0, 127}));
    connect(facade_sw.v_air, v_air_sw) annotation(
      Line(points = {{-60, -44}, {-80, -44}, {-80, -40}, {-140, -40}}, color = {0, 0, 127}));
    connect(window_ne.v_air, v_air_ne) annotation(
      Line(points = {{64, 60}, {80, 60}, {80, 80}, {140, 80}}, color = {0, 0, 127}));
    connect(facade_ne.v_air, v_air_ne) annotation(
      Line(points = {{62, 36}, {80, 36}, {80, 80}, {140, 80}}, color = {0, 0, 127}));
    connect(window_se.v_air, v_air_se) annotation(
      Line(points = {{62, -20}, {80, -20}, {80, -40}, {140, -40}}, color = {0, 0, 127}));
    connect(facade_se.v_air, v_air_se) annotation(
      Line(points = {{62, -44}, {80, -44}, {80, -40}, {140, -40}}, color = {0, 0, 127}));
  end building2;

  model window_test
    Buildings.Fluid.MixingVolumes.MixingVolume vol(nPorts = 1, redeclare package Medium = Buildings.Media.Air "Moist air", m_flow_nominal = 0.001, V = 1) annotation(
      Placement(transformation(origin = {-70, 10}, extent = {{-10, -10}, {10, 10}})));
    Buildings.Fluid.MixingVolumes.MixingVolume vol1(nPorts = 2, redeclare package Medium = Buildings.Media.Air "Moist air", m_flow_nominal = 0.001, V = 1) annotation(
      Placement(transformation(origin = {108, 4}, extent = {{-10, -10}, {10, 10}})));
    Buildings.HeatTransfer.Sources.FixedTemperature preTem(T = 283.15) annotation(
      Placement(transformation(origin = {-110, 10}, extent = {{-10, -10}, {10, 10}})));
    Buildings.HeatTransfer.Sources.PrescribedTemperature preTem1 annotation(
      Placement(transformation(origin = {108, 50}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
    Modelica.Blocks.Sources.Ramp ramp(height = 20, duration = 864000/4, offset = 273.15 + 10, startTime = 864000/2) annotation(
      Placement(transformation(origin = {138, 76}, extent = {{-10, -10}, {10, 10}})));
    Buildings.Fluid.Sensors.Temperature senTem(redeclare package Medium = Buildings.Media.Air "Moist air") annotation(
      Placement(transformation(origin = {54, 18}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
    farm.atmosphere_control.window window1(redeclare package Medium = Buildings.Media.Air "Moist air", window_width = 1, m_flow_nominal = 0.001) annotation(
      Placement(transformation(origin = {0, 10}, extent = {{-10, -10}, {10, 10}})));
  equation
    connect(preTem.port, vol.heatPort) annotation(
      Line(points = {{-100, 10}, {-80, 10}}, color = {191, 0, 0}));
    connect(preTem1.port, vol1.heatPort) annotation(
      Line(points = {{98, 50}, {80, 50}, {80, 4}, {98, 4}}, color = {191, 0, 0}));
    connect(preTem1.T, ramp.y) annotation(
      Line(points = {{120, 50}, {166, 50}, {166, 76}, {150, 76}}, color = {0, 0, 127}));
    connect(vol1.ports[2], senTem.port) annotation(
      Line(points = {{108, -6}, {54, -6}, {54, 8}}, color = {0, 127, 255}));
    connect(senTem.T, window1.T) annotation(
      Line(points = {{48, 18}, {12, 18}}, color = {0, 0, 127}));
    connect(window1.port_b, vol1.ports[2]) annotation(
      Line(points = {{10, 10}, {20, 10}, {20, -6}, {108, -6}}, color = {0, 127, 255}));
    connect(window1.port_a, vol.ports[1]) annotation(
      Line(points = {{-10, 10}, {-20, 10}, {-20, 0}, {-70, 0}}, color = {0, 127, 255}));
    annotation(
      experiment(StartTime = 0, StopTime = 864000, Tolerance = 1e-06, Interval = 1728));
  end window_test;

  model aaaa
    Buildings.HeatTransfer.Convection.Exterior con(A = area, conMod = Buildings.HeatTransfer.Types.ExteriorConvection.TemperatureWind, azi = azimuth, til = 1.5707963267948966) annotation(
      Placement(transformation(origin = {-30, 30}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
    Buildings.HeatTransfer.Convection.Interior con1(A = area, conMod = Buildings.HeatTransfer.Types.InteriorConvection.Temperature, til = 1.5707963267948966) annotation(
      Placement(transformation(origin = {50, 30}, extent = {{-10, -10}, {10, 10}})));
    Buildings.HeatTransfer.Conduction.MultiLayer heaCon(A = area, layers(nLay = 2, material = {concrete, eps})) annotation(
      Placement(transformation(origin = {10, 30}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
    parameter Buildings.HeatTransfer.Data.Solids.Concrete concrete(x = 0.18, k = 2, c = 550, d = 3500) annotation(
      Placement(transformation(origin = {-10, -10}, extent = {{-10, -10}, {10, 10}})));
    parameter Buildings.HeatTransfer.Data.Solids.InsulationBoard eps(x = 0.05, k = 0.035, c = 1200, d = 21) annotation(
      Placement(transformation(origin = {30, -10}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Sources.Constant const(k = 3.14159 - azimuth) annotation(
      Placement(transformation(origin = {10, 70}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a exterior annotation(
      Placement(transformation(extent = {{-110, -10}, {-90, 10}}), iconTransformation(extent = {{-110, -10}, {-90, 10}})));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b interior annotation(
      Placement(transformation(extent = {{90, -10}, {110, 10}}), iconTransformation(extent = {{90, -10}, {110, 10}})));
    parameter Modelica.Units.SI.Area area "facade area";
    parameter Modelica.Units.SI.Angle azimuth "Surface azimuth";
    Modelica.Blocks.Interfaces.RealInput v_air(unit = "m/s") annotation(
      Placement(transformation(origin = {0, -40}, extent = {{-140, 80}, {-100, 120}}), iconTransformation(origin = {10, -40}, extent = {{-140, 80}, {-100, 120}})));
  equation
    connect(con.solid, heaCon.port_b) annotation(
      Line(points = {{-20, 30}, {0, 30}}, color = {191, 0, 0}));
    connect(heaCon.port_a, con1.solid) annotation(
      Line(points = {{20, 30}, {40, 30}}, color = {191, 0, 0}));
    connect(const.y, con.dir) annotation(
      Line(points = {{-1, 70}, {-6, 70}, {-6, 36}, {-18, 36}}, color = {0, 0, 127}));
    connect(v_air, con.v) annotation(
      Line(points = {{-120, 60}, {-12, 60}, {-12, 40}, {-18, 40}}, color = {0, 0, 127}));
    connect(exterior, con.fluid) annotation(
      Line(points = {{-100, 0}, {-60, 0}, {-60, 30}, {-40, 30}}, color = {191, 0, 0}));
    connect(interior, con1.fluid) annotation(
      Line(points = {{100, 0}, {80, 0}, {80, 30}, {60, 30}}, color = {191, 0, 0}));
  end aaaa;

  model bbbb
    Buildings.HeatTransfer.Convection.Exterior con(A = area, azi = azimuth, conMod = Buildings.HeatTransfer.Types.ExteriorConvection.TemperatureWind, til = 1.5707963267948966, roughness = Buildings.HeatTransfer.Types.SurfaceRoughness.VerySmooth) annotation(
      Placement(transformation(origin = {-30, 30}, extent = {{10, -10}, {-10, 10}})));
    Buildings.HeatTransfer.Convection.Interior con1(A = area, conMod = Buildings.HeatTransfer.Types.InteriorConvection.Temperature, til = 1.5707963267948966) annotation(
      Placement(transformation(origin = {50, 30}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Sources.Constant const(k = 3.14159 - azimuth) annotation(
      Placement(transformation(origin = {10, 70}, extent = {{10, -10}, {-10, 10}})));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a exterior annotation(
      Placement(transformation(extent = {{-110, -10}, {-90, 10}}), iconTransformation(extent = {{-110, -10}, {-90, 10}})));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b interior annotation(
      Placement(transformation(extent = {{90, -10}, {110, 10}}), iconTransformation(extent = {{90, -10}, {110, 10}})));
    Modelica.Blocks.Interfaces.RealInput v_air(unit = "m/s") annotation(
      Placement(transformation(origin = {0, -40}, extent = {{-140, 80}, {-100, 120}}), iconTransformation(origin = {10, -40}, extent = {{-140, 80}, {-100, 120}})));
    parameter Modelica.Units.SI.Area area "window area";
    parameter Modelica.Units.SI.Angle azimuth "Surface azimuth";
    Buildings.HeatTransfer.Conduction.SingleLayer lay(A = area, material = glass) annotation(
      Placement(transformation(origin = {10, 30}, extent = {{-10, -10}, {10, 10}})));
    parameter Buildings.HeatTransfer.Data.Solids.Glass glass(x(displayUnit = "mm") = 0.025) annotation(
      Placement(transformation(origin = {10, -22}, extent = {{-10, -10}, {10, 10}})));
  equation
    connect(con.fluid, exterior) annotation(
      Line(points = {{-40, 30}, {-80, 30}, {-80, 0}, {-100, 0}}, color = {191, 0, 0}));
    connect(con1.fluid, interior) annotation(
      Line(points = {{60, 30}, {80, 30}, {80, 0}, {100, 0}}, color = {191, 0, 0}));
    connect(con1.solid, lay.port_b) annotation(
      Line(points = {{40, 30}, {20, 30}}, color = {191, 0, 0}));
    connect(lay.port_a, con.solid) annotation(
      Line(points = {{0, 30}, {-20, 30}}, color = {191, 0, 0}));
    connect(con.dir, const.y) annotation(
      Line(points = {{-18, 36}, {-8, 36}, {-8, 70}, {0, 70}}, color = {0, 0, 127}));
    connect(con.v, v_air) annotation(
      Line(points = {{-18, 40}, {-14, 40}, {-14, 60}, {-120, 60}}, color = {0, 0, 127}));
  end bbbb;

  model cccc
    bbbb bbbb_sw(area = 581.04, azimuth = 0.7330382858376184) annotation(
      Placement(transformation(origin = {-60, -30}, extent = {{-10, -10}, {10, 10}})));
    Buildings.HeatTransfer.Sources.FixedTemperature preTem1(T = 293.15) annotation(
      Placement(transformation(extent = {{10, -10}, {-10, 10}})));
    aaaa aaaa_sw(area = 523.2, azimuth = 0.7330382858376184) annotation(
      Placement(transformation(origin = {-60, -50}, extent = {{-10, -10}, {10, 10}})));
    vertical.bbbb bbbb_sw1(area = 581.04, azimuth = -2.4085543677521746) annotation(
      Placement(transformation(origin = {50, 50}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
    vertical.aaaa aaaa_sw1(area = 523.2, azimuth = -2.4085543677521746) annotation(
      Placement(transformation(origin = {50, 30}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
    vertical.bbbb bbbb_sw11(area = 118.8, azimuth = -0.8377580409572781) annotation(
      Placement(transformation(origin = {70, -30}, extent = {{10, -10}, {-10, 10}})));
    vertical.aaaa aaaa_sw11(area = 369.71, azimuth = -0.8377580409572781) annotation(
      Placement(transformation(origin = {70, -50}, extent = {{10, -10}, {-10, 10}})));
    vertical.bbbb bbbb_sw12(area = 118.8, azimuth = 2.303834612632515) annotation(
      Placement(transformation(origin = {-70, 72}, extent = {{-10, -10}, {10, 10}}, rotation = -0)));
    vertical.aaaa aaaa_sw12(area = 369.71, azimuth = 2.303834612632515) annotation(
      Placement(transformation(origin = {-70, 52}, extent = {{-10, -10}, {10, 10}}, rotation = -0)));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a side_sw annotation(
      Placement(transformation(origin = {-20, -40}, extent = {{-110, -10}, {-90, 10}}), iconTransformation(origin = {0, -40}, extent = {{-110, -10}, {-90, 10}})));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a side_nw annotation(
      Placement(transformation(origin = {-20, 84}, extent = {{-110, -10}, {-90, 10}}), iconTransformation(origin = {0, 80}, extent = {{-110, -10}, {-90, 10}})));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a side_ne annotation(
      Placement(transformation(origin = {20, 80}, extent = {{110, -10}, {90, 10}}), iconTransformation(origin = {0, 80}, extent = {{110, -10}, {90, 10}})));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a side_se annotation(
      Placement(transformation(origin = {22, -20}, extent = {{110, -10}, {90, 10}}), iconTransformation(origin = {0, -40}, extent = {{110, -10}, {90, 10}})));
    Modelica.Blocks.Interfaces.RealInput v_air_sw(unit = "m/s") annotation(
      Placement(transformation(origin = {-20, -180}, extent = {{-140, 80}, {-100, 120}}), iconTransformation(origin = {0, -180}, extent = {{-140, 80}, {-100, 120}})));
    Modelica.Blocks.Interfaces.RealInput v_air_nw(unit = "m/s") annotation(
      Placement(transformation(origin = {-20, -60}, extent = {{-140, 80}, {-100, 120}}), iconTransformation(origin = {0, -60}, extent = {{-140, 80}, {-100, 120}})));
    Modelica.Blocks.Interfaces.RealInput v_air_ne(unit = "m/s") annotation(
      Placement(transformation(origin = {20, -60}, extent = {{140, 80}, {100, 120}}), iconTransformation(origin = {0, -60}, extent = {{140, 80}, {100, 120}})));
    Modelica.Blocks.Interfaces.RealInput v_air_se(unit = "m/s") annotation(
      Placement(transformation(origin = {20, -180}, extent = {{140, 80}, {100, 120}}), iconTransformation(origin = {0, -180}, extent = {{140, 80}, {100, 120}})));
    Modelica.Thermal.HeatTransfer.Sensors.HeatFlowSensor heat_flow_se annotation(
      Placement(transformation(origin = {40, -42}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Thermal.HeatTransfer.Sensors.HeatFlowSensor heat_flow_ne annotation(
      Placement(transformation(origin = {20, 40}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Thermal.HeatTransfer.Sensors.HeatFlowSensor heat_flow_sw annotation(
      Placement(transformation(origin = {-30, -40}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
    Modelica.Thermal.HeatTransfer.Sensors.HeatFlowSensor heat_flow_nw annotation(
      Placement(transformation(origin = {-40, 60}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
    Modelica.Blocks.Interfaces.RealOutput heat_se annotation(
      Placement(transformation(origin = {-40, -110}, extent = {{-10, -10}, {10, 10}}, rotation = -90), iconTransformation(origin = {-60, -110}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Blocks.Interfaces.RealOutput heat_sw annotation(
      Placement(transformation(origin = {-20, -110}, extent = {{-10, -10}, {10, 10}}, rotation = -90), iconTransformation(origin = {-20, -110}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Blocks.Interfaces.RealOutput heat_nw annotation(
      Placement(transformation(origin = {0, -110}, extent = {{-10, -10}, {10, 10}}, rotation = -90), iconTransformation(origin = {20, -110}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Blocks.Interfaces.RealOutput heat_ne annotation(
      Placement(transformation(origin = {20, -110}, extent = {{-10, -10}, {10, 10}}, rotation = -90), iconTransformation(origin = {60, -110}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  equation
    connect(bbbb_sw12.v_air, v_air_nw) annotation(
      Line(points = {{-80, 78}, {-110, 78}, {-110, 40}, {-140, 40}}, color = {0, 0, 127}));
    connect(aaaa_sw12.v_air, v_air_nw) annotation(
      Line(points = {{-80, 58}, {-110, 58}, {-110, 40}, {-140, 40}}, color = {0, 0, 127}));
    connect(bbbb_sw.v_air, v_air_sw) annotation(
      Line(points = {{-71, -24}, {-100, -24}, {-100, -80}, {-140, -80}}, color = {0, 0, 127}));
    connect(aaaa_sw.v_air, v_air_sw) annotation(
      Line(points = {{-71, -44}, {-100, -44}, {-100, -80}, {-140, -80}}, color = {0, 0, 127}));
    connect(bbbb_sw11.v_air, v_air_se) annotation(
      Line(points = {{82, -24}, {100, -24}, {100, -80}, {140, -80}}, color = {0, 0, 127}));
    connect(bbbb_sw1.v_air, v_air_ne) annotation(
      Line(points = {{62, 56}, {100, 56}, {100, 40}, {140, 40}}, color = {0, 0, 127}));
    connect(aaaa_sw1.v_air, v_air_ne) annotation(
      Line(points = {{62, 36}, {100, 36}, {100, 40}, {140, 40}}, color = {0, 0, 127}));
    connect(bbbb_sw12.exterior, side_nw) annotation(
      Line(points = {{-80, 72}, {-100, 72}, {-100, 84}, {-120, 84}}, color = {191, 0, 0}));
    connect(aaaa_sw12.exterior, side_nw) annotation(
      Line(points = {{-80, 52}, {-100, 52}, {-100, 84}, {-120, 84}}, color = {191, 0, 0}));
    connect(bbbb_sw.exterior, side_sw) annotation(
      Line(points = {{-70, -30}, {-96, -30}, {-96, -40}, {-120, -40}}, color = {191, 0, 0}));
    connect(aaaa_sw.exterior, side_sw) annotation(
      Line(points = {{-70, -50}, {-96, -50}, {-96, -40}, {-120, -40}}, color = {191, 0, 0}));
    connect(aaaa_sw1.exterior, side_ne) annotation(
      Line(points = {{60, 30}, {80, 30}, {80, 80}, {120, 80}}, color = {191, 0, 0}));
    connect(bbbb_sw11.exterior, side_se) annotation(
      Line(points = {{80, -30}, {90, -30}, {90, -20}, {122, -20}}, color = {191, 0, 0}));
    connect(aaaa_sw11.exterior, side_se) annotation(
      Line(points = {{80, -50}, {90, -50}, {90, -20}, {122, -20}}, color = {191, 0, 0}));
    connect(bbbb_sw1.exterior, side_ne) annotation(
      Line(points = {{60, 50}, {80, 50}, {80, 80}, {120, 80}}, color = {191, 0, 0}));
    connect(aaaa_sw11.v_air, v_air_se) annotation(
      Line(points = {{82, -44}, {100, -44}, {100, -80}, {140, -80}}, color = {0, 0, 127}));
    connect(preTem1.port, heat_flow_sw.port_a) annotation(
      Line(points = {{-10, 0}, {-20, 0}, {-20, -40}}, color = {191, 0, 0}));
    connect(preTem1.port, heat_flow_nw.port_a) annotation(
      Line(points = {{-10, 0}, {-20, 0}, {-20, 60}, {-30, 60}}, color = {191, 0, 0}));
    connect(preTem1.port, heat_flow_ne.port_a) annotation(
      Line(points = {{-10, 0}, {-20, 0}, {-20, 40}, {10, 40}}, color = {191, 0, 0}));
    connect(preTem1.port, heat_flow_se.port_a) annotation(
      Line(points = {{-10, 0}, {-20, 0}, {-20, -42}, {30, -42}}, color = {191, 0, 0}));
    connect(heat_flow_nw.port_b, bbbb_sw12.interior) annotation(
      Line(points = {{-50, 60}, {-54, 60}, {-54, 72}, {-60, 72}}, color = {191, 0, 0}));
    connect(heat_flow_nw.port_b, aaaa_sw12.interior) annotation(
      Line(points = {{-50, 60}, {-54, 60}, {-54, 52}, {-60, 52}}, color = {191, 0, 0}));
    connect(heat_flow_sw.port_b, bbbb_sw.interior) annotation(
      Line(points = {{-40, -40}, {-46, -40}, {-46, -30}, {-50, -30}}, color = {191, 0, 0}));
    connect(heat_flow_sw.port_b, aaaa_sw.interior) annotation(
      Line(points = {{-40, -40}, {-46, -40}, {-46, -50}, {-50, -50}}, color = {191, 0, 0}));
    connect(heat_flow_ne.port_b, bbbb_sw1.interior) annotation(
      Line(points = {{30, 40}, {36, 40}, {36, 50}, {40, 50}}, color = {191, 0, 0}));
    connect(heat_flow_ne.port_b, aaaa_sw1.interior) annotation(
      Line(points = {{30, 40}, {36, 40}, {36, 30}, {40, 30}}, color = {191, 0, 0}));
    connect(heat_flow_se.port_b, bbbb_sw11.interior) annotation(
      Line(points = {{50, -42}, {54, -42}, {54, -30}, {60, -30}}, color = {191, 0, 0}));
    connect(heat_flow_se.port_b, aaaa_sw11.interior) annotation(
      Line(points = {{50, -42}, {54, -42}, {54, -50}, {60, -50}}, color = {191, 0, 0}));
    connect(heat_flow_sw.Q_flow, heat_sw) annotation(
      Line(points = {{-30, -50}, {-30, -90}, {-20, -90}, {-20, -110}}, color = {0, 0, 127}));
    connect(heat_flow_se.Q_flow, heat_se) annotation(
      Line(points = {{40, -52}, {40, -60}, {-32, -60}, {-32, -90}, {-40, -90}, {-40, -110}}, color = {0, 0, 127}));
    connect(heat_flow_ne.Q_flow, heat_ne) annotation(
      Line(points = {{20, 30}, {20, -110}}, color = {0, 0, 127}));
    connect(heat_flow_nw.Q_flow, heat_nw) annotation(
      Line(points = {{-40, 50}, {-40, 26}, {18, 26}, {18, -90}, {0, -90}, {0, -110}}, color = {0, 0, 127}));
    annotation(
      Icon(graphics = {Bitmap(extent = {{100, -100}, {-100, 100}}, imageSource = "iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAAAXNSR0IArs4c6QAABKpJREFUeF7tnHFS2zwQxdfhIhyFXqQ4MJ3JLQqXaDPTAgkXaY+SHqSodUhwbGQ/rXZlyWTz1zefrJX0ftqntfC0IvtlVaDKOroNTgYg8yYwAAYgswKZh7cMMACZFcg8vGWAAcisQObhLQMMQGYFMg9vGWAA4hW4vl09OUd/nh/Xd/FR8vacZQbU9erSXdATEV0d5NtVf+nTZrPe5ZWTP/rsANQ3X2pXXTTi938752g7t2yYFYDGcoioHttnztH9nCDMAoDHclCuz8aSigdQ36xqV+39nvubhSUVDSDEchCV0i2pSAARloM4FGtJxQE4tRxH5P2T3f7/DzUOoyjSkooCoGE5KBVKs6QiAERZjttnx/KF6LKq6CsSvtdejCVlBxBZ5eyqipabn+vfjbAHgL+I6JIBoghLygog0nI224f1si90A+FlQTU3G3JbUhYAMZbTnLkLR8vN43oztss/36zuuBD+J1E2S5ocgIblIJuZkyVNCqC1nGbYZk+P/Q7PLGiz/fHeckIgxFnS4v758dtk19uTAIixnEbgKsByEIjSLSk5gBDL8bxTdaocJDJq51lSRY5cU+JOUiUlBdCtckJsZy+lt8pBIqP2UqukJABay6musNe30mlYDgLBs6S3TZOsSlIHEGI5rzc8nUNY1XIQBJ4lvUVLYkmqADRfrJCI0vZSLEkFAL/Kec2AKSwHgeJZUpsNWh8BiAF0LCf8inhSy0EQYizJEe1I4SMAMYDr29XezBl39EmqHCQyao+1pO3DWqShqHOzqCMAtECtF6uQcSTPcC1pLgCKshwEiGNJcwBQpOWEQPDeJfXOuSIBnM5ROkEkVOp2ZLHS9SU/A6QTTC0wim8AkEKJ2/cARspr6QazDAAALQMS73AU3gAghRK3+wG0l4lmQVkAtIOWBcBzWEknmFhfGP40A06Xd/xv6frsEI48hA0A3Ls6D9ghrKNjdBQDEC2dTsezA4AWPCZryN+DuIcmmg83Xn/+xR3C6NVfuq+5gp0nAKnKI/0NQHDZF/whFwuXAQgGwNI1+GEDkBVARduH76xzb4Z3QafWob/g4K0+8OCZZYAWAKns8ZdnVgX1tEc7WFsw7XjzeA8QlJHagmnHmwSA5KsI6YKl/fsCDcX7sLehUgGl/UMBHJ9DlohOL1ZJ5gtW2oIH53PYslzBtNc3iQWdDjL1grUF045nAJhfMxsAK0PHjxHtHSKNJ+0/20N46I8hdgaMb+CCqyD/dTQCmjwDyvw8ffjuHgk2lvLtWsO/REsOgHkm2XuAUDBtoFaGWhmatqrS3rHa8c4wA7rnk+RM8m0dbrzJAaBDiNuOFox2LHe8oeeLvg0N+UDKvzD8JYQMwCF+/ATfTRvNBwFP/h6AJjDW7tMJLbiTAZ4A8drHvZeg9RcNIMZzfRbE+GcUkF7lZoB/Z2FL4a44NAPidzpvRmg+KNpZZAASQdJeLoBEWxAteKoq6AgNzQfBnTADJHakdxeEBOG2ZwfAnbA931VAnAEmqEwBAyDTT9zbAIgllAUwADL9xL0NgFhCWQADINNP3NsAiCWUBTAAMv3EvQ2AWEJZAAMg00/c2wCIJZQF+AcXnnmOvw3zCwAAAABJRU5ErkJggg==")}));
  end cccc;

  model zzzz
    cccc cccc1 annotation(
      Placement(transformation(origin = {30, 10}, extent = {{-10, -10}, {10, 10}})));
    Buildings.BoundaryConditions.WeatherData.ReaderTMY3 weaDat(filNam = "/home/marci/dev/ma/irradiance_model/tmy/DEU_BY_Nurnberg.AP.107630_TMYx.mos") annotation(
      Placement(transformation(origin = {-226, 10}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Sources.Constant const(k = 0) annotation(
      Placement(transformation(origin = {70, 50}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
    Buildings.BoundaryConditions.WeatherData.Bus weaBus annotation(
      Placement(transformation(origin = {-196, 10}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {-86, -24}, extent = {{-10, -10}, {10, 10}})));
    envelope envelope1(area = 1236.9, azimuth = 0.7330382858376184) annotation(
      Placement(transformation(origin = {-70, 10}, extent = {{-10, -10}, {10, 10}})));
    Buildings.Fluid.MixingVolumes.MixingVolume vol(redeclare package Medium = Buildings.Media.Air "Moist air", m_flow_nominal = 0.001, V = 594, nPorts = 1) annotation(
      Placement(transformation(origin = {-30, -10}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    plant plant1(redeclare package Medium = Buildings.Media.Air "Moist air", a_plant = 523.2) annotation(
      Placement(transformation(origin = {-30, -50}, extent = {{-10, -10}, {10, 10}})));
    farm.illumination_control.shading shading1(light_setpoint = 200, azimuth = 0.7330382858376184) annotation(
      Placement(transformation(origin = {-150, -30}, extent = {{-10, -10}, {10, 10}})));
    farm.illumination_control.lighting lighting1 annotation(
      Placement(transformation(origin = {-110, -50}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Math.Add total_radiance annotation(
      Placement(transformation(origin = {-76, -38}, extent = {{-10, -10}, {10, 10}})));
    analysis analysis1 annotation(
      Placement(transformation(origin = {150, -50}, extent = {{-10, -10}, {10, 10}})));
    photovoltaik photovoltaik1 annotation(
      Placement(transformation(origin = {-150, -90}, extent = {{-10, -10}, {10, 10}})));
  equation
    connect(const.y, cccc1.v_air_nw) annotation(
      Line(points = {{59, 50}, {10, 50}, {10, 18}, {18, 18}}, color = {0, 0, 127}));
    connect(const.y, cccc1.v_air_sw) annotation(
      Line(points = {{59, 50}, {10, 50}, {10, 6}, {18, 6}}, color = {0, 0, 127}));
    connect(const.y, cccc1.v_air_ne) annotation(
      Line(points = {{60, 50}, {48, 50}, {48, 18}, {42, 18}}, color = {0, 0, 127}));
    connect(const.y, cccc1.v_air_se) annotation(
      Line(points = {{60, 50}, {48, 50}, {48, 6}, {42, 6}}, color = {0, 0, 127}));
    connect(weaDat.weaBus, weaBus) annotation(
      Line(points = {{-216, 10}, {-196, 10}}, color = {255, 204, 51}, thickness = 0.5));
    connect(envelope1.v_air, const.y) annotation(
      Line(points = {{-58, 16}, {-48, 16}, {-48, 50}, {60, 50}}, color = {0, 0, 127}));
    connect(weaBus, envelope1.weaBus) annotation(
      Line(points = {{-196, 10}, {-80, 10}}, thickness = 0.5));
    connect(envelope1.interior, vol.heatPort) annotation(
      Line(points = {{-60, 10}, {-30, 10}, {-30, 0}}, color = {191, 0, 0}));
    connect(vol.heatPort, cccc1.side_nw) annotation(
      Line(points = {{-30, 0}, {-30, 10}, {0, 10}, {0, 14}, {20, 14}}, color = {191, 0, 0}));
    connect(vol.heatPort, cccc1.side_sw) annotation(
      Line(points = {{-30, 0}, {-30, 10}, {0, 10}, {0, 2}, {20, 2}}, color = {191, 0, 0}));
    connect(vol.heatPort, cccc1.side_se) annotation(
      Line(points = {{-30, 0}, {-30, 10}, {0, 10}, {0, -4}, {60, -4}, {60, 2}, {40, 2}}, color = {191, 0, 0}));
    connect(vol.heatPort, cccc1.side_ne) annotation(
      Line(points = {{-30, 0}, {-30, 10}, {0, 10}, {0, -4}, {60, -4}, {60, 14}, {40, 14}}, color = {191, 0, 0}));
    connect(plant1.farm_air, vol.ports[1]) annotation(
      Line(points = {{-40, -50}, {-56, -50}, {-56, -10}, {-40, -10}}, color = {0, 127, 255}));
    connect(plant1.v_air, const.y) annotation(
      Line(points = {{-42, -58}, {-48, -58}, {-48, 50}, {60, 50}}, color = {0, 0, 127}));
    connect(weaBus, shading1.weaBus) annotation(
      Line(points = {{-196, 10}, {-196, -30}, {-160, -30}}, thickness = 0.5));
    connect(shading1.y, total_radiance.u1) annotation(
      Line(points = {{-138, -30}, {-100, -30}, {-100, -32}, {-88, -32}}, color = {0, 0, 127}));
    connect(total_radiance.y, plant1.par) annotation(
      Line(points = {{-64, -38}, {-60, -38}, {-60, -42}, {-42, -42}}, color = {0, 0, 127}));
    connect(lighting1.port, vol.heatPort) annotation(
      Line(points = {{-100, -54}, {-74, -54}, {-74, -76}, {0, -76}, {0, 10}, {-30, 10}, {-30, 0}}, color = {191, 0, 0}));
    connect(lighting1.ppfd_out, total_radiance.u2) annotation(
      Line(points = {{-98, -44}, {-88, -44}}, color = {0, 0, 127}));
    connect(lighting1.p, analysis1.led_power_draw) annotation(
      Line(points = {{-98, -50}, {-60, -50}, {-60, -70}, {120, -70}, {120, -42}, {138, -42}}, color = {0, 0, 127}));
    connect(shading1.y, lighting1.natural_radiation) annotation(
      Line(points = {{-138, -30}, {-130, -30}, {-130, -50}, {-122, -50}}, color = {0, 0, 127}));
    connect(weaBus, photovoltaik1.weaBus) annotation(
      Line(points = {{-196, 10}, {-196, -90}, {-160, -90}}, thickness = 0.5));
    connect(photovoltaik1.pv_power_output, analysis1.pv_power_production) annotation(
      Line(points = {{-138, -90}, {132, -90}, {132, -58}, {138, -58}}, color = {0, 0, 127}));
    annotation(
      experiment(StartTime = 0, StopTime = 31536000, Tolerance = 1e-06, Interval = 315.36),
      Diagram(coordinateSystem(extent = {{-240, 60}, {160, -100}})));
  end zzzz;

  model envelope
    parameter Modelica.Units.SI.Area area "envelope area";
    parameter Modelica.Units.SI.Angle azimuth "Surface azimuth";
    Buildings.HeatTransfer.Convection.Exterior conv_out(A = area, azi = azimuth, conMod = Buildings.HeatTransfer.Types.ExteriorConvection.TemperatureWind, roughness = Buildings.HeatTransfer.Types.SurfaceRoughness.VerySmooth, til = 1.5707963267948966) annotation(
      Placement(transformation(origin = {-30, 30}, extent = {{10, -10}, {-10, 10}})));
    Modelica.Blocks.Sources.Constant air_direction(k = 3.14159 - azimuth) annotation(
      Placement(transformation(origin = {46, 80}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b interior annotation(
      Placement(transformation(extent = {{90, -10}, {110, 10}}), iconTransformation(origin = {0, -60}, extent = {{90, -10}, {110, 10}})));
    Modelica.Blocks.Interfaces.RealInput v_air(unit = "m/s") annotation(
      Placement(transformation(origin = {0, -40}, extent = {{140, 80}, {100, 120}}), iconTransformation(origin = {0, -100}, extent = {{140, 80}, {100, 120}})));
    Buildings.HeatTransfer.Conduction.SingleLayer conduction(A = area, material = glass) annotation(
      Placement(transformation(origin = {10, 30}, extent = {{-10, -10}, {10, 10}})));
    Buildings.HeatTransfer.Convection.Exterior conv_in(A = area, azi = azimuth, conMod = Buildings.HeatTransfer.Types.ExteriorConvection.TemperatureWind, roughness = Buildings.HeatTransfer.Types.SurfaceRoughness.VerySmooth, til = 1.5707963267948966) annotation(
      Placement(transformation(origin = {50, 30}, extent = {{-10, -10}, {10, 10}}, rotation = -0)));
    parameter Buildings.HeatTransfer.Data.Solids.Glass glass(x(displayUnit = "mm") = 0.003) annotation(
      Placement(transformation(origin = {10, -10}, extent = {{-10, -10}, {10, 10}})));
    Buildings.HeatTransfer.Sources.PrescribedTemperature pre_tem annotation(
      Placement(transformation(origin = {-70, 30}, extent = {{-10, -10}, {10, 10}})));
    Buildings.BoundaryConditions.WeatherData.Bus weaBus annotation(
      Placement(transformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {-100, -2}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Math.Gain gaiWin(k = area*0.5) annotation(
      Placement(transformation(origin = {80, 110}, extent = {{10, -10}, {-10, 10}}, rotation = 90)));
    Buildings.HeatTransfer.Sources.PrescribedHeatFlow preHeaFlo annotation(
      Placement(transformation(origin = {80, 80}, extent = {{10, -10}, {-10, 10}}, rotation = 90)));
    Modelica.Blocks.Interfaces.RealInput radiation_shaded(unit = "m/s") annotation(
      Placement(transformation(origin = {0, 30}, extent = {{140, 80}, {100, 120}}), iconTransformation(origin = {0, -40}, extent = {{140, 80}, {100, 120}})));
  equation
    connect(air_direction.y, conv_in.dir) annotation(
      Line(points = {{35, 80}, {30, 80}, {30, 36}, {38, 36}}, color = {0, 0, 127}));
    connect(conv_in.v, v_air) annotation(
      Line(points = {{38, 40}, {32, 40}, {32, 60}, {120, 60}}, color = {0, 0, 127}));
    connect(conv_out.solid, conduction.port_a) annotation(
      Line(points = {{-20, 30}, {0, 30}}, color = {191, 0, 0}));
    connect(conduction.port_b, conv_in.solid) annotation(
      Line(points = {{20, 30}, {40, 30}}, color = {191, 0, 0}));
    connect(conv_in.fluid, interior) annotation(
      Line(points = {{60, 30}, {80, 30}, {80, 0}, {100, 0}}, color = {191, 0, 0}));
    connect(weaBus.TDryBul, pre_tem.T) annotation(
      Line(points = {{-100, 0}, {-90, 0}, {-90, 30}, {-82, 30}}, color = {0, 0, 127}));
    connect(pre_tem.port, conv_out.fluid) annotation(
      Line(points = {{-60, 30}, {-40, 30}}, color = {191, 0, 0}));
    connect(weaBus.winDir, conv_out.dir) annotation(
      Line(points = {{-100, 0}, {-100, 60}, {-6, 60}, {-6, 36}, {-18, 36}}, color = {0, 0, 127}));
    connect(weaBus.winSpe, conv_out.v) annotation(
      Line(points = {{-100, 0}, {-98, 0}, {-98, 58}, {-8, 58}, {-8, 40}, {-18, 40}}, color = {0, 0, 127}));
    connect(radiation_shaded, gaiWin.u) annotation(
      Line(points = {{120, 130}, {80.5, 130}, {80.5, 122}, {80, 122}}, color = {0, 0, 127}));
    connect(gaiWin.y, preHeaFlo.Q_flow) annotation(
      Line(points = {{80, 99}, {80, 89}}, color = {0, 0, 127}));
    connect(preHeaFlo.port, interior) annotation(
      Line(points = {{80, 70}, {80, 0}, {100, 0}}, color = {191, 0, 0}));
    annotation(
      Icon(graphics = {Bitmap(extent = {{-100, 100}, {100, -100}}, imageSource = "iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAAAXNSR0IArs4c6QAABEBJREFUeF7tXdtx2zAQPCiNOJUkbsR6/bCL2FWEH5FMupE4lViNxEyoSJFsEST3eAQIeTXjGX8AuLvdu8VDEOmEn6gIuKjWaVxIQOQkIAEkIDICkc2zAkhAZAQim2cFkIDICEQ2zwogAZERiGyeFZA6AYtFdvM6k4Vz8kVEvo4TT50n1ThD46PuRGT/537LQ1Hk9f/qz6AKuFtl987JN7V1tGPNwSCPUYOd7XdVJeXTNr/vbOlpoA5nvs5+XmT8AaDp4aSFp3e/53KT3/ZufdZQRcBilS0qJ49TAHoKPtR4VpU8aCoBJqDW/GomLxOTAk3ywX06yN45J8viR/6MDAwTEFz3kWiM2vqA7qo2TRXABMzX2ePfFcCiMdYuD00AslsRjeAuPBdoCHgRkZu3WNqBcj7uCACZpMD/Qc4d/Pf/rtzknxEjGgIOC/Jm0MtNDo+JODxfZ60bgtTsw2BdGwAI+XVb6/hJAMgACaAEXZcGgwVACbKWABIAIkACJqDBbfuDMMtQ/74Htc9VUOQKJAEkAEOAc8D7OeCdIKMaiMFvvxOF7a+yqu0oHo3/qiSozoWnsc+iSMB1bQQHVUDTchAtQVgCJrAM9vvspNx8hzCFGo9xGnhdBIigCUgCwAywXoWRAFMCKEGwBID4xz8NvVtnVVvZoBoYG4DY9odL0EfbiBmvwtIlwHMkmloFDifgXQ2nBoBWgnxH4mj8JABkYNrL0MnfpALRbmreEWP0CjAIMekhwhBglelW4/goizB+GAIuAh7nbqi+FAB/jEkKSAAQpB7J5HoGJCABbIyzuzXig60wBIQMLAGez10MQ0BioOzdbUiaMfIoOgGoAyiX1huh2Pa5EwYZ8CZA0DmgxWlWAPYLIVaAqgJ4N9QLW2oVyApQVYC/E5oAJIAEYAhwGWr8nSgG/wQu5xrHTwkCM8C6As0IOG7r0UkIjN/8Xk5s+2YEHM9byi22EYkNQGz7tgQIfjk1NgA6+9yIcSPmQ4BzACbBlCBQg6axCmr5JoMVEKwC+MCmpuJBE5ASlKQEHZzmj/RqIN4qAStg7N8JT+IsiJOw2T6Ec0DKc8DJ95MOohoIxs/DuPk6a3hwKwpjR/uOG1NWF6qsxjmPBk1AjQRdPrbei2eKF3gH+VyUm3yJpCNMgOrh3WOkGhJl37YD/Qzy8O7ux9cfMiiQjPTFNkQ75+R29MfX14EcX+DQO6ggZPSTjoFJ7hfbSpbFNi96Y3JoCEvQ0UDrY+xRL1Jtf2ITfmz9MWQ1AfUAp/mgX/Zd4LwPQNnXlDS9DxrdP3d9EAF7OapfafJp/yal+p0C+tdY+bRhLM1oIrBfQuxfYVVV8mv2KkXU11iZJuEHHWxwBXxQ3MzCJgFmUOoGIgE63Mx6kQAzKHUDkQAdbma9SIAZlLqBSIAON7NeJMAMSt1AJECHm1kvEmAGpW4gEqDDzazXH04uF45wgVW9AAAAAElFTkSuQmCC"), Text(origin = {0, -89}, extent = {{-100, 11}, {100, -11}}, textString = "envelope")}),
      Diagram);
  end envelope;

  model photovoltaik
    Buildings.Electrical.DC.Sources.PVSimple pVSimple(A = 1212, V_nominal = 48, eta = 0.2, fAct = 0.8) annotation(
      Placement(transformation(origin = {24, -12}, extent = {{-10, -10}, {10, 10}})));
    Buildings.BoundaryConditions.SolarIrradiation.DirectTiltedSurface HDirRoo1(azi = Buildings.Types.Azimuth.S, til(displayUnit = "deg") = 0.5235987755982988) annotation(
      Placement(transformation(origin = {-68, 18}, extent = {{20, 20}, {40, 40}})));
    Buildings.BoundaryConditions.SolarIrradiation.DiffuseIsotropic HDifTilIso1(til = 0.5235987755982988) annotation(
      Placement(transformation(origin = {-38, 28}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Math.Add add annotation(
      Placement(transformation(origin = {2, 38}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Electrical.Analog.Basic.Ground ground annotation(
      Placement(transformation(origin = {-8, -2}, extent = {{-92, -40}, {-72, -20}})));
    Buildings.Electrical.DC.Loads.Resistor res(R = 0.5, V_nominal = 12) annotation(
      Placement(transformation(origin = {-36, -14}, extent = {{-2, -10}, {18, 10}}, rotation = -90)));
    Buildings.Electrical.DC.Sources.ConstantVoltage sou(V = 48) annotation(
      Placement(transformation(origin = {4, -12}, extent = {{-82, -10}, {-62, 10}})));
    Buildings.Electrical.DC.Sensors.GeneralizedSensor sen annotation(
      Placement(transformation(origin = {-18, -52}, extent = {{0, 30}, {20, 50}})));
    Buildings.BoundaryConditions.WeatherData.Bus weaBus annotation(
      Placement(transformation(origin = {0, 38}, extent = {{-110, -10}, {-90, 10}}), iconTransformation(extent = {{-120, -18}, {-80, 22}})));
    Modelica.Blocks.Interfaces.RealOutput pv_power_output annotation(
      Placement(transformation(origin = {66, -4}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {110, 0}, extent = {{-10, -10}, {10, 10}})));
  equation
    connect(HDirRoo1.H, add.u1) annotation(
      Line(points = {{-27, 48}, {-21, 48}, {-21, 44}, {-11, 44}}, color = {0, 0, 127}));
    connect(HDifTilIso1.H, add.u2) annotation(
      Line(points = {{-27, 28}, {-21, 28}, {-21, 32}, {-11, 32}}, color = {0, 0, 127}));
    connect(sou.terminal, res.terminal) annotation(
      Line(points = {{-58, -12}, {-36, -12}}, color = {0, 0, 255}));
    connect(sen.terminal_p, pVSimple.terminal) annotation(
      Line(points = {{2, -12}, {14, -12}}, color = {0, 0, 255}));
    connect(sou.n, ground.p) annotation(
      Line(points = {{-78, -12}, {-90, -12}, {-90, -22}}, color = {0, 0, 255}));
    connect(sen.terminal_n, res.terminal) annotation(
      Line(points = {{-18, -12}, {-36, -12}}, color = {0, 0, 255}));
    connect(weaBus, HDirRoo1.weaBus) annotation(
      Line(points = {{-100, 38}, {-60, 38}, {-60, 48}, {-48, 48}}, thickness = 0.5));
    connect(weaBus, HDifTilIso1.weaBus) annotation(
      Line(points = {{-100, 38}, {-60, 38}, {-60, 28}, {-48, 28}}, thickness = 0.5));
  connect(add.y, pVSimple.G) annotation(
      Line(points = {{13, 38}, {24, 38}, {24, 0}}, color = {0, 0, 127}));
  connect(pVSimple.P, pv_power_output) annotation(
      Line(points = {{36, -4}, {66, -4}}, color = {0, 0, 127}));
    annotation(
      Icon(graphics = {Bitmap(origin = {0, -100}, extent = {{-100, 0}, {100, 200}}, imageSource = "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")}));
  end photovoltaik;

  model farmfarm
    parameter Modelica.Units.SI.Length window_width "width of the window array";
    parameter Modelica.Units.SI.Volume volume "volume of the farm";
    parameter Modelica.Units.SI.Area area "area of the farm envelope";
    parameter Modelica.Units.SI.Angle azimuth "azimuth of the farm";
    parameter Modelica.Units.SI.Area plant_area "area for growing plants";
    envelope envelope1(area = area, azimuth = azimuth) annotation(
      Placement(transformation(origin = {-42, -86}, extent = {{-10, -10}, {10, 10}})));
    farm.illumination_control.shading shading1(light_setpoint = 290, azimuth = azimuth) annotation(
      Placement(transformation(origin = {-40, 52}, extent = {{-10, -10}, {10, 10}})));
    farm.illumination_control.lighting lighting1(agricultural_area = plant_area) annotation(
      Placement(transformation(origin = {20, 34}, extent = {{-10, -10}, {10, 10}})));
    plant plant1(redeclare package Medium = Buildings.Media.Air "Moist air", a_plant = plant_area) annotation(
      Placement(transformation(origin = {180, -40}, extent = {{-10, -10}, {10, 10}})));
    Buildings.Fluid.MixingVolumes.MixingVolume vol(nPorts = 4, redeclare package Medium = Buildings.Media.Air "Moist air", m_flow_nominal = 90, V = volume) annotation(
      Placement(transformation(origin = {98, -12}, extent = {{10, -10}, {-10, 10}}, rotation = 90)));
    Buildings.Airflow.Multizone.MediumColumn air_column_bottom(redeclare package Medium = Buildings.Media.Air "Moist air", densitySelection = Buildings.Airflow.Multizone.Types.densitySelection.actual, h = 22.5/2) annotation(
      Placement(transformation(origin = {36, -216}, extent = {{-10, -10}, {10, 10}})));
    Buildings.Airflow.Multizone.MediumColumn air_column_top(redeclare package Medium = Buildings.Media.Air "Moist air", densitySelection = Buildings.Airflow.Multizone.Types.densitySelection.actual, h = 22.5/2) annotation(
      Placement(transformation(origin = {36, -156}, extent = {{-10, -10}, {10, 10}})));
    Buildings.Airflow.Multizone.DoorOperable window_bottom(LClo = 0.001, redeclare package Medium = Buildings.Media.Air "Moist air", hOpe = 0.4, wOpe = window_width) annotation(
      Placement(transformation(origin = {-24, -240}, extent = {{-10, -10}, {10, 10}})));
    Buildings.Airflow.Multizone.DoorOperable window_top(LClo = 0.001, redeclare package Medium = Buildings.Media.Air "Moist air", hOpe = 0.4, wOpe = window_width) annotation(
      Placement(transformation(origin = {-24, -136}, extent = {{-10, -10}, {10, 10}})));
    Buildings.Fluid.Sources.Outside outside_air(redeclare package Medium = Buildings.Media.Air "Moist air", nPorts = 4) annotation(
      Placement(transformation(origin = {-76, -176}, extent = {{-10, -10}, {10, 10}})));
    Buildings.BoundaryConditions.WeatherData.Bus weaBus annotation(
      Placement(transformation(origin = {-20, 0}, extent = {{-110, -10}, {-90, 10}}), iconTransformation(extent = {{-120, -18}, {-80, 22}})));
    Buildings.Controls.Continuous.LimPID conPID(controllerType = Modelica.Blocks.Types.SimpleController.P) annotation(
      Placement(transformation(origin = {-24, -176}, extent = {{10, -10}, {-10, 10}})));
    Modelica.Blocks.Sources.Constant temperature_setpoint(k = 273.15 + 24) annotation(
      Placement(transformation(origin = {-24, -206}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
    Buildings.Fluid.Sensors.Temperature air_temp(redeclare package Medium = Buildings.Media.Air "Moist air") annotation(
      Placement(transformation(origin = {4, -176}, extent = {{10, -10}, {-10, 10}})));
    Modelica.Blocks.Math.Add radiation_total annotation(
      Placement(transformation(origin = {98, 46}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Interfaces.RealOutput led_power_draw annotation(
      Placement(transformation(origin = {260, -210}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {110, -40}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Interfaces.RealOutput air_speed1 annotation(
      Placement(transformation(origin = {258, -106}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {110, 40}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b interior annotation(
      Placement(transformation(origin = {160, 10}, extent = {{90, -10}, {110, 10}}), iconTransformation(origin = {0, 80}, extent = {{90, -10}, {110, 10}})));
    Modelica.Blocks.Sources.SampleTrigger harvest_trigger(period = 3628800) annotation(
      Placement(transformation(origin = {158, -72}, extent = {{-10, -10}, {10, 10}}, rotation = -0)));
    Modelica.Blocks.Math.Gain air_speed(k = 1/(window_width*0.5)) annotation(
      Placement(transformation(origin = {36, -106}, extent = {{-10, -10}, {10, 10}})));
    Buildings.Fluid.Sensors.VolumeFlowRate senVolFlo(redeclare package Medium = Buildings.Media.Air "Moist air", m_flow_nominal = 0.001, tau = 10) annotation(
      Placement(transformation(origin = {10, -136}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Interfaces.RealOutput air_control_power_draw annotation(
      Placement(transformation(origin = {260, -270}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {110, -60}, extent = {{-10, -10}, {10, 10}})));
    window_motors window_motors1 annotation(
      Placement(transformation(origin = {-24, -270}, extent = {{-10, -10}, {10, 10}})));
    pump water_pump annotation(
      Placement(transformation(origin = {-50, -330}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Interfaces.RealOutput pump_power_draw annotation(
      Placement(transformation(origin = {260, -330}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {110, -80}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Interfaces.RealOutput yield annotation(
      Placement(transformation(origin = {260, -150}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {110, 0}, extent = {{-10, -10}, {10, 10}})));
  equation
    connect(temperature_setpoint.y, conPID.u_m) annotation(
      Line(points = {{-24, -195}, {-24, -189}}, color = {0, 0, 127}));
    connect(air_temp.T, conPID.u_s) annotation(
      Line(points = {{-3, -176}, {-13, -176}}, color = {0, 0, 127}));
    connect(air_temp.port, vol.ports[1]) annotation(
      Line(points = {{4, -186}, {108, -186}, {108, -12}}, color = {0, 127, 255}));
    connect(vol.ports[3], air_column_top.port_b) annotation(
      Line(points = {{108, -12}, {108, -166}, {36, -166}}, color = {0, 127, 255}));
    connect(vol.ports[4], air_column_bottom.port_a) annotation(
      Line(points = {{108, -12}, {108, -206}, {36, -206}}, color = {0, 127, 255}));
    connect(window_bottom.port_b1, air_column_bottom.port_b) annotation(
      Line(points = {{-14, -234}, {36, -234}, {36, -226}}, color = {0, 127, 255}));
    connect(window_bottom.port_a2, air_column_bottom.port_b) annotation(
      Line(points = {{-14, -246}, {36, -246}, {36, -226}}, color = {0, 127, 255}));
    connect(window_bottom.port_a1, outside_air.ports[1]) annotation(
      Line(points = {{-34, -234}, {-66, -234}, {-66, -176}}, color = {0, 127, 255}));
    connect(window_bottom.port_b2, outside_air.ports[2]) annotation(
      Line(points = {{-34, -246}, {-66, -246}, {-66, -176}}, color = {0, 127, 255}));
    connect(window_top.port_b2, outside_air.ports[3]) annotation(
      Line(points = {{-34, -142}, {-66, -142}, {-66, -176}}, color = {0, 127, 255}));
    connect(window_top.port_a1, outside_air.ports[4]) annotation(
      Line(points = {{-34, -130}, {-66, -130}, {-66, -176}}, color = {0, 127, 255}));
    connect(weaBus, outside_air.weaBus) annotation(
      Line(points = {{-120, 0}, {-100, 0}, {-100, -176}, {-86, -176}}, thickness = 0.5));
    connect(weaBus, envelope1.weaBus) annotation(
      Line(points = {{-120, 0}, {-99, 0}, {-99, -86}, {-52, -86}}, thickness = 0.5));
    connect(weaBus, shading1.weaBus) annotation(
      Line(points = {{-120, 0}, {-80, 0}, {-80, 52}, {-50, 52}}, thickness = 0.5));
    connect(shading1.y, lighting1.shaded_radiation) annotation(
      Line(points = {{-29, 52}, {-11, 52}, {-11, 34}, {8, 34}}, color = {0, 0, 127}));
    connect(shading1.y, envelope1.radiation_shaded) annotation(
      Line(points = {{-29, 52}, {-21, 52}, {-21, 0}, {-20.5, 0}, {-20.5, -80}, {-30, -80}}, color = {0, 0, 127}));
    connect(lighting1.ppfd_out, radiation_total.u2) annotation(
      Line(points = {{31, 40}, {86, 40}}, color = {0, 0, 127}));
    connect(shading1.y, radiation_total.u1) annotation(
      Line(points = {{-29, 52}, {86, 52}}, color = {0, 0, 127}));
    connect(radiation_total.y, plant1.par) annotation(
      Line(points = {{109, 46}, {152, 46}, {152, -32}, {168, -32}}, color = {0, 0, 127}));
    connect(lighting1.port, vol.heatPort) annotation(
      Line(points = {{30, 28}, {30, 10}, {98, 10}, {98, -2}}, color = {191, 0, 0}));
    connect(envelope1.interior, vol.heatPort) annotation(
      Line(points = {{-32, -92}, {29, -92}, {29, 10}, {97.5, 10}, {97.5, -2}, {98, -2}}, color = {191, 0, 0}));
    connect(lighting1.p, led_power_draw) annotation(
      Line(points = {{31, 34}, {60.5, 34}, {60.5, -56}, {80.25, -56}, {80.25, -210}, {260, -210}}, color = {0, 0, 127}));
    connect(vol.heatPort, interior) annotation(
      Line(points = {{98, -2}, {98, 10}, {260, 10}}, color = {191, 0, 0}));
    connect(harvest_trigger.y, plant1.reset) annotation(
      Line(points = {{169, -72}, {180, -72}, {180, -52}, {185, -52}}, color = {255, 0, 255}));
    connect(air_speed.y, air_speed1) annotation(
      Line(points = {{47, -106}, {258, -106}}, color = {0, 0, 127}));
    connect(air_speed.y, plant1.v_air) annotation(
      Line(points = {{47, -106}, {130, -106}, {130, -48}, {168, -48}}, color = {0, 0, 127}));
    connect(air_speed.y, envelope1.v_air) annotation(
      Line(points = {{47, -106}, {54, -106}, {54, -86}, {-30, -86}}, color = {0, 0, 127}));
    connect(senVolFlo.V_flow, air_speed.u) annotation(
      Line(points = {{10, -125}, {10, -106}, {24, -106}}, color = {0, 0, 127}));
    connect(air_column_top.port_a, senVolFlo.port_b) annotation(
      Line(points = {{36, -146}, {36, -136}, {20, -136}}, color = {0, 127, 255}));
    connect(senVolFlo.port_a, window_top.port_b1) annotation(
      Line(points = {{0, -136}, {-4, -136}, {-4, -130}, {-14, -130}}, color = {0, 127, 255}));
    connect(senVolFlo.port_a, window_top.port_a2) annotation(
      Line(points = {{0, -136}, {-4, -136}, {-4, -142}, {-14, -142}}, color = {0, 127, 255}));
    connect(vol.ports[4], plant1.farm_air) annotation(
      Line(points = {{108, -12}, {108, -40}, {170, -40}}, color = {0, 127, 255}));
    connect(conPID.y, window_bottom.y) annotation(
      Line(points = {{-35, -176}, {-51, -176}, {-51, -240}, {-35, -240}}, color = {0, 0, 127}));
    connect(conPID.y, window_top.y) annotation(
      Line(points = {{-35, -176}, {-51, -176}, {-51, -136}, {-35, -136}}, color = {0, 0, 127}));
    connect(conPID.y, window_motors1.window_control) annotation(
      Line(points = {{-35, -176}, {-51, -176}, {-51, -270}, {-37, -270}}, color = {0, 0, 127}));
    connect(window_motors1.p, air_control_power_draw) annotation(
      Line(points = {{-13, -270}, {260, -270}}, color = {0, 0, 127}));
  connect(weaBus, water_pump.weaBus) annotation(
      Line(points = {{-120, 0}, {-100, 0}, {-100, -330}, {-60, -330}}, thickness = 0.5));
  connect(water_pump.y, pump_power_draw) annotation(
      Line(points = {{-38, -330}, {260, -330}}, color = {0, 0, 127}));
  connect(plant1.y, yield) annotation(
      Line(points = {{192, -40}, {200, -40}, {200, -150}, {260, -150}}, color = {0, 0, 127}));
    annotation(
      Diagram(coordinateSystem(extent = {{-140, 120}, {280, -340}}), graphics = {Rectangle(origin = {4, -178}, extent = {{-64, 106}, {64, -106}}), Bitmap(origin = {-50, -56}, extent = {{-10, 10}, {10, -10}}, imageSource = "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"), Rectangle(origin = {98, 22}, extent = {{-26, 48}, {26, -48}}), Bitmap(origin = {83, 81}, extent = {{-11, 11}, {11, -11}}, imageSource = "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"), Rectangle(origin = {-9, 44}, extent = {{-51, 28}, {51, -28}}), Bitmap(origin = {-47, 83}, extent = {{11, -11}, {-11, 11}}, imageSource = "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"), Text(origin = {35, -170}, extent = {{-13, -2}, {13, 2}}, textString = "air_column_top"), Text(origin = {35, -202}, extent = {{-13, -2}, {13, 2}}, textString = "air_column_bot")}),
      Icon(graphics = {Bitmap(extent = {{100, -100}, {-100, 100}}, imageSource = "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")}));
  end farmfarm;

  model system
    cccc cccc1 annotation(
      Placement(transformation(origin = {18, -12}, extent = {{-10, -10}, {10, 10}})));
    Buildings.BoundaryConditions.WeatherData.ReaderTMY3 weaDat(filNam = "/home/marci/dev/ma/irradiance_model/tmy/DEU_BY_Nurnberg.AP.107630_TMYx.mos") annotation(
      Placement(transformation(origin = {-150, 10}, extent = {{-10, -10}, {10, 10}})));
    Buildings.BoundaryConditions.WeatherData.Bus weaBus annotation(
      Placement(transformation(origin = {-120, 10}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {-86, -24}, extent = {{-10, -10}, {10, 10}})));
    analysis analysis1 annotation(
      Placement(transformation(origin = {70, -88}, extent = {{-10, -10}, {10, 10}})));
    photovoltaik photovoltaik1 annotation(
      Placement(transformation(origin = {10, 70}, extent = {{-10, -10}, {10, 10}})));
    farmfarm farm_sw(window_width = 52.8, volume = 594, area = 1236.9, azimuth = 0.7330382858376184, plant_area = 523.2) annotation(
      Placement(transformation(origin = {-70, -30}, extent = {{-10, -10}, {10, 10}})));
    vertical.farmfarm farm_ne(area = 1236.9, azimuth = -2.4085543677521746, plant_area = 523.2, volume = 594, window_width = 52.8) annotation(
      Placement(transformation(origin = {-70, 30}, extent = {{-10, -10}, {10, 10}})));
    vertical.farmfarm farm_nw(area = 546.9, azimuth = 2.303834612632515, plant_area = 369.72, volume = 256.5, window_width = 22.8) annotation(
      Placement(transformation(origin = {-70, 0}, extent = {{-10, -10}, {10, 10}})));
    vertical.farmfarm farm_se(area = 546.9, azimuth = -0.8377580409572781, plant_area = 369.72, volume = 256.5, window_width = 22.5) annotation(
      Placement(transformation(origin = {-70, -60}, extent = {{-10, -10}, {10, 10}})));
    analysis_insulation analysis_insulation1(temperature_building = 293.15, area_ne = 523.2 + 581.04, area_nw = 369.71 + 118.8, area_sw = 523.2 + 581.04, area_se = 369.71 + 118.8) annotation(
      Placement(transformation(origin = {70, -22}, extent = {{-10, -10}, {10, 10}})));
  analysis_yield analysis_yield1 annotation(
      Placement(transformation(origin = {70, -56}, extent = {{-10, -10}, {10, 10}})));
  equation
    connect(weaDat.weaBus, weaBus) annotation(
      Line(points = {{-140, 10}, {-120, 10}}, color = {255, 204, 51}, thickness = 0.5));
    connect(weaBus, photovoltaik1.weaBus) annotation(
      Line(points = {{-120, 10}, {-120, 70}, {0, 70}}, thickness = 0.5));
    connect(cccc1.side_sw, farm_sw.interior) annotation(
      Line(points = {{8, -16}, {-2, -16}, {-2, -22}, {-60, -22}}, color = {191, 0, 0}));
    connect(farm_sw.air_speed1, cccc1.v_air_sw) annotation(
      Line(points = {{-59, -26}, {2, -26}, {2, -20}, {6, -20}}, color = {0, 0, 127}));
    connect(farm_se.interior, cccc1.side_se) annotation(
      Line(points = {{-60, -52}, {40, -52}, {40, -16}, {28, -16}}, color = {191, 0, 0}));
    connect(farm_se.air_speed1, cccc1.v_air_se) annotation(
      Line(points = {{-59, -56}, {36, -56}, {36, -20}, {30, -20}}, color = {0, 0, 127}));
    connect(farm_nw.interior, cccc1.side_nw) annotation(
      Line(points = {{-60, 8}, {2, 8}, {2, -4}, {8, -4}}, color = {191, 0, 0}));
    connect(farm_nw.air_speed1, cccc1.v_air_nw) annotation(
      Line(points = {{-59, 4}, {-2, 4}, {-2, -8}, {6, -8}}, color = {0, 0, 127}));
    connect(farm_ne.interior, cccc1.side_ne) annotation(
      Line(points = {{-60, 38}, {40, 38}, {40, -4}, {28, -4}}, color = {191, 0, 0}));
    connect(farm_ne.air_speed1, cccc1.v_air_ne) annotation(
      Line(points = {{-59, 34}, {36, 34}, {36, -8}, {30, -8}}, color = {0, 0, 127}));
    connect(weaBus, farm_ne.weaBus) annotation(
      Line(points = {{-120, 10}, {-100, 10}, {-100, 30}, {-80, 30}}, thickness = 0.5));
    connect(weaBus, farm_nw.weaBus) annotation(
      Line(points = {{-120, 10}, {-100, 10}, {-100, 0}, {-80, 0}}, thickness = 0.5));
    connect(weaBus, farm_sw.weaBus) annotation(
      Line(points = {{-120, 10}, {-100, 10}, {-100, -30}, {-80, -30}}, thickness = 0.5));
    connect(weaBus, farm_se.weaBus) annotation(
      Line(points = {{-120, 10}, {-100, 10}, {-100, -60}, {-80, -60}}, thickness = 0.5));
    connect(photovoltaik1.pv_power_output, analysis1.pv_power_production) annotation(
      Line(points = {{21, 70}, {52, 70}, {52, -78}, {59, -78}}, color = {0, 0, 127}));
    connect(farm_se.led_power_draw, analysis1.led_power_draw_se) annotation(
      Line(points = {{-59, -64}, {-26, -64}, {-26, -88}, {59, -88}}, color = {0, 0, 127}));
    connect(farm_sw.led_power_draw, analysis1.led_power_draw_sw) annotation(
      Line(points = {{-59, -34}, {-24, -34}, {-24, -86}, {59, -86}}, color = {0, 0, 127}));
    connect(farm_nw.led_power_draw, analysis1.led_power_draw_nw) annotation(
      Line(points = {{-59, -4}, {-22, -4}, {-22, -84}, {59, -84}}, color = {0, 0, 127}));
    connect(farm_ne.led_power_draw, analysis1.led_power_draw_ne) annotation(
      Line(points = {{-59, 26}, {-20, 26}, {-20, -82}, {59, -82}}, color = {0, 0, 127}));
    connect(cccc1.heat_ne, analysis_insulation1.heat_flow_ne) annotation(
      Line(points = {{24, -23}, {24, -26}, {59, -26}}, color = {0, 0, 127}));
    connect(cccc1.heat_nw, analysis_insulation1.heat_flow_nw) annotation(
      Line(points = {{20, -23}, {20, -28}, {59, -28}}, color = {0, 0, 127}));
    connect(cccc1.heat_sw, analysis_insulation1.heat_flow_sw) annotation(
      Line(points = {{16, -23}, {16, -30}, {59, -30}}, color = {0, 0, 127}));
    connect(cccc1.heat_se, analysis_insulation1.heat_flow_se) annotation(
      Line(points = {{12, -23}, {12, -32}, {59, -32}}, color = {0, 0, 127}));
    connect(weaBus, analysis_insulation1.weaBus) annotation(
      Line(points = {{-120, 10}, {-120, 70}, {-4, 70}, {-4, 82}, {60, 82}, {60, -22}}, thickness = 0.5));
  connect(farm_se.pump_power_draw, analysis1.pump_power_draw_se) annotation(
      Line(points = {{-59, -68}, {-56, -68}, {-56, -111}, {68, -111}, {68, -99}}, color = {0, 0, 127}));
  connect(farm_sw.pump_power_draw, analysis1.pump_power_draw_sw) annotation(
      Line(points = {{-59, -38}, {-54, -38}, {-54, -109}, {66, -109}, {66, -99}}, color = {0, 0, 127}));
  connect(farm_nw.pump_power_draw, analysis1.pump_power_draw_nw) annotation(
      Line(points = {{-59, -8}, {-52, -8}, {-52, -107}, {64, -107}, {64, -99}}, color = {0, 0, 127}));
  connect(farm_ne.pump_power_draw, analysis1.pump_power_draw_ne) annotation(
      Line(points = {{-59, 22}, {-50, 22}, {-50, -105}, {62, -105}, {62, -99}}, color = {0, 0, 127}));
  connect(farm_se.air_control_power_draw, analysis1.air_power_draw_se) annotation(
      Line(points = {{-59, -66}, {-42, -66}, {-42, -98}, {59, -98}}, color = {0, 0, 127}));
  connect(farm_sw.air_control_power_draw, analysis1.air_power_draw_sw) annotation(
      Line(points = {{-59, -36}, {-40, -36}, {-40, -96}, {59, -96}}, color = {0, 0, 127}));
  connect(farm_nw.air_control_power_draw, analysis1.air_power_draw_nw) annotation(
      Line(points = {{-59, -6}, {-38, -6}, {-38, -94}, {59, -94}}, color = {0, 0, 127}));
  connect(farm_ne.air_control_power_draw, analysis1.air_power_draw_ne) annotation(
      Line(points = {{-59, 24}, {-36, 24}, {-36, -92}, {59, -92}}, color = {0, 0, 127}));
  connect(farm_se.yield, analysis_yield1.yield_se) annotation(
      Line(points = {{-58, -60}, {-14, -60}, {-14, -66}, {60, -66}}, color = {0, 0, 127}));
  connect(farm_sw.yield, analysis_yield1.yield_sw) annotation(
      Line(points = {{-58, -30}, {-12, -30}, {-12, -64}, {60, -64}}, color = {0, 0, 127}));
  connect(farm_nw.yield, analysis_yield1.yield_nw) annotation(
      Line(points = {{-58, 0}, {-10, 0}, {-10, -62}, {60, -62}}, color = {0, 0, 127}));
  connect(farm_ne.yield, analysis_yield1.yield_ne) annotation(
      Line(points = {{-58, 30}, {-8, 30}, {-8, -60}, {60, -60}}, color = {0, 0, 127}));
    annotation(
      experiment(StartTime = 0, StopTime = 31536000, Tolerance = 1e-06, Interval = 315.36),
      Diagram(coordinateSystem(extent = {{-240, 60}, {160, -100}})));
  end system;

  model window_motors
    Modelica.Blocks.Continuous.Derivative derivative annotation(
      Placement(transformation(origin = {-50, 0}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Math.Abs absolute annotation(
      Placement(transformation(origin = {-20, 0}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Logical.OnOffController onOffController(bandwidth = 0.00001/4) annotation(
      Placement(transformation(origin = {16, 0}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Sources.Constant threshhold(k = 0.00001) annotation(
      Placement(transformation(origin = {-20, -32}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Math.BooleanToReal booleanToReal(realTrue = 26*2) annotation(
      Placement(transformation(origin = {48, 0}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Interfaces.RealInput window_control annotation(
      Placement(transformation(origin = {-100, 0}, extent = {{-20, -20}, {20, 20}}), iconTransformation(origin = {-120, 0}, extent = {{-20, -20}, {20, 20}})));
    Modelica.Blocks.Interfaces.RealOutput p annotation(
      Placement(transformation(origin = {90, 0}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {110, 0}, extent = {{-10, -10}, {10, 10}})));
  equation
    connect(derivative.y, absolute.u) annotation(
      Line(points = {{-39, 0}, {-32, 0}}, color = {0, 0, 127}));
    connect(threshhold.y, onOffController.u) annotation(
      Line(points = {{-9, -32}, {-2, -32}, {-2, -7}, {4, -7}, {4, -6}}, color = {0, 0, 127}));
    connect(absolute.y, onOffController.reference) annotation(
      Line(points = {{-9, 0}, {-4, 0}, {-4, 6}, {4, 6}}, color = {0, 0, 127}));
    connect(onOffController.y, booleanToReal.u) annotation(
      Line(points = {{27, 0}, {36, 0}}, color = {255, 0, 255}));
    connect(window_control, derivative.u) annotation(
      Line(points = {{-100, 0}, {-62, 0}}, color = {0, 0, 127}));
    connect(booleanToReal.y, p) annotation(
      Line(points = {{60, 0}, {90, 0}}, color = {0, 0, 127}));
    annotation(
      Icon(graphics = {Bitmap(origin = {0, 100}, extent = {{-100, -200}, {100, 0}}, imageSource = "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")}));
  end window_motors;

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
    annotation(
      Icon(graphics = {Bitmap(origin = {0, 100}, extent = {{-100, -200}, {100, 0}}, imageSource = "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")}));
  end pump;

  model analysis
    Modelica.Blocks.Interfaces.RealInput led_power_draw_ne annotation(
      Placement(transformation(origin = {-370, -48}, extent = {{-20, -20}, {20, 20}}), iconTransformation(origin = {-110, 60}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Continuous.Integrator led_daily_power_draw_ne(k = 1/(60*60), use_reset = true) annotation(
      Placement(transformation(origin = {-220, -78}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Sources.SampleTrigger daily_trigger(period = 86400) annotation(
      Placement(transformation(origin = {90, 10}, extent = {{10, -10}, {-10, 10}})));
    Modelica.Blocks.Continuous.Integrator led_yearly_power_draw_ne(k = 1/(60*60), use_reset = false) annotation(
      Placement(transformation(origin = {-220, -48}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Interfaces.RealInput pv_power_production annotation(
      Placement(transformation(origin = {-120, 80}, extent = {{-20, -20}, {20, 20}}), iconTransformation(origin = {-110, 100}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Continuous.Integrator pv_daily_power_production(k = 1/(60*60), use_reset = true) annotation(
      Placement(transformation(origin = {-10, 110}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Continuous.Integrator pv_total_power_production(k = 1/(60*60), use_reset = false) annotation(
      Placement(transformation(origin = {-10, 70}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Continuous.Integrator led_daily_power_draw_nw(k = 1/(60*60), use_reset = true) annotation(
      Placement(transformation(origin = {-220, -158}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Continuous.Integrator led_yearly_power_draw_nw(k = 1/(60*60), use_reset = false) annotation(
      Placement(transformation(origin = {-220, -128}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Continuous.Integrator led_daily_power_draw_sw(k = 1/(60*60), use_reset = true) annotation(
      Placement(transformation(origin = {-220, -238}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Continuous.Integrator led_yearly_power_draw_sw(k = 1/(60*60), use_reset = false) annotation(
      Placement(transformation(origin = {-220, -208}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Continuous.Integrator led_daily_power_draw_se(k = 1/(60*60), use_reset = true) annotation(
      Placement(transformation(origin = {-220, -318}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Continuous.Integrator led_yearly_power_draw_se(k = 1/(60*60), use_reset = false) annotation(
      Placement(transformation(origin = {-220, -288}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Sources.Constant zero(k = 0) annotation(
      Placement(transformation(origin = {-310, -338}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
    Modelica.Blocks.Math.Sum sum_ne(nin = 2) annotation(
      Placement(transformation(origin = {-280, -48}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Math.Sum sum_nw(nin = 2) annotation(
      Placement(transformation(origin = {-280, -128}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Math.Sum sum_sw(nin = 2) annotation(
      Placement(transformation(origin = {-280, -208}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Math.Sum sum_se(nin = 2) annotation(
      Placement(transformation(origin = {-280, -288}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Interfaces.RealInput led_power_draw_nw annotation(
      Placement(transformation(origin = {-370, -128}, extent = {{-20, -20}, {20, 20}}), iconTransformation(origin = {-110, 40}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Interfaces.RealInput led_power_draw_sw annotation(
      Placement(transformation(origin = {-370, -208}, extent = {{-20, -20}, {20, 20}}), iconTransformation(origin = {-110, 20}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Interfaces.RealInput led_power_draw_se annotation(
      Placement(transformation(origin = {-370, -288}, extent = {{-20, -20}, {20, 20}}), iconTransformation(origin = {-110, 0}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Continuous.Integrator led_daily_power_draw_total(k = 1/(60*60), use_reset = true) annotation(
      Placement(transformation(origin = {-158, -418}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Continuous.Integrator led_yearly_power_draw_total(k = 1/(60*60), use_reset = false) annotation(
      Placement(transformation(origin = {-158, -388}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Math.Sum sum_total(nin = 4) annotation(
      Placement(transformation(origin = {-220, -388}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Interfaces.RealInput air_power_draw_ne annotation(
      Placement(transformation(origin = {-60, -50}, extent = {{-20, -20}, {20, 20}}), iconTransformation(origin = {-110, -40}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Interfaces.RealInput air_power_draw_nw annotation(
      Placement(transformation(origin = {-60, -126}, extent = {{-20, -20}, {20, 20}}), iconTransformation(origin = {-110, -60}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Interfaces.RealInput air_power_draw_sw annotation(
      Placement(transformation(origin = {-60, -210}, extent = {{-20, -20}, {20, 20}}), iconTransformation(origin = {-110, -80}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Interfaces.RealInput air_power_draw_se annotation(
      Placement(transformation(origin = {-60, -290}, extent = {{-20, -20}, {20, 20}}), iconTransformation(origin = {-110, -100}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Interfaces.RealInput pump_power_draw_ne annotation(
      Placement(transformation(origin = {240, -52}, extent = {{-20, -20}, {20, 20}}), iconTransformation(origin = {-80, -110}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  Modelica.Blocks.Interfaces.RealInput pump_power_draw_nw annotation(
      Placement(transformation(origin = {240, -130}, extent = {{-20, -20}, {20, 20}}), iconTransformation(origin = {-60, -110}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  Modelica.Blocks.Interfaces.RealInput pump_power_draw_sw annotation(
      Placement(transformation(origin = {238, -210}, extent = {{-20, -20}, {20, 20}}), iconTransformation(origin = {-40, -110}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  Modelica.Blocks.Interfaces.RealInput pump_power_draw_se annotation(
      Placement(transformation(origin = {240, -290}, extent = {{-20, -20}, {20, 20}}), iconTransformation(origin = {-20, -110}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  Modelica.Blocks.Continuous.Integrator air_daily_power_draw_ne(k = 1/(60*60), use_reset = true) annotation(
      Placement(transformation(origin = {90, -80}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Continuous.Integrator air_yearly_power_draw_ne(k = 1/(60*60), use_reset = false) annotation(
      Placement(transformation(origin = {90, -50}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Continuous.Integrator air_daily_power_draw_nw(k = 1/(60*60), use_reset = true) annotation(
      Placement(transformation(origin = {90, -160}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Continuous.Integrator air_yearly_power_draw_nw(k = 1/(60*60), use_reset = false) annotation(
      Placement(transformation(origin = {90, -130}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Continuous.Integrator air_daily_power_draw_sw(k = 1/(60*60), use_reset = true) annotation(
      Placement(transformation(origin = {90, -240}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Continuous.Integrator air_yearly_power_draw_sw(k = 1/(60*60), use_reset = false) annotation(
      Placement(transformation(origin = {90, -210}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Continuous.Integrator air_daily_power_draw_se(k = 1/(60*60), use_reset = true) annotation(
      Placement(transformation(origin = {90, -320}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Continuous.Integrator air_yearly_power_draw_se(k = 1/(60*60), use_reset = false) annotation(
      Placement(transformation(origin = {90, -290}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.Constant zero2(k = 0) annotation(
      Placement(transformation(origin = {0, -340}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  Modelica.Blocks.Math.Sum sum_ne2(nin = 2) annotation(
      Placement(transformation(origin = {30, -50}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Math.Sum sum_nw2(nin = 2) annotation(
      Placement(transformation(origin = {30, -130}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Math.Sum sum_sw2(nin = 2) annotation(
      Placement(transformation(origin = {30, -210}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Math.Sum sum_se2(nin = 2) annotation(
      Placement(transformation(origin = {30, -290}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Continuous.Integrator air_daily_power_draw_total(k = 1/(60*60), use_reset = true) annotation(
      Placement(transformation(origin = {152, -420}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Continuous.Integrator air_yearly_power_draw_total(k = 1/(60*60), use_reset = false) annotation(
      Placement(transformation(origin = {152, -390}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Math.Sum sum_total2(nin = 4) annotation(
      Placement(transformation(origin = {90, -390}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Continuous.Integrator pump_daily_power_draw_ne(k = 1/(60*60), use_reset = true) annotation(
      Placement(transformation(origin = {388, -80}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Continuous.Integrator pump_yearly_power_draw_ne(k = 1/(60*60), use_reset = false) annotation(
      Placement(transformation(origin = {388, -50}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Continuous.Integrator pump_daily_power_draw_nw(k = 1/(60*60), use_reset = true) annotation(
      Placement(transformation(origin = {388, -160}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Continuous.Integrator pump_yearly_power_draw_nw(k = 1/(60*60), use_reset = false) annotation(
      Placement(transformation(origin = {388, -132}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Continuous.Integrator pump_daily_power_draw_sw(k = 1/(60*60), use_reset = true) annotation(
      Placement(transformation(origin = {388, -240}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Continuous.Integrator pump_yearly_power_draw_sw(k = 1/(60*60), use_reset = false) annotation(
      Placement(transformation(origin = {388, -210}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Continuous.Integrator pump_daily_power_draw_se(k = 1/(60*60), use_reset = true) annotation(
      Placement(transformation(origin = {388, -320}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Continuous.Integrator pump_yearly_power_draw_se(k = 1/(60*60), use_reset = false) annotation(
      Placement(transformation(origin = {388, -290}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Sources.Constant zero11(k = 0) annotation(
      Placement(transformation(origin = {298, -340}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  Modelica.Blocks.Math.Sum sum_ne11(nin = 2) annotation(
      Placement(transformation(origin = {328, -50}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Math.Sum sum_nw11(nin = 2) annotation(
      Placement(transformation(origin = {328, -130}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Math.Sum sum_sw11(nin = 2) annotation(
      Placement(transformation(origin = {328, -210}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Math.Sum sum_se11(nin = 2) annotation(
      Placement(transformation(origin = {328, -290}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Continuous.Integrator pump_daily_power_draw_total(k = 1/(60*60), use_reset = true) annotation(
      Placement(transformation(origin = {450, -420}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Continuous.Integrator pump_yearly_power_draw_total(k = 1/(60*60), use_reset = false) annotation(
      Placement(transformation(origin = {450, -390}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Math.Sum sum_total11(nin = 4) annotation(
      Placement(transformation(origin = {388, -390}, extent = {{-10, -10}, {10, 10}})));
  equation
    connect(pv_power_production, pv_daily_power_production.u) annotation(
      Line(points = {{-120, 80}, {-40, 80}, {-40, 110}, {-22, 110}}, color = {0, 0, 127}));
    connect(pv_power_production, pv_total_power_production.u) annotation(
      Line(points = {{-120, 80}, {-40, 80}, {-40, 70}, {-22, 70}}, color = {0, 0, 127}));
    connect(daily_trigger.y, pv_daily_power_production.reset) annotation(
      Line(points = {{79, 10}, {21, 10}, {21, 89.5}, {-4, 89.5}, {-4, 98}}, color = {255, 0, 255}));
    connect(daily_trigger.y, led_daily_power_draw_se.reset) annotation(
      Line(points = {{80, 10}, {-132, 10}, {-132, -330}, {-214, -330}}, color = {255, 0, 255}));
    connect(daily_trigger.y, led_daily_power_draw_sw.reset) annotation(
      Line(points = {{80, 10}, {-132, 10}, {-132, -252}, {-146, -252}, {-146, -250}, {-214, -250}}, color = {255, 0, 255}));
    connect(daily_trigger.y, led_daily_power_draw_nw.reset) annotation(
      Line(points = {{80, 10}, {-132, 10}, {-132, -170}, {-214, -170}}, color = {255, 0, 255}));
    connect(daily_trigger.y, led_daily_power_draw_ne.reset) annotation(
      Line(points = {{80, 10}, {-132, 10}, {-132, -90}, {-214, -90}}, color = {255, 0, 255}));
    connect(led_power_draw_ne, sum_ne.u[1]) annotation(
      Line(points = {{-370, -48}, {-292, -48}}, color = {0, 0, 127}));
    connect(led_power_draw_nw, sum_nw.u[1]) annotation(
      Line(points = {{-370, -128}, {-292, -128}}, color = {0, 0, 127}));
    connect(led_power_draw_sw, sum_sw.u[1]) annotation(
      Line(points = {{-370, -208}, {-292, -208}}, color = {0, 0, 127}));
    connect(led_power_draw_se, sum_se.u[1]) annotation(
      Line(points = {{-370, -288}, {-292, -288}}, color = {0, 0, 127}));
    connect(zero.y, sum_se.u[2]) annotation(
      Line(points = {{-310, -327}, {-310, -291}, {-292, -291}, {-292, -289}}, color = {0, 0, 127}));
    connect(zero.y, sum_sw.u[2]) annotation(
      Line(points = {{-310, -327}, {-310, -211}, {-292, -211}, {-292, -209}}, color = {0, 0, 127}));
    connect(zero.y, sum_nw.u[2]) annotation(
      Line(points = {{-310, -327}, {-310, -131}, {-292, -131}, {-292, -129}}, color = {0, 0, 127}));
    connect(zero.y, sum_ne.u[2]) annotation(
      Line(points = {{-310, -327}, {-310, -51}, {-292, -51}, {-292, -49}}, color = {0, 0, 127}));
    connect(sum_ne.y, led_yearly_power_draw_ne.u) annotation(
      Line(points = {{-269, -48}, {-233, -48}}, color = {0, 0, 127}));
    connect(sum_ne.y, led_daily_power_draw_ne.u) annotation(
      Line(points = {{-269, -48}, {-257, -48}, {-257, -78}, {-233, -78}}, color = {0, 0, 127}));
    connect(sum_nw.y, led_yearly_power_draw_nw.u) annotation(
      Line(points = {{-269, -128}, {-233, -128}}, color = {0, 0, 127}));
    connect(sum_nw.y, led_daily_power_draw_nw.u) annotation(
      Line(points = {{-269, -128}, {-255, -128}, {-255, -158}, {-233, -158}}, color = {0, 0, 127}));
    connect(sum_sw.y, led_yearly_power_draw_sw.u) annotation(
      Line(points = {{-269, -208}, {-233, -208}}, color = {0, 0, 127}));
    connect(sum_sw.y, led_daily_power_draw_sw.u) annotation(
      Line(points = {{-269, -208}, {-253, -208}, {-253, -238}, {-233, -238}}, color = {0, 0, 127}));
    connect(sum_se.y, led_yearly_power_draw_se.u) annotation(
      Line(points = {{-269, -288}, {-233, -288}}, color = {0, 0, 127}));
    connect(sum_se.y, led_daily_power_draw_se.u) annotation(
      Line(points = {{-269, -288}, {-251, -288}, {-251, -318}, {-233, -318}}, color = {0, 0, 127}));
    connect(sum_se.y, sum_total.u[4]) annotation(
      Line(points = {{-269, -288}, {-251, -288}, {-251, -384}, {-233, -384}, {-233, -388}}, color = {0, 0, 127}));
    connect(sum_sw.y, sum_total.u[3]) annotation(
      Line(points = {{-269, -208}, {-253, -208}, {-253, -386}, {-233, -386}, {-233, -388}}, color = {0, 0, 127}));
    connect(sum_nw.y, sum_total.u[2]) annotation(
      Line(points = {{-269, -128}, {-255, -128}, {-255, -388}, {-233, -388}}, color = {0, 0, 127}));
    connect(sum_ne.y, sum_total.u[1]) annotation(
      Line(points = {{-269, -48}, {-257, -48}, {-257, -388}, {-233, -388}}, color = {0, 0, 127}));
    connect(sum_total.y, led_yearly_power_draw_total.u) annotation(
      Line(points = {{-209, -388}, {-171, -388}}, color = {0, 0, 127}));
    connect(sum_total.y, led_daily_power_draw_total.u) annotation(
      Line(points = {{-209, -388}, {-191, -388}, {-191, -418}, {-171, -418}}, color = {0, 0, 127}));
    connect(daily_trigger.y, led_daily_power_draw_total.reset) annotation(
      Line(points = {{80, 10}, {-132, 10}, {-132, -430}, {-152, -430}}, color = {255, 0, 255}));
    connect(zero2.y, sum_se2.u[2]) annotation(
      Line(points = {{0, -329}, {0, -293}, {18, -293}, {18, -291}}, color = {0, 0, 127}));
    connect(zero2.y, sum_sw2.u[2]) annotation(
      Line(points = {{0, -329}, {0, -213}, {18, -213}, {18, -211}}, color = {0, 0, 127}));
    connect(zero2.y, sum_nw2.u[2]) annotation(
      Line(points = {{0, -329}, {0, -133}, {18, -133}, {18, -131}}, color = {0, 0, 127}));
    connect(zero2.y, sum_ne2.u[2]) annotation(
      Line(points = {{0, -329}, {0, -53}, {18, -53}, {18, -51}}, color = {0, 0, 127}));
    connect(sum_ne2.y, air_yearly_power_draw_ne.u) annotation(
      Line(points = {{41, -50}, {77, -50}}, color = {0, 0, 127}));
    connect(sum_ne2.y, air_daily_power_draw_ne.u) annotation(
      Line(points = {{41, -50}, {53, -50}, {53, -80}, {77, -80}}, color = {0, 0, 127}));
    connect(sum_nw2.y, air_yearly_power_draw_nw.u) annotation(
      Line(points = {{41, -130}, {77, -130}}, color = {0, 0, 127}));
    connect(sum_nw2.y, air_daily_power_draw_nw.u) annotation(
      Line(points = {{41, -130}, {55, -130}, {55, -160}, {77, -160}}, color = {0, 0, 127}));
    connect(sum_sw2.y, air_yearly_power_draw_sw.u) annotation(
      Line(points = {{41, -210}, {77, -210}}, color = {0, 0, 127}));
    connect(sum_sw2.y, air_daily_power_draw_sw.u) annotation(
      Line(points = {{41, -210}, {57, -210}, {57, -240}, {77, -240}}, color = {0, 0, 127}));
    connect(sum_se2.y, air_yearly_power_draw_se.u) annotation(
      Line(points = {{41, -290}, {77, -290}}, color = {0, 0, 127}));
    connect(sum_se2.y, air_daily_power_draw_se.u) annotation(
      Line(points = {{41, -290}, {59, -290}, {59, -320}, {77, -320}}, color = {0, 0, 127}));
    connect(sum_se2.y, sum_total2.u[4]) annotation(
      Line(points = {{41, -290}, {59, -290}, {59, -386}, {77, -386}, {77, -390}}, color = {0, 0, 127}));
    connect(sum_sw2.y, sum_total2.u[3]) annotation(
      Line(points = {{41, -210}, {57, -210}, {57, -388}, {77, -388}, {77, -390}}, color = {0, 0, 127}));
    connect(sum_nw2.y, sum_total2.u[2]) annotation(
      Line(points = {{41, -130}, {55, -130}, {55, -390}, {77, -390}}, color = {0, 0, 127}));
    connect(sum_ne2.y, sum_total2.u[1]) annotation(
      Line(points = {{41, -50}, {53, -50}, {53, -390}, {77, -390}}, color = {0, 0, 127}));
    connect(sum_total2.y, air_yearly_power_draw_total.u) annotation(
      Line(points = {{101, -390}, {139, -390}}, color = {0, 0, 127}));
    connect(sum_total2.y, air_daily_power_draw_total.u) annotation(
      Line(points = {{101, -390}, {119, -390}, {119, -420}, {139, -420}}, color = {0, 0, 127}));
    connect(daily_trigger.y, air_daily_power_draw_ne.reset) annotation(
      Line(points = {{80, 10}, {52, 10}, {52, 38}, {140, 38}, {140, -108}, {96, -108}, {96, -92}}, color = {255, 0, 255}));
    connect(daily_trigger.y, air_daily_power_draw_nw.reset) annotation(
      Line(points = {{80, 10}, {52, 10}, {52, 38}, {140, 38}, {140, -182}, {96, -182}, {96, -172}}, color = {255, 0, 255}));
    connect(daily_trigger.y, air_daily_power_draw_sw.reset) annotation(
      Line(points = {{80, 10}, {52, 10}, {52, 38}, {140, 38}, {140, -266}, {96, -266}, {96, -252}}, color = {255, 0, 255}));
    connect(daily_trigger.y, air_daily_power_draw_se.reset) annotation(
      Line(points = {{80, 10}, {52, 10}, {52, 38}, {140, 38}, {140, -348}, {96, -348}, {96, -332}}, color = {255, 0, 255}));
    connect(daily_trigger.y, air_daily_power_draw_total.reset) annotation(
      Line(points = {{80, 10}, {52, 10}, {52, 38}, {140, 38}, {140, -358}, {190, -358}, {190, -450}, {158, -450}, {158, -432}}, color = {255, 0, 255}));
    connect(air_power_draw_ne, sum_ne2.u[1]) annotation(
      Line(points = {{-60, -50}, {18, -50}}, color = {0, 0, 127}));
    connect(air_power_draw_nw, sum_nw2.u[1]) annotation(
      Line(points = {{-60, -126}, {18, -126}, {18, -130}}, color = {0, 0, 127}));
    connect(air_power_draw_sw, sum_sw2.u[1]) annotation(
      Line(points = {{-60, -210}, {18, -210}}, color = {0, 0, 127}));
    connect(air_power_draw_se, sum_se2.u[1]) annotation(
      Line(points = {{-60, -290}, {18, -290}}, color = {0, 0, 127}));
    connect(zero11.y, sum_se11.u[2]) annotation(
      Line(points = {{298, -329}, {298, -293}, {316, -293}, {316, -290}}, color = {0, 0, 127}));
    connect(zero11.y, sum_sw11.u[2]) annotation(
      Line(points = {{298, -329}, {298, -213}, {316, -213}, {316, -210}}, color = {0, 0, 127}));
    connect(zero11.y, sum_nw11.u[2]) annotation(
      Line(points = {{298, -329}, {298, -133}, {316, -133}, {316, -130}}, color = {0, 0, 127}));
    connect(zero11.y, sum_ne11.u[2]) annotation(
      Line(points = {{298, -329}, {298, -53}, {316, -53}, {316, -50}}, color = {0, 0, 127}));
    connect(sum_ne11.y, pump_yearly_power_draw_ne.u) annotation(
      Line(points = {{339, -50}, {376, -50}}, color = {0, 0, 127}));
    connect(sum_ne11.y, pump_daily_power_draw_ne.u) annotation(
      Line(points = {{339, -50}, {351, -50}, {351, -80}, {376, -80}}, color = {0, 0, 127}));
    connect(sum_nw11.y, pump_yearly_power_draw_nw.u) annotation(
      Line(points = {{339, -130}, {357.5, -130}, {357.5, -132}, {376, -132}}, color = {0, 0, 127}));
    connect(sum_nw11.y, pump_daily_power_draw_nw.u) annotation(
      Line(points = {{339, -130}, {353, -130}, {353, -160}, {376, -160}}, color = {0, 0, 127}));
    connect(sum_sw11.y, pump_yearly_power_draw_sw.u) annotation(
      Line(points = {{339, -210}, {376, -210}}, color = {0, 0, 127}));
    connect(sum_sw11.y, pump_daily_power_draw_sw.u) annotation(
      Line(points = {{339, -210}, {355, -210}, {355, -240}, {376, -240}}, color = {0, 0, 127}));
    connect(sum_se11.y, pump_yearly_power_draw_se.u) annotation(
      Line(points = {{339, -290}, {376, -290}}, color = {0, 0, 127}));
    connect(sum_se11.y, pump_daily_power_draw_se.u) annotation(
      Line(points = {{339, -290}, {357, -290}, {357, -320}, {376, -320}}, color = {0, 0, 127}));
    connect(sum_se11.y, sum_total11.u[4]) annotation(
      Line(points = {{339, -290}, {357, -290}, {357, -386}, {376, -386}, {376, -390}}, color = {0, 0, 127}));
    connect(sum_sw11.y, sum_total11.u[3]) annotation(
      Line(points = {{339, -210}, {355, -210}, {355, -388}, {376, -388}, {376, -390}}, color = {0, 0, 127}));
    connect(sum_nw11.y, sum_total11.u[2]) annotation(
      Line(points = {{339, -130}, {353, -130}, {353, -390}, {376, -390}}, color = {0, 0, 127}));
    connect(sum_ne11.y, sum_total11.u[1]) annotation(
      Line(points = {{339, -50}, {351, -50}, {351, -390}, {376, -390}}, color = {0, 0, 127}));
    connect(sum_total11.y, pump_yearly_power_draw_total.u) annotation(
      Line(points = {{399, -390}, {438, -390}}, color = {0, 0, 127}));
    connect(sum_total11.y, pump_daily_power_draw_total.u) annotation(
      Line(points = {{399, -390}, {417, -390}, {417, -420}, {438, -420}}, color = {0, 0, 127}));
  connect(pump_power_draw_ne, sum_ne11.u[1]) annotation(
      Line(points = {{240, -52}, {316, -52}, {316, -50}}, color = {0, 0, 127}));
  connect(pump_power_draw_nw, sum_nw11.u[1]) annotation(
      Line(points = {{240, -130}, {316, -130}}, color = {0, 0, 127}));
  connect(pump_power_draw_sw, sum_sw11.u[1]) annotation(
      Line(points = {{238, -210}, {316, -210}}, color = {0, 0, 127}));
  connect(pump_power_draw_se, sum_se11.u[1]) annotation(
      Line(points = {{240, -290}, {316, -290}}, color = {0, 0, 127}));
  connect(daily_trigger.y, pump_daily_power_draw_ne.reset) annotation(
      Line(points = {{80, 10}, {52, 10}, {52, 36}, {440, 36}, {440, -110}, {394, -110}, {394, -92}}, color = {255, 0, 255}));
  connect(daily_trigger.y, pump_daily_power_draw_nw.reset) annotation(
      Line(points = {{80, 10}, {52, 10}, {52, 36}, {440, 36}, {440, -186}, {394, -186}, {394, -172}}, color = {255, 0, 255}));
  connect(daily_trigger.y, pump_daily_power_draw_sw.reset) annotation(
      Line(points = {{80, 10}, {52, 10}, {52, 36}, {440, 36}, {440, -264}, {394, -264}, {394, -252}}, color = {255, 0, 255}));
  connect(daily_trigger.y, pump_daily_power_draw_se.reset) annotation(
      Line(points = {{80, 10}, {52, 10}, {52, 36}, {440, 36}, {440, -344}, {394, -344}, {394, -332}}, color = {255, 0, 255}));
  connect(daily_trigger.y, pump_daily_power_draw_total.reset) annotation(
      Line(points = {{80, 10}, {52, 10}, {52, 36}, {440, 36}, {440, -354}, {488, -354}, {488, -448}, {456, -448}, {456, -432}}, color = {255, 0, 255}));
    annotation(
      Icon(graphics = {Bitmap(extent = {{-100, -100}, {100, 100}}, imageSource = "iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAAAXNSR0IArs4c6QAABF9JREFUeF7tXW1y2jAQlchF0pO0uUDIBUr5OAjJPQoh9AIhF0h6knIRUMcGBgP+Wj3ZK8xjptP80Hq17+3blWWPbA1/qghYVe90bkiAchKQABKgjICyeyqABCgjoOyeCiAByggou6cCSIAyAsruqYBrJ2D49Hrf2/aG1trvxpgf3vE4Y67kvnxtjEn/bezmZfE+Sv72/kEKGD++PVtrpzvkEgRv7rd2zr3NP349+0buTcCkv/yEMj53xhiRiiL6mq0GDz4keBEwfloOrTOvPg6jtgEYdM69+ChBTEBS8+/c3b+ogRRPDlPe3t3a2s3o9/voS+JeTMCx7kvcoGODAJROAkjyyiB8VCAmYNJfJqVnWDmb2xwg7gU+BCTl514T3yazODeu+g7Xs9XgmwQbHwJK15uz1UB8TcmEtcdO+sug8YvByk4gLzFIgCwBIQKO2XhskiSABDRapaIqQXmRUgEqCjhSQQKUCDg0ZBKgRMBBA7dCQJpwOctAafyBVkHtlaDQTVDasQv978kgAQ3fCIZOACpAKAEScNgKKNifkZYAIf4mDgJKNqdUAMjMpy3/RRBI/bMECSUQhwJKJi3NAGH8aQko2x1uw3/ZnKX+qQBhBlABgffjhfhXNGFrZqufoqQWDU4mGzoDvAG46VWQcg8IWYO9E6DAsPEeMO4vXaFsnDGzD9lmlDYA2v5bL0FoCUPtpYCfj588Ll3ZO6yNKwAFYGdf/J5PVQCof5iAwIsANQX43kmSADADUABReypAmUAS0DkCTvtZVQ87TwC1HlCUiVUBqJegbqyCigtB9ASACqYCwCYQWoFhS1CNO2E0ANQexD/4XpgfAcATMRRA1L4bBACbcSiAp/aXd9RVPYQEKG9FhE2ASzqlCeBXgqJRgBwAEgAu47QBRP0HXIbm72hWSRANIJx9E/O/gkeS4QDMr4OxJ0BABVwnANoJ0DEC5MtQEpBtwh7v12sDWO7/VnvAnsjkv3nF6+naBHasBPE+4AKB2FchVABvxE6S9gq3Iv64suPR2lKg71sd7AGoAnMeSWbJqEqASAjQezGrez3gTItVGaANgLZ/dQWgD0S07VECSQDIIAkAAUTNr34rAgVA2z4qBfDIMvlWCHsAKKGIFOD3SA+MX908IgL8noipIwhOIAYCWju4tf55qQmq0uONpePrMNfOAxnBsfVNBFkHiHpjik69ylrLksAsZqvBqJ733Sjxbqj/4d1SMqTjJWE3M7aVw7uR4+uF2YSjlHWY4zz0fKzdPDR+fH2CSmc/4ABQ7qwZzd8HC+klxCXo4KD4GPts6dj/HTrVpFE2P158bP1hSt4EpEpIPuJj7FTeSQSIRE6eT93PRg8RkFxo3xOm+28K+H/GSsBJOlSPmPQTVs65v9vedqH6GSspZhx/iQCsAIKKIUACMPxgaxIAQ4hdgARg+MHWJACGELsACcDwg61JAAwhdgESgOEHW5MAGELsAiQAww+2JgEwhNgF/gPiE8SOH3BlnAAAAABJRU5ErkJggg==")}, coordinateSystem(extent = {{-100, -100}, {100, 100}})),
      Diagram(graphics = {Rectangle(origin = {60, -120}, extent = {{-470, 424}, {470, -424}})}, coordinateSystem(extent = {{-420, 300}, {540, -560}})));
  end analysis;

  model analysis_insulation
    Modelica.Blocks.Interfaces.RealInput heat_flow_ne annotation(
      Placement(transformation(origin = {-120, -40}, extent = {{-20, -20}, {20, 20}}), iconTransformation(origin = {-110, -40}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Interfaces.RealInput heat_flow_nw annotation(
      Placement(transformation(origin = {-120, -120}, extent = {{-20, -20}, {20, 20}}), iconTransformation(origin = {-110, -60}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Interfaces.RealInput heat_flow_sw annotation(
      Placement(transformation(origin = {-120, -200}, extent = {{-20, -20}, {20, 20}}), iconTransformation(origin = {-110, -80}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Interfaces.RealInput heat_flow_se annotation(
      Placement(transformation(origin = {-120, -280}, extent = {{-20, -20}, {20, 20}}), iconTransformation(origin = {-110, -100}, extent = {{-10, -10}, {10, 10}})));
    Buildings.BoundaryConditions.WeatherData.Bus weaBus annotation(
      Placement(transformation(origin = {-20, 0}, extent = {{-110, -10}, {-90, 10}}), iconTransformation(extent = {{-120, -18}, {-80, 22}})));
    parameter Modelica.Units.SI.Area area_ne "total area of the north-east building face";
    parameter Modelica.Units.SI.Area area_nw "total area of the north-west building face";
    parameter Modelica.Units.SI.Area area_sw "total area of the south-west building face";
    parameter Modelica.Units.SI.Area area_se "total area of the south-east building face";
    parameter Modelica.Units.SI.Temperature temperature_building "inside temperature of the building";
    Modelica.Units.SI.TemperatureDifference temperature_diff;
    Real u_ne;
    Real u_nw;
    Real u_sw;
    Real u_se;
    Buildings.Controls.OBC.CDL.Reals.MovingAverage u_ne_avg(delta = 86400) annotation(
      Placement(transformation(origin = {10, -40}, extent = {{-10, -10}, {10, 10}})));
    Buildings.Controls.OBC.CDL.Reals.MovingAverage u_nw_avg(delta = 86400) annotation(
      Placement(transformation(origin = {10, -120}, extent = {{-10, -10}, {10, 10}})));
    Buildings.Controls.OBC.CDL.Reals.MovingAverage u_sw_avg(delta = 86400) annotation(
      Placement(transformation(origin = {10, -200}, extent = {{-10, -10}, {10, 10}})));
    Buildings.Controls.OBC.CDL.Reals.MovingAverage u_se_avg(delta = 86400) annotation(
      Placement(transformation(origin = {10, -280}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Sources.RealExpression realExpression(y = u_ne) annotation(
      Placement(transformation(origin = {-40, -40}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Sources.RealExpression realExpression1(y = u_nw) annotation(
      Placement(transformation(origin = {-40, -120}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Sources.RealExpression realExpression2(y = u_sw) annotation(
      Placement(transformation(origin = {-42, -200}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Sources.RealExpression realExpression3(y = u_se) annotation(
      Placement(transformation(origin = {-40, -280}, extent = {{-10, -10}, {10, 10}})));
  equation
    temperature_diff = abs(temperature_building - weaBus.TDryBul);
    u_ne = max(heat_flow_ne/max(area_ne*temperature_diff, Modelica.Constants.small), 0);
    u_nw = max(heat_flow_nw/max(area_nw*temperature_diff, Modelica.Constants.small), 0);
    u_sw = max(heat_flow_sw/max(area_sw*temperature_diff, Modelica.Constants.small), 0);
    u_se = max(heat_flow_se/max(area_se*temperature_diff, Modelica.Constants.small), 0);
    connect(realExpression.y, u_ne_avg.u) annotation(
      Line(points = {{-28, -40}, {-2, -40}}, color = {0, 0, 127}));
    connect(realExpression1.y, u_nw_avg.u) annotation(
      Line(points = {{-28, -120}, {-2, -120}}, color = {0, 0, 127}));
    connect(realExpression2.y, u_sw_avg.u) annotation(
      Line(points = {{-30, -200}, {-2, -200}}, color = {0, 0, 127}));
    connect(realExpression3.y, u_se_avg.u) annotation(
      Line(points = {{-28, -280}, {-2, -280}}, color = {0, 0, 127}));
    annotation(
      Icon(graphics = {Bitmap(extent = {{-100, -100}, {100, 100}}, imageSource = "iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAAAXNSR0IArs4c6QAABI9JREFUeF7tnUty00AQhltODmJ2ugKYFOQGFDF7+yRJThKzxqG4AaHAcAXt4oMkFiXLTvyQZtzzz6hHUnutefT/9WNGkkcJ6U9UgUR0dB2cFICwEygABSCsgPDwGgEKQFgB4eE1AhSAsALCw2sEKABhBYSH1whoO4D0078hnT9NiM4+EOUfTfYUtHNhg+Hhc1pSkiwpz5f0PLjNfrxdIn1CEZCOf98QDa6RCcTdtnAXg0RrGKuv2fzixtUOZwDp+O9Pm8e7Tip4O++hmDxk83eXLvN2ApB+XkxoQHe7A3q3ac8aiyfWWb6elGNbtnGrW5dIYAMoc/7qsdLmsBRcHKy5NmU6mmbziwfOoHwAnc/7ZdZ3Wyzwo4AP4GpxRwlNyikeNHefOcdp3K5tZG78WuAC4JESGrqpEFOrutoA1Iycltn96A3HSj6A8cIYndl8xO6TM2Hpa1Oj/Tll8/cs+1kXF8abJ0DUbwB8+xUAM6R8O6BfAAlR9q3PKUgjgOnP/MvjjgDiewBfAtkW6ZdFbtokcGug3xTUBwCeV4EKgBlQvU9BvgVg6u99Gd66CKgHUO5guTlYATAV0AjYKUKv97de75+0zQOZ/GNKQdU3rRQAbyPaoRpQ+nLbHEABMHOQ7xqkALgA+r4T9u2BTP03Rbj+oQ03BbYyAkxPF7kCuAGob8Udv5UATKJxBYAAVHgCd3wFwCTgOwU2DqDWgI032TzItwBM/WPaiFVPPbSABQCtAQa3aQJANDWgYiI2+w+bxJOCNjOzGaApyPRE6ISH8qiAaHtuzj+83vf4shHgsIzzLQAXiO/xZQE45FDfArQYwGY7fuDF1hwO3ktRAOBbAaiAaHuux0dcAxwjQBhghwDIbcR0H6AbsVoFbDWw9xsxtIag7T0CcHsojxrQ9vYeAbSzBkgD7DCA096MgwGA+5gOAyhNsxVBGAC4jA4GYLshjl2AzgLYklUAkb8Zt+eBAndDfUVA3VM5mwMGS0HtiwCZZbQCAIsoGkHBAGgRPm0VJg8AXEejHijdPhiAVtUAw3sttiKKAlQA3awBp/9FCfWgdPwnNx2k5+rBTdWwaCLgKAswXk1E3ozDHaA4rmd3Cbu/nLU5gDgA9JGgdPtjgPsw9LygwITQCPITAcAqIrA+wbuPAwDwTDi4QoEHqAKw+3/pZlJQRQScuooIrE/w7ksAOSWUbE6t0SIcXPTdATQFNSr38WAKIAoAkn9TvVrUH9xq2iF5E67e+EaGt9jRwEasxcfWH4jnHVhOs+x+NOX4msP/A7r+0QaOfIfXNnF4t+n4+s18vHsWokmQtnVpcHUZ/Pj6wp6qDzgc2dl9Cvsmr2iafR/NuLzZKWg7QPpyjP3xkH60B04x56pQdRQ/qw/+sfXb7p0BrCOhBx9zsHPg5/3dPiEAawhFTThbXVOSDF0/6uMnYuxSvVzhOOC62fYzVvT8i57OZ6KfsWKYrJfWKABHgCqLKaAAMP3g1goAlhDrQAFg+sGtFQAsIdaBAsD0g1srAFhCrAMFgOkHt1YAsIRYBwoA0w9urQBgCbEO/gOX7M2O9c+W1gAAAABJRU5ErkJggg==")}),
      Diagram(graphics = {Rectangle(origin = {17, -68}, extent = {{-599, 566}, {599, -566}})}, coordinateSystem(extent = {{-580, 500}, {620, -640}})));
  end analysis_insulation;

  model analysis_yield
  Modelica.Blocks.Interfaces.RealInput yield_ne annotation(
      Placement(transformation(origin = {-120, -40}, extent = {{-20, -20}, {20, 20}}), iconTransformation(origin = {-110, -40}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Interfaces.RealInput yield_nw annotation(
      Placement(transformation(origin = {-120, -72}, extent = {{-20, -20}, {20, 20}}), iconTransformation(origin = {-110, -60}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Interfaces.RealInput yield_sw annotation(
      Placement(transformation(origin = {-120, -100}, extent = {{-20, -20}, {20, 20}}), iconTransformation(origin = {-110, -80}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Blocks.Interfaces.RealInput yield_se annotation(
      Placement(transformation(origin = {-120, -128}, extent = {{-20, -20}, {20, 20}}), iconTransformation(origin = {-110, -100}, extent = {{-10, -10}, {10, 10}})));
  equation

  annotation(
      Icon(graphics = {Bitmap(origin = {0, 100}, extent = {{-100, -200}, {100, 0}}, imageSource = "iVBORw0KGgoAAAANSUhEUgAAAGAAAABgCAYAAADimHc4AAAAAXNSR0IArs4c6QAABHNJREFUeF7tXUFy00AQHJE8JNwoYt4AuZlfOH/AVdxIblQ5f0h4BboR3oBDcYsfkiCQ4sSyJVma7bFmZXdupHZ2drqnZ0fKskqEP64IJK7e6VxIgHMSkAAS4IyAs3sqgAQ4I+DsngogAc4IOLunAkiAMwLO7qmAoRMw/vzmRB6OJ5LIexH54BxPD+6zhUiyEMkWcvR4mX79s0CcQgoYfxpdSCJfpHijkb2sY/1fq+U1/b4cQGVMPm2y6QEJ2dI2W0iWfEuv5hehswYTMJ6OfhxGxneC9jadzc86jdwYFETAePp2IvLqemvmhqxmyDaZXIYoQU1AUfMfj++bsepSaLogbTVPF1+KMcuSWC2sRTk6T6/mt4rZ9K+jV3Vf42bb2EiBrlly60oDVKBXwPT0+n8XMLGCP455WqHtukz1XhBCwL1IcrIfOV0fRTgd2SKd3b3uylY+LoCA0arfrPGUzubqOTUL9h47ntrGrwbLegHegGr9P8Xf/NyjTUASoGTAOgFJwNAI+DgdZds2Aa0ElfG7D49KAXXdAgnQNSEsQUpNVRWwTMPlE7I2AUkATMD6BHtPgHUNVuIv1v4HpwBrAAZNgMcmTAKMH8W9M9DbP1SCqIBN+hJJZ79UmKoG5+7iKgHVFNB2IYNWQN3ihwaALQEHp4AqfENLgIGXIBLwdByoROPQMtC2BIlo46cClAxYNyEk4NAIQDMItVfiXRlu7b93BaABoPYkAHyVQQKcASQBJMD0XBT3AOWmYK1AXwIqJ43bH2SsAVDib/4y0peAgKONJMBkD2g+/tr2KE8CTAhoFj4JaCmKaAY+2zdpgAT0RECTGxIwRAJK3dSuCUQrQPWvyMo+DF1A3fn68hKCAEQIUP49A40/EgL8NmEUQNR+rwnIN/bvLf9FCgUQbSJ2S0Amkl5tP55tBYDZJs4StN6QBu0BJTZit9+tAmT373KiU9AGom0JQAKcn+TjIgB+G6o/muitoLgIMH4b2mcXFNoE2BPQcxdhmcEhp7tR//YEKDchNABLexJgXILy6dq6kAqBsIJ1bTQVwC4Iu23EsgQV2QgrYD2n2xS4Bwp4l5VvaCwHdJhdUIHAqg5qM0D5Ntx9OKpAdwW4IwguoErAwDZhMH5386gUENJHuyMILoAKAAFEzRsVwNtSUGi72UdVguqWfFBdUMDbXHZB3RK9cVQECjhtvbgVjLFf85osRhagrQAhp6P119YbB4kAtAvbVTeY3aSzu3ONDz0BLx9t6OKm+RR0F+vOY2IhuJfLu1uvr+8M2/AHbhKfydnOr6/PUav7gEMvaMaS6bXB/j1PZ79vtDioS9Czg3HoNfa9g9hLGVRfW/+MYzABhRJU+4E2N5rHlyHtBd5tSw+o++XpIAIKEoo94Sj/klL+TYHwz1j1rozQhFh+xiqTn3L8cOP6GavQEGi3QgBWAMHEECABGH6wNQmAIcQmIAEYfrA1CYAhxCYgARh+sDUJgCHEJiABGH6wNQmAIcQmIAEYfrA1CYAhxCb4Byn9vY5cOlskAAAAAElFTkSuQmCC")}));
end analysis_yield;
  annotation(
    uses(Buildings(version = "11.0.0"), Modelica(version = "4.0.0")));
end vertical;