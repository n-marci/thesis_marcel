package greenshell
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
      if (par == 0) then
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
      Modelica.Blocks.Interfaces.RealOutput y annotation(
        Placement(transformation(origin = {110, 0}, extent = {{-10, -10}, {10, 10}}), iconTransformation(origin = {110, 0}, extent = {{-10, -10}, {10, 10}})));
      parameter Modelica.Units.SI.Area a_plant "planting area";
      Real x_nsdw(start = 0.001, final unit = "g/m2") "state variable non-structural dry weight";
      Real x_sdw(start = 0.001, final unit = "g/m2") "state variable structural dry weight";
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
      Real g_bnd = /*(final unit = "m/s")*/0.007 "Boundary layer conductance";
      Real g_stm = /*(final unit = "m/s")*/0.005 "Stomatal conductance";
      Real g_car "Carboxylation conductance";
      /*(final unit = "m/s")*/
      Real r_gr(final unit = "1/s") "specific growth rate";
      Real c_gr_max(final unit = "1/s") = 5e-6 "saturation growth rate at 20°C";
      Real c_gamma(final unit = "1") = 1.0 "growth rate coefficient";
      Real c_q_10_gr(final unit = "1") = 1.6 "growth rate sensitivity to the canopy temperature";
    equation
      der(x_nsdw) = c_alpha*f_phot - r_gr*x_sdw - f_resp - ((1 - c_beta)/c_beta)*r_gr*x_sdw;
      f_resp = (c_resp_sht*(1 - c_tau)*x_sdw + c_resp_rt*c_tau*x_sdw)*c_q_10_resp^(((air_temperature - 273.15) - 25)/10);
      f_phot = (1 - exp(-c_k*c_lar*(1 - c_tau)*x_sdw))*f_phot_max;
      f_phot_max = (epsilon*par*g_co2*c_omega*(co2_concentration - capital_gamma))/(epsilon*par + g_co2*c_omega*(co2_concentration - capital_gamma));
      capital_gamma = c_cap_gamma*c_q_10_cap_gamma^(((air_temperature - 273.15) - 20)/10);
      epsilon = c_epsilon*(co2_concentration - capital_gamma)/(co2_concentration + 2*capital_gamma);
      g_co2 = 1/((1/g_bnd) + (1/g_stm) + (1/g_car));
      g_car = -1.32e-5*(air_temperature - 273.15)^2 + 5.94e-4*(air_temperature - 273.15) - 2.64e-3;
      der(x_sdw) = r_gr*x_sdw;
      r_gr = c_gr_max*(x_nsdw/(c_gamma*x_sdw + x_nsdw))*c_q_10_gr^(((air_temperature - 273.15) - 20)/10);
      dw = a_plant*(x_nsdw + x_sdw)/*/ 1000*/;
      y = dw;
    end plant_yield;

    plant_transpiration transpiration(a_plant = a_plant) annotation(
      Placement(transformation(origin = {10, -30}, extent = {{-10, -10}, {10, 10}})));
    plant_yield yield(x_nsdw(start = 0.25*0.72), x_sdw(start = 0.75*0.72), a_plant = a_plant) annotation(
      Placement(transformation(origin = {10, 50}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Fluid.Sensors.Pressure pressure(redeclare package Medium = Buildings.Media.Air "Moist air") annotation(
      Placement(transformation(origin = {-70, 10}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Fluid.Sensors.Temperature temperature(redeclare package Medium = Buildings.Media.Air "Moist air") annotation(
      Placement(transformation(origin = {-90, 30}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Sources.Constant co2_concentration(k = 400) annotation(
      Placement(transformation(origin = {-20, 80}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Blocks.Sources.Constant basal_crop_coefficient(k = 0.1) annotation(
      Placement(transformation(origin = {-20, 2}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Blocks.Sources.Constant relative_humidity(k = 0.8) annotation(
      Placement(transformation(origin = {-20, -82}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  equation
    connect(pressure.port, farm_air) annotation(
      Line(points = {{-70, 0}, {-100, 0}}, color = {0, 127, 255}));
    connect(temperature.port, farm_air) annotation(
      Line(points = {{-90, 20}, {-90, 0}, {-100, 0}}, color = {0, 127, 255}));
    connect(yield.air_temperature, temperature.T) annotation(
      Line(points = {{-2, 58}, {-60, 58}, {-60, 30}, {-82, 30}}, color = {0, 0, 127}));
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
    connect(yield.par, par) annotation(
      Line(points = {{-2, 42}, {-40, 42}, {-40, 80}, {-120, 80}}, color = {0, 0, 127}));
    connect(transpiration.v_air, v_air) annotation(
      Line(points = {{-2, -20}, {-10, -20}, {-10, -60}, {-120, -60}}, color = {0, 0, 127}));
    connect(yield.y, y) annotation(
      Line(points = {{22, 50}, {60, 50}, {60, 0}, {110, 0}}, color = {0, 0, 127}));
    connect(relative_humidity.y, transpiration.relative_humidity) annotation(
      Line(points = {{-20, -70}, {-20, -34}, {0, -34}}, color = {0, 0, 127}));
    annotation(
      Icon(graphics = {Bitmap(extent = {{-100, -100}, {100, 100}}, imageSource = "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")}));
  end plant;

  model led
  equation

  end led;

  model pump
  equation

  end pump;

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
        parameter Real light_setpoint "setpoint for natural lighting";
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
          Placement(transformation(origin = {10, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
        Modelica.Blocks.Math.Product shade annotation(
          Placement(transformation(origin = {10, 0}, extent = {{-10, -10}, {10, 10}})));
      equation
        connect(shading_control.y, shade.u2) annotation(
          Line(points = {{0, -44}, {-10, -44}, {-10, -6}, {-2, -6}}, color = {0, 0, 127}));
        connect(add.y, shade.u1) annotation(
          Line(points = {{-18, 0}, {-14, 0}, {-14, 6}, {-2, 6}}, color = {0, 0, 127}));
        connect(HDifTilIso.H, add.u1) annotation(
          Line(points = {{-58, 10}, {-50, 10}, {-50, 6}, {-42, 6}}, color = {0, 0, 127}));
        connect(HDirTil.H, add.u2) annotation(
          Line(points = {{-58, -10}, {-50, -10}, {-50, -6}, {-42, -6}}, color = {0, 0, 127}));
        connect(shade.y, y) annotation(
          Line(points = {{22, 0}, {110, 0}}, color = {0, 0, 127}));
        connect(weaBus, HDifTilIso.weaBus) annotation(
          Line(points = {{-100, 0}, {-88, 0}, {-88, 10}, {-80, 10}}, thickness = 0.5));
        connect(weaBus, HDirTil.weaBus) annotation(
          Line(points = {{-100, 0}, {-88, 0}, {-88, -10}, {-80, -10}}, thickness = 0.5));
  connect(shading_control.u_m, shade.y) annotation(
          Line(points = {{10, -56}, {10, -60}, {40, -60}, {40, 0}, {22, 0}}, color = {0, 0, 127}));
  connect(shading_control.u_s, light_set.y) annotation(
          Line(points = {{22, -44}, {30, -44}, {30, -66}, {10, -66}, {10, -68}}, color = {0, 0, 127}));
        annotation(
          Icon(graphics = {Bitmap(extent = {{-100, 100}, {100, -100}}, imageSource = "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")}));
      end shading;
    equation

    end illumination_control;

    model atmosphere_control
      model window
        
        extends Buildings.Fluid.Interfaces.PartialTwoPortTransport;
        replaceable package Medium = Modelica.Media.Interfaces.PartialMedium "Medium in the component" annotation(
          choices(choice(redeclare package Medium = Buildings.Media.Air "Moist air")));
      //  parameter Integer nPorts=4 "Number of ports"
      //    annotation(Evaluate=true, Dialog(connectorSizing=true, tab="General",group="Ports"));
        Modelica.Blocks.Interfaces.RealInput T annotation (Placement(transformation(origin = {0, 60},extent = {{140, -20}, {100, 20}}, rotation = -0), iconTransformation(origin = {0, 80},extent = {{140, -20}, {100, 20}}, rotation = -0)));
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

    Buildings.BoundaryConditions.WeatherData.ReaderTMY3 weaDat1(filNam = "/home/marci/Downloads/DEU_BY_Nurnberg.AP.107630_TMYx.2007-2021/DEU_BY_Nurnberg.AP.107630_TMYx.2007-2021.mos") annotation(
      Placement(transformation(origin = {0, -20}, extent = {{-60, 20}, {-40, 40}})));
    illumination_control.shading shading1(light_setpoint = 400, azimuth = 0.7330382858376184) annotation(
      Placement(transformation(origin = {-10, 10}, extent = {{-10, -10}, {10, 10}})));
    Buildings.HeatTransfer.Radiosity.OpaqueSurface sur(A = 1, absIR = 0.9) annotation(
      Placement(transformation(origin = {50, 10}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
    Buildings.HeatTransfer.Sources.FixedTemperature preTem(T = 293.15) annotation(
      Placement(transformation(origin = {10, -30}, extent = {{-10, -10}, {10, 10}})));
  equation
    connect(weaDat1.weaBus, shading1.weaBus) annotation(
      Line(points = {{-40, 10}, {-20, 10}}, color = {255, 204, 51}, thickness = 0.5));
    connect(sur.JIn, shading1.y) annotation(
      Line(points = {{40, 6}, {8, 6}, {8, 10}, {2, 10}}));
    connect(preTem.port, sur.heatPort) annotation(
      Line(points = {{20, -30}, {50, -30}, {50, 0}}, color = {191, 0, 0}));
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
      Buildings.HeatTransfer.Conduction.MultiLayer conduction(A = area, layers (nLay = 2, material = {concrete, eps})) annotation(
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
  Modelica.Thermal.HeatTransfer.Components.ThermalConductor thermalConductor(G = 346.4900662251656)  annotation(
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

  model system
    farm.illumination_control.shading shading_sw(light_setpoint = 400, azimuth = 0.7330382858376184) annotation(
      Placement(transformation(origin = {-70, 10}, extent = {{-10, -10}, {10, 10}})));
    plant lettuce_sw(redeclare package Medium = Buildings.Media.Air "Moist air", a_plant = 523.2) annotation(
      Placement(transformation(origin = {0, -20}, extent = {{-10, -10}, {10, 10}})));
    Buildings.Fluid.MixingVolumes.MixingVolume farm_air_sw(redeclare package Medium = Buildings.Media.Air "Moist air", m_flow_nominal = 0.001, V = 594, nPorts = 2) annotation(
      Placement(transformation(origin = {0, 20}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    farm.atmosphere_control.window window_sw(redeclare package Medium = Buildings.Media.Air "Moist air", window_width = 52.8) annotation(
      Placement(transformation(origin = {-70, 30}, extent = {{-10, -10}, {10, 10}})));
    Buildings.BoundaryConditions.WeatherData.ReaderTMY3 weaDat1(filNam = "/home/marci/Downloads/DEU_BY_Nurnberg.AP.107630_TMYx.2007-2021/DEU_BY_Nurnberg.AP.107630_TMYx.2007-2021.mos") annotation(
      Placement(transformation(origin = {-100, -20}, extent = {{-60, 20}, {-40, 40}})));
    Buildings.Fluid.Sources.Outside outside_air(redeclare package Medium = Buildings.Media.Air "Moist air", nPorts = 1) annotation(
      Placement(transformation(origin = {-110, 30}, extent = {{-10, -10}, {10, 10}})));
    Modelica.Blocks.Sources.Constant air_speed(k = 0) annotation(
      Placement(transformation(origin = {0, -70}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  equation
    connect(air_speed.y, lettuce_sw.v_air) annotation(
      Line(points = {{0, -58}, {0, -46}, {-40, -46}, {-40, -28}, {-12, -28}}, color = {0, 0, 127}));
    connect(lettuce_sw.par, shading_sw.y) annotation(
      Line(points = {{-12, -12}, {-32, -12}, {-32, 10}, {-58, 10}}, color = {0, 0, 127}));
    connect(lettuce_sw.farm_air, farm_air_sw.ports[1]) annotation(
      Line(points = {{-10, -20}, {-20, -20}, {-20, 20}, {-10, 20}}, color = {0, 127, 255}));
    connect(farm_air_sw.ports[2], window_sw.farm_air) annotation(
      Line(points = {{-10, 20}, {-20, 20}, {-20, 30}, {-60, 30}}, color = {0, 127, 255}));
    connect(window_sw.city_air, outside_air.ports[1]) annotation(
      Line(points = {{-80, 30}, {-100, 30}}, color = {0, 127, 255}));
    connect(weaDat1.weaBus, outside_air.weaBus) annotation(
      Line(points = {{-140, 10}, {-130, 10}, {-130, 30}, {-120, 30}}, color = {255, 204, 51}, thickness = 0.5));
    connect(weaDat1.weaBus, shading_sw.weaBus) annotation(
      Line(points = {{-140, 10}, {-80, 10}}, color = {255, 204, 51}, thickness = 0.5));
    annotation(
      experiment(StartTime = 0, StopTime = 3.07584e+07, Tolerance = 1e-06, Interval = 307.584));
  end system;

  model building2
  building.facade2 facade_sw(area = 523.2, azimuth = 0.7330382858376184)  annotation(
      Placement(transformation(origin = {-50, -50}, extent = {{-10, -10}, {10, 10}})));
  Buildings.HeatTransfer.Sources.FixedTemperature air_temperature(T = 293.15) annotation(
      Placement(transformation(extent = {{-10, -10}, {10, 10}})));
  building.facade2 facade_ne(area = 523.2, azimuth = 3.8746309394274117)  annotation(
      Placement(transformation(origin = {50, 30}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
  building.facade2 facade_nw(area = 369.7, azimuth = 2.303834612632515)  annotation(
      Placement(transformation(origin = {-50, 50}, extent = {{-10, -10}, {10, 10}})));
  building.facade2 facade_se(area = 369.7, azimuth = 5.445427266222308)  annotation(
      Placement(transformation(origin = {50, -50}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
  building.window2 window_sw(area = 581.04, azimuth = 0.7330382858376184)  annotation(
      Placement(transformation(origin = {-50, -28}, extent = {{-10, -10}, {10, 10}})));
  building.window2 window_nw(area = 118.8, azimuth = 2.303834612632515)  annotation(
      Placement(transformation(origin = {-50, 72}, extent = {{-10, -10}, {10, 10}})));
  building.window2 window_ne(area = 581.04, azimuth = 3.8746309394274117)  annotation(
      Placement(transformation(origin = {52, 54}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
  building.window2 window_se(area = 118.8, azimuth = 5.445427266222308)  annotation(
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
  Buildings.Fluid.MixingVolumes.MixingVolume vol(nPorts = 1, redeclare package Medium = Buildings.Media.Air "Moist air", m_flow_nominal = 0.001, V = 1)  annotation(
      Placement(transformation(origin = {-70, 10}, extent = {{-10, -10}, {10, 10}})));
  Buildings.Fluid.MixingVolumes.MixingVolume vol1(nPorts = 2, redeclare package Medium = Buildings.Media.Air "Moist air", m_flow_nominal = 0.001, V = 1)  annotation(
      Placement(transformation(origin = {108, 4}, extent = {{-10, -10}, {10, 10}})));
  Buildings.HeatTransfer.Sources.FixedTemperature preTem(T = 283.15)  annotation(
      Placement(transformation(origin = {-110, 10}, extent = {{-10, -10}, {10, 10}})));
  Buildings.HeatTransfer.Sources.PrescribedTemperature preTem1 annotation(
      Placement(transformation(origin = {108, 50}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
  Modelica.Blocks.Sources.Ramp ramp(height = 20, duration = 864000/4, offset = 273.15 +10, startTime = 864000/2)  annotation(
      Placement(transformation(origin = {138, 76}, extent = {{-10, -10}, {10, 10}})));
  Buildings.Fluid.Sensors.Temperature senTem(redeclare package Medium = Buildings.Media.Air "Moist air")  annotation(
      Placement(transformation(origin = {54, 18}, extent = {{10, -10}, {-10, 10}}, rotation = -0)));
  farm.atmosphere_control.window window1(redeclare package Medium = Buildings.Media.Air "Moist air", window_width = 1, m_flow_nominal = 0.001)  annotation(
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
  annotation(
    uses(Buildings(version = "11.0.0"), Modelica(version = "4.0.0")));
end greenshell;