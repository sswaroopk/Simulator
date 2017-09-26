package unitoperations
  model parameters
    constant Integer noc = 4;
    constant Chemsep_Database.Toluene comp1;
    constant Chemsep_Database.Hydrogen comp2;
    constant Chemsep_Database.Benzene comp3;
    constant Chemsep_Database.Methane comp4;
    constant Chemsep_Database.General_Properties components[noc] = {comp1, comp2, comp3, comp4};
  end parameters;

  model compounds
    protected
    parameter Integer NOC = 4;
    parameter Chemsep_Database.Toluene comp1 annotation(Dialog(group = "compounds"));
    parameter Chemsep_Database.Hydrogen comp2 annotation(Dialog(group = "compounds"));
    parameter Chemsep_Database.Benzene comp3 annotation(Dialog(group = "compounds"));
    parameter Chemsep_Database.Methane comp4 annotation(Dialog(group = "compounds"));
    parameter Chemsep_Database.General_Properties comp[NOC] = {comp1, comp2, comp3, comp4} annotation(Dialog(group = "compounds"));
  end compounds;

  connector sensor
    Real var;
    annotation(Diagram(graphics = {Ellipse(origin = {-1, 1}, fillColor = {65, 252, 255}, fillPattern = FillPattern.Solid, extent = {{-67, 65}, {67, -65}}, endAngle = 360)}), Icon(graphics = {Ellipse(origin = {-4, -1}, fillColor = {255, 85, 0}, fillPattern = FillPattern.Solid, extent = {{-60, 61}, {60, -61}}, endAngle = 360)}, coordinateSystem(initialScale = 0.1)));
  end sensor;

  connector port
    parameter Integer NOC = 4;
    Real moleflow;
    Real pressure;
    Real temperature;
    Real molefrac[NOC];
    Real liquidmoleflow, vapormoleflow;
    Real liquidmolefrac[NOC], vapormolefrac[NOC];
    Real enthalpy;
    annotation(Diagram(graphics = {Ellipse(origin = {-1, 1}, fillColor = {65, 252, 255}, fillPattern = FillPattern.Solid, extent = {{-67, 65}, {67, -65}}, endAngle = 360)}), Icon(graphics = {Ellipse(origin = {-4, -1}, fillColor = {81, 253, 248}, fillPattern = FillPattern.Solid, extent = {{-60, 61}, {60, -61}}, endAngle = 360)}));
  end port;

  model valve
    parameter Real coeff(fixed = false) "Coeff for valve", valveCv = 0.4 "valve Cv if not control valve";
    parameter Boolean Control = false;
    parameter Boolean OutletPfixed = false;
    parameter Real OutletPressure = 1e5 "used only when OutletPfixed is true";
    Real flowrate, Cv, outletP, delP;
    unitoperations.sensor sensor1 annotation(Placement(visible = true, transformation(origin = {1, 77}, extent = {{-19, -19}, {19, 19}}, rotation = 0), iconTransformation(origin = {2, 56}, extent = {{-18, -18}, {18, 18}}, rotation = 0)));
    unitoperations.port port1 annotation(Placement(visible = true, transformation(origin = {-80, -4}, extent = {{-18, -18}, {18, 18}}, rotation = 0), iconTransformation(origin = {-82, -2}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
    unitoperations.port port2 annotation(Placement(visible = true, transformation(origin = {79, 1}, extent = {{-17, -17}, {17, 17}}, rotation = 0), iconTransformation(origin = {82, 2}, extent = {{-18, -18}, {18, 18}}, rotation = 0)));
  equation
    if Control == false then
      sensor1.var = 0;
    end if;
    if Control == true then
      Cv = sensor1.var;
    else
      Cv = valveCv;
    end if;
    flowrate = coeff * Cv * delP ^ 0.5;
    port1.moleflow = flowrate;
    port1.moleflow = port2.moleflow;
    if OutletPfixed == true then
      outletP = OutletPressure;
      port2.pressure = outletP;
    end if;
    if OutletPfixed == false then
      outletP = 0;
    end if;
    delP = port1.pressure - port2.pressure;
    port2.temperature = port1.temperature;
    port2.molefrac[:] = port1.molefrac[:];
    port2.liquidmoleflow = port1.liquidmoleflow;
    port2.vapormoleflow = port1.vapormoleflow;
    port2.liquidmolefrac[:] = port1.liquidmolefrac[:];
    port2.vapormolefrac[:] = port1.vapormolefrac[:];
    port2.enthalpy = port1.enthalpy;
    
    annotation(Icon(graphics = {Polygon(origin = {-48.26, -2.99}, fillColor = {85, 170, 255}, fillPattern = FillPattern.Solid, points = {{-47.7381, 48.9887}, {-47.7381, -49.0113}, {48.2619, 2.98874}, {48.2619, 2.98874}, {-47.7381, 48.9887}}), Polygon(origin = {49.25, -4.98}, fillColor = {85, 170, 255}, fillPattern = FillPattern.Solid, points = {{-47.2509, 4.98071}, {46.7491, 48.9807}, {46.7491, -49.0193}, {-47.2509, 4.98071}}), Rectangle(origin = {1, 35}, fillColor = {85, 170, 255}, fillPattern = FillPattern.Solid, extent = {{-15, 35}, {15, -35}})}, coordinateSystem(initialScale = 0.1)));
  end valve;

  model MaterialStream
    extends compounds;
    parameter Real Flowrate = 100, Pressure = 1e5, Temperature = 300, molefraction[NOC] = zeros(NOC);
    parameter Boolean unspecified = true;
    parameter Boolean stepchange = false;
    parameter Real stepchangetime = 0.01;
    parameter Real step_value = 1;
    Real kf[NOC], zl[NOC](each min = 0, each max = 1, start = {0.5, 1e-18, 0.5, 0}), zv[NOC](each min = 0, each max = 1, start = {0, 0.25, 0, 0.75}), Fl(min = 0, start = 100), Fv(min = 0, start = 140), Tbf(start = 62), Tdf(start = 495.5), Psat_Tdf[NOC], Psat_Tbf[NOC], Psat_Tf[NOC], Pf, Tf, z[NOC], F;
    Real Hvf[NOC], Hlf[NOC], H;
    unitoperations.port port2 annotation(Placement(visible = true, transformation(origin = {80, -4}, extent = {{-18, -18}, {18, 18}}, rotation = 0), iconTransformation(origin = {85, 1}, extent = {{-21, -21}, {21, 21}}, rotation = 0)));
    unitoperations.port port1 annotation(Placement(visible = true, transformation(origin = {-80, -4}, extent = {{-18, -18}, {18, 18}}, rotation = 0), iconTransformation(origin = {-82, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  equation
    if unspecified == false then
      if stepchange == true then
        if time < stepchangetime then
          port1.moleflow = Flowrate;
        else
          port1.moleflow = Flowrate + step_value;
        end if;
      else
        port1.moleflow = Flowrate;
      end if;
      port1.pressure = Pressure;
      port1.temperature = Temperature;
      port1.molefrac[:] = molefraction[:];
      port1.liquidmolefrac[:] = zl[:];
      port1.vapormolefrac[:] = zv[:];
      port1.liquidmoleflow = Fl;
      port1.vapormoleflow = Fv;
      port1.enthalpy = H;
    end if;
    if unspecified == true then
      port1.liquidmolefrac[:] = zl[:];
      port1.vapormolefrac[:] = zv[:];
      port1.liquidmoleflow = Fl;
      port1.vapormoleflow = Fv;
      port1.enthalpy = H;
    end if;
    if unspecified == false then
      Pf = Pressure;
      z[:] = molefraction[:];
      Tf = Temperature;
      F = port1.moleflow;
    else
      Pf = port1.pressure;
      z[:] = port1.molefrac[:];
      Tf = port1.temperature;
      F = port1.moleflow;
    end if;
//flash calculations
    sum(Pf * z[:] ./ Psat_Tdf[:]) = 1;
    sum(Psat_Tbf[:] ./ Pf .* z[:]) = 1;
    if Tf < Tbf then
      zl[:] = z[:];
      zv = zeros(NOC);
      Fl = F;
      Fv = 0;
      kf = zeros(NOC);
    elseif Tf > Tdf then
      zv[:] = z[:];
      zl = zeros(NOC);
      Fl = 0;
      Fv = F;
      kf = zeros(NOC);
    else
      sum(zl[:]) = 1;
      sum(zv[:]) = 1;
      zv[:] = kf[:] .* zl[:];
      kf[:] = Psat_Tf[:] ./ Pf;
      F = Fl + Fv;
      for i in 1:NOC - 1 loop
        F * z[i] = Fl * zl[i] + Fv * zv[i];
      end for;
    end if;
    H = Fl * sum(zl[:] .* Hlf[:]) + Fv * sum(zv[:] .* Hvf[:]);
    for i in 1:NOC loop
      Psat_Tbf[i] = Functions.Psat(comp[i].VP, Tbf);
      Psat_Tdf[i] = Functions.Psat(comp[i].VP, Tdf);
      Psat_Tf[i] = Functions.Psat(comp[i].VP, Tf);
      Hlf[i] = Functions.HLiqId(comp[i].VapCp, comp[i].HOV, comp[i].Tc, Tf);
      Hvf[i] = Functions.HVapId(comp[i].VapCp, comp[i].HOV, comp[i].Tc, Tf);
    end for;
    port2.moleflow = F;
    port2.pressure = port1.pressure;
    port2.temperature = port1.temperature;
    port2.molefrac[:] = port1.molefrac[:];
    port2.liquidmolefrac[:] = zl[:];
    port2.vapormolefrac[:] = zv[:];
    port2.liquidmoleflow = Fl;
    port2.vapormoleflow = Fv;
    port2.enthalpy = H;
    annotation(Icon(graphics = {Rectangle(origin = {-10, 3}, fillColor = {255, 0, 0}, fillPattern = FillPattern.Solid, extent = {{-90, 45}, {90, -45}}), Text(origin = {-17, -1}, extent = {{-51, 27}, {51, -27}}, textString = "MS1")}, coordinateSystem(initialScale = 0.1)));
  end MaterialStream;

  model materialtest
    MaterialStream materialStream1(Flowrate = 100, Name = "MS2", Pressure = 5e5, Temperature = 350, molefraction = {0.5, 0, 0.5, 0}, unspecified = false) annotation(Placement(visible = true, transformation(origin = {-24, 22}, extent = {{-24, -24}, {24, 24}}, rotation = 0)));
  end materialtest;

  model CSTR
    extends unitoperations.compounds;
    parameter Real s[NOC] = {-1, -1, 1, 1} "stoichiometric coefficients";
    parameter Integer NOIS = 2 "No of input streams";
    parameter Integer n = 1 "base compound identity in comp[n]";
    parameter Real V_Total = 2.321 "Volume of reactor", P_init = 25e5 "Pressure at t = 0";
    constant Real R = 8.314;
    parameter Real Af "frequency factor for forward reaction" annotation(Dialog(tab = "Reactions", group = "forward reaction rate constants"));
    parameter Real order_f[4] "order wrt to components for forward reaction" annotation(Dialog(tab = "Reactions", group = "forward reaction rate constants"));
    parameter Real order_b[4] "order wrt to components for backward reaction" annotation(Dialog(tab = "Reactions", group = "backward reaction rate constants"));
    parameter Real Ab "frequency factor for backward reaction" annotation(Dialog(tab = "Reactions", group = "backward reaction rate constants"));
    parameter Real Eaf "Activation energy for forward reaction" annotation(Dialog(tab = "Reactions", group = "forward reaction rate constants"));
    parameter Real Eab "Activation energy for backward reaction" annotation(Dialog(tab = "Reactions", group = "backward reaction rate constants"));
    parameter Real delH_r = 12.6e3 "Heat of reaction" annotation(Dialog(tab = "Reactions", group = "Reaction rate constants"));
    type operation_type = enumeration(Isothermal, Adiabatic);
    parameter operation_type operation_mode; 
    parameter Real  T_iso = 900;
    Real F_in[NOIS], z[NOIS, NOC], Hin[NOIS] "Feed values";
    Real M_Total(start = 1.5), M[NOC], x[NOC], F_out, densityi[NOC], P;
    Real r, kf, kb, c[NOC];
    Real H_r, Hout, Q, T;
    unitoperations.port port1 annotation(Placement(visible = true, transformation(origin = {-82, 2}, extent = {{-18, -18}, {18, 18}}, rotation = 0), iconTransformation(origin = {-85, 1}, extent = {{-17, -17}, {17, 17}}, rotation = 0)));
    unitoperations.port port2 annotation(Placement(visible = true, transformation(origin = {2, 84}, extent = {{-18, -18}, {18, 18}}, rotation = 0), iconTransformation(origin = {3, 83}, extent = {{-17, -17}, {17, 17}}, rotation = 0)));
  unitoperations.port port3 annotation(Placement(visible = true, transformation(origin = {89, -75}, extent = {{-17, -17}, {17, 17}}, rotation = 0), iconTransformation(origin = {82, -68}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  initial equation
  //calculates steady state solution at t=0
    P = P_init;
    der(M_Total) = 0;
    for i in 1:NOC - 1 loop
      der(M[i]) = 0;
    end for;
  equation
    F_in[1] = port1.moleflow;
    F_in[2] = port2.moleflow;
    z[1, :] = port1.molefrac[:];
    z[2, :] = port2.molefrac[:];
    Hin[1] = port1.enthalpy;
    Hin[2] = port2.enthalpy;
    for i in 1:NOC loop
      densityi[i] = Functions.Density(comp[i].LiqDen, comp[i].Tc, T, P);
    end for;
    for i in 1:NOC loop
      x[i] = M[i] / M_Total;
      //c[i] = M[i]*1000 / (V_Total*1000);
      c[i] = x[i] * P/(R*T);
    end for;
  //reaction rate
    kf = Af * exp(-Eaf / (R * T));
    kb = Ab * exp(-Eab / (R * T)); 
    r = kf * product(c[:].^order_f[:]) - kb * product(c[:].^order_b[:]);
  //unsteady state mass balance
    for i in 1:NOC - 1 loop
      der(M[i]) = sum(z[:, i] .* F_in[:]) - x[i] * F_out + s[i] * r * V_Total*1000 / abs(s[n]);
    end for;
    der(M_Total) = sum(F_in[:]) - F_out;
  //Pressure
    M_Total = sum(M[:]);
    P * V_Total = M_Total * R * T * 1000;
  //Energy balance
    if Integer(operation_mode) == 1 then
      T = T_iso;
    else
      Q = 0;
    end if;   
    H_r = r * V_Total * delH_r;
    sum(Hin[:]) + Q + H_r = Hout;
    Hout = port3.enthalpy;
    port3.moleflow = F_out;
    port3.temperature = T;
    port3.pressure = P;
    port3.molefrac[:] = x[:];
    annotation(Icon(graphics = {Rectangle(origin = {-1, 1}, fillColor = {255, 170, 0}, fillPattern = FillPattern.Solid, extent = {{-97, 95}, {97, -95}}), Text(origin = {2, 0}, extent = {{-80, 62}, {80, -62}}, textString = "CSTR")}, coordinateSystem(initialScale = 0.1)));
  end CSTR;

  model CSTR_DynamicTest
    MaterialStream materialStream1(Flowrate = 100, Pressure = 24e5, Temperature = 600, molefraction = {0.9, 0, 0.1, 0}, step_value = 20, stepchange = false, stepchangetime = 0.01, unspecified = false) annotation(Placement(visible = true, transformation(origin = {-87, -1}, extent = {{-19, -19}, {19, 19}}, rotation = 0)));
    MaterialStream materialStream2(Flowrate = 100, Pressure = 24e5, Temperature = 350, molefraction = {0, 0.9, 0, 0.1}, unspecified = false) annotation(Placement(visible = true, transformation(origin = {-45, 65}, extent = {{-17, -17}, {17, 17}}, rotation = 0)));
    valve valve1(Control = false, OutletPfixed = true, OutletPressure = 5e5, valveCv = 0.4) annotation(Placement(visible = true, transformation(origin = {56, -20}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
    MaterialStream materialStream3 annotation(Placement(visible = true, transformation(origin = {110, -22}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    CSTR cSTR1(Ab = 0, Af = 5.1e11, Eab = 0, Eaf = 230e3, V_Total = 4, operation_mode = modelicatest.CSTR.operation_type.Isothermal, order_b = {0, 0, 0, 0}, order_f = {1, 0.5, 0, 0})  annotation(Placement(visible = true, transformation(origin = {-14, -10}, extent = {{-22, -22}, {22, 22}}, rotation = 0)));
  equation
    connect(cSTR1.port3, valve1.port1) annotation(Line(points = {{4, -26}, {38, -26}, {38, -20}, {40, -20}}));
    connect(materialStream2.port2, cSTR1.port2) annotation(Line(points = {{-30, 66}, {-14, 66}, {-14, 8}, {-14, 8}}));
    connect(materialStream1.port2, cSTR1.port1) annotation(Line(points = {{-70, 0}, {-34, 0}, {-34, -10}, {-32, -10}}));
    connect(valve1.port2, materialStream3.port1) annotation(Line(points = {{72, -20}, {102, -20}, {102, -22}, {102, -22}}));
  end CSTR_DynamicTest;

  model Flash
    parameter Real R = 8.314, A = 1.5, hset = 3.7, Pset = 5e5, V_Total = 8.11;
    extends compounds;
    parameter Boolean connectedToInput = false;
    parameter Real Ti = 310;
    Real z[NOC], Tf;
    Real y[NOC], x[NOC](start = {0.7, 1e-18, 0.3, 0}), k[NOC], L(start = 100, min = 0), V(start = 140, min = 0), Psat_T[NOC], M[NOC], M_Total, ML(start = 50), MG(start = 0.5), VL, VG, Q, hv[NOC], hl[NOC], Hf, Hv, Hl, H_M_Total, F, densityi[NOC], P, h;
    unitoperations.sensor sensor1 annotation(Placement(visible = true, transformation(origin = {2, 82}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {8.88178e-16, 82}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
    unitoperations.sensor sensor3 annotation(Placement(visible = true, transformation(origin = {82, -32}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {77, -31}, extent = {{-19, -19}, {19, 19}}, rotation = 0)));
    unitoperations.port port1 annotation(Placement(visible = true, transformation(origin = {1, -83}, extent = {{-17, -17}, {17, 17}}, rotation = 0), iconTransformation(origin = {-8.88178e-16, -74}, extent = {{-18, -18}, {18, 18}}, rotation = 0)));
    unitoperations.port port2 annotation(Placement(visible = true, transformation(origin = {80, 48}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {76, 52}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
    unitoperations.port port3 annotation(Placement(visible = true, transformation(origin = {-84, 0}, extent = {{-18, -18}, {18, 18}}, rotation = 0), iconTransformation(origin = {-80, 4}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  initial equation
    h = hset;
    P = Pset;
    for i in 1:NOC - 1 loop
      der(M[i]) = 0;
    end for;
    der(M_Total) = 0;
//der(H_M_Total) = 0;
  equation
    F = port3.moleflow;
    z[:] = port3.molefrac[:];
    Tf = port3.temperature;
    if connectedToInput == false then
     port3.pressure = P; 
    end if;
    for i in 1:NOC loop
      Psat_T[i] = Functions.Psat(comp[i].VP, Ti);
      hv[i] = Functions.HVapId(comp[i].VapCp, comp[i].HOV, comp[i].Tc, Ti);
      hl[i] = Functions.HLiqId(comp[i].VapCp, comp[i].HOV, comp[i].Tc, Ti);
      densityi[i] = Functions.Density(comp[i].LiqDen, comp[i].Tc, Ti, P);
    end for;
/* Mass Blanace equations */
    der(M_Total) = F - L - V;
//F - L - V = 0;
    for i in 1:NOC - 1 loop
      der(M[i]) = F * z[i] - L * x[i] - V * y[i];
//F * z[i] - L * x[i] - V * y[i] = 0;
      M[i] = ML * x[i] + MG * y[i];
    end for;
    sum(M[:]) = M_Total;
    M_Total = MG + ML;
    VL = ML / sum(x[:] .* densityi[:]);
    VG = V_Total - VL;
    P * VG = MG * R * Ti * 1000;
//ideal gas law for gas phase
/*energy balance */
    Hv = V * sum(y[:] .* hv[:]);
    Hf = port3.enthalpy;
    Hl = L * sum(x[:] .* hl[:]);
    H_M_Total = ML * sum(x[:] .* hl[:]) + MG * sum(y[:] .* hv[:]);
//Hf - Hv - Hl = der(H_M_Total);
    Hf - Hv - Hl + Q = 0;
/*Thermodynamic equations */
    sum(x[:]) = 1;
    sum(y[:]) = 1;
    y[:] = k[:] .* x[:];
    k[:] = Psat_T[:] / P;
/* Control */
    A * h = VL;
//connector equations
    sensor1.var = P;
    sensor3.var = h;
    port1.moleflow = L;
    port1.pressure = P;
    port1.temperature = Ti;
    port1.molefrac[:] = x[:];
    port2.moleflow = V;
    port2.pressure = P;
    port2.temperature = Ti;
    port2.molefrac[:] = y[:];

    annotation(Icon(graphics = {Text(origin = {-30, 86}, extent = {{-22, 32}, {22, -32}}, textString = "Pressure"), Text(origin = {56, 46}, extent = {{-10, 24}, {2, -2}}, textString = "V"), Text(origin = {60, -33}, extent = {{-12, 25}, {4, -3}}, textString = "h"), Text(origin = {-16, -82}, extent = {{-14, 26}, {2, -6}}, textString = "L"), Text(origin = {0, 15}, extent = {{-46, 41}, {46, -41}}, textString = "PT flash"), Rectangle(origin = {1, 0}, extent = {{-87, 96}, {87, -96}}), Text(origin = {-64, -21}, extent = {{-10, 17}, {10, -17}}, textString = "F")}, coordinateSystem(initialScale = 0.1)));
  end Flash;

  model dynFlashtest
    MaterialStream materialStream1(Flowrate = 100, Pressure = 10e5, Temperature = 350, molefraction = {0.25, 0.25, 0.25, 0.25}, stepchange = true, unspecified = false) annotation(Placement(visible = true, transformation(origin = {-78, -2}, extent = {{-18, -18}, {18, 18}}, rotation = 0)));
    Flash flash1(connectedToInput = true)  annotation(Placement(visible = true, transformation(origin = {-5, -1}, extent = {{-29, -29}, {29, 29}}, rotation = 0)));
    MaterialStream materialStream2 annotation(Placement(visible = true, transformation(origin = {101, 27}, extent = {{-19, -19}, {19, 19}}, rotation = 0)));
    MaterialStream materialStream3 annotation(Placement(visible = true, transformation(origin = {84, -38}, extent = {{-24, -24}, {24, 24}}, rotation = 0)));
    valve valve1(OutletPfixed = true, OutletPressure = 1e5) annotation(Placement(visible = true, transformation(origin = {51, 25}, extent = {{-17, -17}, {17, 17}}, rotation = 0)));
    valve valve2(OutletPfixed = true, OutletPressure = 1e5) annotation(Placement(visible = true, transformation(origin = {26, -48}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(materialStream1.port2, flash1.port3) annotation(Line(points = {{-62, -2}, {-28, -2}, {-28, 0}, {-28, 0}}));
    connect(valve2.port2, materialStream3.port1) annotation(Line(points = {{34, -48}, {64, -48}, {64, -38}, {64, -38}}));
    connect(flash1.port1, valve2.port1) annotation(Line(points = {{-4, -22}, {-4, -22}, {-4, -48}, {18, -48}, {18, -48}}));
    connect(valve1.port2, materialStream2.port1) annotation(Line(points = {{64, 26}, {84, 26}, {84, 28}, {86, 28}}));
    connect(flash1.port2, valve1.port1) annotation(Line(points = {{18, 14}, {36, 14}, {36, 24}, {38, 24}}));
  end dynFlashtest;

  model HeatExchanger
    parameter Integer NOC = 4;
    parameter Chemsep_Database.Toluene comp1;
    parameter Chemsep_Database.Hydrogen comp2;
    parameter Chemsep_Database.Benzene comp3;
    parameter Chemsep_Database.Methane comp4;
    parameter Chemsep_Database.General_Properties comp[NOC] = {comp1, comp2, comp3, comp4};
    parameter Real Tout = 38 + 273;
    parameter Real pressure_drop = 0.5e5;
    Real Hin, Hout, Q "heat exchanger duty";
    port port1 annotation(Placement(visible = true, transformation(origin = {-86, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-86, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    port port2 annotation(Placement(visible = true, transformation(origin = {88, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {88, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    Hin = port1.enthalpy;
    Hout = port2.enthalpy;
    Tout = port2.temperature;
    port2.pressure = port1.pressure - pressure_drop;
    port2.moleflow = port1.moleflow;
    port2.molefrac[:] = port1.molefrac[:];
    Hin = Hout + Q;
    annotation(Icon(graphics = {Rectangle(origin = {1, 1}, fillColor = {255, 255, 0}, fillPattern = FillPattern.HorizontalCylinder, extent = {{-89, 47}, {89, -47}}), Text(origin = {6, 0}, extent = {{-58, 30}, {58, -30}}, textString = "HX")}));
  end HeatExchanger;

  model reactor_Flash_test
    unitoperations.MaterialStream materialStream1(Flowrate = 100, Pressure = 25e5, Temperature = 600, molefraction = {0.9, 0, 0.1, 0}, stepchange = true, unspecified = false) annotation(Placement(visible = true, transformation(origin = {-445, -30}, extent = {{-40, -40}, {40, 40}}, rotation = 0)));
    unitoperations.CSTR cSTR1(Operation = true) annotation(Placement(visible = true, transformation(origin = {-310, -30}, extent = {{-40, -40}, {40, 40}}, rotation = 0)));
    unitoperations.MaterialStream materialStream2(Flowrate = 100, Pressure = 25e5, Temperature = 400, molefraction = {0, 0.9, 0, 0.1}, stepchange = false, unspecified = false) annotation(Placement(visible = true, transformation(origin = {-447.5, 77.5}, extent = {{-47.5, -47.5}, {47.5, 47.5}}, rotation = 0)));
    unitoperations.valve valve1(OutletPfixed = false) annotation(Placement(visible = true, transformation(origin = {-200, -60}, extent = {{-30, -30}, {30, 30}}, rotation = 0)));
    unitoperations.MaterialStream materialStream3 annotation(Placement(visible = true, transformation(origin = {-90, -60}, extent = {{-30, -30}, {30, 30}}, rotation = 0)));
    unitoperations.Flash flash1 annotation(Placement(visible = true, transformation(origin = {132.5, -62.5}, extent = {{-52.5, -52.5}, {52.5, 52.5}}, rotation = 0)));
    unitoperations.valve valve2(OutletPfixed = true, OutletPressure = 1e5) annotation(Placement(visible = true, transformation(origin = {262.5, 47.5}, extent = {{-37.5, -37.5}, {37.5, 37.5}}, rotation = 0)));
    unitoperations.valve valve3(OutletPfixed = true, OutletPressure = 1e5) annotation(Placement(visible = true, transformation(origin = {215, -155}, extent = {{-35, -35}, {35, 35}}, rotation = 0)));
    unitoperations.MaterialStream materialStream4 annotation(Placement(visible = true, transformation(origin = {367.5, 47.5}, extent = {{-47.5, -47.5}, {47.5, 47.5}}, rotation = 0)));
    unitoperations.MaterialStream materialStream5 annotation(Placement(visible = true, transformation(origin = {375, -155}, extent = {{-30, -30}, {30, 30}}, rotation = 0)));
    HeatExchanger heatExchanger1 annotation(Placement(visible = true, transformation(origin = {-22.5, -57.5}, extent = {{-22.5, -22.5}, {22.5, 22.5}}, rotation = 0)));
    MaterialStream materialStream6(zl(start = {0.442, 7e-17, 0.5547, 0.0033}), zv(start = {0.0062, 0.424, 0.0245, 0.5453})) annotation(Placement(visible = true, transformation(origin = {45, -60}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  equation
    connect(valve3.port2, materialStream5.port1) annotation(Line(points = {{244, -154}, {278.5, -154}, {278.5, -155}, {350, -155}}));
    connect(materialStream6.port2, flash1.port3) annotation(Line(points = {{60, -60}, {95, -60}, {95, -60}, {90, -60}}));
    connect(heatExchanger1.port2, materialStream6.port1) annotation(Line(points = {{-5, -55}, {25, -55}, {25, -60}, {30, -60}}));
    connect(flash1.port1, valve3.port1) annotation(Line(points = {{132.5, -101}, {111.75, -101}, {111.75, -156}, {186, -156}}));
    connect(flash1.port2, valve2.port1) annotation(Line(points = {{172, -35}, {178.5, -35}, {178.5, 47}, {232, 47}}));
    connect(materialStream3.port2, heatExchanger1.port1) annotation(Line(points = {{-65, -60}, {-42, -60}, {-42, -57}}));
    connect(valve1.port2, materialStream3.port1) annotation(Line(points = {{-175, -59}, {-178, -59}, {-178, -60}, {-115, -60}}));
    connect(cSTR1.port3, valve1.port1) annotation(Line(points = {{-279, -60}, {-252.5, -60}, {-252.5, -61}, {-225, -61}}));
    connect(materialStream2.port2, cSTR1.port2) annotation(Line(points = {{-407, 78}, {-309, 78}, {-309, 3}}));
    connect(materialStream1.port2, cSTR1.port1) annotation(Line(points = {{-411, -30}, {-344, -30}}));
    connect(valve2.port2, materialStream4.port1) annotation(Line(points = {{293, 48}, {243.5, 48}, {243.5, 47.5}, {329, 47.5}}));
    annotation(Diagram(coordinateSystem(extent = {{-500, -200}, {500, 200}}, grid = {5, 5})), Icon(coordinateSystem(extent = {{-500, -200}, {500, 200}}, grid = {5, 5})), version = "", uses);
  end reactor_Flash_test;

  model CompoundSeperator
    port port1 annotation(Placement(visible = true, transformation(origin = {-88, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-88, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    port port2 annotation(Placement(visible = true, transformation(origin = {90, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {90, 70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    port port3 annotation(Placement(visible = true, transformation(origin = {92, -78}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {92, -78}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    parameter Real P = 3e5;
    Real M;
  initial equation
    der(M) = 0;
  equation
    der(M) = port1.moleflow - port2.moleflow - port3.moleflow;
    (port1.molefrac[1] + port1.molefrac[3]) * port1.moleflow = port2.moleflow;
    port2.molefrac[1] = port1.molefrac[1] * port1.moleflow / port2.moleflow;
    port2.molefrac[3] = 1 - port2.molefrac[1];
    port2.molefrac[2] = 0;
    port2.molefrac[4] = 0;
    port3.molefrac[2] = port1.molefrac[2] * port1.moleflow / port3.moleflow;
    port3.molefrac[4] = 1 - port3.molefrac[2];
    port3.molefrac[1] = 0;
    port3.molefrac[3] = 0;
    port1.temperature = port2.temperature;
    port1.temperature = port3.temperature;
    port1.pressure = port2.pressure;
    port1.pressure = port3.pressure;
//port1.pressure = P;
    annotation(Icon(graphics = {Polygon(origin = {3.29, -3.99}, fillColor = {255, 170, 0}, fillPattern = FillPattern.CrossDiag, points = {{-91.2946, 7.98643}, {88.7054, 73.9864}, {90.7054, -74.0136}, {-91.2946, 7.98643}})}));
  end CompoundSeperator;

  model Distillation
    parameter Chemsep_Database.Benzene comp1;
    parameter Chemsep_Database.Toluene comp2;
    parameter Chemsep_Database.General_Properties comp[2] = {comp1, comp2};
    parameter Integer N_Trays = 20, NOC = 2, N_Feed = 10;
    parameter Real M_eff[N_Trays, NOC] = fill(0.99, N_Trays, NOC);
    parameter Real Pressure_drop = 2432, P_condenser = 101350, R = 8.314;
    parameter Real A_active = 1, h_weir = 0.01, d_weir = 0.8, taul = 0.062, Tray_volume = 0.1;
    parameter Real pi1 = -12.55, pi2 = 0.91;
    Real y[N_Trays, NOC](start = fill(0.5, N_Trays, NOC)), x[N_Trays, NOC](start = fill(0.5, N_Trays, NOC)), y_eq[N_Trays, NOC], Tf[N_Trays];
    Real V[N_Trays](each start = 70), L[N_Trays](each start = 100);
    Real T[N_Trays](start = linspace(381, 369, N_Trays)), TC(start = 368), TB(start = 377), L0(start = 50), VNT, D, B, xc[NOC], xr[NOC], QC, QB;
    Real Keq[N_Trays, NOC], Psat[N_Trays, NOC], PsatC[NOC], PsatB[NOC];
    Real yNT[NOC], M[N_Trays](each start = 1), den[N_Trays, NOC];
    Real hv[N_Trays, NOC], hl[N_Trays, NOC], hf[N_Trays, NOC], hv_B[NOC], hl_B[NOC], hl_C[NOC];
    Real P[N_Trays], Ks[N_Trays];
    Real dx[N_Trays, NOC], dM[N_Trays], dhl[N_Trays, NOC];
    Real F[N_Trays](each start = 0), z[N_Trays, NOC](start = fill(0.5, N_Trays, NOC));
    port port1 annotation(Placement(visible = true, transformation(origin = {-88, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-88, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    port port2 annotation(Placement(visible = true, transformation(origin = {92, 80}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {92, 80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    port port3 annotation(Placement(visible = true, transformation(origin = {94, -74}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {94, -74}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  initial equation
    for i in 1:N_Trays loop
//M[i] = M0[i];
      der(M[i]) = 0;
//der(x[i,1]) = 0;
//dM[i] = 0;
      dhl[i, 1] = 0;
//dx[i,1] = 0;
    end for;
  equation
    port1.pressure = P[N_Feed];
    for i in 1:N_Trays loop
      if i == N_Feed then
        F[i] = port1.moleflow;
      else
        F[i] = 0;
      end if;
      Tf[i] = port1.temperature;
      z[i, 1] = port1.molefrac[3];
      z[i, 2] = port1.molefrac[1];
    end for;
//functions required
    for i in 1:NOC loop
      for j in 1:N_Trays loop
        Psat[j, i] = Functions.Psat(comp[i].VP, T[j]);
        hv[j, i] = Functions.HVapId(comp[i].VapCp, comp[i].HOV, comp[i].Tc, T[j]);
        hl[j, i] = Functions.HLiqId(comp[i].VapCp, comp[i].HOV, comp[i].Tc, T[j]);
        hf[j, i] = Functions.HLiqId(comp[i].VapCp, comp[i].HOV, comp[i].Tc, Tf[j]);
        den[j, i] = Functions.Density(comp[i].LiqDen, comp[i].Tc, T[j], P[j]);
      end for;
      hv_B[i] = Functions.HVapId(comp[i].VapCp, comp[i].HOV, comp[i].Tc, TB);
      hl_B[i] = Functions.HLiqId(comp[i].VapCp, comp[i].HOV, comp[i].Tc, TB);
      hl_C[i] = Functions.HLiqId(comp[i].VapCp, comp[i].HOV, comp[i].Tc, TC);
      PsatC[i] = Functions.Psat(comp[i].VP, TC);
      PsatB[i] = Functions.Psat(comp[i].VP, TB);
    end for;
    for i in 1:N_Trays loop
      Ks[i] = V[i] * R * T[i] / (A_active * P[i] * 10 ^ (-3)) * (P[i] * 10 ^ (-3) / (R * T[i] * (sum(x[i, :] .* den[i, :]) - P[i] * 10 ^ (-3) / (R * T[i]))));
    end for;
//defining state variables!
    for i in 1:N_Trays loop
      dx[i, :] = der(x[i, :]);
      dM[i] = der(M[i]);
      dhl[i, :] = der(hl[i, :]);
    end for;
//tray mass balance
    xr[:] = x[1, :];
    yNT[:] = xr[:];
    L[1] - VNT - B = 0;
//when time>0 then
    for i in 1:N_Trays loop
      M[i] = A_active * sum(x[i, :] .* den[i, :]) * (exp(pi1 * Ks[i] ^ pi2) * h_weir + 44300 * 10 ^ (-3) * (L[i] / (sum(x[i, :] .* den[i, :]) * d_weir * 1000 * 3600)) ^ 0.704);
    end for;
//end when;
    M[1] * dx[1, :] + x[1, :] * dM[1] = VNT .* yNT[:] + L[2] .* x[2, :] - V[1] .* y[1, :] - L[1] .* x[1, :] + F[1] .* z[1, :];
//M0[1] = Tray_volume * sum(x[1,:].*den[1,:]);
//L[1] = L0 + (M[1]-M0[1])/taul;
//M[1] = A_active * sum(x[1,:].*den[1,:])*(exp(pi1*Ks[1]^pi2)*h_weir + 44300*(L[1]/(sum(x[1,:].*den[1,:]) * 0.5* d_weir))^0.704);
    for i in 2:N_Trays - 1 loop
      dM[i] * x[i, :] + dx[i, :] * M[i] = V[i - 1] .* y[i - 1, :] + L[i + 1] .* x[i + 1, :] - V[i] .* y[i, :] - L[i] .* x[i, :] + F[i] .* z[i, :];
//L[i] = L0 + (M[i]-M0[i])/taul;
//M0[i] = Tray_volume * sum(x[i,:].*den[i,:]);
    end for;
    M[N_Trays] * dx[N_Trays, :] + x[N_Trays, :] * dM[N_Trays] = V[N_Trays - 1] .* y[N_Trays - 1, :] + L0 .* xc[:] - V[N_Trays] .* y[N_Trays, :] - L[N_Trays] .* x[N_Trays, :] + F[N_Trays] .* z[N_Trays, :];
//M0[N_Trays] = Tray_volume * sum(x[N_Trays,:].*den[N_Trays,:]);
//L[N_Trays] = L0 + (M[N_Trays]-M0[N_Trays])/taul;
//M[N_Trays] = A_active * sum(x[N_Trays,:].*den[N_Trays,:])*(exp(pi1*Ks[N_Trays]^pi2)*h_weir + 44300*(L[N_Trays]/(sum(x[N_Trays,:].*den[N_Trays,:]) * 0.5* d_weir))^0.704);
    V[N_Trays] - L0 - D = 0;
    y[N_Trays, :] = xc[:];
//energy balance
    VNT * sum(yNT[:] .* hv_B[:]) - V[1] * sum(y[1, :] .* hv[1, :]) + L[2] * sum(x[2, :] .* hl[2, :]) - L[1] * sum(x[1, :] .* hl[1, :]) + F[1] * sum(z[1, :] .* hf[1, :]) = M[1] * x[1, :] * dhl[1, :] + M[1] * dx[1, :] * hl[1, :] + dM[1] * x[1, :] * hl[1, :];
    for i in 2:N_Trays - 1 loop
      V[i - 1] * sum(y[i - 1, :] .* hv[i - 1, :]) - V[i] * sum(y[i, :] .* hv[i, :]) + L[i + 1] * sum(x[i + 1, :] .* hl[i + 1, :]) - L[i] * sum(x[i, :] .* hl[i, :]) + F[i] * sum(z[i, :] .* hf[i, :]) = M[i] * x[i, :] * dhl[i, :] + M[i] * dx[i, :] * hl[i, :] + dM[i] * x[i, :] * hl[i, :];
    end for;
    V[N_Trays - 1] * sum(y[N_Trays - 1, :] .* hv[N_Trays - 1, :]) - V[N_Trays] * sum(y[N_Trays, :] .* hv[N_Trays, :]) + L0 * sum(xc[:] .* hl_C[:]) - L[N_Trays] * sum(x[N_Trays, :] .* hl[N_Trays, :]) + F[N_Trays] * sum(z[N_Trays, :] .* hf[N_Trays, :]) = dM[N_Trays] * x[N_Trays, :] * hl[N_Trays, :] + M[N_Trays] * dx[N_Trays, :] * hl[N_Trays, :] + M[N_Trays] * x[N_Trays, :] * dhl[N_Trays, :];
    V[N_Trays] * sum(y[N_Trays, :] .* hv[N_Trays, :]) - (L0 + D) * sum(xc[:] .* hl_C[:]) = QC;
    L[1] * sum(x[1, :] .* hl[1, :]) - B * sum(xr[:] .* hl_B[:]) - VNT * sum(xr[:] .* hv_B[:]) = QB;
//pressure
    for i in 1:N_Trays loop
      P[i] = P_condenser + (N_Trays - i + 1) * Pressure_drop;
    end for;
//Equilibrium
    for i in 1:N_Trays loop
      Keq[i, :] = Psat[i, :] ./ P[i];
      y_eq[i, :] = Keq[i, :] .* x[i, :];
      M_eff[i, :] = (y[i, :] - y[i - 1, :]) ./ (y_eq[i, :] - y[i - 1, :]);
      sum(x[i, :]) = 1;
      sum(y_eq[i, :]) = 1;
    end for;
    sum(y[N_Trays, :] .* PsatC[:] / P_condenser) = 1;
//sum(xc[:]) =1
    sum(x[1, :] .* (P[1] + Pressure_drop) ./ PsatB[:]) = 1;
//sum(yNT[:]) =1;
    D = 0.5 * L0;
    B = 61.1;
    port2.moleflow = D;
    port2.molefrac = {xc[2], 0, xc[1], 0};
    port2.temperature = TC;
    port2.pressure = P_condenser;
    port3.moleflow = B;
    port3.molefrac = {xr[2], 0, xr[1], 0};
    port3.temperature = TB;
    port3.pressure = P_condenser + Pressure_drop * (N_Trays + 1);
    annotation(Icon(graphics = {Rectangle(origin = {-2, -1}, fillColor = {0, 85, 127}, fillPattern = FillPattern.VerticalCylinder, extent = {{-94, 95}, {94, -95}})}));
  end Distillation;

  model reactor_Flash_Seperator
    unitoperations.MaterialStream materialStream1(Flowrate = 100, Pressure = 25e5, Temperature = 600, molefraction = {0.9, 0, 0.1, 0}, stepchange = true, unspecified = false) annotation(Placement(visible = true, transformation(origin = {-445, -30}, extent = {{-40, -40}, {40, 40}}, rotation = 0)));
    unitoperations.CSTR cSTR1(Operation = true) annotation(Placement(visible = true, transformation(origin = {-310, -30}, extent = {{-40, -40}, {40, 40}}, rotation = 0)));
    unitoperations.MaterialStream materialStream2(Flowrate = 140, Pressure = 25e5, Temperature = 400, molefraction = {0, 0.9, 0, 0.1}, stepchange = false, unspecified = false) annotation(Placement(visible = true, transformation(origin = {-447.5, 77.5}, extent = {{-47.5, -47.5}, {47.5, 47.5}}, rotation = 0)));
    unitoperations.valve valve1(Control = false, OutletPfixed = false, valveCv = 0.4) annotation(Placement(visible = true, transformation(origin = {-200, -60}, extent = {{-30, -30}, {30, 30}}, rotation = 0)));
    unitoperations.MaterialStream materialStream3 annotation(Placement(visible = true, transformation(origin = {-90, -60}, extent = {{-30, -30}, {30, 30}}, rotation = 0)));
    unitoperations.Flash flash1 annotation(Placement(visible = true, transformation(origin = {132.5, -62.5}, extent = {{-52.5, -52.5}, {52.5, 52.5}}, rotation = 0)));
    unitoperations.valve valve2(Control = false, OutletPfixed = true, OutletPressure = 1e5) annotation(Placement(visible = true, transformation(origin = {262.5, 47.5}, extent = {{-37.5, -37.5}, {37.5, 37.5}}, rotation = 0)));
    unitoperations.MaterialStream materialStream4 annotation(Placement(visible = true, transformation(origin = {367.5, 47.5}, extent = {{-47.5, -47.5}, {47.5, 47.5}}, rotation = 0)));
    unitoperations.MaterialStream materialStream5 annotation(Placement(visible = true, transformation(origin = {212.5, -157.5}, extent = {{-32.5, -32.5}, {32.5, 32.5}}, rotation = 0)));
    HeatExchanger heatExchanger1 annotation(Placement(visible = true, transformation(origin = {-22.5, -57.5}, extent = {{-22.5, -22.5}, {22.5, 22.5}}, rotation = 0)));
    MaterialStream materialStream6(zl(start = {0.7, 1e-16, 0.3, 0.001}), zv(start = {0.01, 0.75, 0.0145, 0.23})) annotation(Placement(visible = true, transformation(origin = {45, -60}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
    unitoperations.CompoundSeperator compoundSeperator1 annotation(Placement(visible = true, transformation(origin = {315, -160}, extent = {{-35, -35}, {35, 35}}, rotation = 0)));
    unitoperations.MaterialStream materialStream7 annotation(Placement(visible = true, transformation(origin = {502.5, -47.5}, extent = {{-27.5, -27.5}, {27.5, 27.5}}, rotation = 0)));
    unitoperations.MaterialStream materialStream8 annotation(Placement(visible = true, transformation(origin = {530, -180}, extent = {{-30, -30}, {30, 30}}, rotation = 0)));
    unitoperations.valve valve4(Control = false, OutletPfixed = true, OutletPressure = 1e5) annotation(Placement(visible = true, transformation(origin = {435, -135}, extent = {{-30, -30}, {30, 30}}, rotation = 0)));
    unitoperations.valve valve5(Control = false, OutletPfixed = true, OutletPressure = 1e5) annotation(Placement(visible = true, transformation(origin = {437.5, -182.5}, extent = {{-27.5, -27.5}, {27.5, 27.5}}, rotation = 0)));
  equation
    connect(valve5.port2, materialStream8.port1) annotation(Line(points = {{460, -182}, {500, -182}, {500, -180}, {505, -180}}));
    connect(compoundSeperator1.port3, valve5.port1) annotation(Line(points = {{347, -187}, {360.875, -187}, {360.875, -182}, {387.5, -182}, {387.5, -183}, {415, -183}}));
    connect(valve4.port2, materialStream7.port1) annotation(Line(points = {{460, -134}, {460, -47.5}, {480, -47.5}}));
    connect(valve4.port1, compoundSeperator1.port2) annotation(Line(points = {{410, -136}, {410, -135.5}, {346.5, -135.5}}));
    connect(materialStream5.port2, compoundSeperator1.port1) annotation(Line(points = {{240, -157}, {270.25, -157}, {270.25, -159}, {284, -159}}));
    connect(flash1.port1, materialStream5.port1) annotation(Line(points = {{135, -100}, {135, -157.5}, {186, -157.5}}));
    connect(materialStream6.port2, flash1.port3) annotation(Line(points = {{60, -60}, {95, -60}, {95, -60}, {90, -60}}));
    connect(heatExchanger1.port2, materialStream6.port1) annotation(Line(points = {{-5, -55}, {25, -55}, {25, -60}, {30, -60}}));
    connect(flash1.port2, valve2.port1) annotation(Line(points = {{172, -35}, {178.5, -35}, {178.5, 47}, {232, 47}}));
    connect(materialStream3.port2, heatExchanger1.port1) annotation(Line(points = {{-65, -60}, {-42, -60}, {-42, -57}}));
    connect(valve1.port2, materialStream3.port1) annotation(Line(points = {{-175, -59}, {-178, -59}, {-178, -60}, {-115, -60}}));
    connect(cSTR1.port3, valve1.port1) annotation(Line(points = {{-279, -60}, {-252.5, -60}, {-252.5, -61}, {-225, -61}}));
    connect(materialStream2.port2, cSTR1.port2) annotation(Line(points = {{-407, 78}, {-309, 78}, {-309, 3}}));
    connect(materialStream1.port2, cSTR1.port1) annotation(Line(points = {{-411, -30}, {-344, -30}}));
    connect(valve2.port2, materialStream4.port1) annotation(Line(points = {{293, 48}, {243.5, 48}, {243.5, 47.5}, {329, 47.5}}));
    annotation(Diagram(coordinateSystem(extent = {{-500, -200}, {600, 200}}, grid = {5, 5})), Icon(coordinateSystem(extent = {{-500, -200}, {600, 200}}, grid = {5, 5})), version = "", uses);
  end reactor_Flash_Seperator;

  model separator_test
    MaterialStream materialStream1(Flowrate = 100, Pressure = 5e5, Temperature = 310, molefraction = {0.5, 0.001, 0.4, 0.099}, stepchange = false, unspecified = false) annotation(Placement(visible = true, transformation(origin = {-68, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    CompoundSeperator compoundSeperator1 annotation(Placement(visible = true, transformation(origin = {-26, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    MaterialStream materialStream2 annotation(Placement(visible = true, transformation(origin = {54, 24}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    MaterialStream materialStream3 annotation(Placement(visible = true, transformation(origin = {16, -18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    valve valve1(Control = false, OutletPfixed = true, OutletPressure = 2e5) annotation(Placement(visible = true, transformation(origin = {10, 24}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(valve1.port2, materialStream2.port1) annotation(Line(points = {{18, 24}, {44, 24}, {44, 24}, {46, 24}}));
    connect(compoundSeperator1.port2, valve1.port1) annotation(Line(points = {{-16, 8}, {-14, 8}, {-14, 22}, {2, 22}, {2, 24}}));
    connect(materialStream1.port2, compoundSeperator1.port1) annotation(Line(points = {{-60, 0}, {-34, 0}, {-34, 0}, {-34, 0}}));
    connect(compoundSeperator1.port3, materialStream3.port1) annotation(Line(points = {{-16, -8}, {-16, -8}, {-16, -18}, {8, -18}, {8, -18}}));
  end separator_test;

  model Distillationtest
    MaterialStream materialStream1(Flowrate = 98, Pressure = 1e5, Temperature = 365, molefraction = {0.45, 0, 0.55, 0}, unspecified = false) annotation(Placement(visible = true, transformation(origin = {-70, 4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  MaterialStream materialStream3 annotation(Placement(visible = true, transformation(origin = {54, -34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  MaterialStream materialStream2 annotation(Placement(visible = true, transformation(origin = {41, 27}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
  Distillation distillation1 annotation(Placement(visible = true, transformation(origin = {-14, 6}, extent = {{-22, -22}, {22, 22}}, rotation = 0)));
  equation
    connect(materialStream1.port2, distillation1.port1) annotation(Line(points = {{-62, 4}, {-34, 4}, {-34, 6}, {-34, 6}}));
    connect(distillation1.port2, materialStream2.port1) annotation(Line(points = {{6, 24}, {30, 24}, {30, 28}, {32, 28}}));
    connect(materialStream3.port1, distillation1.port3) annotation(Line(points = {{46, -34}, {8, -34}, {8, -10}, {6, -10}}));
  end Distillationtest;

  model flowsheet
    unitoperations.MaterialStream materialStream1(Flowrate = 100, Pressure = 25e5, Temperature = 600, molefraction = {0.9, 0, 0.1, 0}, step_value = 2, stepchange = true, stepchangetime = 0.01, unspecified = false) annotation(Placement(visible = true, transformation(origin = {-445, -30}, extent = {{-40, -40}, {40, 40}}, rotation = 0)));
    unitoperations.CSTR cSTR1(Operation = true) annotation(Placement(visible = true, transformation(origin = {-310, -30}, extent = {{-40, -40}, {40, 40}}, rotation = 0)));
    unitoperations.MaterialStream materialStream2(Flowrate = 140, Pressure = 25e5, Temperature = 400, molefraction = {0, 0.9, 0, 0.1}, stepchange = false, unspecified = false) annotation(Placement(visible = true, transformation(origin = {-447.5, 77.5}, extent = {{-47.5, -47.5}, {47.5, 47.5}}, rotation = 0)));
    unitoperations.valve valve1(Control = false, OutletPfixed = false, valveCv = 0.4) annotation(Placement(visible = true, transformation(origin = {-200, -60}, extent = {{-30, -30}, {30, 30}}, rotation = 0)));
    unitoperations.MaterialStream materialStream3 annotation(Placement(visible = true, transformation(origin = {-90, -60}, extent = {{-30, -30}, {30, 30}}, rotation = 0)));
    unitoperations.Flash flash1 annotation(Placement(visible = true, transformation(origin = {132.5, -62.5}, extent = {{-52.5, -52.5}, {52.5, 52.5}}, rotation = 0)));
    unitoperations.valve valve2(Control = false, OutletPfixed = true, OutletPressure = 1e5) annotation(Placement(visible = true, transformation(origin = {262.5, 47.5}, extent = {{-37.5, -37.5}, {37.5, 37.5}}, rotation = 0)));
    unitoperations.MaterialStream materialStream4 annotation(Placement(visible = true, transformation(origin = {367.5, 47.5}, extent = {{-47.5, -47.5}, {47.5, 47.5}}, rotation = 0)));
    unitoperations.MaterialStream materialStream5 annotation(Placement(visible = true, transformation(origin = {212.5, -157.5}, extent = {{-32.5, -32.5}, {32.5, 32.5}}, rotation = 0)));
    HeatExchanger heatExchanger1 annotation(Placement(visible = true, transformation(origin = {-22.5, -57.5}, extent = {{-22.5, -22.5}, {22.5, 22.5}}, rotation = 0)));
    MaterialStream materialStream6(zl(start = {0.7, 1e-16, 0.3, 0.001}), zv(start = {0.01, 0.75, 0.0145, 0.23})) annotation(Placement(visible = true, transformation(origin = {45, -60}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
    unitoperations.CompoundSeperator compoundSeperator1 annotation(Placement(visible = true, transformation(origin = {315, -160}, extent = {{-35, -35}, {35, 35}}, rotation = 0)));
    unitoperations.MaterialStream materialStream7 annotation(Placement(visible = true, transformation(origin = {502.5, -47.5}, extent = {{-27.5, -27.5}, {27.5, 27.5}}, rotation = 0)));
    unitoperations.MaterialStream materialStream8 annotation(Placement(visible = true, transformation(origin = {530, -180}, extent = {{-30, -30}, {30, 30}}, rotation = 0)));
    unitoperations.valve valve4(Control = false, OutletPfixed = false, OutletPressure = 1e5) annotation(Placement(visible = true, transformation(origin = {435, -135}, extent = {{-30, -30}, {30, 30}}, rotation = 0)));
    unitoperations.valve valve5(Control = false, OutletPfixed = true, OutletPressure = 1e5) annotation(Placement(visible = true, transformation(origin = {437.5, -182.5}, extent = {{-27.5, -27.5}, {27.5, 27.5}}, rotation = 0)));
  Distillation distillation1 annotation(Placement(visible = true, transformation(origin = {600, -40}, extent = {{-40, -40}, {40, 40}}, rotation = 0)));
  MaterialStream materialStream9 annotation(Placement(visible = true, transformation(origin = {712.5, 7.5}, extent = {{-27.5, -27.5}, {27.5, 27.5}}, rotation = 0)));
  MaterialStream materialStream10(Fl(start = 0), Fv(start = 50),Tbf(start = 392), Tdf(start = 395),zl(start = {0.9, 0, 0.1, 0}), zv(start = {0.86, 0, 0.14, 0}))  annotation(Placement(visible = true, transformation(origin = {715, -85}, extent = {{-30, -30}, {30, 30}}, rotation = 0)));
  equation
    connect(distillation1.port2, materialStream9.port1) annotation(Line(points = {{635, -10}, {650, -10}, {650, 10}, {690, 10}, {690, 5}}));
    connect(distillation1.port3, materialStream10.port1) annotation(Line(points = {{640, -70}, {660, -70}, {660, -85}, {690, -85}, {690, -85}}));
    connect(materialStream7.port2, distillation1.port1) annotation(Line(points = {{525, -45}, {565, -45}, {565, -40}, {565, -40}}));
    connect(valve5.port2, materialStream8.port1) annotation(Line(points = {{460, -182}, {500, -182}, {500, -180}, {505, -180}}));
    connect(compoundSeperator1.port3, valve5.port1) annotation(Line(points = {{347, -187}, {360.875, -187}, {360.875, -182}, {387.5, -182}, {387.5, -183}, {415, -183}}));
    connect(valve4.port2, materialStream7.port1) annotation(Line(points = {{460, -134}, {460, -47.5}, {480, -47.5}}));
    connect(valve4.port1, compoundSeperator1.port2) annotation(Line(points = {{410, -136}, {410, -135.5}, {346.5, -135.5}}));
    connect(materialStream5.port2, compoundSeperator1.port1) annotation(Line(points = {{240, -157}, {270.25, -157}, {270.25, -159}, {284, -159}}));
    connect(flash1.port1, materialStream5.port1) annotation(Line(points = {{135, -100}, {135, -157.5}, {186, -157.5}}));
    connect(materialStream6.port2, flash1.port3) annotation(Line(points = {{60, -60}, {95, -60}, {95, -60}, {90, -60}}));
    connect(heatExchanger1.port2, materialStream6.port1) annotation(Line(points = {{-5, -55}, {25, -55}, {25, -60}, {30, -60}}));
    connect(flash1.port2, valve2.port1) annotation(Line(points = {{172, -35}, {178.5, -35}, {178.5, 47}, {232, 47}}));
    connect(materialStream3.port2, heatExchanger1.port1) annotation(Line(points = {{-65, -60}, {-42, -60}, {-42, -57}}));
    connect(valve1.port2, materialStream3.port1) annotation(Line(points = {{-175, -59}, {-178, -59}, {-178, -60}, {-115, -60}}));
    connect(cSTR1.port3, valve1.port1) annotation(Line(points = {{-279, -60}, {-252.5, -60}, {-252.5, -61}, {-225, -61}}));
    connect(materialStream2.port2, cSTR1.port2) annotation(Line(points = {{-407, 78}, {-309, 78}, {-309, 3}}));
    connect(materialStream1.port2, cSTR1.port1) annotation(Line(points = {{-411, -30}, {-344, -30}}));
    connect(valve2.port2, materialStream4.port1) annotation(Line(points = {{293, 48}, {243.5, 48}, {243.5, 47.5}, {329, 47.5}}));
    annotation(Diagram(coordinateSystem(extent = {{-500, -200}, {600, 200}}, grid = {5, 5})), Icon(coordinateSystem(extent = {{-500, -200}, {600, 200}}, grid = {5, 5})), version = "", uses);
  end flowsheet;

model Distillation1
  parameter Chemsep_Database.Benzene comp1;
  parameter Chemsep_Database.Toluene comp2;
  parameter Chemsep_Database.General_Properties comp[2] = {comp1, comp2};
  parameter Integer N_Trays = 20, NOC = 2, N_Feed = 10;
  parameter Real M_eff[N_Trays, NOC] = fill(0.99, N_Trays, NOC);
  parameter Real Pressure_drop = 2432, P_condenser = 101350, R = 8.314;
  parameter Real A_active = 1, h_weir = 0.01, d_weir = 0.8, taul = 0.062, Tray_volume = 0.1;
  parameter Real pi1 = -12.55, pi2 = 0.91;
  Real y[N_Trays, NOC](start = fill(0.5, N_Trays, NOC)), x[N_Trays, NOC](start = fill(0.5, N_Trays, NOC)), y_eq[N_Trays, NOC], Tf[N_Trays];
  Real V[N_Trays](each start = 70), L[N_Trays](each start = 100);
  Real T[N_Trays](start = linspace(381, 369, N_Trays)), TC(start = 368), TB(start = 377), L0(start = 50), VNT, D, B, xc[NOC], xr[NOC], QC, QB;
  Real Keq[N_Trays, NOC], Psat[N_Trays, NOC], PsatC[NOC], PsatB[NOC];
  Real yNT[NOC], M[N_Trays](each start = 1), den[N_Trays, NOC];
  Real hv[N_Trays, NOC], hl[N_Trays, NOC], hf[N_Trays, NOC], hv_B[NOC], hl_B[NOC], hl_C[NOC];
  Real P[N_Trays], Ks[N_Trays];
  Real dx[N_Trays, NOC], dM[N_Trays], dhl[N_Trays, NOC];
  Real F[N_Trays](each start = 0), z[N_Trays, NOC](start = fill(0.5, N_Trays, NOC));
  port port1 annotation(Placement(visible = true, transformation(origin = {-88, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-88, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  port port2 annotation(Placement(visible = true, transformation(origin = {92, 80}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {92, 80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  port port3 annotation(Placement(visible = true, transformation(origin = {94, -74}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {94, -74}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
initial equation
  for i in 1:N_Trays loop
//M[i] = M0[i];
    der(M[i]) = 0;
//der(x[i,1]) = 0;
//dM[i] = 0;
    dhl[i, 1] = 0;
//dx[i,1] = 0;
  end for;
equation
  port1.pressure = P[N_Feed];
  for i in 1:N_Trays loop
    if i == N_Feed then
      F[i] = port1.moleflow;
    else
      F[i] = 0;
    end if;
    Tf[i] = port1.temperature;
    z[i, 1] = port1.molefrac[3];
    z[i, 2] = port1.molefrac[1];
  end for;
//functions required
  for i in 1:NOC loop
    for j in 1:N_Trays loop
      Psat[j, i] = Functions.Psat(comp[i].VP, T[j]);
      hv[j, i] = Functions.HVapId(comp[i].VapCp, comp[i].HOV, comp[i].Tc, T[j]);
      hl[j, i] = Functions.HLiqId(comp[i].VapCp, comp[i].HOV, comp[i].Tc, T[j]);
      hf[j, i] = Functions.HLiqId(comp[i].VapCp, comp[i].HOV, comp[i].Tc, Tf[j]);
      den[j, i] = Functions.Density(comp[i].LiqDen, comp[i].Tc, T[j], P[j]);
    end for;
    hv_B[i] = Functions.HVapId(comp[i].VapCp, comp[i].HOV, comp[i].Tc, TB);
    hl_B[i] = Functions.HLiqId(comp[i].VapCp, comp[i].HOV, comp[i].Tc, TB);
    hl_C[i] = Functions.HLiqId(comp[i].VapCp, comp[i].HOV, comp[i].Tc, TC);
    PsatC[i] = Functions.Psat(comp[i].VP, TC);
    PsatB[i] = Functions.Psat(comp[i].VP, TB);
  end for;
  for i in 1:N_Trays loop
    Ks[i] = V[i] * R * T[i] / (A_active * P[i] * 10 ^ (-3)) * (P[i] * 10 ^ (-3) / (R * T[i] * (sum(x[i, :] .* den[i, :]) - P[i] * 10 ^ (-3) / (R * T[i]))));
  end for;
//defining state variables!
  for i in 1:N_Trays loop
    dx[i, :] = der(x[i, :]);
    dM[i] = der(M[i]);
    dhl[i, :] = der(hl[i, :]);
  end for;
//tray mass balance
  xr[:] = x[1, :];
  yNT[:] = xr[:];
  L[1] - VNT - B = 0;
//when time>0 then
  for i in 1:N_Trays loop
    M[i] = A_active * sum(x[i, :] .* den[i, :]) * (exp(pi1 * Ks[i] ^ pi2) * h_weir + 44300 * 10 ^ (-3) * (L[i] / (sum(x[i, :] .* den[i, :]) * d_weir * 1000 * 3600)) ^ 0.704);
  end for;
//end when;
  M[1] * dx[1, :] + x[1, :] * dM[1] = VNT .* yNT[:] + L[2] .* x[2, :] - V[1] .* y[1, :] - L[1] .* x[1, :] + F[1] .* z[1, :];
//M0[1] = Tray_volume * sum(x[1,:].*den[1,:]);
//L[1] = L0 + (M[1]-M0[1])/taul;
//M[1] = A_active * sum(x[1,:].*den[1,:])*(exp(pi1*Ks[1]^pi2)*h_weir + 44300*(L[1]/(sum(x[1,:].*den[1,:]) * 0.5* d_weir))^0.704);
  for i in 2:N_Trays - 1 loop
    dM[i] * x[i, :] + dx[i, :] * M[i] = V[i - 1] .* y[i - 1, :] + L[i + 1] .* x[i + 1, :] - V[i] .* y[i, :] - L[i] .* x[i, :] + F[i] .* z[i, :];
//L[i] = L0 + (M[i]-M0[i])/taul;
//M0[i] = Tray_volume * sum(x[i,:].*den[i,:]);
  end for;
  M[N_Trays] * dx[N_Trays, :] + x[N_Trays, :] * dM[N_Trays] = V[N_Trays - 1] .* y[N_Trays - 1, :] + L0 .* xc[:] - V[N_Trays] .* y[N_Trays, :] - L[N_Trays] .* x[N_Trays, :] + F[N_Trays] .* z[N_Trays, :];
//M0[N_Trays] = Tray_volume * sum(x[N_Trays,:].*den[N_Trays,:]);
//L[N_Trays] = L0 + (M[N_Trays]-M0[N_Trays])/taul;
//M[N_Trays] = A_active * sum(x[N_Trays,:].*den[N_Trays,:])*(exp(pi1*Ks[N_Trays]^pi2)*h_weir + 44300*(L[N_Trays]/(sum(x[N_Trays,:].*den[N_Trays,:]) * 0.5* d_weir))^0.704);
  V[N_Trays] - L0 - D = 0;
  y[N_Trays, :] = xc[:];
//energy balance
  VNT * sum(yNT[:] .* hv_B[:]) - V[1] * sum(y[1, :] .* hv[1, :]) + L[2] * sum(x[2, :] .* hl[2, :]) - L[1] * sum(x[1, :] .* hl[1, :]) + F[1] * sum(z[1, :] .* hf[1, :]) = M[1] * x[1, :] * dhl[1, :] + M[1] * dx[1, :] * hl[1, :] + dM[1] * x[1, :] * hl[1, :];
  for i in 2:N_Trays - 1 loop
    V[i - 1] * sum(y[i - 1, :] .* hv[i - 1, :]) - V[i] * sum(y[i, :] .* hv[i, :]) + L[i + 1] * sum(x[i + 1, :] .* hl[i + 1, :]) - L[i] * sum(x[i, :] .* hl[i, :]) + F[i] * sum(z[i, :] .* hf[i, :]) = M[i] * x[i, :] * dhl[i, :] + M[i] * dx[i, :] * hl[i, :] + dM[i] * x[i, :] * hl[i, :];
  end for;
  V[N_Trays - 1] * sum(y[N_Trays - 1, :] .* hv[N_Trays - 1, :]) - V[N_Trays] * sum(y[N_Trays, :] .* hv[N_Trays, :]) + L0 * sum(xc[:] .* hl_C[:]) - L[N_Trays] * sum(x[N_Trays, :] .* hl[N_Trays, :]) + F[N_Trays] * sum(z[N_Trays, :] .* hf[N_Trays, :]) = dM[N_Trays] * x[N_Trays, :] * hl[N_Trays, :] + M[N_Trays] * dx[N_Trays, :] * hl[N_Trays, :] + M[N_Trays] * x[N_Trays, :] * dhl[N_Trays, :];
  V[N_Trays] * sum(y[N_Trays, :] .* hv[N_Trays, :]) - (L0 + D) * sum(xc[:] .* hl_C[:]) = QC;
  L[1] * sum(x[1, :] .* hl[1, :]) - B * sum(xr[:] .* hl_B[:]) - VNT * sum(xr[:] .* hv_B[:]) = QB;
//pressure
  for i in 1:N_Trays loop
    P[i] = P_condenser + (N_Trays - i + 1) * Pressure_drop;
  end for;
//Equilibrium
  for i in 1:N_Trays loop
    Keq[i, :] = Psat[i, :] ./ P[i];
    y_eq[i, :] = Keq[i, :] .* x[i, :];
    M_eff[i, :] = (y[i, :] - y[i - 1, :]) ./ (y_eq[i, :] - y[i - 1, :]);
    sum(x[i, :]) = 1;
    sum(y_eq[i, :]) = 1;
  end for;
  sum(y[N_Trays, :] .* PsatC[:] / P_condenser) = 1;
//sum(xc[:]) =1
  sum(x[1, :] .* (P[1] + Pressure_drop) ./ PsatB[:]) = 1;
//sum(yNT[:]) =1;
  D = 0.5 * L0;
  B = 61.1;
  port2.moleflow = D;
  port2.molefrac = {xc[2], 0, xc[1], 0};
  port2.temperature = TC;
  port2.pressure = P_condenser;
  port3.moleflow = B;
  port3.molefrac = {xr[2], 0, xr[1], 0};
  port3.temperature = TB;
  port3.pressure = P_condenser + Pressure_drop * (N_Trays + 1);
  /*
  Real liquidmoleflow, vapormoleflow;
  Real liquidmolefrac[NOC], vapormolefrac[NOC];
  Real enthalpy; */
  port2.liquidmoleflow = 0;
  port2.vapormoleflow = 0;
  port2.liquidmolefrac = zeros(4);
  port2.vapormolefrac = zeros(4);
  port2.enthalpy = 0;
  port3.liquidmoleflow = 0;
  port3.vapormoleflow = 0;
  port3.liquidmolefrac = zeros(4);
  port3.vapormolefrac = zeros(4);
  port3.enthalpy = 0;
  annotation(Icon(graphics = {Rectangle(origin = {-2, -1}, fillColor = {0, 85, 127}, fillPattern = FillPattern.VerticalCylinder, extent = {{-94, 95}, {94, -95}})}));
end Distillation1;

model flowsheet1
  unitoperations.MaterialStream materialStream1(Flowrate = 100, Pressure = 25e5, Temperature = 600, molefraction = {0.9, 0, 0.1, 0}, step_value = 2, stepchange = true, stepchangetime = 0.01, unspecified = false) annotation(Placement(visible = true, transformation(origin = {-445, -30}, extent = {{-40, -40}, {40, 40}}, rotation = 0)));
  unitoperations.CSTR cSTR1(Operation = true) annotation(Placement(visible = true, transformation(origin = {-310, -30}, extent = {{-40, -40}, {40, 40}}, rotation = 0)));
  unitoperations.MaterialStream materialStream2(Flowrate = 140, Pressure = 25e5, Temperature = 400, molefraction = {0, 0.9, 0, 0.1}, stepchange = false, unspecified = false) annotation(Placement(visible = true, transformation(origin = {-447.5, 77.5}, extent = {{-47.5, -47.5}, {47.5, 47.5}}, rotation = 0)));
  unitoperations.valve valve1(Control = false, OutletPfixed = false, valveCv = 0.4) annotation(Placement(visible = true, transformation(origin = {-200, -60}, extent = {{-30, -30}, {30, 30}}, rotation = 0)));
  unitoperations.MaterialStream materialStream3 annotation(Placement(visible = true, transformation(origin = {-90, -60}, extent = {{-30, -30}, {30, 30}}, rotation = 0)));
  unitoperations.Flash flash1 annotation(Placement(visible = true, transformation(origin = {132.5, -62.5}, extent = {{-52.5, -52.5}, {52.5, 52.5}}, rotation = 0)));
  unitoperations.valve valve2(Control = false, OutletPfixed = true, OutletPressure = 1e5) annotation(Placement(visible = true, transformation(origin = {262.5, 47.5}, extent = {{-37.5, -37.5}, {37.5, 37.5}}, rotation = 0)));
  unitoperations.MaterialStream materialStream4 annotation(Placement(visible = true, transformation(origin = {367.5, 47.5}, extent = {{-47.5, -47.5}, {47.5, 47.5}}, rotation = 0)));
  unitoperations.MaterialStream materialStream5 annotation(Placement(visible = true, transformation(origin = {212.5, -157.5}, extent = {{-32.5, -32.5}, {32.5, 32.5}}, rotation = 0)));
  HeatExchanger heatExchanger1 annotation(Placement(visible = true, transformation(origin = {-22.5, -57.5}, extent = {{-22.5, -22.5}, {22.5, 22.5}}, rotation = 0)));
  MaterialStream materialStream6(zl(start = {0.7, 1e-16, 0.3, 0.001}), zv(start = {0.01, 0.75, 0.0145, 0.23})) annotation(Placement(visible = true, transformation(origin = {45, -60}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  unitoperations.CompoundSeperator compoundSeperator1 annotation(Placement(visible = true, transformation(origin = {315, -160}, extent = {{-35, -35}, {35, 35}}, rotation = 0)));
  unitoperations.MaterialStream materialStream7 annotation(Placement(visible = true, transformation(origin = {502.5, -47.5}, extent = {{-27.5, -27.5}, {27.5, 27.5}}, rotation = 0)));
  unitoperations.MaterialStream materialStream8 annotation(Placement(visible = true, transformation(origin = {530, -180}, extent = {{-30, -30}, {30, 30}}, rotation = 0)));
  unitoperations.valve valve4(Control = false, OutletPfixed = false, OutletPressure = 1e5) annotation(Placement(visible = true, transformation(origin = {435, -135}, extent = {{-30, -30}, {30, 30}}, rotation = 0)));
  unitoperations.valve valve5(Control = false, OutletPfixed = true, OutletPressure = 1e5) annotation(Placement(visible = true, transformation(origin = {437.5, -182.5}, extent = {{-27.5, -27.5}, {27.5, 27.5}}, rotation = 0)));
  Distillation1 distillation11 annotation(Placement(visible = true, transformation(origin = {630, -55}, extent = {{-35, -35}, {35, 35}}, rotation = 0)));
  equation
    connect(materialStream7.port2, distillation11.port1) annotation(Line(points = {{525, -45}, {600, -45}, {600, -55}, {600, -55}}));
    connect(valve5.port2, materialStream8.port1) annotation(Line(points = {{460, -182}, {500, -182}, {500, -180}, {505, -180}}));
    connect(compoundSeperator1.port3, valve5.port1) annotation(Line(points = {{347, -187}, {360.875, -187}, {360.875, -182}, {387.5, -182}, {387.5, -183}, {415, -183}}));
    connect(valve4.port2, materialStream7.port1) annotation(Line(points = {{460, -134}, {460, -47.5}, {480, -47.5}}));
    connect(valve4.port1, compoundSeperator1.port2) annotation(Line(points = {{410, -136}, {410, -135.5}, {346.5, -135.5}}));
    connect(materialStream5.port2, compoundSeperator1.port1) annotation(Line(points = {{240, -157}, {270.25, -157}, {270.25, -159}, {284, -159}}));
    connect(flash1.port1, materialStream5.port1) annotation(Line(points = {{135, -100}, {135, -157.5}, {186, -157.5}}));
    connect(materialStream6.port2, flash1.port3) annotation(Line(points = {{60, -60}, {95, -60}, {95, -60}, {90, -60}}));
    connect(heatExchanger1.port2, materialStream6.port1) annotation(Line(points = {{-5, -55}, {25, -55}, {25, -60}, {30, -60}}));
    connect(flash1.port2, valve2.port1) annotation(Line(points = {{172, -35}, {178.5, -35}, {178.5, 47}, {232, 47}}));
    connect(materialStream3.port2, heatExchanger1.port1) annotation(Line(points = {{-65, -60}, {-42, -60}, {-42, -57}}));
    connect(valve1.port2, materialStream3.port1) annotation(Line(points = {{-175, -59}, {-178, -59}, {-178, -60}, {-115, -60}}));
    connect(cSTR1.port3, valve1.port1) annotation(Line(points = {{-279, -60}, {-252.5, -60}, {-252.5, -61}, {-225, -61}}));
    connect(materialStream2.port2, cSTR1.port2) annotation(Line(points = {{-407, 78}, {-309, 78}, {-309, 3}}));
    connect(materialStream1.port2, cSTR1.port1) annotation(Line(points = {{-411, -30}, {-344, -30}}));
    connect(valve2.port2, materialStream4.port1) annotation(Line(points = {{293, 48}, {243.5, 48}, {243.5, 47.5}, {329, 47.5}}));
    annotation(Diagram(coordinateSystem(extent = {{-500, -200}, {600, 200}}, grid = {5, 5})), Icon(coordinateSystem(extent = {{-500, -200}, {600, 200}}, grid = {5, 5})), version = "", uses);
end flowsheet1;

  model DistillationWithSizing
    parameter Chemsep_Database.Benzene comp1 annotation(Dialog(tab = "General", group = "compounds"));
    parameter Chemsep_Database.Toluene comp2 annotation(Dialog(tab = "General", group = "compounds"));
    parameter Chemsep_Database.General_Properties comp[2] = {comp1, comp2} annotation(Dialog(tab = "General", group = "compounds"));
    parameter Integer N_Trays = 20 "No. of trays without condensor and reboiler" annotation(Dialog(tab = "General", group = "Trays"));
    parameter Integer NOC = 2 "No. of compounds" annotation(Dialog(tab = "General", group = "compounds"));
    parameter Integer N_Feed = 10 "Feed tray location" annotation(Dialog(tab = "General", group = "Trays"));
    parameter Real M_eff[N_Trays, NOC] = fill(0.99, N_Trays, NOC) "Murfy's efficiency of trays" annotation(Dialog(tab = "General", group = "Efficiency"));
    parameter Real Pressure_drop = 2432 "Pressure Drop per tray, Pa" annotation(Dialog(tab = "General", group = "Pressure Profile"));
    parameter Real P_condenser = 101350 "Pressure in Condensor" annotation(Dialog(tab = "General", group = "Pressure Profile"));
    constant Real R = 8.314 "Gas Constant";
    parameter Real A_active = 1 "Active area of Tray, sq.m" annotation(Dialog(tab = "dynamic", group = "Tray data"));
    parameter Real h_weir = 0.01 "Height of weir, m" annotation(Dialog(tab = "dynamic", group = "Tray data"));
    parameter Real d_weir = 0.8 "Diameter of weir, m" annotation(Dialog(tab = "dynamic", group = "Tray data"));
    type spec1 = enumeration(CondensorHeatLoad, ProductMolarFlow, CompoundMolarFlow, CompoundFractionInStream, RefluxRatio, Temperature) "condensor";
    type spec2 = enumeration(ReboilerHeatLoad, ProductMolarFlow, CompoundMolarFlow, CompoundFractionInStream, BoilUpRatio, Temperature) "reboiler";
    parameter spec1 specification1;
    parameter Real specification1_value = 1;
    parameter spec2 specification2;
    parameter Real specification2_value = 1;
    parameter Real Tray_volume = 0.1 annotation(Dialog(tab = "dynamic", group = "Tray data"));
    parameter Real pi1 = -12.55 "constant for calculationg froth density, phi" annotation(Dialog(tab = "dynamic", group = "vapor flow"));
    parameter Real pi2 = 0.91 "constant for calculationg froth density, phi" annotation(Dialog(tab = "dynamic", group = "vapor flow"));
    Real y[N_Trays, NOC](start = fill(0.5, N_Trays, NOC)), x[N_Trays, NOC](start = fill(0.5, N_Trays, NOC)), y_eq[N_Trays, NOC], Tf[N_Trays];
    Real V[N_Trays](each start = 70), L[N_Trays](each start = 100);
    Real T[N_Trays](start = linspace(386, 354, N_Trays)), TC(start = 368), TB(start = 377), L0(start = 50), VNT, D, B, xc[NOC], xr[NOC], QC, QB;
    Real Keq[N_Trays, NOC], Psat[N_Trays, NOC], PsatC[NOC], PsatB[NOC];
    Real yNT[NOC], M[N_Trays](each start = 1), den[N_Trays, NOC];
    Real hv[N_Trays, NOC] annotation(each HideResult = true), hl[N_Trays, NOC] annotation(each HideResult = true), hf[N_Trays, NOC] annotation(each HideResult = true), hv_B[NOC] annotation(each HideResult = true), hl_B[NOC] annotation(each HideResult = true), hl_C[NOC] annotation(each HideResult = true);
    Real P[N_Trays], Ks[N_Trays] annotation(each HideResult = true);
    Real dx[N_Trays, NOC], dM[N_Trays], dhl[N_Trays, NOC];
    Real F[N_Trays](each start = 0) annotation(each HideResult = true), z[N_Trays, NOC](start = fill(0.5, N_Trays, NOC)) annotation(each HideResult = true);
    port port1 annotation(Placement(visible = true, transformation(origin = {-88, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-88, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    port port2 annotation(Placement(visible = true, transformation(origin = {92, 80}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {92, 80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  initial equation
    for i in 1:N_Trays loop
  //M[i] = M0[i];
      der(M[i]) = 0;
  //der(x[i,1]) = 0;
  //dM[i] = 0;
      dhl[i, 1] = 0;
  //dx[i,1] = 0;
    end for;
  equation
  //  port1.pressure = P[N_Feed]; "should be automatic"
    for i in 1:N_Trays loop
      if i == N_Feed then
        F[i] = port1.moleflow;
      else
        F[i] = 0;
      end if;
      Tf[i] = port1.temperature;
      z[i, 1] = port1.molefrac[3];
      z[i, 2] = port1.molefrac[1];
    end for;
  //functions required
    for i in 1:NOC loop
      for j in 1:N_Trays loop
        Psat[j, i] = Functions.Psat(comp[i].VP, T[j]);
        hv[j, i] = Functions.HVapId(comp[i].VapCp, comp[i].HOV, comp[i].Tc, T[j]);
        hl[j, i] = Functions.HLiqId(comp[i].VapCp, comp[i].HOV, comp[i].Tc, T[j]);
        hf[j, i] = Functions.HLiqId(comp[i].VapCp, comp[i].HOV, comp[i].Tc, Tf[j]);
        den[j, i] = Functions.Density(comp[i].LiqDen, comp[i].Tc, T[j], P[j]);
      end for;
      hv_B[i] = Functions.HVapId(comp[i].VapCp, comp[i].HOV, comp[i].Tc, TB);
      hl_B[i] = Functions.HLiqId(comp[i].VapCp, comp[i].HOV, comp[i].Tc, TB);
      hl_C[i] = Functions.HLiqId(comp[i].VapCp, comp[i].HOV, comp[i].Tc, TC);
      PsatC[i] = Functions.Psat(comp[i].VP, TC);
      PsatB[i] = Functions.Psat(comp[i].VP, TB);
    end for;
    for i in 1:N_Trays loop
      Ks[i] = V[i] * R * T[i] / (A_active * P[i] * 10 ^ (-3)) * (P[i] * 10 ^ (-3) / (R * T[i] * (sum(x[i, :] .* den[i, :]) - P[i] * 10 ^ (-3) / (R * T[i]))));
    end for;
  //defining state variables!
    for i in 1:N_Trays loop
      dx[i, :] = der(x[i, :]);
      dM[i] = der(M[i]);
      dhl[i, :] = der(hl[i, :]);
    end for;
  //tray mass balance
    xr[:] = x[1, :];
    yNT[:] = xr[:];
    L[1] - VNT - B = 0;
  //when time>0 then
    for i in 1:N_Trays loop
      M[i] = A_active * sum(x[i, :] .* den[i, :]) * (exp(pi1 * Ks[i] ^ pi2) * h_weir + 44300 * 10 ^ (-3) * (L[i] / (sum(x[i, :] .* den[i, :]) * d_weir * 1000 * 3600)) ^ 0.704);
    end for;
  //end when;
    M[1] * dx[1, :] + x[1, :] * dM[1] = VNT .* yNT[:] + L[2] .* x[2, :] - V[1] .* y[1, :] - L[1] .* x[1, :] + F[1] .* z[1, :];
  //M0[1] = Tray_volume * sum(x[1,:].*den[1,:]);
  //L[1] = L0 + (M[1]-M0[1])/taul;
  //M[1] = A_active * sum(x[1,:].*den[1,:])*(exp(pi1*Ks[1]^pi2)*h_weir + 44300*(L[1]/(sum(x[1,:].*den[1,:]) * 0.5* d_weir))^0.704);
    for i in 2:N_Trays - 1 loop
      dM[i] * x[i, :] + dx[i, :] * M[i] = V[i - 1] .* y[i - 1, :] + L[i + 1] .* x[i + 1, :] - V[i] .* y[i, :] - L[i] .* x[i, :] + F[i] .* z[i, :];
  //L[i] = L0 + (M[i]-M0[i])/taul;
  //M0[i] = Tray_volume * sum(x[i,:].*den[i,:]);
    end for;
    M[N_Trays] * dx[N_Trays, :] + x[N_Trays, :] * dM[N_Trays] = V[N_Trays - 1] .* y[N_Trays - 1, :] + L0 .* xc[:] - V[N_Trays] .* y[N_Trays, :] - L[N_Trays] .* x[N_Trays, :] + F[N_Trays] .* z[N_Trays, :];
  //M0[N_Trays] = Tray_volume * sum(x[N_Trays,:].*den[N_Trays,:]);
  //L[N_Trays] = L0 + (M[N_Trays]-M0[N_Trays])/taul;
  //M[N_Trays] = A_active * sum(x[N_Trays,:].*den[N_Trays,:])*(exp(pi1*Ks[N_Trays]^pi2)*h_weir + 44300*(L[N_Trays]/(sum(x[N_Trays,:].*den[N_Trays,:]) * 0.5* d_weir))^0.704);
    V[N_Trays] - L0 - D = 0;
    y[N_Trays, :] = xc[:];
  //energy balance
    VNT * sum(yNT[:] .* hv_B[:]) - V[1] * sum(y[1, :] .* hv[1, :]) + L[2] * sum(x[2, :] .* hl[2, :]) - L[1] * sum(x[1, :] .* hl[1, :]) + F[1] * sum(z[1, :] .* hf[1, :]) = M[1] * x[1, :] * dhl[1, :] + M[1] * dx[1, :] * hl[1, :] + dM[1] * x[1, :] * hl[1, :];
    for i in 2:N_Trays - 1 loop
      V[i - 1] * sum(y[i - 1, :] .* hv[i - 1, :]) - V[i] * sum(y[i, :] .* hv[i, :]) + L[i + 1] * sum(x[i + 1, :] .* hl[i + 1, :]) - L[i] * sum(x[i, :] .* hl[i, :]) + F[i] * sum(z[i, :] .* hf[i, :]) = M[i] * x[i, :] * dhl[i, :] + M[i] * dx[i, :] * hl[i, :] + dM[i] * x[i, :] * hl[i, :];
    end for;
    V[N_Trays - 1] * sum(y[N_Trays - 1, :] .* hv[N_Trays - 1, :]) - V[N_Trays] * sum(y[N_Trays, :] .* hv[N_Trays, :]) + L0 * sum(xc[:] .* hl_C[:]) - L[N_Trays] * sum(x[N_Trays, :] .* hl[N_Trays, :]) + F[N_Trays] * sum(z[N_Trays, :] .* hf[N_Trays, :]) = dM[N_Trays] * x[N_Trays, :] * hl[N_Trays, :] + M[N_Trays] * dx[N_Trays, :] * hl[N_Trays, :] + M[N_Trays] * x[N_Trays, :] * dhl[N_Trays, :];
    V[N_Trays] * sum(y[N_Trays, :] .* hv[N_Trays, :]) - (L0 + D) * sum(xc[:] .* hl_C[:]) = QC;
    L[1] * sum(x[1, :] .* hl[1, :]) - B * sum(xr[:] .* hl_B[:]) - VNT * sum(xr[:] .* hv_B[:]) = QB;
  //pressure
    for i in 1:N_Trays loop
      P[i] = P_condenser + (N_Trays - i + 1) * Pressure_drop;
    end for;
  //Equilibrium
    for i in 1:N_Trays loop
      Keq[i, :] = Psat[i, :] ./ P[i];
      y_eq[i, :] = Keq[i, :] .* x[i, :];
      M_eff[i, :] = (y[i, :] - y[i - 1, :]) ./ (y_eq[i, :] - y[i - 1, :]);
      sum(x[i, :]) = 1;
      sum(y_eq[i, :]) = 1;
    end for;
    sum(y[N_Trays, :] .* PsatC[:] / P_condenser) = 1;
  //sum(xc[:]) =1
    sum(x[1, :] .* (P[1] + Pressure_drop) ./ PsatB[:]) = 1;
  //sum(yNT[:]) =1;
  //  D = 0.5 * L0;
  //  B = 61.1;
    port2.moleflow = D;
    port2.molefrac = {xc[2], 0, xc[1], 0};
    port2.temperature = TC;
    port2.pressure = P_condenser;
   /* port3.moleflow = B;
    port3.molefrac = {xr[2], 0, xr[1], 0};
    port3.temperature = TB;
    port3.pressure = P_condenser + Pressure_drop * (N_Trays + 1); */
  //Equations for Specification
  if Integer(specification1) == 1 then
  QC = specification1_value;
  elseif Integer(specification1) == 2 then 
  D = specification1_value;
  elseif Integer(specification1) == 3 then
  xc[1] * D = specification1_value; //yet to modify
  elseif Integer(specification1) == 4 then
  xc[1] = specification1_value; //yet to modify
  elseif Integer(specification1) == 5 then
  L0 = specification1_value * D;
  else TC = specification1_value;
  end if;
  
  if Integer(specification2) == 1 then
  QB = specification2_value;
  elseif Integer(specification2) == 2 then 
  B = specification2_value;
  elseif Integer(specification2) == 3 then
  xr[1] * D = specification2_value; //yet to modify
  elseif Integer(specification2) == 4 then
  xc[1] = specification2_value; //yet to modify
  elseif Integer(specification2) == 5 then
  VNT = specification2_value * B;
  else TB = specification2_value;
  end if;
    annotation(Icon(graphics = {Rectangle(origin = {-2, -1}, fillColor = {0, 85, 127}, fillPattern = FillPattern.VerticalCylinder, extent = {{-94, 95}, {94, -95}})}), Documentation(info = "<HTML> <p> This is a generalized model for distilation column </p> </HTML>"));
  end DistillationWithSizing;

  model FlashWithSizing
    parameter Real hset = 3.7 "units = m" annotation(Dialog(group = "Operating conditions")), Pset = 5e5 "units = Pa" annotation(Dialog(group = "Operating conditions"));
    extends compounds;
    parameter Boolean connectedToInput = false;
    parameter Boolean OverrideSizeCalculations(start = false) annotation(Dialog(tab = "Sizing"));
    parameter Real k_drum = 0.3 "units = ft/s" annotation(Dialog(tab = "Sizing"));
    parameter Real Area = 4 "units = m2" annotation(Dialog(tab = "Sizing")), Volume = 8 "units = m3" annotation(Dialog(tab = "Sizing"));
    parameter Real Ti = 310 "units = K" annotation(Dialog(group = "Operating conditions"));
    parameter Real R = 8.314 "units = kJ/kmol.K", A(fixed = false), V_Total(fixed = false);
    Real z[NOC], Tf;
    Real y[NOC], x[NOC](start = {0.7, 1e-18, 0.3, 0}), k[NOC], L(start = 100, min = 0), V(start = 140, min = 0), Psat_T[NOC], M[NOC], M_Total, ML(start = 50), MG(start = 0.5), VL, VG, Q, hv[NOC], hl[NOC], Hf, Hv, Hl, H_M_Total, F, densityi[NOC], P, h;
    unitoperations.sensor sensor1 annotation(Placement(visible = true, transformation(origin = {2, 82}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {8.88178e-16, 82}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
    unitoperations.sensor sensor3 annotation(Placement(visible = true, transformation(origin = {82, -32}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {77, -31}, extent = {{-19, -19}, {19, 19}}, rotation = 0)));
    unitoperations.port port1 annotation(Placement(visible = true, transformation(origin = {1, -83}, extent = {{-17, -17}, {17, 17}}, rotation = 0), iconTransformation(origin = {-8.88178e-16, -74}, extent = {{-18, -18}, {18, 18}}, rotation = 0)));
    unitoperations.port port2 annotation(Placement(visible = true, transformation(origin = {80, 48}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {76, 52}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
    unitoperations.port port3 annotation(Placement(visible = true, transformation(origin = {-84, 0}, extent = {{-18, -18}, {18, 18}}, rotation = 0), iconTransformation(origin = {-80, 4}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  initial equation
    h = hset;
    P = Pset;
    for i in 1:NOC - 1 loop
      der(M[i]) = 0;
    end for;
    der(M_Total) = 0;
  //der(H_M_Total) = 0;
    if OverrideSizeCalculations == false then
      k_drum *0.3048* ((sum(x[:] .* densityi[:]) - P/(R*Ti*1000))*P/(R*Ti*1000))^0.5 * A = V;
      V_Total = A * 4 * (4 * A / 3.14)^0.5;
    else
      A = Area;
      V_Total = Volume;
    end if;
  equation
    F = port3.moleflow;
    z[:] = port3.molefrac[:];
    Tf = port3.temperature;
    if connectedToInput == false then
     port3.pressure = P; 
    end if;
    for i in 1:NOC loop
      Psat_T[i] = Functions.Psat(comp[i].VP, Ti);
      hv[i] = Functions.HVapId(comp[i].VapCp, comp[i].HOV, comp[i].Tc, Ti);
      hl[i] = Functions.HLiqId(comp[i].VapCp, comp[i].HOV, comp[i].Tc, Ti);
      densityi[i] = Functions.Density(comp[i].LiqDen, comp[i].Tc, Ti, P);
    end for;
  /* Mass Blanace equations */
    der(M_Total) = F - L - V;
  //F - L - V = 0;
    for i in 1:NOC - 1 loop
      der(M[i]) = F * z[i] - L * x[i] - V * y[i];
  //F * z[i] - L * x[i] - V * y[i] = 0;
      M[i] = ML * x[i] + MG * y[i];
    end for;
    sum(M[:]) = M_Total;
    M_Total = MG + ML;
    VL = ML / sum(x[:] .* densityi[:]);
    VG = V_Total - VL;
    P * VG = MG * R * Ti * 1000;
  //ideal gas law for gas phase
  /*energy balance */
    Hv = V * sum(y[:] .* hv[:]);
    Hf = port3.enthalpy;
    Hl = L * sum(x[:] .* hl[:]);
    H_M_Total = ML * sum(x[:] .* hl[:]) + MG * sum(y[:] .* hv[:]);
  //Hf - Hv - Hl = der(H_M_Total);
    Hf - Hv - Hl + Q = 0;
  /*Thermodynamic equations */
    sum(x[:]) = 1;
    sum(y[:]) = 1;
    y[:] = k[:] .* x[:];
    k[:] = Psat_T[:] / P;
  /* Control */
    A * h = VL;
  //connector equations
    sensor1.var = P;
    sensor3.var = h;
    port1.moleflow = L;
    port1.pressure = P;
    port1.temperature = Ti;
    port1.molefrac[:] = x[:];
    port2.moleflow = V;
    port2.pressure = P;
    port2.temperature = Ti;
    port2.molefrac[:] = y[:];
  //port3.pressure = P;
    annotation(Icon(graphics = {Text(origin = {-30, 86}, extent = {{-22, 32}, {22, -32}}, textString = "Pressure"), Text(origin = {56, 46}, extent = {{-10, 24}, {2, -2}}, textString = "V"), Text(origin = {60, -33}, extent = {{-12, 25}, {4, -3}}, textString = "h"), Text(origin = {-16, -82}, extent = {{-14, 26}, {2, -6}}, textString = "L"), Text(origin = {0, 15}, extent = {{-46, 41}, {46, -41}}, textString = "PT flash"), Rectangle(origin = {1, 0}, extent = {{-87, 96}, {87, -96}}), Text(origin = {-64, -21}, extent = {{-10, 17}, {10, -17}}, textString = "F")}, coordinateSystem(initialScale = 0.1)));
  
  end FlashWithSizing;

  model Pump
  end Pump;

  model FlashWithSizingTest
    MaterialStream materialStream1(Flowrate = 1, Pressure = 10e5, Temperature = 300, molefraction = {0.25, 0.25, 0.25, 0.25}, step_value = 0.001, stepchange = true, unspecified = false) annotation(Placement(visible = true, transformation(origin = {-78, -2}, extent = {{-18, -18}, {18, 18}}, rotation = 0)));
    FlashWithSizing flash1(OverrideSizeCalculations = false, connectedToInput = true)  annotation(Placement(visible = true, transformation(origin = {-5, -1}, extent = {{-29, -29}, {29, 29}}, rotation = 0)));
    MaterialStream materialStream2(Tdf(start = 279))  annotation(Placement(visible = true, transformation(origin = {101, 27}, extent = {{-19, -19}, {19, 19}}, rotation = 0)));
    MaterialStream materialStream3(Tdf(start = 370))  annotation(Placement(visible = true, transformation(origin = {84, -38}, extent = {{-24, -24}, {24, 24}}, rotation = 0)));
    valve valve1(OutletPfixed = true, OutletPressure = 1e5) annotation(Placement(visible = true, transformation(origin = {51, 25}, extent = {{-17, -17}, {17, 17}}, rotation = 0)));
    valve valve2(OutletPfixed = true, OutletPressure = 1e5) annotation(Placement(visible = true, transformation(origin = {26, -48}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(materialStream1.port2, flash1.port3) annotation(Line(points = {{-62, -2}, {-28, -2}, {-28, 0}, {-28, 0}}));
    connect(valve2.port2, materialStream3.port1) annotation(Line(points = {{34, -48}, {64, -48}, {64, -38}, {64, -38}}));
    connect(flash1.port1, valve2.port1) annotation(Line(points = {{-4, -22}, {-4, -22}, {-4, -48}, {18, -48}, {18, -48}}));
    connect(valve1.port2, materialStream2.port1) annotation(Line(points = {{64, 26}, {84, 26}, {84, 28}, {86, 28}}));
    connect(flash1.port2, valve1.port1) annotation(Line(points = {{18, 14}, {36, 14}, {36, 24}, {38, 24}}));
  end FlashWithSizingTest;

  model PhFlashWithSizing
    parameter Real hset = 3.7 "units = m" annotation(Dialog(group = "Operating conditions")), Pset = 5e5 "units = Pa" annotation(Dialog(group = "Operating conditions"));
    extends compounds;
    parameter Boolean connectedToInput = false;
    parameter Boolean OverrideSizeCalculations(start = false) annotation(Dialog(tab = "Sizing"));
    parameter Real k_drum = 0.3 "units = ft/s" annotation(Dialog(tab = "Sizing"));
    parameter Real Area = 4 "units = m2" annotation(Dialog(tab = "Sizing")), Volume = 8 "units = m3" annotation(Dialog(tab = "Sizing"));
  //  parameter Real Ti = 310 "units = K" annotation(Dialog(group = "Operating conditions"));
    protected parameter Real R = 8.314 "units = kJ/kmol.K", A(fixed = false), V_Total(fixed = false);
    Real z[NOC];
    Real y[NOC], x[NOC](start = {0.5, 1e-15, 0.5, 0}, each min = 0), k[NOC], L(start = 0.5, min = 0), V(start = 0.5, min = 0), Psat_T[NOC], M[NOC], M_Total, ML(start = 50), MG(start = 0.5), VL, VG, Q, hv[NOC], hl[NOC], Hf, Hv, Hl, H_M_Total, F, densityi[NOC], P, h, Ti(start = 290);
    unitoperations.sensor sensor1 annotation(Placement(visible = true, transformation(origin = {2, 82}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {8.88178e-16, 82}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
    unitoperations.sensor sensor3 annotation(Placement(visible = true, transformation(origin = {82, -32}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {77, -31}, extent = {{-19, -19}, {19, 19}}, rotation = 0)));
    unitoperations.port port1 annotation(Placement(visible = true, transformation(origin = {1, -83}, extent = {{-17, -17}, {17, 17}}, rotation = 0), iconTransformation(origin = {-8.88178e-16, -74}, extent = {{-18, -18}, {18, 18}}, rotation = 0)));
    unitoperations.port port2 annotation(Placement(visible = true, transformation(origin = {80, 48}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {76, 52}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
    unitoperations.port port3 annotation(Placement(visible = true, transformation(origin = {-84, 0}, extent = {{-18, -18}, {18, 18}}, rotation = 0), iconTransformation(origin = {-80, 4}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  initial equation
    h = hset;
    P = Pset;
    for i in 1:NOC - 1 loop
      der(M[i]) = 0;
    end for;
    der(M_Total) = 0;
  //der(H_M_Total) = 0;
    if OverrideSizeCalculations == false then
      k_drum *0.3048* ((sum(x[:] .* densityi[:]) - P/(R*Ti*1000))*P/(R*Ti*1000))^0.5 * A = V;
      V_Total = A * 4 * (4 * A / 3.14)^0.5;
    else
      A = Area;
      V_Total = Volume;
    end if;
  equation
    F = port3.moleflow;
    z[:] = port3.molefrac[:];
    Q = 0;
    if connectedToInput == false then
     port3.pressure = P; 
    end if;
    for i in 1:NOC loop
      Psat_T[i] = Functions.Psat(comp[i].VP, Ti);
      hv[i] = Functions.HVapId(comp[i].VapCp, comp[i].HOV, comp[i].Tc, Ti);
      hl[i] = Functions.HLiqId(comp[i].VapCp, comp[i].HOV, comp[i].Tc, Ti);
      densityi[i] = Functions.Density(comp[i].LiqDen, comp[i].Tc, Ti, P);
    end for;
  /* Mass Blanace equations */
    der(M_Total) = F - L - V;
  //F - L - V = 0;
    for i in 1:NOC - 1 loop
      der(M[i]) = F * z[i] - L * x[i] - V * y[i];
  //F * z[i] - L * x[i] - V * y[i] = 0;
      M[i] = ML * x[i] + MG * y[i];
    end for;
    sum(M[:]) = M_Total;
    M_Total = MG + ML;
    VL = ML / sum(x[:] .* densityi[:]);
    VG = V_Total - VL;
    P * VG = MG * R * Ti * 1000;
  //ideal gas law for gas phase
  /*energy balance */
    Hv = V * sum(y[:] .* hv[:]);
    Hf = port3.enthalpy;
    Hl = L * sum(x[:] .* hl[:]);
    H_M_Total = ML * sum(x[:] .* hl[:]) + MG * sum(y[:] .* hv[:]);
  //Hf - Hv - Hl = der(H_M_Total);
    Hf - Hv - Hl + Q = 0;
  /*Thermodynamic equations */
    sum(x[:]) = 1;
    sum(y[:]) = 1;
    y[:] = k[:] .* x[:];
    k[:] = Psat_T[:] / P;
  /* Control */
    A * h = VL;
  //connector equations
    sensor1.var = P;
    sensor3.var = h;
    port1.moleflow = L;
    port1.pressure = P;
    port1.temperature = Ti;
    port1.molefrac[:] = x[:];
    port2.moleflow = V;
    port2.pressure = P;
    port2.temperature = Ti;
    port2.molefrac[:] = y[:];
    annotation(Icon(graphics = {Text(origin = {-30, 86}, extent = {{-22, 32}, {22, -32}}, textString = "Pressure"), Text(origin = {56, 46}, extent = {{-10, 24}, {2, -2}}, textString = "V"), Text(origin = {60, -33}, extent = {{-12, 25}, {4, -3}}, textString = "h"), Text(origin = {-16, -82}, extent = {{-14, 26}, {2, -6}}, textString = "L"), Text(origin = {0, 15}, extent = {{-46, 41}, {46, -41}}, textString = "PT flash"), Rectangle(origin = {1, 0}, extent = {{-87, 96}, {87, -96}}), Text(origin = {-64, -21}, extent = {{-10, 17}, {10, -17}}, textString = "F")}, coordinateSystem(initialScale = 0.1)));
  
  end PhFlashWithSizing;

  model PhFlashWithSizingTest
    MaterialStream materialStream1(Flowrate = 1, Pressure = 10e5, Temperature = 300, molefraction = {0.25, 0.25, 0.25, 0.25}, step_value = 0.001, stepchange = true, unspecified = false) annotation(Placement(visible = true, transformation(origin = {-78, -2}, extent = {{-18, -18}, {18, 18}}, rotation = 0)));
    PhFlashWithSizing flash1(OverrideSizeCalculations = false, connectedToInput = true)  annotation(Placement(visible = true, transformation(origin = {-5, -1}, extent = {{-29, -29}, {29, 29}}, rotation = 0)));
    MaterialStream materialStream2 annotation(Placement(visible = true, transformation(origin = {101, 27}, extent = {{-19, -19}, {19, 19}}, rotation = 0)));
    MaterialStream materialStream3(Tbf(start = 350))  annotation(Placement(visible = true, transformation(origin = {84, -38}, extent = {{-24, -24}, {24, 24}}, rotation = 0)));
    valve valve1(OutletPfixed = true, OutletPressure = 1e5) annotation(Placement(visible = true, transformation(origin = {51, 25}, extent = {{-17, -17}, {17, 17}}, rotation = 0)));
    valve valve2(OutletPfixed = true, OutletPressure = 1e5) annotation(Placement(visible = true, transformation(origin = {26, -48}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(materialStream1.port2, flash1.port3) annotation(Line(points = {{-62, -2}, {-28, -2}, {-28, 0}, {-28, 0}}));
    connect(valve2.port2, materialStream3.port1) annotation(Line(points = {{34, -48}, {64, -48}, {64, -38}, {64, -38}}));
    connect(flash1.port1, valve2.port1) annotation(Line(points = {{-4, -22}, {-4, -22}, {-4, -48}, {18, -48}, {18, -48}}));
    connect(valve1.port2, materialStream2.port1) annotation(Line(points = {{64, 26}, {84, 26}, {84, 28}, {86, 28}}));
    connect(flash1.port2, valve1.port1) annotation(Line(points = {{18, 14}, {36, 14}, {36, 24}, {38, 24}}));
  end PhFlashWithSizingTest;
end unitoperations;
