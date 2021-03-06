package unitoperations

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
    annotation(Diagram(graphics = {Ellipse(origin = {-1, 1}, fillColor = {0, 85, 255}, fillPattern = FillPattern.Solid, extent = {{-67, 65}, {67, -65}}, endAngle = 360)}, coordinateSystem(initialScale = 0.1)), Icon(graphics = {Ellipse(origin = {-4, -1}, fillColor = {0, 85, 255}, fillPattern = FillPattern.Solid, extent = {{-60, 61}, {60, -61}}, endAngle = 360)}, coordinateSystem(initialScale = 0.1)));
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
    annotation(Diagram(graphics = {Ellipse(origin = {-1, 1}, fillColor = {65, 252, 255}, fillPattern = FillPattern.Solid, extent = {{-67, 65}, {67, -65}}, endAngle = 360)}), Icon(graphics = {Ellipse(origin = {-4, -1}, fillColor = {170, 170, 255}, fillPattern = FillPattern.Sphere, extent = {{-60, 61}, {60, -61}}, endAngle = 360)}, coordinateSystem(initialScale = 0.1)));
  end port;

  model valve
    parameter Real coeff(fixed = false) "Coeff for valve", valveCv = 0.4 "valve Cv if not control valve";
    parameter Boolean Control = false;
    parameter Boolean OutletPfixed = false;
    parameter Real OutletPressure(unit = "atm") = 1 "used only when OutletPfixed is true";
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
      outletP = OutletPressure * 101325;
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
    annotation(Icon(graphics = {Polygon(origin = {-48.26, -2.99}, fillColor = {170, 170, 255}, fillPattern = FillPattern.Solid, points = {{-47.7381, 48.9887}, {-47.7381, -49.0113}, {48.2619, 2.98874}, {48.2619, 2.98874}, {-47.7381, 48.9887}}), Polygon(origin = {49.25, -4.98}, fillColor = {170, 170, 255}, fillPattern = FillPattern.Solid, points = {{-47.2509, 4.98071}, {46.7491, 48.9807}, {46.7491, -49.0193}, {-47.2509, 4.98071}}), Rectangle(origin = {1, 35}, fillColor = {170, 170, 255}, fillPattern = FillPattern.Solid, extent = {{-15, 35}, {15, -35}}), Text(origin = {0, -73}, extent = {{-52, 25}, {52, -25}}, textString = "Valve")}, coordinateSystem(initialScale = 0.1)));
  end valve;

  model MaterialStream
    extends compounds;
    parameter Boolean specified_stream = false;
    parameter Real Flowrate(unit = "mol/s") = 100 annotation(HideResult = true), Pressure(unit = "atm") = 1 annotation(HideResult = true), Temperature(unit = "K") = 300 annotation(HideResult = true), molefraction[NOC] = zeros(NOC) annotation(each HideResult = true);
    
    parameter Boolean stepchange = false annotation(Dialog(group = "Dynamic"));
    parameter Real stepchangetime = 0.01 annotation(Dialog(group = "Dynamic"));
    parameter Real step_value = 1 annotation(Dialog(group = "Dynamic"));
    
    Real kf[NOC] annotation(each HideResult = true);
    
    Real zl[NOC](each min = 0, each max = 1, start = {0.5, 1e-18, 0.5, 0}) "liquid composition", zv[NOC](each min = 0, each max = 1, start = {0, 0.25, 0, 0.75}) "Vapor composition";
    
    Real Fl(unit = "mol/s", min = 0, start = 100) "Liquid flowrate", Fv(unit = "mol/s", min = 0, start = 140) "Vapor flowrate";
    
    Real Tbf(unit = "K", start = 62) "Bubble point of stream", Tdf(unit = "K", start = 495.5) "Dew Point of stream";
    Real Psat_Tdf[NOC](each unit = "Pa") annotation(each HideResult = true), Psat_Tbf[NOC](each unit = "Pa") annotation(each HideResult = true), Psat_Tf[NOC](each unit = "Pa") annotation(each HideResult = true);
    
    Real Pf(unit = "Pa") "Pressure of stream", Tf(unit = "K") "Temperature of stream", z[NOC] "Composition of stream", F(unit = "mol/s") "Flowrate of stream";
    
    Real Hvf[NOC](each unit = "J/mol") annotation(each HideResult = true), Hlf[NOC](each unit = "J/mol") annotation(each HideResult = true), H(unit = "J/s");
    
    unitoperations.port port2 annotation(Placement(visible = true, transformation(origin = {80, -4}, extent = {{-18, -18}, {18, 18}}, rotation = 0), iconTransformation(origin = {85, 1}, extent = {{-21, -21}, {21, 21}}, rotation = 0)));
    unitoperations.port port1 annotation(Placement(visible = true, transformation(origin = {-80, -4}, extent = {{-18, -18}, {18, 18}}, rotation = 0), iconTransformation(origin = {-82, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  
  equation
  
    if specified_stream == true then
  
      if stepchange == true then  
  
        if time < stepchangetime then
          port1.moleflow = Flowrate;
        else
          port1.moleflow = Flowrate + step_value;
        end if;
  
      else
        port1.moleflow = Flowrate;
      end if;
  
      port1.pressure = Pressure * 101325;
      port1.temperature = Temperature;
      port1.molefrac[:] = molefraction[:];
      port1.liquidmolefrac[:] = zl[:];
      port1.vapormolefrac[:] = zv[:];
      port1.liquidmoleflow = Fl;
      port1.vapormoleflow = Fv;
      port1.enthalpy = H;
    end if;
    
    if specified_stream == false then
      port1.liquidmolefrac[:] = zl[:];
      port1.vapormolefrac[:] = zv[:];
      port1.liquidmoleflow = Fl;
      port1.vapormoleflow = Fv;
      port1.enthalpy = H;
    end if;
    
    if specified_stream == true then
      Pf = Pressure * 101325;
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
      
      if Tf <= Tbf + 0.01 then
        zl[:] = z[:];
        zv = zeros(NOC);
        Fl = F;
        Fv = 0;
        kf = zeros(NOC);
      elseif Tf >= Tdf - 0.01 then
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
    
    annotation(Icon(graphics = {Text(origin = {-17, -1}, extent = {{-51, 27}, {51, -27}}, textString = ""),  Text(origin = {0, -70}, extent = {{-86, 20}, {86, -20}}, textString = "Material"), Rectangle(origin = {-22, 0}, fillColor = {170, 170, 255}, fillPattern = FillPattern.Solid, extent = {{-72, 20}, {72, -20}}), Polygon(origin = {74.81, 1.98}, fillColor = {170, 170, 255}, fillPattern = FillPattern.Solid, points = {{-24.8123, 38.0153}, {25.1877, -1.98467}, {-24.8123, -37.9847}, {-24.8123, 18.0153}, {-24.8123, 38.0153}})}, coordinateSystem(initialScale = 0.1)));
  end MaterialStream;




































  model CSTR
    extends unitoperations.compounds;
    parameter Real s[NOC] = {-1, -1, 1, 1} "array for A + B --> C is {-1,-1,1}" annotation(Dialog(tab = "Reactions", group = "Stoichiometric constants"));
    constant Integer NOIS = 2 "No of input streams" ;
    parameter Integer n = 1 "base compound for reaction rate, ri then n =i" annotation(Dialog(tab = "Reactions", group = "Stoichiometric constants"));
    parameter Real V_Total(unit = "m3") = 2.321 "Volume of reactor", P_init(unit = "atm") = 25 "Pressure at t = 0";
    constant Real R(unit = "J/mol.K") = 8.314;
    parameter Real Af "frequency factor for forward reaction" annotation(Dialog(tab = "Reactions", group = "forward reaction rate constants"));
    parameter Real order_f[4] "order wrt to components for forward reaction" annotation(Dialog(tab = "Reactions", group = "forward reaction rate constants"));
    parameter Real order_b[4] "order wrt to components for backward reaction" annotation(Dialog(tab = "Reactions", group = "backward reaction rate constants"));
    parameter Real Ab "frequency factor for backward reaction" annotation(Dialog(tab = "Reactions", group = "backward reaction rate constants"));
    parameter Real Eaf "Activation energy for forward reaction" annotation(Dialog(tab = "Reactions", group = "forward reaction rate constants"));
    parameter Real Eab "Activation energy for backward reaction" annotation(Dialog(tab = "Reactions", group = "backward reaction rate constants"));
    parameter Real delH_r(unit = "J/mol") = 12.6e3 "Heat of reaction" annotation(Dialog(tab = "Reactions", group = "Reaction rate constants"));
    type operation_type = enumeration(Isothermal, Adiabatic);
    parameter operation_type operation_mode;
    parameter Real T_iso(unit = "K") = 900;
    
    Real F_in[NOIS](each unit = "mol/s") annotation(each HideResult = true), z[NOIS, NOC] annotation(each HideResult = true), Hin[NOIS](each unit = "J/s") annotation(each HideResult = true);
    Real M_Total(unit = "mol", start = 1.5), M[NOC](each unit = "mol"), x[NOC], F_out(unit = "mol/s"), densityi[NOC](each unit = "kmol/m3") annotation(each HideResult = true), P(unit = "Pa");
    Real r(unit = "mol/s"), kf, kb, c[NOC](each unit = "mol/m3");
    Real H_r(unit = "J/s"), Hout(unit = "J/s"), Q(unit = "J/s"), T(unit = "K");
    
    unitoperations.port port1 annotation(Placement(visible = true, transformation(origin = {-82, 2}, extent = {{-18, -18}, {18, 18}}, rotation = 0), iconTransformation(origin = {-87, 1}, extent = {{-17, -17}, {17, 17}}, rotation = 0)));
    unitoperations.port port2 annotation(Placement(visible = true, transformation(origin = {-84, 52}, extent = {{-18, -18}, {18, 18}}, rotation = 0), iconTransformation(origin = {-87, 37}, extent = {{-17, -17}, {17, 17}}, rotation = 0)));
    unitoperations.port port3 annotation(Placement(visible = true, transformation(origin = {89, -75}, extent = {{-17, -17}, {17, 17}}, rotation = 0), iconTransformation(origin = {87, -47}, extent = {{-17, -17}, {17, 17}}, rotation = 0)));
  
  initial equation
//calculates steady state solution at t=0
    P = P_init * 101325;
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
      c[i] = x[i] * P / (R * T);
    end for;

  //reaction rate
    kf = Af * exp(-Eaf / (R * T));
    kb = Ab * exp(-Eab / (R * T));
    r = kf * product(c[:] .^ order_f[:]) - kb * product(c[:] .^ order_b[:]);

  //unsteady state mass balance
    for i in 1:NOC - 1 loop
      der(M[i]) = sum(z[:, i] .* F_in[:]) - x[i] * F_out + s[i] * r * V_Total * 1000 / abs(s[n]);
    end for;
    
    der(M_Total) = sum(F_in[:]) - F_out;

  //Pressure
    M_Total = sum(M[:]);
    P * V_Total = M_Total * R * T;

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
    
    annotation(Icon(coordinateSystem(extent = {{-100, -120}, {100, 100}}), graphics = {Polygon(origin = {1, -1}, fillColor = {170, 170, 255}, fillPattern = FillPattern.Solid, points = {{-95, 65}, {-95, -59}, {-81, -73}, {-55, -79}, {-33, -83}, {-7, -83}, {23, -83}, {53, -79}, {81, -69}, {95, -59}, {95, 65}, {75, 77}, {43, 83}, {-1, 83}, {-47, 83}, {-81, 75}, {-95, 65}}), Line(origin = {4.06, 27.01}, points = {{-2.06154, 70.9903}, {-2.06154, -53.0097}, {-16.0615, -43.0097}, {-30.0615, -41.0097}, {-44.0615, -43.0097}, {-52.0615, -49.0097}, {-58.0615, -59.0097}, {-52.0615, -69.0097}, {-32.0615, -73.0097}, {-18.0615, -67.0097}, {-2.06154, -53.0097}, {7.93846, -45.0097}, {13.9385, -41.0097}, {27.9385, -41.0097}, {45.9385, -43.0097}, {53.9385, -51.0097}, {55.9385, -61.0097}, {41.9385, -69.0097}, {27.9385, -71.0097}, {17.9385, -69.0097}, {7.93846, -61.0097}, {-2.06154, -55.0097}}, thickness = 1), Text(origin = {-2, -98}, extent = {{-70, 16}, {76, -24}}, textString = "CSTR")}), Diagram(coordinateSystem(extent = {{-100, -120}, {100, 100}})), version = "", uses);
  end CSTR;
























  model HeatExchanger
    extends compounds;
    parameter Real Tout = 38 + 273;
    parameter Real pressure_drop = 0.5;
    
    Real Hin, Hout, Q "heat exchanger duty";
    
    port port1 annotation(Placement(visible = true, transformation(origin = {-86, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-86, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    port port2 annotation(Placement(visible = true, transformation(origin = {88, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {88, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  
  equation
    Hin = port1.enthalpy;
    Hout = port2.enthalpy;
    Tout = port2.temperature;
    port2.pressure = port1.pressure - pressure_drop * 101325;
    port2.moleflow = port1.moleflow;
    port2.molefrac[:] = port1.molefrac[:];
    Hin = Hout + Q;
    
    annotation(Icon(coordinateSystem(extent = {{-100, -140}, {100, 100}}), graphics = {Ellipse(fillColor = {170, 170, 255}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}}, endAngle = 360), Line(origin = {-1, 7.01}, points = {{-87, -7.00598}, {-37, -7.00598}, {-9, 30.994}, {13, -39.006}, {31, -3.00598}, {87, -3.00598}}, thickness = 1), Text(origin = {-25, -106}, extent = {{-49, 16}, {93, -48}}, textString = "HX")}), Diagram(coordinateSystem(extent = {{-100, -140}, {100, 100}})), version = "", uses);
  end HeatExchanger;





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
    annotation(Icon(graphics = {Polygon(origin = {3.29, -3.99}, fillColor = {170, 170, 255}, fillPattern = FillPattern.CrossDiag, points = {{-91.2946, 7.98643}, {88.7054, 73.9864}, {90.7054, -74.0136}, {-91.2946, 7.98643}})}, coordinateSystem(initialScale = 0.1)));
  end CompoundSeperator;

  model Distillation
    parameter Chemsep_Database.Benzene comp1 annotation(Dialog(tab = "General", group = "compounds"));
    parameter Chemsep_Database.Toluene comp2 annotation(Dialog(tab = "General", group = "compounds"));
    parameter Chemsep_Database.General_Properties comp[2] = {comp1, comp2} annotation(Dialog(tab = "General", group = "compounds"));
    parameter Integer N_Trays = 20 "No. of trays without condensor and reboiler" annotation(Dialog(tab = "General", group = "Trays"));
    parameter Integer NOC = 2 "No. of compounds" annotation(Dialog(tab = "General", group = "compounds"));
    parameter Integer N_Feed = 10 "Feed tray location" annotation(Dialog(tab = "General", group = "Trays"));
    parameter Real M_eff = 0.99 "Murphree efficiency of trays" annotation(Dialog(tab = "General", group = "Efficiency"), each HideResult = true);
    parameter Real Pressure_drop(unit = "atm") = 0.007 "Pressure Drop per tray" annotation(Dialog(tab = "General", group = "Pressure Profile"));
    parameter Real P_condenser(unit = "atm") = 1 "Pressure in Condensor" annotation(Dialog(tab = "General", group = "Pressure Profile"));
    constant Real R = 8.314 "Gas Constant";
    parameter Real Active_area(unit = "m2") = 1 "Active area of Tray" annotation(Dialog(tab = "dynamic", group = "Tray data"));
    parameter Real Weir_diameter(unit = "m") = 0.8 "Diameter of weir" annotation(Dialog(tab = "dynamic", group = "Tray data"));
    parameter Real A_active(unit = "m2", fixed = false, start = 0.6);
    parameter Real d_weir(unit = "m", fixed = false);
    parameter Boolean Override_Sizing_Calculations;
    parameter Real Kv(unit = "ft/s") = 0.3 "constant for calculating max velocity permissible, ft/s" annotation(Dialog(tab = "dynamic", group = "Tray data"));
    parameter Real h_weir(unit = "m") = 0.1 "Height of weir" annotation(Dialog(tab = "dynamic", group = "Tray data"));
    type spec1 = enumeration(CondensorHeatLoad, ProductMolarFlow, CompoundMolarFlow, CompoundFractionInStream, RefluxRatio, Temperature) "condensor";
    type spec2 = enumeration(ReboilerHeatLoad, ProductMolarFlow, CompoundMolarFlow, CompoundFractionInStream, BoilUpRatio, Temperature) "reboiler";
    parameter spec1 specification1 "condensor";
    parameter Real specification1_value = 1;
    parameter spec2 specification2 "reboiler";
    parameter Real specification2_value = 1;
    parameter Real QC(unit = "J/s", fixed = false, start = 1.35e6, nominal = 1e6), QB(unit = "J/s", fixed = false, start = 1.5e6, nominal = 1e6);
    parameter Real Tray_volume(unit = "m3") = 0.1 annotation(Dialog(tab = "dynamic", group = "Tray data"));
    parameter Real pi1 = -12.55 "constant for calculationg froth density, phi" annotation(Dialog(tab = "dynamic", group = "vapor flow"));
    parameter Real pi2 = 0.91 "constant for calculationg froth density, phi" annotation(Dialog(tab = "dynamic", group = "vapor flow"));
  
    Real y[N_Trays, NOC](start = fill(0.5, N_Trays, NOC)), x[N_Trays, NOC](start = fill(0.5, N_Trays, NOC), each nominal = 1e-1), y_eq[N_Trays, NOC] annotation(each HideResult = true), Tf[N_Trays](each unit = "K") annotation(each HideResult = true);
  
    Real V[N_Trays](each unit = "mol/s", each start = 70), L[N_Trays](each unit = "mol/s", each start = 100);
  
    Real T[N_Trays](each unit = "K", start = linspace(386, 354, N_Trays), each nominal = 1e2), TC(unit = "K", start = 368), TB(unit = "K", start = 377), L0(unit = "mol/s", start = 50), VNT(unit = "mol/s"), D(unit = "mol/s"), B(unit = "mol/s"), xc[NOC], xr[NOC];
  
    Real Keq[N_Trays, NOC] annotation(each HideResult = true), Psat[N_Trays, NOC](each unit = "Pa") annotation(each HideResult = true), PsatC[NOC](each unit = "Pa") annotation(each HideResult = true), PsatB[NOC](each unit = "Pa") annotation(each HideResult = true);
  
    Real yNT[NOC], M[N_Trays](each unit = "mol", each start = 2000), den[N_Trays, NOC](each unit = "kmol/m3") annotation(each HideResult = true);
  
    Real hv[N_Trays, NOC](each unit = "J/mol", each nominal = 1e3) annotation(each HideResult = true), hl[N_Trays, NOC](each unit = "J/mol", each nominal = 1e4) annotation(each HideResult = true), hf[N_Trays, NOC](each unit = "J/mol") annotation(each HideResult = true), hv_B[NOC](each unit = "J/mol") annotation(each HideResult = true), hl_B[NOC](each unit = "J/mol") annotation(each HideResult = true), hl_C[NOC](each unit = "J/mol") annotation(each HideResult = true), Cpl[N_Trays,NOC](each unit = "J/mol.K") annotation(each HideResult = true);
  
    Real P[N_Trays](each unit = "Pa"), Ks[N_Trays] annotation(each HideResult = true);
  
    Real F[N_Trays](each start = 0) annotation(each HideResult = true), z[N_Trays, NOC](start = fill(0.5, N_Trays, NOC)) annotation(each HideResult = true);
  
    unitoperations.port port1 annotation(Placement(visible = true, transformation(origin = {-54, 0}, extent = {{-18, -18}, {18, 18}}, rotation = 0), iconTransformation(origin = {-56, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    unitoperations.port port2 annotation(Placement(visible = true, transformation(origin = {60, 82}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {50, 64}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    unitoperations.port port3 annotation(Placement(visible = true, transformation(origin = {62, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {50, -62}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  
  initial equation
/*sizing*/
    if Override_Sizing_Calculations == false then
      Kv * 0.3048 * ((den[1, :] * x[1, :] - P[1] / (R * T[1] * 1000)) * (P[1] / (R * T[1] * 1000))) ^ 0.5 * 1000 * A_active = max(V);
      d_weir = (A_active * 4 / 3.14) ^ 0.5;
    else
      A_active = Active_area;
      d_weir = Weir_diameter;
    end if;
    
    for i in 1:N_Trays loop
      der(M[i]) = 0;
      der(T[i]) = 0;
    end for;
    
    V[N_Trays] * sum(y[N_Trays, :] .* hv[N_Trays, :]) - (L0 + D) * sum(xc[:] .* hl_C[:]) = QC;
    L[1] * sum(x[1, :] .* hl[1, :]) - B * sum(xr[:] .* hl_B[:]) - VNT * sum(xr[:] .* hv_B[:]) = QB;
  
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
      Psat[:, i] = Functions.Psat(comp[i].VP, T[:]);
      hv[:, i] = Functions.HVapId(comp[i].VapCp, comp[i].HOV, comp[i].Tc, T[:]);
      hf[:, i] = Functions.HLiqId(comp[i].VapCp, comp[i].HOV, comp[i].Tc, Tf[:]);
      den[:, i] = Functions.Density(comp[i].LiqDen, comp[i].Tc, T[:], P[:]);
      Cpl[:,i] =Functions.LiqCpId(comp[i].LiqCp, T[:]);
      hl[:, i] = Functions.HLiqId(comp[i].VapCp, comp[i].HOV, comp[i].Tc, T[:])"have to change this"; 
      hv_B[i] = Functions.HVapId(comp[i].VapCp, comp[i].HOV, comp[i].Tc, TB);
      hl_B[i] = Functions.HLiqId(comp[i].VapCp, comp[i].HOV, comp[i].Tc, TB);
      hl_C[i] = Functions.HLiqId(comp[i].VapCp, comp[i].HOV, comp[i].Tc, TC);
      PsatC[i] = Functions.Psat(comp[i].VP, TC);
      PsatB[i] = Functions.Psat(comp[i].VP, TB);
    end for;
    
    for i in 1:N_Trays loop
      Ks[i] = V[i] * R * T[i] / (A_active * P[i]) * (P[i] / (R * T[i] * (1000 * sum(x[i, :] .* den[i, :]) - P[i] / (R * T[i])))) ^ 0.5;
    end for;

  //tray mass balance
    xr[:] = x[1, :];
    yNT[:] = xr[:];
    L[1] - VNT - B = 0;
    
    for i in 1:N_Trays loop
      M[i] = A_active * sum(x[i, :] .* den[i, :]) * 1000 * (exp(pi1 * Ks[i] ^ pi2) * h_weir + 44300 * 10 ^ (-3) * (L[i] / (sum(x[i, :] .* den[i, :]) * 1000 * d_weir * 1000)) ^ 0.704);
    end for;
    
    der(x[1, :] * M[1]) = VNT .* yNT[:] + L[2] .* x[2, :] - V[1] .* y[1, :] - L[1] .* x[1, :] + F[1] .* z[1, :];
    
    for i in 2:N_Trays - 1 loop
      der(x[i, :] * M[i]) = V[i - 1] .* y[i - 1, :] + L[i + 1] .* x[i + 1, :] - V[i] .* y[i, :] - L[i] .* x[i, :] + F[i] .* z[i, :];
    end for;
    
    der(x[N_Trays, :] * M[N_Trays]) = V[N_Trays - 1] .* y[N_Trays - 1, :] + L0 .* xc[:] - V[N_Trays] .* y[N_Trays, :] - L[N_Trays] .* x[N_Trays, :] + F[N_Trays] .* z[N_Trays, :];
    V[N_Trays] - L0 - D = 0;
    y[N_Trays, :] = xc[:];

  //energy balance
    VNT * sum(yNT[:] .* hv_B[:]) - V[1] * sum(y[1, :] .* hv[1, :]) + L[2] * sum(x[2, :] .* hl[2, :]) - L[1] * sum(x[1, :] .* hl[1, :]) + F[1] * sum(z[1, :] .* hf[1, :]) = der(x[1, :] * M[1]) * hl[1, :] + x[1, :] * M[1] * Cpl[1,:] * der(T[1]);
  
    for i in 2:N_Trays - 1 loop
      V[i - 1] * sum(y[i - 1, :] .* hv[i - 1, :]) - V[i] * sum(y[i, :] .* hv[i, :]) + L[i + 1] * sum(x[i + 1, :] .* hl[i + 1, :]) - L[i] * sum(x[i, :] .* hl[i, :]) + F[i] * sum(z[i, :] .* hf[i, :]) = der(x[i, :] * M[i]) * hl[i, :] + x[i, :] * M[i] * Cpl[i,:] * der(T[i]);
    end for;
  
    V[N_Trays - 1] * sum(y[N_Trays - 1, :] .* hv[N_Trays - 1, :]) - V[N_Trays] * sum(y[N_Trays, :] .* hv[N_Trays, :]) + L0 * sum(xc[:] .* hl_C[:]) - L[N_Trays] * sum(x[N_Trays, :] .* hl[N_Trays, :]) + F[N_Trays] * sum(z[N_Trays, :] .* hf[N_Trays, :]) = der(x[N_Trays, :] * M[N_Trays]) * hl[N_Trays, :] + x[N_Trays,:] * M[N_Trays] * Cpl[N_Trays,:] * der(T[N_Trays]);
    
//pressure
    for i in 1:N_Trays loop
      P[i] = (P_condenser + (N_Trays - i + 1) * Pressure_drop) * 101325;
    end for;

  //Equilibrium
    for i in 1:N_Trays loop
      Keq[i, :] = Psat[i, :] ./ P[i];
      y_eq[i, :] = Keq[i, :] .* x[i, :];
      {M_eff, M_eff} = (y[i, :] - y[i - 1, :]) ./ (y_eq[i, :] - y[i - 1, :]);
      sum(x[i, :]) = 1;
      sum(y_eq[i, :]) = 1;
    end for;
    
    sum(y[N_Trays, :] .* PsatC[:] / (P_condenser * 101325)) = 1;
    sum(x[1, :] .* (P[1] + Pressure_drop * 101325) ./ PsatB[:]) = 1;
    port2.moleflow = D;
    port2.molefrac = {xc[2], 0, xc[1], 0};
    port2.temperature = TC;
    port2.pressure = P_condenser * 101325;
    port3.moleflow = B;
    port3.molefrac = {xr[2], 0, xr[1], 0};
    port3.temperature = TB;
    port3.pressure = P[1] + Pressure_drop * 101325;

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
    else
      TC = specification1_value;
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
    else
      TB = specification2_value;
    end if;
  
    annotation(Icon(coordinateSystem(extent = {{-70, -140}, {70, 100}}, preserveAspectRatio = false), graphics = {Polygon(origin = {-1, 2}, fillColor = {170, 170, 255}, fillPattern = FillPattern.Solid, points = {{-63, -62.0013}, {-51, -78.0013}, {-33, -86.0013}, {-13, -90.0013}, {15, -90.0013}, {39, -84.0013}, {55, -76.0013}, {63, -62.0013}, {63, 71.9987}, {45, 81.9987}, {29, 87.9987}, {17, 89.9987}, {-15, 89.9987}, {-33, 85.9987}, {-49, 79.9987}, {-63, 69.9987}, {-63, -62.0013}}), Line(origin = {8, 55}, points = {{-54, 5}, {-54, -5}, {54, -5}}), Line(origin = {-10, 36}, points = {{-54, -4}, {54, -4}, {54, 4}}), Line(origin = {3, 14}, points = {{-59, 6}, {-59, -6}, {59, -6}}), Line(origin = {-9, -9}, points = {{-55, -5}, {53, -5}, {53, 5}}), Line(origin = {6, -35}, points = {{-58, 5}, {-58, -5}, {56, -5}}), Text(origin = {1, -102}, extent = {{-49, 18}, {57, -28}}, textString = "Distillation")}), Diagram(coordinateSystem(extent = {{-70, -140}, {70, 100}}, preserveAspectRatio = false)), version = "", uses, Documentation(info = "<HTML> <p> This is a generalized model for distilation column </p> </HTML>"));
  end Distillation;



















































  model PTFlash
    parameter Real hset(unit = "m") = 0.7 annotation(Dialog(group = "Operating conditions")), Pset(unit = "atm") = 5 annotation(Dialog(group = "Operating conditions"));
    extends compounds;
    parameter Boolean connectedToInput = false;
    parameter Boolean OverrideSizeCalculations = false annotation(Dialog(tab = "Sizing"));
    parameter Real k_drum(unit = "ft/s") = 0.3 annotation(Dialog(tab = "Sizing"));
    parameter Real Area(unit = "m2") = 4 annotation(Dialog(tab = "Sizing")), Volume(unit = "m3") = 8 annotation(Dialog(tab = "Sizing"));
    parameter Real Ti(unit = "K") = 310 annotation(Dialog(group = "Operating conditions"));
    constant Real R(unit = "J/mol.K") = 8.314;
    parameter Real A(unit = "m2", fixed = false), V_Total(unit = "m3", fixed = false);
    
    Real z[NOC] annotation(each HideResult = true), Tf(unit = "K") "Temperature of feed", F(unit = "mol") "Feed flowrate";
  
    Real P(unit = "Pa") "Pressure inside column", h(unit = "m") "Liquid level inside column";
      
    Real y[NOC](each min = 0) "molefraction in vapor phase", x[NOC](each min = 0, start = {0.7, 1e-18, 0.3, 0}) "molefraction in liquid phase";
    
    Real k[NOC] annotation(each HideResult = true);
    
    Real L(unit = "mol/s", start = 100, min = 0) "Liquid outlet flowrate", V(start = 140, min = 0) "Vapor outlet flowrate";
    
    Real Psat_T[NOC](each unit = "Pa") annotation(each HideResult = true), M[NOC](each unit = "mol") annotation(each HideResult = true);
    
    Real M_Total(unit = "mol") "Total number of moles in column", ML(unit = "mol", start = 50) "Amount of liquid in column", MG(unit = "mol", start = 0.5) "Amount of vapor in column";
    
    Real VL(unit = "m3") "Liquid volume in column", VG(unit = "m3") "vapor volume in column";
    
    Real Q(unit = "J/s") "Energy required/removed", hv[NOC](each unit = "J/mol") annotation(each HideResult = true), hl[NOC](each unit = "J/mol") annotation(each HideResult = true), Hf(unit = "J/s") "Enthalpy of feed", Hv(unit = "J/s") "Enthalpy of vapor", Hl(unit = "J/s") "Enthaply of liquid", H_M_Total(unit = "J");
  
    Real densityi[NOC](each unit = "kmol/m3") annotation(each HideResult = true);
    
    unitoperations.sensor sensor1 annotation(Placement(visible = true, transformation(origin = {2, 82}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {8.88178e-16, 82}, extent = {{-14, -14}, {14, 14}}, rotation = 0)));
    unitoperations.sensor sensor3 annotation(Placement(visible = true, transformation(origin = {54, -2}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {57, 1}, extent = {{-15, -15}, {15, 15}}, rotation = 0)));
    unitoperations.port port1 annotation(Placement(visible = true, transformation(origin = {47, -69}, extent = {{-17, -17}, {17, 17}}, rotation = 0), iconTransformation(origin = {51, -67}, extent = {{-15, -15}, {15, 15}}, rotation = 0)));
    unitoperations.port port2 annotation(Placement(visible = true, transformation(origin = {50, 62}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {54, 66}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
    unitoperations.port port3 annotation(Placement(visible = true, transformation(origin = {-54, 4}, extent = {{-18, -18}, {18, 18}}, rotation = 0), iconTransformation(origin = {-54, 2}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  
  initial equation
    h = hset;
    P = Pset * 101325;
    
    for i in 1:NOC - 1 loop
      der(M[i]) = 0;
    end for;
    
    der(M_Total) = 0;
//der(H_M_Total) = 0;
    if OverrideSizeCalculations == false then
      k_drum * 0.3048 * ((sum(x[:] .* densityi[:]) - P / (R * Ti * 1000)) * P / (R * Ti * 1000)) ^ 0.5 * 1000 * A = V;
      V_Total = A * 4 * (4 * A / 3.14) ^ 0.5;
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
    VL = ML / (sum(x[:] .* densityi[:]) * 1000);
    VG = V_Total - VL;
    P * VG = MG * R * Ti; //ideal gas law for gas phase

/*energy balance */
    Hv = V * sum(y[:] .* hv[:]);
    Hf = port3.enthalpy;
    Hl = L * sum(x[:] .* hl[:]);
    H_M_Total = ML * sum(x[:] .* hl[:]) + MG * sum(y[:] .* hv[:]);
//Hf - Hv - Hl + Q  = der(H_M_Total);
    Hf - Hv - Hl + Q = 0;

  /*Thermodynamic equations */
    sum(x[:]) = 1;
    sum(y[:]) = 1;
    y[:] = k[:] .* x[:];
    k[:] = Psat_T[:] / P;
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
  
    annotation(Icon(coordinateSystem(extent = {{-70, -140}, {70, 100}}, preserveAspectRatio = false, initialScale = 0.1), graphics = {Polygon(origin = {-1, 2}, fillColor = {170, 170, 255}, fillPattern = FillPattern.Solid, points = {{-63, -62.0013}, {-51, -78.0013}, {-33, -86.0013}, {-13, -90.0013}, {15, -90.0013}, {39, -84.0013}, {55, -76.0013}, {63, -62.0013}, {63, 71.9987}, {45, 81.9987}, {29, 87.9987}, {17, 89.9987}, {-15, 89.9987}, {-33, 85.9987}, {-49, 79.9987}, {-63, 69.9987}, {-63, -62.0013}}), Line(origin = {-1.01, -25.28}, points = {{-62.9906, -4.72122}, {-52.9906, 3.27878}, {-42.9906, -4.72122}, {-32.9906, 3.27878}, {-24.9906, -4.72122}, {-8.99059, 5.27878}, {5.00941, -4.72122}, {21.0094, 5.27878}, {33.0094, -4.72122}, {45.0094, 5.27878}, {59.0094, -4.72122}, {63.0094, 1.27878}}, smooth = Smooth.Bezier), Text(origin = {1, -113}, extent = {{-59, 21}, {59, -21}}, textString = "PT Flash")}), Diagram(coordinateSystem(extent = {{-70, -140}, {70, 100}}, preserveAspectRatio = false)), version = "", uses);
  end PTFlash;




















  model Pump
    extends compounds;
    parameter Real Final_Pressure(unit = "atm") = 25;
    parameter Real Efficiency = 0.7;
    
    Real power(unit = "kW"), den[NOC] annotation(each HideResult = true);
    
    port port1 annotation(Placement(visible = true, transformation(origin = {-81, -1}, extent = {{-17, -17}, {17, 17}}, rotation = 0), iconTransformation(origin = {-88, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    port port2 annotation(Placement(visible = true, transformation(origin = {86, 44}, extent = {{-18, -18}, {18, 18}}, rotation = 0), iconTransformation(origin = {92, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  
  equation
    for i in 1:NOC loop
      den[i] = Functions.Density(comp[i].LiqDen, comp[i].Tc, port2.temperature, port2.pressure);
    end for;
  
    port2.pressure = Final_Pressure * 101325;
    port1.moleflow = port2.moleflow;
    port1.molefrac = port2.molefrac;
    port1.temperature = port2.temperature;
    power = (port2.pressure - port1.pressure) * port2.moleflow / (sum(port1.molefrac[:] .* den[:]) * Efficiency);
  
    annotation(Icon(graphics = {Ellipse(origin = {-3, 3}, fillColor = {170, 170, 255}, fillPattern = FillPattern.Solid, extent = {{-57, 61}, {57, -61}}, endAngle = 360), Polygon(origin = {3.92, -58}, fillColor = {170, 170, 255}, fillPattern = FillPattern.Solid, points = {{-45.9243, 16}, {-71.9243, -16}, {58.0757, -16}, {32.0757, 16}, {-45.9243, 16}}), Rectangle(origin = {-51, 0}, fillColor = {170, 170, 255}, fillPattern = FillPattern.Solid, extent = {{-49, 16}, {49, -16}}), Polygon(origin = {59.21, 48}, fillColor = {170, 170, 255}, fillPattern = FillPattern.Solid, points = {{-39.2094, 12}, {-25.2094, 2}, {-13.2094, -12}, {38.7906, -12}, {38.7906, 12}, {-39.2094, 12}}), Text(origin = {3, -102}, extent = {{-75, 18}, {77, -26}}, textString = "Pump")}, coordinateSystem(extent = {{-100, -140}, {100, 100}})), Diagram(coordinateSystem(extent = {{-100, -140}, {100, 100}})), version = "", uses);
  end Pump;










  model PHFlash
    parameter Real hset(unit = "m") = 3.7 annotation(Dialog(group = "Operating conditions")), Pset(unit = "atm") = 5 annotation(Dialog(group = "Operating conditions"));
    extends compounds;
    parameter Boolean connectedToInput = false;
    parameter Boolean OverrideSizeCalculations = false annotation(Dialog(tab = "Sizing"));
    parameter Real k_drum(unit = "ft/s") = 0.3 annotation(Dialog(tab = "Sizing"));
    parameter Real Area(unit = "m2") = 4 annotation(Dialog(tab = "Sizing")), Volume(unit = "m3") = 8 annotation(Dialog(tab = "Sizing"));
    //  parameter Real Ti = 310 "K" annotation(Dialog(group = "Operating conditions"));
    constant Real R(unit = "J/mol.K") = 8.314;
    parameter Real A(unit = "m2", fixed = false), V_Total(unit = "m3", fixed = false);
    
    Real z[NOC] annotation(each HideResult = true);
    
    Real F(unit = "mol/s") "Feed flowrate", P(unit = "Pa") "Pressure inside column", h(unit = "m") "Liuid level inside column", Ti(unit = "K", start = 290) "Temperature inside column";
    
    Real y[NOC] "molefraction in vapor phase", x[NOC](start = {0.5, 1e-15, 0.5, 0}, each min = 0) "molefraction in liquid phase";
    
    Real k[NOC] annotation(each HideResult = true);
    
    Real L(unit = "mol/s", start = 0.5, min = 0), V(unit = "mol/s", start = 0.5, min = 0);
    
    Real Psat_T[NOC](each unit = "Pa") annotation(each HideResult = true);
    
    Real M[NOC](each unit = "mol") annotation(each HideResult = true), M_Total(unit = "mol") "Total number of moles in column", ML(unit = "mol", start = 50) "number of moles in liquid phase", MG(unit = "mol", start = 0.5) "number of moles in vapor phase";
    
    Real VL(unit = "m3") "Volume of liquid in column", VG(unit = "m3") "volume of vapor in column";
    
    Real Q(unit = "J/s") annotation(HideResult = true), hv[NOC](each unit = "J/mol") annotation(each HideResult = true), hl[NOC](each unit = "J/mol") annotation(each HideResult = true), Hf(unit = "J/s"), Hv(unit = "J/s"), Hl(unit = "J/s"), H_M_Total(unit = "J");
    
    Real densityi[NOC](each unit = "kmol/m3") annotation(each HideResult = true);
    
    unitoperations.sensor sensor1 annotation(Placement(visible = true, transformation(origin = {2, 82}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {8.88178e-16, 82}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
    unitoperations.sensor sensor3 annotation(Placement(visible = true, transformation(origin = {58, -2}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {48, 4}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
    unitoperations.port port1 annotation(Placement(visible = true, transformation(origin = {45, -59}, extent = {{-17, -17}, {17, 17}}, rotation = 0), iconTransformation(origin = {48, -62}, extent = {{-18, -18}, {18, 18}}, rotation = 0)));
    unitoperations.port port2 annotation(Placement(visible = true, transformation(origin = {56, 64}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {52, 68}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
    unitoperations.port port3 annotation(Placement(visible = true, transformation(origin = {-50, 0}, extent = {{-18, -18}, {18, 18}}, rotation = 0), iconTransformation(origin = {-55, 7}, extent = {{-17, -17}, {17, 17}}, rotation = 0)));
  
  initial equation
    h = hset;
    P = Pset * 101325;
  
    for i in 1:NOC - 1 loop
      der(M[i]) = 0;
    end for;
  
    der(M_Total) = 0;
    der(H_M_Total) = 0;
  
    if OverrideSizeCalculations == false then
      k_drum * 0.3048 * ((sum(x[:] .* densityi[:]) - P / (R * Ti * 1000)) * P / (R * Ti * 1000)) ^ 0.5 * 1000 * A = V;
      V_Total = A * 4 * (4 * A / 3.14) ^ 0.5;
    else
      A = Area;
      V_Total = Volume;
    end if;
  
  equation
    F = port3.moleflow;
    z[:] = port3.molefrac[:];
  
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
    VL = ML / (sum(x[:] .* densityi[:]) * 1000);
    VG = V_Total - VL;
    P * VG = MG * R * Ti; //ideal gas law for gas phase

/*energy balance */
    Hv = V * sum(y[:] .* hv[:]);
    Hf = port3.enthalpy;
    Hl = L * sum(x[:] .* hl[:]);
    H_M_Total = ML * sum(x[:] .* hl[:]) + MG * sum(y[:] .* hv[:]);
    Hf - Hv - Hl = der(H_M_Total);
    //Hf - Hv - Hl + Q = 0;
    Q = 0;

  /*Thermodynamic equations */
    sum(x[:]) = 1;
    sum(y[:]) = 1;
    y[:] = k[:] .* x[:];
    k[:] = Psat_T[:] / P;
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
  
    annotation(Icon(coordinateSystem(extent = {{-70, -140}, {70, 100}}, preserveAspectRatio = false, initialScale = 0.1), graphics = {Polygon(origin = {-1, 2}, fillColor = {170, 170, 255}, fillPattern = FillPattern.Solid, points = {{-63, -62.0013}, {-51, -78.0013}, {-33, -86.0013}, {-13, -90.0013}, {15, -90.0013}, {39, -84.0013}, {55, -76.0013}, {63, -62.0013}, {63, 71.9987}, {45, 81.9987}, {29, 87.9987}, {17, 89.9987}, {-15, 89.9987}, {-33, 85.9987}, {-49, 79.9987}, {-63, 69.9987}, {-63, -62.0013}}), Line(origin = {-1.01, -25.28}, points = {{-62.9906, -4.72122}, {-52.9906, 3.27878}, {-42.9906, -4.72122}, {-32.9906, 3.27878}, {-24.9906, -4.72122}, {-8.99059, 5.27878}, {5.00941, -4.72122}, {21.0094, 5.27878}, {33.0094, -4.72122}, {45.0094, 5.27878}, {59.0094, -4.72122}, {63.0094, 1.27878}}, smooth = Smooth.Bezier), Text(origin = {-4, -113}, extent = {{-60, 23}, {60, -23}}, textString = "PH Flash")}), Diagram(coordinateSystem(extent = {{-70, -140}, {70, 100}}, preserveAspectRatio = false)), version = "", uses);
  end PHFlash;





















  model Mixer
    extends compounds;
    
    Real F, z[NOC], Fl(start = 37), zl[NOC](start = {1,0,1e-3,0}, each nominal = 1e-2), Fv, zv[NOC](each nominal = 1e-2);
    Real H, Hl[NOC], Hv[NOC], Hl_Tb[NOC], Hv_Td[NOC], H_Tb, H_Td;
    Real T(start = 315, nominal = 1e2), Tb(start = 350), Td(start = 350);
    Real P, Psat_T[NOC], Psat_Tb[NOC], Psat_Td[NOC], k[NOC];
    
    unitoperations.port port1 annotation(Placement(visible = true, transformation(origin = {-90, 90}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-87, 79}, extent = {{-21, -21}, {21, 21}}, rotation = 0)));
    unitoperations.port port2 annotation(Placement(visible = true, transformation(origin = {-90, -92}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-87,-81}, extent = {{-19, -19}, {19, 19}}, rotation = 0)));
    unitoperations.port port3 annotation(Placement(visible = true, transformation(origin = {90, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {88, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  
  equation
  /*  port3.moleflow = port1.moleflow + port2.moleflow;
    port3.pressure = min({port1.pressure, port2.pressure});
    port3.molefrac[:] = (port1.molefrac[:] * port1.moleflow + port2.molefrac[:] * port2.moleflow) / port3.moleflow;
    port3.temperature = min({port1.temperature, port2.temperature});
    port3.liquidmoleflow = 0;
    port3.vapormoleflow = 0;
    port3.liquidmolefrac = zeros(NOC);
    port3.vapormolefrac = zeros(NOC);
    port3.enthalpy = 0; */
    port3.moleflow = port1.moleflow + port2.moleflow;
    port3.pressure = min({port1.pressure, port2.pressure});
    port3.molefrac[:] = (port1.molefrac[:] * port1.moleflow + port2.molefrac[:] * port2.moleflow) / port3.moleflow;
    port3.temperature = T;
    port3.molefrac = z;
    port3.pressure = P;
    port3.moleflow = F;
    H = port1.enthalpy +  port2.enthalpy;
  
    for i in 1:NOC loop
     Psat_Tb[i] = Functions.Psat(comp[i].VP, Tb);
     Psat_Td[i] = Functions.Psat(comp[i].VP, Td);
     Hl_Tb[i] = Functions.HLiqId(comp[i].VapCp, comp[i].HOV, comp[i].Tc, Tb);
     Hv_Td[i] = Functions.HVapId(comp[i].VapCp, comp[i].HOV, comp[i].Tc, Td);
    end for;
  
    sum(P * z[:] ./ Psat_Td[:]) = 1;
    sum(Psat_Tb[:] ./ P .* z[:]) = 1;
    H_Tb = F * sum(z[:] .* Hl_Tb[:]);
    H_Td = F * sum(z[:] .* Hv_Td[:]);
  
    if H <= H_Tb then
    k = zeros(NOC);
    zl = z;
    zv = zeros(NOC);
    F = Fl;
    Fv = 0;
    H = Fl * sum(zl[:] .* Hl[:]);
    elseif H >= H_Td then
    k = zeros(NOC);
    zl = zeros(NOC);
    zv = z;
    Fl = 0;
    Fv = F;
    H = Fv * sum(zv[:] .* Hv[:]);
    else 
    sum(zl[:]) = 1;
    sum(zv[:]) = 1;
    zv[:] = k[:] .* zl[:];
    k[:] = Psat_T[:] ./ P;
    F = Fl + Fv;
  
    for i in 1:NOC - 1 loop
      F * z[i] = Fl * zl[i] + Fv * zv[i];
    end for;
  
    H = Fl * sum(zl[:] .* Hl[:]) + Fv * sum(zv[:] .* Hv[:]);
    end if;
  
    for i in 1:NOC loop
     Psat_T[i] = Functions.Psat(comp[i].VP, T);
     Hl[i] = Functions.HLiqId(comp[i].VapCp, comp[i].HOV, comp[i].Tc, T);
     Hv[i] = Functions.HVapId(comp[i].VapCp, comp[i].HOV, comp[i].Tc, T);
    end for;
  
    annotation(Icon(graphics = {Polygon(origin = {-0.28, 0}, fillColor = {170, 170, 255}, fillPattern = FillPattern.Solid, points = {{-99.7236, 100}, {100.276, 0}, {-99.7236, -100}, {-99.7236, 100}}), Text(origin = {-1, -107}, extent = {{-81, 19}, {85, -27}}, textString = "Mixer")}, coordinateSystem(extent = {{-100, -140}, {100, 100}})), Diagram(coordinateSystem(extent = {{-100, -140}, {100, 100}})), version = "", uses);
  end Mixer;



end unitoperations;
