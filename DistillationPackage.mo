package DistillationPackage
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
    sensor sensor1 annotation(Placement(visible = true, transformation(origin = {1, 77}, extent = {{-19, -19}, {19, 19}}, rotation = 0), iconTransformation(origin = {2, 56}, extent = {{-18, -18}, {18, 18}}, rotation = 0)));
    port port1 annotation(Placement(visible = true, transformation(origin = {-80, -4}, extent = {{-18, -18}, {18, 18}}, rotation = 0), iconTransformation(origin = {-82, -2}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
    port port2 annotation(Placement(visible = true, transformation(origin = {79, 1}, extent = {{-17, -17}, {17, 17}}, rotation = 0), iconTransformation(origin = {82, 2}, extent = {{-18, -18}, {18, 18}}, rotation = 0)));
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
    parameter Integer NOC = 4;
    parameter Chemsep_Database.Toluene comp1;
    parameter Chemsep_Database.Hydrogen comp2;
    parameter Chemsep_Database.Benzene comp3;
    parameter Chemsep_Database.Methane comp4;
    parameter String Name = "MS1";
    parameter Chemsep_Database.General_Properties comp[NOC] = {comp1, comp2, comp3, comp4};
    parameter Real Flowrate = 100, Pressure = 1e5, Temperature = 300, molefraction[NOC] = zeros(NOC);
    parameter Boolean unspecified = true;
    parameter Boolean stepchange = false;
    parameter Real stepchangetime = 0.01;
    parameter Real step_value = 1;
    Real kf[NOC], zl[NOC](each min = 0, each max = 1, start = {0.5, 1e-18, 0.5, 0}), zv[NOC](each min = 0, each max = 1, start = {0, 0.25, 0, 0.75}), Fl(min = 0, start = 100), Fv(min = 0, start = 140), Tbf(start = 62), Tdf(start = 495.5), Psat_Tdf[NOC], Psat_Tbf[NOC], Psat_Tf[NOC], Pf, Tf, z[NOC], F;
    Real Hvf[NOC], Hlf[NOC], H;
    port port2 annotation(Placement(visible = true, transformation(origin = {80, -4}, extent = {{-18, -18}, {18, 18}}, rotation = 0), iconTransformation(origin = {85, 1}, extent = {{-21, -21}, {21, 21}}, rotation = 0)));
    port port1 annotation(Placement(visible = true, transformation(origin = {-80, -4}, extent = {{-18, -18}, {18, 18}}, rotation = 0), iconTransformation(origin = {-82, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
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

  model Distillation
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
      xc[1] * D = specification1_value;
//yet to modify
    elseif Integer(specification1) == 4 then
      xc[1] = specification1_value;
//yet to modify
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
      xr[1] * D = specification2_value;
//yet to modify
    elseif Integer(specification2) == 4 then
      xc[1] = specification2_value;
//yet to modify
    elseif Integer(specification2) == 5 then
      VNT = specification2_value * B;
    else
      TB = specification2_value;
    end if;
    annotation(Icon(graphics = {Rectangle(origin = {-2, -1}, fillColor = {0, 85, 127}, fillPattern = FillPattern.VerticalCylinder, extent = {{-94, 95}, {94, -95}})}), Documentation(info = "<HTML> <p> This is a generalized model for distilation column </p> </HTML>"));
  end Distillation;

  model test1
    Distillation distillation1(specification1 = DistillationPackage.Distillation.spec1.RefluxRatio, specification1_value = 2, specification2 = DistillationPackage.Distillation.spec2.ProductMolarFlow, specification2_value = 61.5) annotation(Placement(visible = true, transformation(origin = {-3, -7}, extent = {{-19, -19}, {19, 19}}, rotation = 0)));
    MaterialStream materialStream1(Flowrate = 97.5, Pressure = 1e5, Temperature = 360, molefraction = {0.45, 0, 0.55, 0}, stepchange = false, unspecified = false) annotation(Placement(visible = true, transformation(origin = {-70, -8}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    MaterialStream materialStream2 annotation(Placement(visible = true, transformation(origin = {50, 16}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(distillation1.port2, materialStream2.port1) annotation(Line(points = {{14, 8}, {24, 8}, {24, 18}, {42, 18}, {42, 16}}));
    connect(materialStream1.port2, distillation1.port1) annotation(Line(points = {{-62, -8}, {-20, -8}, {-20, -6}, {-20, -6}}));
  end test1;

  model DistillationC
    parameter Chemsep_Database.Benzene comp1 annotation(Dialog(tab = "General", group = "compounds"));
    parameter Chemsep_Database.Toluene comp2 annotation(Dialog(tab = "General", group = "compounds"));
    parameter Chemsep_Database.General_Properties comp[2] = {comp1, comp2} annotation(Dialog(tab = "General", group = "compounds"));
    parameter Integer N_Trays = 20 "No. of trays without condensor and reboiler" annotation(Dialog(tab = "General", group = "Trays"));
    parameter Integer NOC = 2 "No. of compounds" annotation(Dialog(tab = "General", group = "compounds"));
    parameter Integer N_Feed = 10 "Feed tray location" annotation(Dialog(tab = "General", group = "Trays"));
    parameter Real M_eff[N_Trays, NOC] = fill(1, N_Trays, NOC) "Murfy's efficiency of trays" annotation(Dialog(tab = "General", group = "Efficiency"));
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
    /*variables*/
    Real y[N_Trays, NOC](start = fill(0.5, N_Trays, NOC)) "vapor molefraction", x[N_Trays, NOC](start = fill(0.5, N_Trays, NOC)) "liquid molefraction", y_eq[N_Trays, NOC] "equilibrium vapor molefraction", Tf[N_Trays] "temperature of feed";
    Real V[N_Trays](each start = 70) "vapor flowrates", L[N_Trays](each start = 100) "Liquid flowrates";
    Real T[N_Trays](start = linspace(386, 354, N_Trays)) "Temperature of trays", TC(start = 368) "Temperature of condenser", TB(start = 377) "Temeprature of reboiler", L0(start = 50) "", VNT "Boil up", D "Distillate", B "Bottoms", xc[NOC] "composition in condenser", xr[NOC] "composition in reboiler", QC "Condenser heat load", QB "Reboiler heat load";
    Real Keq[N_Trays, NOC], Psat[N_Trays, NOC], PsatC[NOC], PsatB[NOC];
    Real yNT[NOC], M[N_Trays](each start = 1) "moles of liquid in tray", den[N_Trays, NOC] "density of component in tray";
    Real hv[N_Trays, NOC] annotation(each HideResult = true), hl[N_Trays, NOC] annotation(each HideResult = true), hf[N_Trays, NOC] annotation(each HideResult = true), hv_B[NOC] annotation(each HideResult = true), hl_B[NOC] annotation(each HideResult = true), hl_C[NOC] annotation(each HideResult = true);
    Real P[N_Trays], Ks[N_Trays] annotation(each HideResult = true);
    Real dx[N_Trays, NOC], dM[N_Trays], dhl[N_Trays, NOC];
    Real F[N_Trays](each start = 0) annotation(each HideResult = true), z[N_Trays, NOC](start = fill(0.5, N_Trays, NOC)) annotation(each HideResult = true);
    /*ports*/
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
//  port1.pressure = P[N_Feed];
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
//equation for moles on trays
    for i in 1:N_Trays loop
      M[i] = A_active * sum(x[i, :] .* den[i, :]) * (exp(pi1 * Ks[i] ^ pi2) * h_weir + 44300 * 10 ^ (-3) * (L[i] / (sum(x[i, :] .* den[i, :]) * d_weir * 1000 * 3600)) ^ 0.704);
    end for;
//mass balance for tray 1
    M[1] * dx[1, :] + x[1, :] * dM[1] = VNT .* yNT[:] + L[2] .* x[2, :] - V[1] .* y[1, :] - L[1] .* x[1, :] + F[1] .* z[1, :];
//mass balance for trays in between
    for i in 2:N_Trays - 1 loop
      dM[i] * x[i, :] + dx[i, :] * M[i] = V[i - 1] .* y[i - 1, :] + L[i + 1] .* x[i + 1, :] - V[i] .* y[i, :] - L[i] .* x[i, :] + F[i] .* z[i, :];
    end for;
//mass balance for top tray
    M[N_Trays] * dx[N_Trays, :] + x[N_Trays, :] * dM[N_Trays] = V[N_Trays - 1] .* y[N_Trays - 1, :] + L0 .* xc[:] - V[N_Trays] .* y[N_Trays, :] - L[N_Trays] .* x[N_Trays, :] + F[N_Trays] .* z[N_Trays, :];
//mass balance for condenser
    V[N_Trays] - L0 - D = 0;
    y[N_Trays, :] = xc[:];
//mass balance for reboiler
    xr[:] = x[1, :];
    yNT[:] = xr[:];
    L[1] - VNT - B = 0;
//energy balance
//at tray 1 (from bottom)
    VNT * sum(yNT[:] .* hv_B[:]) - V[1] * sum(y[1, :] .* hv[1, :]) + L[2] * sum(x[2, :] .* hl[2, :]) - L[1] * sum(x[1, :] .* hl[1, :]) + F[1] * sum(z[1, :] .* hf[1, :]) = M[1] * x[1, :] * dhl[1, :] + M[1] * dx[1, :] * hl[1, :] + dM[1] * x[1, :] * hl[1, :];
// for trays in between
    for i in 2:N_Trays - 1 loop
      V[i - 1] * sum(y[i - 1, :] .* hv[i - 1, :]) - V[i] * sum(y[i, :] .* hv[i, :]) + L[i + 1] * sum(x[i + 1, :] .* hl[i + 1, :]) - L[i] * sum(x[i, :] .* hl[i, :]) + F[i] * sum(z[i, :] .* hf[i, :]) = M[i] * x[i, :] * dhl[i, :] + M[i] * dx[i, :] * hl[i, :] + dM[i] * x[i, :] * hl[i, :];
    end for;
//at top tray
    V[N_Trays - 1] * sum(y[N_Trays - 1, :] .* hv[N_Trays - 1, :]) - V[N_Trays] * sum(y[N_Trays, :] .* hv[N_Trays, :]) + L0 * sum(xc[:] .* hl_C[:]) - L[N_Trays] * sum(x[N_Trays, :] .* hl[N_Trays, :]) + F[N_Trays] * sum(z[N_Trays, :] .* hf[N_Trays, :]) = dM[N_Trays] * x[N_Trays, :] * hl[N_Trays, :] + M[N_Trays] * dx[N_Trays, :] * hl[N_Trays, :] + M[N_Trays] * x[N_Trays, :] * dhl[N_Trays, :];
//condenser
    V[N_Trays] * sum(y[N_Trays, :] .* hv[N_Trays, :]) - (L0 + D) * sum(xc[:] .* hl_C[:]) = QC;
//reboiler
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
/*sending data to ports*/
    port2.moleflow = D;
    port2.molefrac = {xc[2], 0, xc[1], 0};
    port2.temperature = TC;
    port2.pressure = P_condenser;
    port3.moleflow = B;
    port3.molefrac = {xr[2], 0, xr[1], 0};
    port3.temperature = TB;
    port3.pressure = P_condenser + Pressure_drop * (N_Trays + 1);
//Equations for Specification
    if Integer(specification1) == 1 then
      QC = specification1_value;
    elseif Integer(specification1) == 2 then
      D = specification1_value;
    elseif Integer(specification1) == 3 then
      xc[1] * D = specification1_value;
//yet to modify
    elseif Integer(specification1) == 4 then
      xc[1] = specification1_value;
//yet to modify
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
      xr[1] * D = specification2_value;
//yet to modify
    elseif Integer(specification2) == 4 then
      xc[1] = specification2_value;
//yet to modify
    elseif Integer(specification2) == 5 then
      VNT = specification2_value * B;
    else
      TB = specification2_value;
    end if;
    annotation(Icon(graphics = {Rectangle(origin = {-2, -1}, fillColor = {0, 85, 127}, fillPattern = FillPattern.VerticalCylinder, extent = {{-94, 95}, {94, -95}})}));
  end DistillationC;

  model Distillationwithsizing
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
    parameter Real Active_area = 1 "Active area of Tray, sq.m" annotation(Dialog(tab = "dynamic", group = "Tray data"));
    parameter Real A_active(fixed = false);
    parameter Real d_weir(fixed = false);
    parameter Boolean Override_Sizing_Calculations(start = false);
    parameter Real Kv = 0.3 "constant for calculating max velocity permissible, ft/s" annotation(Dialog(tab = "dynamic", group = "Tray data"));
    parameter Real h_weir = 0.01 "Height of weir, m" annotation(Dialog(tab = "dynamic", group = "Tray data"));
    parameter Real Weir_diameter = 0.8 "Diameter of weir, m" annotation(Dialog(tab = "dynamic", group = "Tray data"));
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
    unitoperations.port port2 annotation(Placement(visible = true, transformation(origin = {92, 86}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {92, 86}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  initial equation
/*sizing*/
    if Override_Sizing_Calculations == false then
      Kv * 0.3048 * ((den[1, :] * x[1, :] - P[1] / (R * T[1] * 1000)) * (P[1] / (R * T[1] * 1000))) ^ 0.5 * A_active = max(V) / 3600;
      d_weir = (A_active * 4 / 3.14) ^ 0.5;
    else
      A_active = Active_area;
      d_weir = Weir_diameter;
    end if;
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
    for i in 1:N_Trays loop
      M[i] = A_active * sum(x[i, :] .* den[i, :]) * (exp(pi1 * Ks[i] ^ pi2) * h_weir + 44300 * 10 ^ (-3) * (L[i] / (sum(x[i, :] .* den[i, :]) * d_weir * 1000 * 3600)) ^ 0.704);
    end for;
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
      xc[1] * D = specification1_value;
//yet to modify
    elseif Integer(specification1) == 4 then
      xc[1] = specification1_value;
//yet to modify
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
      xr[1] * D = specification2_value;
//yet to modify
    elseif Integer(specification2) == 4 then
      xc[1] = specification2_value;
//yet to modify
    elseif Integer(specification2) == 5 then
      VNT = specification2_value * B;
    else
      TB = specification2_value;
    end if;
    annotation(Icon(graphics = {Rectangle(origin = {-2, -1}, fillColor = {0, 85, 127}, fillPattern = FillPattern.VerticalCylinder, extent = {{-94, 95}, {94, -95}})}), Documentation(info = "<HTML> <p> This is a generalized model for distilation column </p> </HTML>"));
  end Distillationwithsizing;

  model Distillationwithsizingtest
    Distillationwithsizing distillation1(specification1 = DistillationPackage.Distillation.spec1.RefluxRatio, specification1_value = 2, specification2 = DistillationPackage.Distillation.spec2.ProductMolarFlow, specification2_value = 61.1) annotation(Placement(visible = true, transformation(origin = {-3, -7}, extent = {{-19, -19}, {19, 19}}, rotation = 0)));
    MaterialStream materialStream1(Flowrate = 95.7, Pressure = 1e5, Tdf(start = 370), Temperature = 360, molefraction = {0.45, 0, 0.55, 0}, stepchange = false, unspecified = false) annotation(Placement(visible = true, transformation(origin = {-70, -8}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    unitoperations.MaterialStream materialStream2(Fl(start = 34.6), Fv(start = 0), Tbf(start = 353), Tdf(start = 353.3), zl(start = {0.001, 0, 0.9986, 0}), zv(start = {0.25, 0.25, 0.25, 0.25})) annotation(Placement(visible = true, transformation(origin = {54, 18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(distillation1.port2, materialStream2.port1) annotation(Line(points = {{14, 10}, {44, 10}, {44, 18}, {46, 18}}));
    connect(materialStream1.port2, distillation1.port1) annotation(Line(points = {{-62, -8}, {-20, -8}, {-20, -6}, {-20, -6}}));
    annotation(experiment(StartTime = 0, StopTime = 0, Tolerance = 1e-06, Interval = 0));
  end Distillationwithsizingtest;

  model materialtest
    MaterialStream materialStream1(Flowrate = 61.1, Pressure = 152422, Temperature = 391.5, molefraction = {0.704, 0, 0.296, 0}, unspecified = false) annotation(Placement(visible = true, transformation(origin = {-22, 8}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    annotation(experiment(StartTime = 0, StopTime = 0, Tolerance = 1e-06, Interval = 0));
  end materialtest;
end DistillationPackage;