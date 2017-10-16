package modelicatest
  connector port
    parameter Integer NOC = 4;
    Real moleflow;
    parameter Boolean datagiven(fixed = false);
    annotation(
      Diagram(graphics = {Ellipse(origin = {-1, 1}, fillColor = {65, 252, 255}, fillPattern = FillPattern.Solid, extent = {{-67, 65}, {67, -65}}, endAngle = 360)}),
      Icon(graphics = {Ellipse(origin = {-4, -1}, fillColor = {81, 253, 248}, fillPattern = FillPattern.Solid, extent = {{-60, 61}, {60, -61}}, endAngle = 360)}));
  end port;

  model portparametertest
    parameter Boolean unspecified = false;
    parameter Boolean check(fixed = false);
    port port1 annotation(
      Placement(visible = true, transformation(origin = {-90, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-90, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    port port2 annotation(
      Placement(visible = true, transformation(origin = {90, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {90, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Real a;
  initial equation
    port1.datagiven = if unspecified == true then false else true;
    port1.datagiven = port2.datagiven;
    check = port1.datagiven;
  equation
    connect(port1, port2) annotation(
      Line(points = {{-90, 0}, {92, 0}, {92, 2}, {90, 2}}));
    port1.moleflow = 20;
    if check == true then
      a = 10;
    end if;
  end portparametertest;

  model phFlash
    extends unitoperations.compounds;
    parameter Real V_Total = 8.11, F(fixed = false), z[NOC](each fixed = false), Hf(fixed = false), P = 1e5;
    parameter Real R = 8.314, A = 2, Q = 0;
    Real h "liquid level";
    unitoperations.MaterialStream materialStream1(Flowrate = 1, Pressure = 10e5, molefraction = {0.25, 0.25, 0.25, 0.25}, unspecified = false) annotation(
      Placement(visible = true, transformation(origin = {-52, 8}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Real Psat_T[NOC], k[NOC], hv[NOC], hl[NOC], densityi[NOC], L(start = 0.50), V(start = 0.50), x[NOC], y[NOC], M[NOC], ML, MG, M_Total;
    Real Ti(start = 290);
    Real VL, VG, Hv, Hl, H_M_Total;
  initial equation
    F = materialStream1.port2.moleflow;
    z = materialStream1.port2.molefrac;
    Hf = materialStream1.port2.enthalpy;
/*  h = hset;
    P = Pset;
    for i in 1:NOC - 1 loop
      der(M[i]) = 0;
    end for;
    der(M_Total) = 0;*/
//der(H_M_Total) = 0;
  equation
    for i in 1:NOC loop
      Psat_T[i] = Functions.Psat(comp[i].VP, Ti);
      hv[i] = Functions.HVapId1(comp[i].VapCp, comp[i].HOV, comp[i].Tc, Ti);
      hl[i] = Functions.HLiqId(comp[i].VapCp, comp[i].HOV, comp[i].Tc, Ti);
      densityi[i] = Functions.Density(comp[i].LiqDen, comp[i].Tc, Ti, P);
    end for;
/* Mass Blanace equations */
//der(M_Total) = F - L - V;
    F - L - V = 0;
    for i in 1:NOC - 1 loop
// der(M[i]) = F * z[i] - L * x[i] - V * y[i];
      F * z[i] - L * x[i] - V * y[i] = 0;
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
//  Hf = materialStream1.port2.enthalpy;
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
    h = 2;
  end phFlash;

  model ms
    unitoperations.MaterialStream materialStream1(Flowrate = 43.0402, Pressure = 1, Temperature = 310, molefraction = {0.5059, 1e-16, 0.4911, 0.003}, specified_stream = true) annotation(
      Placement(visible = true, transformation(origin = {-12, 4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  end ms;

  model phFlashDynamic
    extends unitoperations.compounds;
    parameter Real V_Total = 8.11;
    parameter Real hset = 2, Pset = 5e5, Pout = 1e5;
    parameter Real R = 8.314, A = 2, Q = 0;
    Real Cv1, Cv2;
    Real h "liquid level";
    unitoperations.MaterialStream materialStream1(Pressure = 10e5, molefraction = {0.3, 0.3, 0.2, 0.2}, unspecified = false) annotation(
      Placement(visible = true, transformation(origin = {-52, 8}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Real F, z[NOC], Hf;
    Real Psat_T[NOC], k[NOC], hv[NOC], hl[NOC], densityi[NOC], L(start = 50), V(start = 50), x[NOC](each start = 0.25), y[NOC], M[NOC], ML, MG, M_Total;
    Real Ti(start = 290), P;
    Real VL, VG, Hv, Hl, H_M_Total;
  initial equation
    h = hset;
    P = Pset;
    for i in 1:NOC - 1 loop
      der(M[i]) = 0;
    end for;
    der(M_Total) = 0;
    der(H_M_Total) = 0;
  equation
    F = materialStream1.port2.moleflow;
    z = materialStream1.port2.molefrac;
    Hf = materialStream1.port2.enthalpy;
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
    Hl = L * sum(x[:] .* hl[:]);
    H_M_Total = ML * sum(x[:] .* hl[:]) + MG * sum(y[:] .* hv[:]);
    Hf - Hv - Hl = der(H_M_Total);
//  Hf - Hv - Hl + Q = 0;
/*Thermodynamic equations */
    sum(x[:]) = 1;
    sum(y[:]) = 1;
    y[:] = k[:] .* x[:];
    k[:] = Psat_T[:] / P;
/* Control */
    A * h = VL;
// h = 2;
    L = Cv1 * 0.4 * (P - Pout) ^ 0.5;
    V = Cv2 * 0.4 * (P - Pout) ^ 0.5;
    der(Cv2) = 0;
    der(Cv1) = 0;
  end phFlashDynamic;

  model parameterTest
    parameter Real a(fixed = false);
    Real b, c;
  initial equation
    a = 2 * b;
  equation
    b = 4;
    c = a * b;
  end parameterTest;

  model matrixtest
    parameter Real a[4, 2] = ones(4, 2), b[4, 2] = [1, 2; 3, 4; 5, 6; 7, 8];
    Real c;
  equation
    c = a[1, :] * b[1, :];
  end matrixtest;

  type compounds = enumeration(Air, Argon, Bromine, Carbontetrachloride, Carbonmonoxide, Carbondioxide, Carbondisulfide, Phosgene, Trichloroacetylchloride, Hydrogenchloride, Chlorine, Hydrogeniodide, Hydrogen, Water, Hydrogensulfide, Ammonia, Neon, Nitricacid, Nitricoxide, Nitrogendioxide, Nitrogen, Nitrousoxide, Oxygen, Sulfurdioxide, Sulfurtrioxide, Chloroform, Hydrogencyanide, Formaldehyde, Methylchloride, Methyliodide, Methane, Methanol, Methylamine, Trichloroethylene, Dichloroacetylchloride, Trichloroacetaldehyde, Acetylene, Dichloroacetaldehyde, Vinylchloride, Acetylchloride, OneOneTwotrichloroethane, Acetonitrile, Ethylene, OneOnedichloroethane, OneTwodichloroethane, Acetaldehyde, Ethyleneoxide, Aceticacid, Methylformate, Ethylchloride, Ethane, Ethanol, Dimethylether, Ethyleneglycol, Dimethylsulfide, Ethylmercaptan, Ethylamine, Acrylonitrile, Methylacetylene, Propadiene, Propylene, Acetone, Ethylformate, Methylacetate, Propionicacid, Nndimethylformamide, Propane, Isopropanol, Onepropanol, Trimethylamine, Vinylacetylene, Thiophene, Methacrylonitrile, Dimethylacetylene, Ethylacetylene, OneTwobutadiene, OneThreebutadiene, Onebutene, CisTwobutene, TransTwobutene, Isobutene, Twomethylpropanal, Methylethylketone, Tetrahydrofuran, OneFourdioxane, Nbutyricacid, Ethylacetate, Methylpropionate, Npropylformate, Sulfolane, Nndimethylacetamide, Nbutane, Isobutane, Onebutanol, TwomethylOnepropanol, Twobutanol, TwomethylTwopropanol, Diethylether, Diethyleneglycol, Diethylamine, Furfural, Pyridine, Isoprene, Cyclopentane, TwomethylOnebutene, ThreemethylOnebutene, TwomethylTwobutene, Onepentene, CisTwopentene, TransTwopentene, Threepentanone, Methylisopropylketone, Npropylacetate, Isopentane, Npentane, Neopentane, OneTwoFourtrichlorobenzene, Mdichlorobenzene, Odichlorobenzene, Pdichlorobenzene, Bromobenzene, Monochlorobenzene, Iodobenzene, Nitrobenzene, Benzene, Phenol, Aniline, Cyclohexanone, Cyclohexane, Onehexene, Methylcyclopentane, Cyclohexanol, TwoTwodimethylbutane, TwoThreedimethylbutane, Nhexane, Twomethylpentane, Threemethylpentane, Triethyleneglycol, Triethylamine, Toluene, Mcresol, Ocresol, Pcresol, Methylcyclohexane, Ethylcyclopentane, Oneheptene, Nheptane, Styrene, Ethylbenzene, Mxylene, Oxylene, Pxylene, Ethylcyclohexane, Npropylcyclopentane, Noctane, TwoTwoThreetrimethylpentane, TwoTwoFourtrimethylpentane, TwoThreeThreetrimethylpentane, TwoThreeFourtrimethylpentane, Tetraethyleneglycol, Indene, Indane, Cumene, Npropylbenzene, Npropylcyclohexane, Nnonane, Naphthalene, Onemethylindene, Twomethylindene, Dicyclopentadiene, Nbutylbenzene, Nbutylcyclohexane, Ndecane, Onemethylnaphthalene, Twomethylnaphthalene, Nundecane, Acenaphthene, Biphenyl, Ndodecane, Fluorene, Ntridecane, Phenanthrene, Ntetradecane, Npentadecane, Fluoranthene, Pyrene, Onephenylnaphthalene, Nhexadecane, Chrysene, Cisdecahydronaphthalene, Transdecahydronaphthalene, Methyltertbutylether, Methyltertpentylether, TwomethylTwobutanol, Nitrogentrioxide, Nitrogentetroxide, HeliumFour, Fluorine, Krypton, Xenon, Ozone, Carbonylsulfide, Sulfurhexafluoride, Dimethylsulfoxide, Nheptadecane, Noctadecane, Nnonadecane, Nheneicosane, Ndocosane, Ntricosane, Ntetracosane, Npentacosane, Nhexacosane, Nheptacosane, Noctacosane, Nnonacosane, Squalane);

  model compoundsloadingfile
    parameter Integer NOC = 4;
    parameter compounds comp1;
    parameter Chemsep_Database.Toluene comp2;
    annotation(
      Diagram(graphics = {Rectangle(origin = {0, -3}, fillColor = {255, 85, 255}, fillPattern = FillPattern.Solid, extent = {{-88, 79}, {88, -79}})}),
      Icon(graphics = {Rectangle(origin = {-2, -7}, fillColor = {255, 85, 255}, fillPattern = FillPattern.Solid, extent = {{-94, 81}, {94, -81}})}));
  end compoundsloadingfile;

  model compoundstest
    compoundsloadingfile compoundsloadingfile1(comp1 = modelicatest.compounds.Air) annotation(
      Placement(visible = true, transformation(origin = {0, -8}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  end compoundstest;

  model CSTR
    extends unitoperations.compounds;
    parameter Real s[NOC] = {-1, -1, 1, 1} "stoichiometric coefficients";
    parameter Integer NOIS = 2 "No of input streams";
    parameter Integer n = 1 "base compound identity in comp[n]";
    parameter Real V_Total = 2.321 "Volume of reactor", P_init = 25e5 "Pressure at t = 0";
    constant Real R = 8.314;
    parameter Real Af "frequency factor for forward reaction" annotation(
      Dialog(tab = "Reactions", group = "forward reaction rate constants"));
    parameter Real order_f[4] "order wrt to components for forward reaction" annotation(
      Dialog(tab = "Reactions", group = "forward reaction rate constants"));
    parameter Real order_b[4] "order wrt to components for backward reaction" annotation(
      Dialog(tab = "Reactions", group = "backward reaction rate constants"));
    parameter Real Ab "frequency factor for backward reaction" annotation(
      Dialog(tab = "Reactions", group = "backward reaction rate constants"));
    parameter Real Eaf "Activation energy for forward reaction" annotation(
      Dialog(tab = "Reactions", group = "forward reaction rate constants"));
    parameter Real Eab "Activation energy for backward reaction" annotation(
      Dialog(tab = "Reactions", group = "backward reaction rate constants"));
    parameter Real delH_r = 12.6e3 "Heat of reaction" annotation(
      Dialog(tab = "Reactions", group = "Reaction rate constants"));
    type operation_type = enumeration(Isothermal, Adiabatic);
    parameter operation_type operation_mode;
    parameter Real T_iso = 900;
    Real F_in[NOIS], z[NOIS, NOC], Hin[NOIS] "Feed values";
    Real M_Total(start = 1.5), M[NOC], x[NOC], F_out, densityi[NOC], P;
    Real r, kf, kb, c[NOC];
    Real H_r, Hout, Q, T(start = 650);
    unitoperations.port port1 annotation(
      Placement(visible = true, transformation(origin = {-82, 2}, extent = {{-18, -18}, {18, 18}}, rotation = 0), iconTransformation(origin = {-85, 1}, extent = {{-17, -17}, {17, 17}}, rotation = 0)));
    unitoperations.port port2 annotation(
      Placement(visible = true, transformation(origin = {2, 84}, extent = {{-18, -18}, {18, 18}}, rotation = 0), iconTransformation(origin = {3, 83}, extent = {{-17, -17}, {17, 17}}, rotation = 0)));
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
      c[i] = M[i] / V_Total;
    end for;
//reaction rate
    kf = Af * exp(-Eaf / (R * T));
    kb = Ab * exp(-Eab / (R * T));
    r = kf * product(c[:] .^ order_f[:]) - kb * product(c[:] .^ order_b[:]);
//unsteady state mass balance
    for i in 1:NOC - 1 loop
      der(M[i]) = sum(z[:, i] .* F_in[:]) - x[i] * F_out + s[i] * r * V_Total / abs(s[n]);
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
//Hout = port3.enthalpy;
    Hout = 0;
    annotation(
      Icon(graphics = {Rectangle(origin = {-1, 1}, fillColor = {255, 170, 0}, fillPattern = FillPattern.Solid, extent = {{-97, 95}, {97, -95}}), Text(origin = {2, 0}, extent = {{-80, 62}, {80, -62}}, textString = "CSTR")}, coordinateSystem(initialScale = 0.1)));
  end CSTR;

  model PFR
  end PFR;

  model CSTRtest
    MaterialStream materialStream1(Flowrate = 66, Pressure = 22e5, Temperature = 800, molefraction = {0, 1, 0, 0}, unspecified = false) annotation(
      Placement(visible = true, transformation(origin = {-78, 4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    MaterialStream materialStream2(Flowrate = 22.2, Pressure = 22e5, Temperature = 800, molefraction = {1, 0, 0, 0}, unspecified = false) annotation(
      Placement(visible = true, transformation(origin = {-72, -44}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    CSTR cSTR1(Ab = 0, Eab = 0, V_Total = 6, operation_mode = modelicatest.CSTR.operation_type.Isothermal, order_b = {0, 0, 0, 0}, order_f = {1, 0.5, 0, 0}) annotation(
      Placement(visible = true, transformation(origin = {6, -20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(materialStream2.port2, cSTR1.port1) annotation(
      Line(points = {{-64, -44}, {-26, -44}, {-26, -20}, {-2, -20}, {-2, -20}}));
    connect(materialStream1.port2, cSTR1.port2) annotation(
      Line(points = {{-70, 4}, {6, 4}, {6, -12}, {6, -12}}));
  end CSTRtest;

  model MaterialStream
    extends unitoperations.compounds;
    parameter Real Flowrate = 100, Pressure = 1e5, Temperature = 300, molefraction[NOC] = zeros(NOC);
    parameter Boolean unspecified = true;
    parameter Boolean stepchange = false;
    parameter Real stepchangetime = 0.01;
    parameter Real step_value = 1;
    Real kf[NOC], zl[NOC](each min = 0, each max = 1, start = {0.5, 1e-18, 0.5, 0}), zv[NOC](each min = 0, each max = 1, start = {0, 0.25, 0, 0.75}), Fl(min = 0, start = 100), Fv(min = 0, start = 140), Tbf(start = 62), Tdf(start = 495.5), Psat_Tdf[NOC], Psat_Tbf[NOC], Psat_Tf[NOC], Pf, Tf, z[NOC], F;
    Real Hvf[NOC], Hlf[NOC], H;
    unitoperations.port port2 annotation(
      Placement(visible = true, transformation(origin = {80, -4}, extent = {{-18, -18}, {18, 18}}, rotation = 0), iconTransformation(origin = {85, 1}, extent = {{-21, -21}, {21, 21}}, rotation = 0)));
    unitoperations.port port1 annotation(
      Placement(visible = true, transformation(origin = {-80, -4}, extent = {{-18, -18}, {18, 18}}, rotation = 0), iconTransformation(origin = {-82, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
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
    annotation(
      Icon(graphics = {Rectangle(origin = {-10, 3}, fillColor = {255, 0, 0}, fillPattern = FillPattern.Solid, extent = {{-90, 45}, {90, -45}}), Text(origin = {-17, -1}, extent = {{-51, 27}, {51, -27}}, textString = "MS1")}, coordinateSystem(initialScale = 0.1)));
  end MaterialStream;

  model PTFlash
    extends unitoperationsModified.compounds;
    import Modelica.Constants.R;
    import Modelica.SIunits.{Temperature,Volume,Pressure,MolarFlowRate,MoleFraction,Time,Enthalpy,MolarEnthalpy,EnthalpyFlowRate,Area,Length};
    parameter Boolean Dynamic = true annotation(
      choices(checkBox = true));
    parameter Length hset = 0.7 annotation(
      Dialog(group = "Operating conditions"));
    parameter Pressure Pset = 5e5 annotation(
      Dialog(group = "Operating conditions"));
    parameter Boolean connectedToInput = false annotation(
      choices(checkBox = true));
    parameter Boolean OverrideSizeCalculations = false annotation(
      choices(checkBox = true),
      Dialog(tab = "Sizing"));
    parameter Real k_drum(unit = "ft/s") = 0.3 annotation(
      Dialog(tab = "Sizing"));
    parameter Area area = 4 annotation(
      Dialog(tab = "Sizing")), A(fixed = false, start = 1) annotation(
      Dialog(enable = false));
    parameter Volume volume = 8 annotation(
      Dialog(tab = "Sizing")), V_Total(fixed = false) annotation(
      Dialog(enable = false));
    parameter Temperature Ti = 310 annotation(
      Dialog(group = "Operating conditions"));
    MoleFraction z[NOC] annotation(
      each HideResult = true), y[NOC](each min = 0) "molefraction in vapor phase", x[NOC](each min = 0, start = {0.7, 1e-18, 0.3, 0}) "molefraction in liquid phase";
    MolarFlowRate F "Feed flowrate", L(start = 100, min = 0) "Liquid outlet flowrate", V(start = 140, min = 0) "Vapor outlet flowrate";
    Pressure P "Pressure inside column", Psat_T[NOC] annotation(
      each HideResult = true);
    Length h "Liquid level inside column";
    Volume VL "Liquid volume in column", VG "vapor volume in column";
    MolarEnthalpy hl[NOC] annotation(
      each HideResult = true), hv[NOC] annotation(
      each HideResult = true);
    Enthalpy H_M_Total;
    EnthalpyFlowRate Hl "Enthaply of liquid", Q "Energy required/removed", Hf "Enthalpy of feed", Hv "Enthalpy of vapor";
    Real densityi[NOC](each unit = "kmol/m3") annotation(
      each HideResult = true), k[NOC] annotation(
      each HideResult = true), M[NOC](each unit = "mol") annotation(
      each HideResult = true), M_Total(unit = "mol") "Total number of moles in column", ML(unit = "mol", start = 50) "Amount of liquid in column", MG(unit = "mol", start = 0.5) "Amount of vapor in column";
    unitoperations.sensor sensor1 annotation(
      Placement(visible = true, transformation(origin = {2, 82}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {8.88178e-16, 82}, extent = {{-14, -14}, {14, 14}}, rotation = 0)));
    unitoperations.sensor sensor3 annotation(
      Placement(visible = true, transformation(origin = {54, -2}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {57, 1}, extent = {{-15, -15}, {15, 15}}, rotation = 0)));
    unitoperationsModified.port port1 annotation(
      Placement(visible = true, transformation(origin = {47, -69}, extent = {{-17, -17}, {17, 17}}, rotation = 0), iconTransformation(origin = {51, -67}, extent = {{-15, -15}, {15, 15}}, rotation = 0)));
    unitoperationsModified.port port2 annotation(
      Placement(visible = true, transformation(origin = {50, 62}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {54, 66}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
    unitoperationsModified.port port3 annotation(
      Placement(visible = true, transformation(origin = {-54, 4}, extent = {{-18, -18}, {18, 18}}, rotation = 0), iconTransformation(origin = {-54, 2}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
    //parameter Real Cv(fixed = false);
  initial equation
    if Dynamic == true then
      h = hset;
      P = Pset;
      der(M_Total) = 0;
      
      for i in 1:NOC - 1 loop
        der(M[i]) = 0;
      end for;
    
    end if;
    
    if OverrideSizeCalculations == false then
      k_drum * 0.3048 * ((sum(x[:] .* densityi[:]) - P / (R * Ti * 1000)) * P / (R * Ti * 1000)) ^ 0.5 * 1000 * A = V;
      V_Total = A * 4 * (4 * A / 3.14) ^ 0.5;
    else
      A = area;
      V_Total = volume;
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
    F - L - V = if Dynamic == true then der(M_Total) else 0;
    for i in 1:NOC - 1 loop
      F * z[i] - L * x[i] - V * y[i] = if Dynamic == true then der(M[i]) else 0;
      M[i] = ML * x[i] + MG * y[i];
    end for;
    sum(M[:]) = M_Total;
    M_Total = MG + ML;
    VL = ML / (sum(x[:] .* densityi[:]) * 1000);
    VG = V_Total - VL;
//ideal gas law for gas phase
    P * VG = MG * R * Ti;
/*energy balance */
    Hv = V * sum(y[:] .* hv[:]);
    Hf = port3.enthalpy;
    Hl = L * sum(x[:] .* hl[:]);
    H_M_Total = ML * sum(x[:] .* hl[:]) + MG * sum(y[:] .* hv[:]);
    Hf - Hv - Hl + Q = 0;
/*Thermodynamic equations */
    sum(x[:]) = 1;
    sum(y[:]) = 1;
    y[:] = k[:] .* x[:];
    k[:] = Psat_T[:] / P;
    A * h = VL;
    if Dynamic == false then
      h = hset;
      P = Pset;
    end if;
//connector equations
    sensor1.var = P;
    sensor3.var = h;
// L = Cv * (P - 1e5)^0.5;
    port1.moleflow = L;
    port1.pressure = P;
    port1.temperature = Ti;
    port1.molefrac[:] = x[:];
    port2.moleflow = V;
    port2.pressure = P;
    port2.temperature = Ti;
    port2.molefrac[:] = y[:];
    annotation(
      Icon(coordinateSystem(extent = {{-70, -140}, {70, 100}}, preserveAspectRatio = false, initialScale = 0.1), graphics = {Polygon(origin = {-1, 2}, fillColor = {170, 170, 255}, fillPattern = FillPattern.Solid, points = {{-63, -62.0013}, {-51, -78.0013}, {-33, -86.0013}, {-13, -90.0013}, {15, -90.0013}, {39, -84.0013}, {55, -76.0013}, {63, -62.0013}, {63, 71.9987}, {45, 81.9987}, {29, 87.9987}, {17, 89.9987}, {-15, 89.9987}, {-33, 85.9987}, {-49, 79.9987}, {-63, 69.9987}, {-63, -62.0013}}), Line(origin = {-1.01, -25.28}, points = {{-62.9906, -4.72122}, {-52.9906, 3.27878}, {-42.9906, -4.72122}, {-32.9906, 3.27878}, {-24.9906, -4.72122}, {-8.99059, 5.27878}, {5.00941, -4.72122}, {21.0094, 5.27878}, {33.0094, -4.72122}, {45.0094, 5.27878}, {59.0094, -4.72122}, {63.0094, 1.27878}}, smooth = Smooth.Bezier), Text(origin = {1, -113}, extent = {{-59, 21}, {59, -21}}, textString = "PT Flash")}),
      Diagram(coordinateSystem(extent = {{-70, -140}, {70, 100}}, preserveAspectRatio = false)),
      version = "",
      uses);
  end PTFlash;



  model pttest
    unitoperationsModified.MaterialStream materialStream1(Flowrate = 88, molefraction = {0.25, 0.25, 0.25, 0.25}, pressure = 600000, specified_stream = true, stepchange = true) annotation(
      Placement(visible = true, transformation(origin = {-84, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    modelicatest.PTFlash pTFlash1(Dynamic = false, OverrideSizeCalculations = true, area = 0.3, connectedToInput = true, volume = 1) annotation(
      Placement(visible = true, transformation(origin = {-34, -4.17014}, extent = {{-23, -31.0382}, {23, 22.1701}}, rotation = 0)));
    unitoperationsModified.MaterialStream materialStream2(Tbf(displayUnit = "K", start = 20), Tdf(displayUnit = "K", start = 278)) annotation(
      Placement(visible = true, transformation(origin = {56, 18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    unitoperationsModified.MaterialStream materialStream3(Tdf(displayUnit = "K", start = 370)) annotation(
      Placement(visible = true, transformation(origin = {72, -34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    valve valve1(Dynamic = false, OutletPfixed = true) annotation(
      Placement(visible = true, transformation(origin = {22, 12}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    valve valve2(Dynamic = false, OutletPfixed = true) annotation(
      Placement(visible = true, transformation(origin = {30, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(valve2.port2, materialStream3.port1) annotation(
      Line(points = {{38, -40}, {53, -40}, {53, -34}, {64, -34}}));
    connect(valve1.port2, materialStream2.port1) annotation(
      Line(points = {{30, 12}, {43, 12}, {43, 18}, {48, 18}}));
    connect(pTFlash1.port1, valve2.port1) annotation(
      Line(points = {{-18, -22}, {-8, -22}, {-8, -40}, {22, -40}, {22, -40}}));
    connect(pTFlash1.port2, valve1.port1) annotation(
      Line(points = {{-16, 8}, {12, 8}, {12, 12}, {14, 12}}));
    connect(materialStream1.port2, pTFlash1.port3) annotation(
      Line(points = {{-76, -10}, {-52, -10}, {-52, -6}, {-52, -6}}));
  end pttest;

  model valve
    import Modelica.SIunits.{MolarFlowRate,Pressure};
    parameter Boolean Dynamic = true;
    parameter Real coeff(fixed = false, start = 0.5) "Coeff for valve", valveCv = 0.4 "valve Cv if not control valve";
    parameter Boolean Control = false;
    parameter Boolean OutletPfixed = false;
    parameter Pressure OutletPressure = 1e5 "used only when OutletPfixed is true";
    MolarFlowRate flowrate;
    Real Cv;
    Pressure outletP, delP;
    unitoperations.sensor sensor1 annotation(
      Placement(visible = true, transformation(origin = {1, 77}, extent = {{-19, -19}, {19, 19}}, rotation = 0), iconTransformation(origin = {2, 56}, extent = {{-18, -18}, {18, 18}}, rotation = 0)));
    unitoperations.port port1 annotation(
      Placement(visible = true, transformation(origin = {-80, -4}, extent = {{-18, -18}, {18, 18}}, rotation = 0), iconTransformation(origin = {-82, -2}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
    unitoperations.port port2 annotation(
      Placement(visible = true, transformation(origin = {79, 1}, extent = {{-17, -17}, {17, 17}}, rotation = 0), iconTransformation(origin = {82, 2}, extent = {{-18, -18}, {18, 18}}, rotation = 0)));
  equation
    if Control == false then
      sensor1.var = 0;
    end if;
    if Control == true and Dynamic == true then
      Cv = sensor1.var;
    elseif Control == false and Dynamic == true then
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
    annotation(
      Icon(graphics = {Polygon(origin = {-48.26, -2.99}, fillColor = {170, 170, 255}, fillPattern = FillPattern.Solid, points = {{-47.7381, 48.9887}, {-47.7381, -49.0113}, {48.2619, 2.98874}, {48.2619, 2.98874}, {-47.7381, 48.9887}}), Polygon(origin = {49.25, -4.98}, fillColor = {170, 170, 255}, fillPattern = FillPattern.Solid, points = {{-47.2509, 4.98071}, {46.7491, 48.9807}, {46.7491, -49.0193}, {-47.2509, 4.98071}}), Rectangle(origin = {1, 35}, fillColor = {170, 170, 255}, fillPattern = FillPattern.Solid, extent = {{-15, 35}, {15, -35}}), Text(origin = {0, -73}, extent = {{-52, 25}, {52, -25}}, textString = "Valve")}, coordinateSystem(initialScale = 0.1)));
  end valve;
end modelicatest;
