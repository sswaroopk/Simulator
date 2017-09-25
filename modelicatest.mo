package modelicatest
















  connector port
    parameter Integer NOC = 4;
    Real moleflow;
    parameter Boolean datagiven(fixed = false);
    annotation(Diagram(graphics = {Ellipse(origin = {-1, 1}, fillColor = {65, 252, 255}, fillPattern = FillPattern.Solid, extent = {{-67, 65}, {67, -65}}, endAngle = 360)}), Icon(graphics = {Ellipse(origin = {-4, -1}, fillColor = {81, 253, 248}, fillPattern = FillPattern.Solid, extent = {{-60, 61}, {60, -61}}, endAngle = 360)}));
  end port;

  model portparametertest
    parameter Boolean unspecified = false;
    parameter Boolean check(fixed = false);
    port port1 annotation(Placement(visible = true, transformation(origin = {-90, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-90, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    port port2 annotation(Placement(visible = true, transformation(origin = {90, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {90, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
      Real a;
      
    initial equation
    port1.datagiven = if unspecified == true then false else true;
    port1.datagiven = port2.datagiven;
    check = port1.datagiven;
  equation
    connect(port1, port2) annotation(Line(points = {{-90, 0}, {92, 0}, {92, 2}, {90, 2}}));
    port1.moleflow =20;
    if check == true then a =10; end if;
  end portparametertest;

  model phFlash
  extends unitoperations.compounds;
  parameter Real V_Total =8.11, F(fixed =false), z[NOC](each fixed =false), Hf(fixed =false),P = 1e5;
  parameter Real R = 8.314, A = 2, Q=0;
  Real h "liquid level";
  unitoperations.MaterialStream materialStream1(Flowrate = 1,Pressure = 10e5, molefraction = {0.25, 0.25, 0.25, 0.25}, unspecified = false)  annotation(Placement(visible = true, transformation(origin = {-52, 8}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Real Psat_T[NOC], k[NOC], hv[NOC], hl[NOC], densityi[NOC],L(start = 0.50),V(start = 0.50), x[NOC], y[NOC],M[NOC], ML, MG, M_Total; 
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
    unitoperations.MaterialStream materialStream1(unspecified = false)  annotation(Placement(visible = true, transformation(origin = {-12, 4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  end ms;

   model phFlashDynamic
  extends unitoperations.compounds;
  parameter Real V_Total =8.11;
  parameter Real hset = 2, Pset = 5e5, Pout = 1e5;
  parameter Real R = 8.314, A = 2, Q=0;
  Real Cv1, Cv2;
  Real h "liquid level";
  unitoperations.MaterialStream materialStream1(Pressure = 10e5, molefraction = {0.3, 0.3, 0.2, 0.2}, unspecified = false)  annotation(Placement(visible = true, transformation(origin = {-52, 8}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Real F, z[NOC], Hf;
  Real Psat_T[NOC], k[NOC], hv[NOC], hl[NOC], densityi[NOC],L(start = 50),V(start = 50), x[NOC](each start = 0.25), y[NOC],M[NOC], ML, MG, M_Total; 
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
   L = Cv1*0.4*(P-Pout)^0.5;
   V = Cv2*0.4*(P-Pout)^0.5;
   der(Cv2) = 0;
   der(Cv1) = 0;
  end phFlashDynamic;

  model parameterTest
    parameter Real a(fixed = false);
    Real b, c;
    initial equation
    a = 2*b;
    equation
      b = 4;
      c = a * b;
  end parameterTest;

  model matrixtest
  parameter Real a[4,2] = ones(4,2), b[4,2] = [1,2;3,4;5,6;7,8];
  Real c;
  equation
  c = a[1,:]*b[1,:];
  end matrixtest;

type compounds = enumeration(Air, Argon, Bromine, Carbontetrachloride, Carbonmonoxide, Carbondioxide, Carbondisulfide, Phosgene, Trichloroacetylchloride, Hydrogenchloride, Chlorine, Hydrogeniodide, Hydrogen, Water, Hydrogensulfide, Ammonia, Neon, Nitricacid, Nitricoxide, Nitrogendioxide, Nitrogen, Nitrousoxide, Oxygen, Sulfurdioxide, Sulfurtrioxide, Chloroform, Hydrogencyanide, Formaldehyde, Methylchloride, Methyliodide, Methane, Methanol, Methylamine, Trichloroethylene, Dichloroacetylchloride, Trichloroacetaldehyde, Acetylene, Dichloroacetaldehyde, Vinylchloride, Acetylchloride, OneOneTwotrichloroethane, Acetonitrile, Ethylene, OneOnedichloroethane, OneTwodichloroethane, Acetaldehyde, Ethyleneoxide, Aceticacid, Methylformate, Ethylchloride, Ethane, Ethanol, Dimethylether, Ethyleneglycol, Dimethylsulfide, Ethylmercaptan, Ethylamine, Acrylonitrile, Methylacetylene, Propadiene, Propylene, Acetone, Ethylformate, Methylacetate, Propionicacid, Nndimethylformamide, Propane, Isopropanol, Onepropanol, Trimethylamine, Vinylacetylene, Thiophene, Methacrylonitrile, Dimethylacetylene, Ethylacetylene, OneTwobutadiene, OneThreebutadiene, Onebutene, CisTwobutene, TransTwobutene, Isobutene, Twomethylpropanal, Methylethylketone, Tetrahydrofuran, OneFourdioxane, Nbutyricacid, Ethylacetate, Methylpropionate, Npropylformate, Sulfolane, Nndimethylacetamide, Nbutane, Isobutane, Onebutanol, TwomethylOnepropanol, Twobutanol, TwomethylTwopropanol, Diethylether, Diethyleneglycol, Diethylamine, Furfural, Pyridine, Isoprene, Cyclopentane, TwomethylOnebutene, ThreemethylOnebutene, TwomethylTwobutene, Onepentene, CisTwopentene, TransTwopentene, Threepentanone, Methylisopropylketone, Npropylacetate, Isopentane, Npentane, Neopentane, OneTwoFourtrichlorobenzene, Mdichlorobenzene, Odichlorobenzene, Pdichlorobenzene, Bromobenzene, Monochlorobenzene, Iodobenzene, Nitrobenzene, Benzene, Phenol, Aniline, Cyclohexanone, Cyclohexane, Onehexene, Methylcyclopentane, Cyclohexanol, TwoTwodimethylbutane, TwoThreedimethylbutane, Nhexane, Twomethylpentane, Threemethylpentane, Triethyleneglycol, Triethylamine, Toluene, Mcresol, Ocresol, Pcresol, Methylcyclohexane, Ethylcyclopentane, Oneheptene, Nheptane, Styrene, Ethylbenzene, Mxylene, Oxylene, Pxylene, Ethylcyclohexane, Npropylcyclopentane, Noctane, TwoTwoThreetrimethylpentane, TwoTwoFourtrimethylpentane, TwoThreeThreetrimethylpentane, TwoThreeFourtrimethylpentane, Tetraethyleneglycol, Indene, Indane, Cumene, Npropylbenzene, Npropylcyclohexane, Nnonane, Naphthalene, Onemethylindene, Twomethylindene, Dicyclopentadiene, Nbutylbenzene, Nbutylcyclohexane, Ndecane, Onemethylnaphthalene, Twomethylnaphthalene, Nundecane, Acenaphthene, Biphenyl, Ndodecane, Fluorene, Ntridecane, Phenanthrene, Ntetradecane, Npentadecane, Fluoranthene, Pyrene, Onephenylnaphthalene, Nhexadecane, Chrysene, Cisdecahydronaphthalene, Transdecahydronaphthalene, Methyltertbutylether, Methyltertpentylether, TwomethylTwobutanol, Nitrogentrioxide, Nitrogentetroxide, HeliumFour, Fluorine, Krypton, Xenon, Ozone, Carbonylsulfide, Sulfurhexafluoride, Dimethylsulfoxide, Nheptadecane, Noctadecane, Nnonadecane, Nheneicosane, Ndocosane, Ntricosane, Ntetracosane, Npentacosane, Nhexacosane, Nheptacosane, Noctacosane, Nnonacosane, Squalane);

  model compoundsloadingfile
  parameter Integer NOC = 4;
  parameter compounds comp1;
  parameter Chemsep_Database.Toluene comp2;
  annotation(Diagram(graphics = {Rectangle(origin = {0, -3}, fillColor = {255, 85, 255}, fillPattern = FillPattern.Solid, extent = {{-88, 79}, {88, -79}})}), Icon(graphics = {Rectangle(origin = {-2, -7}, fillColor = {255, 85, 255}, fillPattern = FillPattern.Solid, extent = {{-94, 81}, {94, -81}})}));
  end compoundsloadingfile;

  model compoundstest
    compoundsloadingfile compoundsloadingfile1(comp1 = modelicatest.compounds.Air)  annotation(Placement(visible = true, transformation(origin = {0, -8}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  end compoundstest;

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
    Real H_r, Hout, Q, T(start = 650);
    unitoperations.port port1 annotation(Placement(visible = true, transformation(origin = {-82, 2}, extent = {{-18, -18}, {18, 18}}, rotation = 0), iconTransformation(origin = {-85, 1}, extent = {{-17, -17}, {17, 17}}, rotation = 0)));
    unitoperations.port port2 annotation(Placement(visible = true, transformation(origin = {2, 84}, extent = {{-18, -18}, {18, 18}}, rotation = 0), iconTransformation(origin = {3, 83}, extent = {{-17, -17}, {17, 17}}, rotation = 0)));
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
    r = kf * product(c[:].^order_f[:]) - kb * product(c[:].^order_b[:]);
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
    annotation(Icon(graphics = {Rectangle(origin = {-1, 1}, fillColor = {255, 170, 0}, fillPattern = FillPattern.Solid, extent = {{-97, 95}, {97, -95}}), Text(origin = {2, 0}, extent = {{-80, 62}, {80, -62}}, textString = "CSTR")}, coordinateSystem(initialScale = 0.1)));
  end CSTR;

  model PFR
  end PFR;

  model CSTRtest
  
    MaterialStream materialStream1(Flowrate = 66, Pressure = 22e5, Temperature = 800, molefraction = {0, 1, 0, 0}, unspecified = false)  annotation(Placement(visible = true, transformation(origin = {-78, 4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    MaterialStream materialStream2(Flowrate = 22.2, Pressure = 22e5, Temperature = 800, molefraction = {1, 0, 0, 0}, unspecified = false)  annotation(Placement(visible = true, transformation(origin = {-72, -44}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  CSTR cSTR1(Ab = 0, Eab = 0, V_Total = 6, operation_mode = modelicatest.CSTR.operation_type.Isothermal, order_b = {0, 0, 0, 0}, order_f = {1, 0.5, 0, 0})  annotation(Placement(visible = true, transformation(origin = {6, -20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(materialStream2.port2, cSTR1.port1) annotation(Line(points = {{-64, -44}, {-26, -44}, {-26, -20}, {-2, -20}, {-2, -20}}));
    connect(materialStream1.port2, cSTR1.port2) annotation(Line(points = {{-70, 4}, {6, 4}, {6, -12}, {6, -12}}));
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
end modelicatest;