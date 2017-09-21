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
  unitoperations.MaterialStream materialStream1(Pressure = 10e5, molefraction = {0.3, 0.3, 0.2, 0.2}, unspecified = false)  annotation(Placement(visible = true, transformation(origin = {-52, 8}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Real Psat_T[NOC], k[NOC], hv[NOC], hl[NOC], densityi[NOC],L(start = 50),V(start = 50), x[NOC], y[NOC],M[NOC], ML, MG, M_Total; 
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
  c = max(a*transpose(b));
  end matrixtest;
end modelicatest;
