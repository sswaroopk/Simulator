package BenzeneProduction
  model Reactor
    unitoperations.MaterialStream materialStream1(Flowrate = 22, Pressure = 24, Temperature = 600, molefraction = {0.9, 0, 0.1, 0}, step_value = 2, stepchange = true, unspecified = false) annotation(Placement(visible = true, transformation(origin = {-82, -22}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
  unitoperations.MaterialStream materialStream2(Flowrate = 66, Pressure = 24, Temperature = 600, molefraction = {0, 0.9, 0, 0.1}, unspecified = false) annotation(Placement(visible = true, transformation(origin = {-84, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  unitoperations.CSTR cSTR1(Ab = 0, Af = 5.1e11, Eab = 0, Eaf = 230e3, T_iso = 700, V_Total = 1, operation_mode = unitoperations.CSTR.operation_type.Isothermal, order_b = {0, 0, 0, 0}, order_f = {1, 0.5, 0, 0}) annotation(Placement(visible = true, transformation(origin = {12, 4}, extent = {{-18, -18}, {18, 18}}, rotation = 0)));
  unitoperations.MaterialStream materialStream3 annotation(Placement(visible = true, transformation(origin = {102, -4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  unitoperations.valve valve1(OutletPfixed = true) annotation(Placement(visible = true, transformation(origin = {56, -4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(valve1.port2, materialStream3.port1) annotation(Line(points = {{64, -4}, {94, -4}}));
    connect(cSTR1.port3, valve1.port1) annotation(Line(points = {{28, -3}, {38.5, -3}, {38.5, -4}, {48, -4}}));
    connect(materialStream2.port2, cSTR1.port2) annotation(Line(points = {{-75.5, 10}, {-4, 10}, {-4, 11}}));
    connect(materialStream1.port2, cSTR1.port1) annotation(Line(points = {{-72, -22}, {-40.5, -22}, {-40.5, 5}, {-4, 5}}));
  end Reactor;

  model Flashafterreactor
  unitoperations.MaterialStream materialStream1(Flowrate = 88, Pressure = 24, Temperature = 300, molefraction = {0.175, 0.625, 0.75, 0.125}, unspecified = false) annotation(Placement(visible = true, transformation(origin = {-80, 8}, extent = {{-14, -14}, {14, 14}}, rotation = 0)));
  unitoperations.PHFlash PHFlash1(L(start = 22), OverrideSizeCalculations = false, Ti(start = 280), V(start = 66), connectedToInput = true) annotation(Placement(visible = true, transformation(origin = {-22, 10.1631}, extent = {{-30, -44.4964}, {30, 31.7832}}, rotation = 0)));
  unitoperations.MaterialStream materialStream2 annotation(Placement(visible = true, transformation(origin = {86, 28}, extent = {{-14, -14}, {14, 14}}, rotation = 0)));
  unitoperations.valve valve1(OutletPfixed = true) annotation(Placement(visible = true, transformation(origin = {48, 28}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
  unitoperations.valve valve2(OutletPfixed = true) annotation(Placement(visible = true, transformation(origin = {50, -20}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
  unitoperations.MaterialStream materialStream3(Tdf(start = 353)) annotation(Placement(visible = true, transformation(origin = {88, -20}, extent = {{-14, -14}, {14, 14}}, rotation = 0)));
  equation
    connect(materialStream1.port2, PHFlash1.port3) annotation(Line(points = {{-68, 8}, {-46, 8}}));
    connect(PHFlash1.port2, valve1.port1) annotation(Line(points = {{0, 27}, {-17.5, 27}, {-17.5, 28}, {38, 28}}));
    connect(PHFlash1.port1, valve2.port1) annotation(Line(points = {{-1, -14}, {23.5, -14}, {23.5, -20}, {40, -20}}));
    connect(valve2.port2, materialStream3.port1) annotation(Line(points = {{60, -20}, {77, -20}}));
    connect(valve1.port2, materialStream2.port1) annotation(Line(points = {{58, 28}, {75, 28}}));
  end Flashafterreactor;

  model PTflashAfterReactor
  unitoperations.MaterialStream materialStream1(Flowrate = 88, Pressure = 24, Temperature = 300, molefraction = {0.175, 0.625, 0.75, 0.125}, unspecified = false) annotation(Placement(visible = true, transformation(origin = {-75, 11}, extent = {{-15, -15}, {15, 15}}, rotation = 0)));
  unitoperations.MaterialStream materialStream2 annotation(Placement(visible = true, transformation(origin = {85, 27}, extent = {{-13, -13}, {13, 13}}, rotation = 0)));
  unitoperations.valve valve1(OutletPfixed = true) annotation(Placement(visible = true, transformation(origin = {43, 27}, extent = {{-17, -17}, {17, 17}}, rotation = 0)));
  unitoperations.valve valve2(OutletPfixed = true) annotation(Placement(visible = true, transformation(origin = {41, -37}, extent = {{-13, -13}, {13, 13}}, rotation = 0)));
  unitoperations.MaterialStream materialStream3(Tbf(start = 300), Tdf(start = 350)) annotation(Placement(visible = true, transformation(origin = {86, -38}, extent = {{-14, -14}, {14, 14}}, rotation = 0)));
  unitoperations.PTFlash PTFlash1(Area = 1, L(start = 22), MG(start = 400), ML(start = 700), OverrideSizeCalculations = true, Pset = 5, Ti = 300, V(start = 60), Volume = 2, connectedToInput = true, hset = 0.2, x(start = {0.7, 1e-16, 0.28, 0.001})) annotation(Placement(visible = true, transformation(origin = {-23, 13.3364}, extent = {{-25, -34.9956}, {25, 24.9969}}, rotation = 0)));
  equation
    connect(materialStream1.port2, PTFlash1.port3) annotation(Line(points = {{-62, 11}, {-42, 11}}));
    connect(PTFlash1.port1, valve2.port1) annotation(Line(points = {{-5, -6}, {-5, -37}, {30, -37}}));
    connect(PTFlash1.port2, valve1.port1) annotation(Line(points = {{-4, 27}, {29, 27}}));
    connect(valve2.port2, materialStream3.port1) annotation(Line(points = {{52, -37}, {52, -38}, {75, -38}}));
    connect(valve1.port2, materialStream2.port1) annotation(Line(points = {{57, 27}, {74, 27}}));
  end PTflashAfterReactor;

  model reactorAndFlash
    unitoperations.MaterialStream materialStream1(Flowrate = 22, Pressure = 24, Temperature = 600, molefraction = {0.9, 0, 0.1, 0}, unspecified = false) annotation(Placement(visible = true, transformation(origin = {-82, -22}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
  unitoperations.MaterialStream materialStream2(Flowrate = 66, Pressure = 24, Temperature = 600, molefraction = {0, 0.9, 0, 0.1}, unspecified = false) annotation(Placement(visible = true, transformation(origin = {-86, 2}, extent = {{10, 10}, {-10, -10}}, rotation = 0)));
  unitoperations.CSTR cSTR1(Ab = 0, Af = 5.1e11, Eab = 0, Eaf = 230e3, T_iso = 700, V_Total = 1, operation_mode = unitoperations.CSTR.operation_type.Isothermal, order_b = {0, 0, 0, 0}, order_f = {1, 0.5, 0, 0}) annotation(Placement(visible = true, transformation(origin = {-22, -4}, extent = {{-18, -18}, {18, 18}}, rotation = 0)));
  unitoperations.MaterialStream materialStream3 annotation(Placement(visible = true, transformation(origin = {56, -14}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  unitoperations.valve valve1(OutletPfixed = false) annotation(Placement(visible = true, transformation(origin = {20, -12}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  unitoperations.PTFlash PTFlash1 annotation(Placement(visible = true, transformation(origin = {97, -12.0255}, extent = {{-19, -28.3078}, {19, 20.2199}}, rotation = 0)));
  unitoperations.valve valve2(OutletPfixed = true) annotation(Placement(visible = true, transformation(origin = {150, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  unitoperations.valve valve3(OutletPfixed = true) annotation(Placement(visible = true, transformation(origin = {144, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  unitoperations.MaterialStream materialStream4 annotation(Placement(visible = true, transformation(origin = {188, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  unitoperations.MaterialStream materialStream5 annotation(Placement(visible = true, transformation(origin = {184, -50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(valve3.port2, materialStream5.port1) annotation(Line(points = {{152, -50}, {176, -50}}));
    connect(PTFlash1.port1, valve3.port1) annotation(Line(points = {{111, -28}, {94, -28}, {94, -50}, {136, -50}}));
    connect(valve2.port2, materialStream4.port1) annotation(Line(points = {{158, -2}, {180, -2}}));
    connect(PTFlash1.port2, valve2.port1) annotation(Line(points = {{112, -1}, {142, -1}, {142, -2}}));
    connect(materialStream3.port2, PTFlash1.port3) annotation(Line(points = {{64.5, -14}, {82, -14}}));
    connect(valve1.port2, materialStream3.port1) annotation(Line(points = {{28, -12}, {38, -12}, {38, -14}, {48, -14}}));
    connect(cSTR1.port3, valve1.port1) annotation(Line(points = {{-6, -11}, {-4.5, -11}, {-4.5, -12}, {12, -12}}));
    connect(materialStream2.port2, cSTR1.port2) annotation(Line(points = {{-94.5, 2}, {-38, 2}, {-38, 3}}));
    connect(materialStream1.port2, cSTR1.port1) annotation(Line(points = {{-72, -22}, {-58.5, -22}, {-58.5, -3}, {-38, -3}}));
    annotation(Diagram(coordinateSystem(extent = {{-100, -100}, {200, 100}})), Icon(coordinateSystem(extent = {{-100, -100}, {200, 100}})), version = "", uses);end reactorAndFlash;

  model WithoutRecycle
  unitoperations.MaterialStream materialStream1(Flowrate = 22, Pressure = 24, Temperature = 600, molefraction = {0.9, 0, 0.1, 0}, step_value = 2, stepchange = true, unspecified = false) annotation(Placement(visible = true, transformation(origin = {-136, -22}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
  unitoperations.MaterialStream materialStream2(Flowrate = 66, Pressure = 24, Temperature = 600, molefraction = {0, 0.9, 0, 0.1}, unspecified = false) annotation(Placement(visible = true, transformation(origin = {-138, 16}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  unitoperations.CSTR cSTR1(Ab = 0, Af = 5.1e11, Eab = 0, Eaf = 230e3, T_iso = 700, V_Total = 1, operation_mode = unitoperations.CSTR.operation_type.Isothermal, order_b = {0, 0, 0, 0}, order_f = {1, 0.5, 0, 0}) annotation(Placement(visible = true, transformation(origin = {-70, 10}, extent = {{-18, -18}, {18, 18}}, rotation = 0)));
  unitoperations.MaterialStream materialStream3 annotation(Placement(visible = true, transformation(origin = {3, 3}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
  unitoperations.valve valve1(OutletPfixed = false) annotation(Placement(visible = true, transformation(origin = {-29, 3}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
  unitoperations.PTFlash PTFlash1 annotation(Placement(visible = true, transformation(origin = {46, 4.36213}, extent = {{-20, -28.6955}, {20, 20.4968}}, rotation = 0)));
  unitoperations.valve valve2(OutletPfixed = true) annotation(Placement(visible = true, transformation(origin = {88, 16}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  unitoperations.valve valve3(OutletPfixed = true) annotation(Placement(visible = true, transformation(origin = {80, -46}, extent = {{14, -14}, {-14, 14}}, rotation = 0)));
  unitoperations.MaterialStream materialStream4 annotation(Placement(visible = true, transformation(origin = {122, 16}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  unitoperations.MaterialStream materialStream5 annotation(Placement(visible = true, transformation(origin = {122, -46}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  unitoperations.Distillation distillation1(Override_Sizing_Calculations = false, h_weir = 0.1, specification1 = unitoperations.Distillation.spec1.RefluxRatio, specification1_value = 2, specification2 = unitoperations.Distillation.spec2.ProductMolarFlow, specification2_value = 7) annotation(Placement(visible = true, transformation(origin = {163, -44.441}, extent = {{-21, -25.059}, {21, 17.8993}}, rotation = 0)));
  unitoperations.MaterialStream materialStream6 annotation(Placement(visible = true, transformation(origin = {223, -17}, extent = {{-13, -13}, {13, 13}}, rotation = 0)));
  unitoperations.MaterialStream materialStream7(Tbf(start = 380), Tdf(start = 380), flashCalculations = false) annotation(Placement(visible = true, transformation(origin = {223, -83}, extent = {{-13, -13}, {13, 13}}, rotation = 0)));
  equation
    connect(distillation1.port2, materialStream6.port1) annotation(Line(points = {{178, -35}, {186, -35}, {186, -17}, {212, -17}}));
    connect(distillation1.port3, materialStream7.port1) annotation(Line(points = {{178, -57}, {203, -57}, {203, -83}, {212, -83}}));
    connect(materialStream5.port2, distillation1.port1) annotation(Line(points = {{130.5, -46}, {146, -46}}));
    connect(valve3.port2, materialStream5.port1) annotation(Line(points = {{69, -46}, {114, -46}}));
    connect(PTFlash1.port1, valve3.port1) annotation(Line(points = {{61, -11}, {60, -11}, {60, -46}, {91, -46}}));
    connect(valve2.port2, materialStream4.port1) annotation(Line(points = {{96, 16}, {114, 16}}));
    connect(PTFlash1.port2, valve2.port1) annotation(Line(points = {{61, 16}, {80, 16}}));
    connect(materialStream3.port2, PTFlash1.port3) annotation(Line(points = {{12, 3}, {31, 3}}));
    connect(valve1.port2, materialStream3.port1) annotation(Line(points = {{-20, 3}, {-6, 3}}));
    connect(cSTR1.port3, valve1.port1) annotation(Line(points = {{-54, 3}, {-38, 3}}));
    connect(materialStream1.port2, cSTR1.port1) annotation(Line(points = {{-126, -22}, {-99.5, -22}, {-99.5, 11}, {-86, 11}}));
    connect(materialStream2.port2, cSTR1.port2) annotation(Line(points = {{-129.5, 16}, {-86, 16}, {-86, 17}}));
    annotation(Diagram(coordinateSystem(extent = {{-200, -100}, {200, 100}})), Icon(coordinateSystem(extent = {{-200, -100}, {200, 100}})), version = "", uses);end WithoutRecycle;

  model WithoutrecyclePhflash
  unitoperations.MaterialStream materialStream1(Flowrate = 22, Pressure = 24, Temperature = 600, molefraction = {0.9, 0, 0.1, 0}, step_value = 2, stepchange = true, unspecified = false) annotation(Placement(visible = true, transformation(origin = {-284, -4}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
  unitoperations.MaterialStream materialStream2(Flowrate = 66, Pressure = 24, Temperature = 600, molefraction = {0, 0.9, 0, 0.1}, unspecified = false) annotation(Placement(visible = true, transformation(origin = {-284, 46}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  unitoperations.CSTR cSTR1(Ab = 0, Af = 5.1e11, Eab = 0, Eaf = 230e3, T_iso = 700, V_Total = 1, operation_mode = unitoperations.CSTR.operation_type.Isothermal, order_b = {0, 0, 0, 0}, order_f = {1, 0.5, 0, 0}) annotation(Placement(visible = true, transformation(origin = {-202.92, 28.1767}, extent = {{-25.0804, -36.4167}, {25.0804, 30.3473}}, rotation = 0)));
  unitoperations.MaterialStream materialStream3 annotation(Placement(visible = true, transformation(origin = {-66, -2}, extent = {{-14, -14}, {14, 14}}, rotation = 0)));
  unitoperations.valve valve1(OutletPfixed = false) annotation(Placement(visible = true, transformation(origin = {-112, -2}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
  unitoperations.valve valve2(OutletPfixed = true) annotation(Placement(visible = true, transformation(origin = {116, 16}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  unitoperations.valve valve3(OutletPfixed = true) annotation(Placement(visible = true, transformation(origin = {120, -46}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
  unitoperations.MaterialStream materialStream4 annotation(Placement(visible = true, transformation(origin = {152, 16}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  unitoperations.MaterialStream materialStream5 annotation(Placement(visible = true, transformation(origin = {160, -46}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
  unitoperations.Distillation distillation1(Override_Sizing_Calculations = false, h_weir = 0.1, specification1 = unitoperations.Distillation.spec1.RefluxRatio, specification1_value = 2, specification2 = unitoperations.Distillation.spec2.ProductMolarFlow, specification2_value = 7) annotation(Placement(visible = true, transformation(origin = {211.129, -42.7404}, extent = {{-21.9441, -45.7842}, {21.9441, 32.703}}, rotation = 0)));
  unitoperations.MaterialStream materialStream6 annotation(Placement(visible = true, transformation(origin = {272, -22}, extent = {{-14, -14}, {14, 14}}, rotation = 0)));
  unitoperations.MaterialStream materialStream7(Tbf(start = 380), Tdf(start = 380), flashCalculations = false) annotation(Placement(visible = true, transformation(origin = {270, -72}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
  unitoperations.PHFlash PHFlash1(hset = 2) annotation(Placement(visible = true, transformation(origin = {72.7743, -0.576747}, extent = {{-21.2869, -42.4196}, {21.2869, 30.2997}}, rotation = 0)));
  unitoperations.HeatExchanger heatExchanger1(pressure_drop = 0.5) annotation(Placement(visible = true, transformation(origin = {-21, -1}, extent = {{-15, -21}, {15, 15}}, rotation = 0)));
  unitoperations.MaterialStream materialStream8 annotation(Placement(visible = true, transformation(origin = {24, -2}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
  equation
    connect(distillation1.port2, materialStream6.port1) annotation(Line(points = {{222, -22}, {261, -22}}));
    connect(distillation1.port3, materialStream7.port1) annotation(Line(points = {{222, -72}, {260, -72}}));
    connect(materialStream5.port2, distillation1.port1) annotation(Line(points = {{170, -46}, {199, -46}}));
    connect(valve3.port2, materialStream5.port1) annotation(Line(points = {{130, -46}, {150, -46}}));
    connect(PHFlash1.port1, valve3.port1) annotation(Line(points = {{83, -27}, {83, -47}, {100, -47}, {100, -46.5}, {110, -46.5}, {110, -46}}));
    connect(PHFlash1.port2, valve2.port1) annotation(Line(points = {{84, 20}, {104, 20}, {104, 16}, {108, 16}}));
    connect(materialStream8.port2, PHFlash1.port3) annotation(Line(points = {{34, -2}, {61, -2}}));
    connect(materialStream2.port2, cSTR1.port2) annotation(Line(points = {{-275.5, 46}, {-247.25, 46}, {-247.25, 38}, {-225, 38}}));
    connect(materialStream1.port2, cSTR1.port1) annotation(Line(points = {{-274, -4}, {-235.5, -4}, {-235.5, 26}, {-225, 26}}));
    connect(cSTR1.port3, valve1.port1) annotation(Line(points = {{-181, 10}, {-181, 9}, {-125, 9}, {-125, -2}}));
    connect(valve2.port2, materialStream4.port1) annotation(Line(points = {{124, 16}, {144, 16}}));
    connect(heatExchanger1.port2, materialStream8.port1) annotation(Line(points = {{-8, -2}, {14, -2}}));
    connect(materialStream3.port2, heatExchanger1.port1) annotation(Line(points = {{-54, -2}, {-34, -2}}));
    connect(valve1.port2, materialStream3.port1) annotation(Line(points = {{-99, -2}, {-77, -2}}));
    annotation(Diagram(coordinateSystem(extent = {{-300, -100}, {300, 100}}, preserveAspectRatio = false)), Icon(coordinateSystem(extent = {{-300, -100}, {300, 100}}, preserveAspectRatio = false)), version = "", uses);end WithoutrecyclePhflash;

  model WithRecycle
  unitoperations.MaterialStream materialStream1(Flowrate = 22, Pressure = 24, Temperature = 450, molefraction = {0.9, 0, 0.1, 0}, step_value = 2, stepchange = true, unspecified = false) annotation(Placement(visible = true, transformation(origin = {-280, 18}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
  unitoperations.MaterialStream materialStream2(Flowrate = 66, Pressure = 24, Temperature = 300, molefraction = {0, 0.9, 0, 0.1}, unspecified = false) annotation(Placement(visible = true, transformation(origin = {-280, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  unitoperations.CSTR cSTR1(Ab = 0, Af = 5.1e11, Eab = 0, Eaf = 230e3, T_iso = 700, V_Total = 1, operation_mode = unitoperations.CSTR.operation_type.Isothermal, order_b = {0, 0, 0, 0}, order_f = {1, 0.5, 0, 0}) annotation(Placement(visible = true, transformation(origin = {-142.216, 18.9457}, extent = {{-19.7835, -21.582}, {19.7835, 17.985}}, rotation = 0)));
  unitoperations.MaterialStream materialStream3 annotation(Placement(visible = true, transformation(origin = {-58, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  unitoperations.valve valve1(OutletPfixed = false) annotation(Placement(visible = true, transformation(origin = {-94, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  unitoperations.PTFlash PTFlash1 annotation(Placement(visible = true, transformation(origin = {1, 13.5154}, extent = {{-23, -42.0154}, {23, 30.011}}, rotation = 0)));
  unitoperations.valve valve2(OutletPfixed = true) annotation(Placement(visible = true, transformation(origin = {52, 28}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  unitoperations.valve valve3(OutletPfixed = true) annotation(Placement(visible = true, transformation(origin = {60, -26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  unitoperations.MaterialStream materialStream4 annotation(Placement(visible = true, transformation(origin = {92, 56}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  unitoperations.MaterialStream materialStream5 annotation(Placement(visible = true, transformation(origin = {92, -26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  unitoperations.Distillation distillation1(Override_Sizing_Calculations = false, h_weir = 0.1, specification1 = unitoperations.Distillation.spec1.RefluxRatio, specification1_value = 2, specification2 = unitoperations.Distillation.spec2.ProductMolarFlow, specification2_value = 15) annotation(Placement(visible = true, transformation(origin = {147, -22.8324}, extent = {{-31, -38.6205}, {31, 27.5861}}, rotation = 0)));
  unitoperations.MaterialStream materialStream6(flashCalculations = false) annotation(Placement(visible = true, transformation(origin = {213, 11}, extent = {{-17, -17}, {17, 17}}, rotation = 0)));
  unitoperations.MaterialStream materialStream7(flashCalculations = false) annotation(Placement(visible = true, transformation(origin = {215, -43}, extent = {{-17, -17}, {17, 17}}, rotation = 0)));
  unitoperations.Mixer mixer1(T(start = 530))  annotation(Placement(visible = true, transformation(origin = {-234, 4}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
  unitoperations.Pump pump1 annotation(Placement(visible = true, transformation(origin = {-57, -65}, extent = {{19, -26.6}, {-19, 19}}, rotation = 0)));
  unitoperations.MaterialStream materialStream8 annotation(Placement(visible = true, transformation(origin = {-163, -57}, extent = {{15, -15}, {-15, 15}}, rotation = 0)));
  unitoperations.MaterialStream materialStream9 annotation(Placement(visible = true, transformation(origin = {-194, 6}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(materialStream9.port2, cSTR1.port1) annotation(Line(points = {{-186, 6}, {-176, 6}, {-176, 18}, {-160, 18}, {-160, 18}}));
    connect(mixer1.port3, materialStream9.port1) annotation(Line(points = {{-220, 6}, {-202, 6}, {-202, 6}, {-202, 6}}));
    connect(distillation1.port3, materialStream7.port1) annotation(Line(points = {{169, -43}, {201, -43}}));
    connect(materialStream7.port2, pump1.port1) annotation(Line(points = {{229, -43}, {244, -43}, {244, -66}, {-40, -66}}));
    connect(distillation1.port2, materialStream6.port1) annotation(Line(points = {{169, -9}, {190, -9}, {190, 11}, {199, 11}}));
    connect(materialStream8.port2, mixer1.port2) annotation(Line(points = {{-176, -56}, {-278, -56}, {-278, -10}, {-248, -10}, {-248, -8}, {-248, -8}}));
    connect(materialStream8.port1, pump1.port2) annotation(Line(points = {{-150, -56}, {-74, -56}, {-74, -58}, {-74, -58}}));
    connect(materialStream5.port2, distillation1.port1) annotation(Line(points = {{100.5, -26}, {122, -26}}));
    connect(valve3.port2, materialStream5.port1) annotation(Line(points = {{68, -26}, {84, -26}}));
    connect(PTFlash1.port1, valve3.port1) annotation(Line(points = {{18, -11}, {18, -26}, {52, -26}}));
    connect(valve2.port2, materialStream4.port1) annotation(Line(points = {{60, 28}, {72, 28}, {72, 56}, {84, 56}}));
    connect(PTFlash1.port2, valve2.port1) annotation(Line(points = {{19, 29}, {44, 29}, {44, 28}}));
    connect(materialStream3.port2, PTFlash1.port3) annotation(Line(points = {{-49.5, 10}, {-17, 10}}));
    connect(valve1.port2, materialStream3.port1) annotation(Line(points = {{-86, 10}, {-66, 10}}));
    connect(cSTR1.port3, valve1.port1) annotation(Line(points = {{-125, 10}, {-102, 10}}));
    connect(materialStream2.port2, cSTR1.port2) annotation(Line(points = {{-271.5, 60}, {-203, 60}, {-203, 23.25}, {-159, 23.25}, {-159, 25}}));
    connect(materialStream1.port2, mixer1.port1) annotation(Line(points = {{-270, 18}, {-248, 18}}));
    annotation(Diagram(coordinateSystem(extent = {{-300, -100}, {300, 100}}, preserveAspectRatio = false)), Icon(coordinateSystem(extent = {{-300, -100}, {300, 100}}, preserveAspectRatio = false)), version = "", uses);end WithRecycle;

  model Withrecycle1
  unitoperations.MaterialStream materialStream1(Flowrate = 22, Pressure = 24, Temperature = 600, molefraction = {0.9, 0, 0.1, 0}, step_value = 2, stepchange = true, unspecified = false) annotation(Placement(visible = true, transformation(origin = {-280, 18}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
  unitoperations.MaterialStream materialStream2(Flowrate = 66, Pressure = 24, Temperature = 600, molefraction = {0, 0.9, 0, 0.1}, unspecified = false) annotation(Placement(visible = true, transformation(origin = {-280, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  unitoperations.CSTR cSTR1(Ab = 0, Af = 5.1e11, Eab = 0, Eaf = 230e3, T_iso = 700, V_Total = 1, operation_mode = unitoperations.CSTR.operation_type.Isothermal, order_b = {0, 0, 0, 0}, order_f = {1, 0.5, 0, 0}) annotation(Placement(visible = true, transformation(origin = {-142.216, 18.9457}, extent = {{-19.7835, -21.582}, {19.7835, 17.985}}, rotation = 0)));
  unitoperations.MaterialStream materialStream3 annotation(Placement(visible = true, transformation(origin = {-58, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  unitoperations.valve valve1(OutletPfixed = false) annotation(Placement(visible = true, transformation(origin = {-94, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  unitoperations.PTFlash PTFlash1 annotation(Placement(visible = true, transformation(origin = {1, 13.5154}, extent = {{-23, -42.0154}, {23, 30.011}}, rotation = 0)));
  unitoperations.valve valve2(OutletPfixed = true) annotation(Placement(visible = true, transformation(origin = {52, 28}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  unitoperations.valve valve3(OutletPfixed = true) annotation(Placement(visible = true, transformation(origin = {60, -26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  unitoperations.MaterialStream materialStream4 annotation(Placement(visible = true, transformation(origin = {92, 56}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  unitoperations.MaterialStream materialStream5 annotation(Placement(visible = true, transformation(origin = {92, -26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  unitoperations.Distillation distillation1(Override_Sizing_Calculations = false, h_weir = 0.01, specification1 = unitoperations.Distillation.spec1.RefluxRatio, specification1_value = 2, specification2 = unitoperations.Distillation.spec2.ProductMolarFlow, specification2_value = 15) annotation(Placement(visible = true, transformation(origin = {147, -22.8324}, extent = {{-31, -38.6205}, {31, 27.5861}}, rotation = 0)));
  unitoperations.MaterialStream materialStream6(flashCalculations = false) annotation(Placement(visible = true, transformation(origin = {225, -9}, extent = {{-17, -17}, {17, 17}}, rotation = 0)));
  unitoperations.MaterialStream materialStream7(flashCalculations = false) annotation(Placement(visible = true, transformation(origin = {207, -43}, extent = {{-17, -17}, {17, 17}}, rotation = 0)));
  unitoperations.Mixer mixer1 annotation(Placement(visible = true, transformation(origin = {-234, 4}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
  unitoperations.MaterialStream materialStream8 annotation(Placement(visible = true, transformation(origin = {-188, 4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(materialStream8.port2, cSTR1.port1) annotation(Line(points = {{-180, 4}, {-172, 4}, {-172, 18}, {-160, 18}, {-160, 18}}));
    connect(mixer1.port3, materialStream8.port1) annotation(Line(points = {{-220, 6}, {-196, 6}, {-196, 4}, {-196, 4}}));
    connect(materialStream7.port2, mixer1.port2) annotation(Line(points = {{222, -42}, {238, -42}, {238, -84}, {-276, -84}, {-276, -6}, {-248, -6}, {-248, -4}}));
    connect(distillation1.port3, materialStream7.port1) annotation(Line(points = {{169, -43}, {193, -43}}));
    connect(distillation1.port2, materialStream6.port1) annotation(Line(points = {{169, -9}, {211, -9}}));
    connect(materialStream5.port2, distillation1.port1) annotation(Line(points = {{100.5, -26}, {122, -26}}));
    connect(valve3.port2, materialStream5.port1) annotation(Line(points = {{68, -26}, {84, -26}}));
    connect(PTFlash1.port1, valve3.port1) annotation(Line(points = {{18, -11}, {18, -26}, {52, -26}}));
    connect(valve2.port2, materialStream4.port1) annotation(Line(points = {{60, 28}, {72, 28}, {72, 56}, {84, 56}}));
    connect(PTFlash1.port2, valve2.port1) annotation(Line(points = {{19, 29}, {44, 29}, {44, 28}}));
    connect(materialStream3.port2, PTFlash1.port3) annotation(Line(points = {{-49.5, 10}, {-17, 10}}));
    connect(valve1.port2, materialStream3.port1) annotation(Line(points = {{-86, 10}, {-66, 10}}));
    connect(cSTR1.port3, valve1.port1) annotation(Line(points = {{-125, 10}, {-102, 10}}));
    connect(materialStream2.port2, cSTR1.port2) annotation(Line(points = {{-271.5, 60}, {-203, 60}, {-203, 23.25}, {-159, 23.25}, {-159, 25}}));
    connect(materialStream1.port2, mixer1.port1) annotation(Line(points = {{-270, 18}, {-248, 18}}));
    annotation(Diagram(coordinateSystem(extent = {{-300, -100}, {300, 100}}, preserveAspectRatio = false)), Icon(coordinateSystem(extent = {{-300, -100}, {300, 100}}, preserveAspectRatio = false)), version = "", uses);end Withrecycle1;
end BenzeneProduction;
