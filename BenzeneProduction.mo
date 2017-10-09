package BenzeneProduction
  model Reactor
    unitoperations.MaterialStream materialStream1(Flowrate = 22, Pressure = 24, Temperature = 600, molefraction = {0.9, 0, 0.1, 0}, step_value = 2, stepchange = true, unspecified = false) annotation(Placement(visible = true, transformation(origin = {-82, -22}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
    unitoperations.MaterialStream materialStream2(Flowrate = 66, Pressure = 24, Temperature = 600, molefraction = {0, 0.9, 0, 0.1}, unspecified = false) annotation(Placement(visible = true, transformation(origin = {-80, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    unitoperations.CSTR cSTR1(Ab = 0, Af = 5.1e11, Eab = 0, Eaf = 230e3, T_iso = 700, V_Total = 1, operation_mode = unitoperations.CSTR.operation_type.Isothermal, order_b = {0, 0, 0, 0}, order_f = {1, 0.5, 0, 0})  annotation(Placement(visible = true, transformation(origin = {6, -22}, extent = {{-18, -18}, {18, 18}}, rotation = 0)));
    unitoperations.MaterialStream materialStream3 annotation(Placement(visible = true, transformation(origin = {116, -34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    unitoperations.valve valve1(OutletPfixed = true) annotation(Placement(visible = true, transformation(origin = {64, -34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(valve1.port2, materialStream3.port1) annotation(Line(points = {{72, -34}, {108, -34}}));
    connect(cSTR1.port3, valve1.port1) annotation(Line(points = {{21, -34}, {56, -34}}));
    connect(materialStream1.port2, cSTR1.port1) annotation(Line(points = {{-72, -22}, {-9, -22}}));
    connect(materialStream2.port2, cSTR1.port2) annotation(Line(points = {{-71, 42}, {7, 42}, {7, -7}}));
  end Reactor;

  model Flashafterreactor
  
    unitoperations.MaterialStream materialStream1(Flowrate = 88, Pressure = 24, Temperature = 300, molefraction = {0.175, 0.625, 0.75, 0.125}, unspecified = false)  annotation(Placement(visible = true, transformation(origin = {-131, 5}, extent = {{-23, -23}, {23, 23}}, rotation = 0)));
    unitoperations.PhFlashWithSizing phFlashWithSizing1(L(start = 22), OverrideSizeCalculations = false, Ti(start = 280), V(start = 66), connectedToInput = true)  annotation(Placement(visible = true, transformation(origin = {-56, 6}, extent = {{-26, -26}, {26, 26}}, rotation = 0)));
  unitoperations.MaterialStream materialStream2 annotation(Placement(visible = true, transformation(origin = {73, 27}, extent = {{-19, -19}, {19, 19}}, rotation = 0)));
  unitoperations.valve valve1(OutletPfixed = true)  annotation(Placement(visible = true, transformation(origin = {11, 25}, extent = {{-19, -19}, {19, 19}}, rotation = 0)));
  unitoperations.valve valve2(OutletPfixed = true)  annotation(Placement(visible = true, transformation(origin = {15, -37}, extent = {{-17, -17}, {17, 17}}, rotation = 0)));
  unitoperations.MaterialStream materialStream3(Tdf(start = 353))  annotation(Placement(visible = true, transformation(origin = {82, -38}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  equation
    connect(valve2.port2, materialStream3.port1) annotation(Line(points = {{28, -36}, {64, -36}, {64, -38}, {66, -38}}));
    connect(phFlashWithSizing1.port1, valve2.port1) annotation(Line(points = {{-56, -16}, {2, -16}, {2, -38}, {2, -38}}));
    connect(valve1.port2, materialStream2.port1) annotation(Line(points = {{26, 26}, {56, 26}, {56, 28}, {58, 28}}));
    connect(phFlashWithSizing1.port2, valve1.port1) annotation(Line(points = {{-36, 20}, {-4, 20}, {-4, 24}, {-4, 24}}));
    connect(materialStream1.port2, phFlashWithSizing1.port3) annotation(Line(points = {{-112, 6}, {-76, 6}, {-76, 8}, {-76, 8}}));
  end Flashafterreactor;

  model PTflashAfterReactor
  
    unitoperations.MaterialStream materialStream1(Flowrate = 88, Pressure = 24, Temperature = 300, molefraction = {0.175, 0.625, 0.75, 0.125}, unspecified = false)  annotation(Placement(visible = true, transformation(origin = {-131, 5}, extent = {{-23, -23}, {23, 23}}, rotation = 0)));
  unitoperations.MaterialStream materialStream2 annotation(Placement(visible = true, transformation(origin = {73, 27}, extent = {{-19, -19}, {19, 19}}, rotation = 0)));
  unitoperations.valve valve1(OutletPfixed = true)  annotation(Placement(visible = true, transformation(origin = {11, 25}, extent = {{-19, -19}, {19, 19}}, rotation = 0)));
  unitoperations.valve valve2(OutletPfixed = true) annotation(Placement(visible = true, transformation(origin = {17, -39}, extent = {{-17, -17}, {17, 17}}, rotation = 0)));
  unitoperations.MaterialStream materialStream3(Tbf(start = 300), Tdf(start = 350))  annotation(Placement(visible = true, transformation(origin = {82, -38}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  unitoperations.FlashWithSizing flashWithSizing1( Area = 1,L(start = 22), MG(start = 400), ML(start = 700), OverrideSizeCalculations = true, Pset = 5, Ti = 300, V(start = 60), Volume = 2, connectedToInput = true, hset = 0.2, x(start = {0.7, 1e-16, 0.28, 0.001})) annotation(Placement(visible = true, transformation(origin = {-52, 6}, extent = {{-22, -22}, {22, 22}}, rotation = 0)));
  equation
    connect(materialStream1.port2, flashWithSizing1.port3) annotation(Line(points = {{-112, 6}, {-94, 6}, {-94, 7}, {-70, 7}}));
    connect(flashWithSizing1.port2, valve1.port1) annotation(Line(points = {{-35, 17}, {-6, 17}, {-6, 24}, {-4, 24}}));
    connect(flashWithSizing1.port1, valve2.port1) annotation(Line(points = {{-52, -10}, {-52, -39}, {3, -39}}));
    connect(valve2.port2, materialStream3.port1) annotation(Line(points = {{31, -39}, {64, -39}, {64, -38}, {66, -38}}));
    connect(valve1.port2, materialStream2.port1) annotation(Line(points = {{26, 26}, {56, 26}, {56, 28}, {58, 28}}));
  end PTflashAfterReactor;

  model reactorAndFlash
    unitoperations.MaterialStream materialStream1(Flowrate = 22, Pressure = 24, Temperature = 600, molefraction = {0.9, 0, 0.1, 0}, unspecified = false) annotation(Placement(visible = true, transformation(origin = {-82, -22}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
    unitoperations.MaterialStream materialStream2(Flowrate = 66, Pressure = 24, Temperature = 600, molefraction = {0, 0.9, 0, 0.1}, unspecified = false) annotation(Placement(visible = true, transformation(origin = {-80, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  unitoperations.CSTR cSTR1(Ab = 0, Af = 5.1e11, Eab = 0, Eaf = 230e3, T_iso = 700, V_Total = 1, operation_mode = unitoperations.CSTR.operation_type.Isothermal, order_b = {0, 0, 0, 0}, order_f = {1, 0.5, 0, 0}) annotation(Placement(visible = true, transformation(origin = {-30, -22}, extent = {{-18, -18}, {18, 18}}, rotation = 0)));
  unitoperations.MaterialStream materialStream3 annotation(Placement(visible = true, transformation(origin = {46, -34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  unitoperations.valve valve1(OutletPfixed = false) annotation(Placement(visible = true, transformation(origin = {14, -34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  unitoperations.FlashWithSizing flashWithSizing1 annotation(Placement(visible = true, transformation(origin = {92, -34}, extent = {{-14, -14}, {14, 14}}, rotation = 0)));
  unitoperations.valve valve2(OutletPfixed = true)  annotation(Placement(visible = true, transformation(origin = {146, -12}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  unitoperations.valve valve3(OutletPfixed = true)  annotation(Placement(visible = true, transformation(origin = {144, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  unitoperations.MaterialStream materialStream4 annotation(Placement(visible = true, transformation(origin = {198, -14}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  unitoperations.MaterialStream materialStream5 annotation(Placement(visible = true, transformation(origin = {198, -68}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(valve3.port2, materialStream5.port1) annotation(Line(points = {{152, -60}, {188, -60}, {188, -68}, {190, -68}}));
    connect(flashWithSizing1.port1, valve3.port1) annotation(Line(points = {{92, -44}, {94, -44}, {94, -60}, {136, -60}, {136, -60}}));
    connect(valve2.port2, materialStream4.port1) annotation(Line(points = {{154, -12}, {190, -12}, {190, -14}, {190, -14}}));
    connect(flashWithSizing1.port2, valve2.port1) annotation(Line(points = {{102, -26}, {138, -26}, {138, -12}, {138, -12}}));
    connect(materialStream3.port2, flashWithSizing1.port3) annotation(Line(points = {{54, -34}, {80, -34}, {80, -34}, {80, -34}}));
    connect(valve1.port2, materialStream3.port1) annotation(Line(points = {{22, -34}, {38, -34}}));
    connect(cSTR1.port3, valve1.port1) annotation(Line(points = {{-15, -34}, {6, -34}}));
    connect(materialStream2.port2, cSTR1.port2) annotation(Line(points = {{-71, 42}, {-29, 42}, {-29, -7}}));
    connect(materialStream1.port2, cSTR1.port1) annotation(Line(points = {{-72, -22}, {-45, -22}}));
  end reactorAndFlash;

  model WithoutRecycle
    unitoperations.MaterialStream materialStream1(Flowrate = 22, Pressure = 24, Temperature = 600, molefraction = {0.9, 0, 0.1, 0}, step_value = 2, stepchange = true, unspecified = false) annotation(Placement(visible = true, transformation(origin = {-82, -22}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
    unitoperations.MaterialStream materialStream2(Flowrate = 66, Pressure = 24, Temperature = 600, molefraction = {0, 0.9, 0, 0.1}, unspecified = false) annotation(Placement(visible = true, transformation(origin = {-80, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  unitoperations.CSTR cSTR1(Ab = 0, Af = 5.1e11, Eab = 0, Eaf = 230e3, T_iso = 700, V_Total = 1, operation_mode = unitoperations.CSTR.operation_type.Isothermal, order_b = {0, 0, 0, 0}, order_f = {1, 0.5, 0, 0}) annotation(Placement(visible = true, transformation(origin = {-30, -22}, extent = {{-18, -18}, {18, 18}}, rotation = 0)));
  unitoperations.MaterialStream materialStream3 annotation(Placement(visible = true, transformation(origin = {46, -34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  unitoperations.valve valve1(OutletPfixed = false) annotation(Placement(visible = true, transformation(origin = {14, -34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  unitoperations.FlashWithSizing flashWithSizing1 annotation(Placement(visible = true, transformation(origin = {92, -34}, extent = {{-14, -14}, {14, 14}}, rotation = 0)));
  unitoperations.valve valve2(OutletPfixed = true)  annotation(Placement(visible = true, transformation(origin = {146, -12}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  unitoperations.valve valve3(OutletPfixed = true)  annotation(Placement(visible = true, transformation(origin = {144, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  unitoperations.MaterialStream materialStream4 annotation(Placement(visible = true, transformation(origin = {198, -14}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  unitoperations.MaterialStream materialStream5 annotation(Placement(visible = true, transformation(origin = {198, -68}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  unitoperations.Distillation distillation1(Override_Sizing_Calculations = false, h_weir = 0.1, specification1 = unitoperations.Distillation.spec1.RefluxRatio, specification1_value = 2, specification2 = unitoperations.Distillation.spec2.ProductMolarFlow, specification2_value = 7)  annotation(Placement(visible = true, transformation(origin = {247, -69}, extent = {{-21, -21}, {21, 21}}, rotation = 0)));
  unitoperations.MaterialStream materialStream6 annotation(Placement(visible = true, transformation(origin = {308, -42}, extent = {{-14, -14}, {14, 14}}, rotation = 0)));
  unitoperations.MaterialStream materialStream7(Tbf(start = 380), Tdf(start = 380), flashCalculations = false)  annotation(Placement(visible = true, transformation(origin = {305, -89}, extent = {{-13, -13}, {13, 13}}, rotation = 0)));
  equation
    connect(distillation1.port3, materialStream7.port1) annotation(Line(points = {{266, -88}, {294, -88}, {294, -88}, {294, -88}}));
    connect(distillation1.port2, materialStream6.port1) annotation(Line(points = {{266, -50}, {270, -50}, {270, -42}, {296, -42}, {296, -42}}));
    connect(materialStream5.port2, distillation1.port1) annotation(Line(points = {{206, -68}, {228, -68}, {228, -68}, {228, -68}}));
    connect(valve3.port2, materialStream5.port1) annotation(Line(points = {{152, -60}, {188, -60}, {188, -68}, {190, -68}}));
    connect(flashWithSizing1.port1, valve3.port1) annotation(Line(points = {{92, -44}, {94, -44}, {94, -60}, {136, -60}, {136, -60}}));
    connect(valve2.port2, materialStream4.port1) annotation(Line(points = {{154, -12}, {190, -12}, {190, -14}, {190, -14}}));
    connect(flashWithSizing1.port2, valve2.port1) annotation(Line(points = {{102, -26}, {138, -26}, {138, -12}, {138, -12}}));
    connect(materialStream3.port2, flashWithSizing1.port3) annotation(Line(points = {{54, -34}, {80, -34}, {80, -34}, {80, -34}}));
    connect(valve1.port2, materialStream3.port1) annotation(Line(points = {{22, -34}, {38, -34}}));
    connect(cSTR1.port3, valve1.port1) annotation(Line(points = {{-15, -34}, {6, -34}}));
    connect(materialStream2.port2, cSTR1.port2) annotation(Line(points = {{-71, 42}, {-29, 42}, {-29, -7}}));
    connect(materialStream1.port2, cSTR1.port1) annotation(Line(points = {{-72, -22}, {-45, -22}}));
  end WithoutRecycle;

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
  unitoperations.PhFlashWithSizing phFlashWithSizing1(hset = 2) annotation(Placement(visible = true, transformation(origin = {72.7743, -0.576747}, extent = {{-21.2869, -42.4196}, {21.2869, 30.2997}}, rotation = 0)));
  unitoperations.HeatExchanger heatExchanger1(pressure_drop = 0.5) annotation(Placement(visible = true, transformation(origin = {-21, -1}, extent = {{-15, -21}, {15, 15}}, rotation = 0)));
  unitoperations.MaterialStream materialStream8 annotation(Placement(visible = true, transformation(origin = {24, -2}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
  equation
    connect(distillation1.port2, materialStream6.port1) annotation(Line(points = {{222, -22}, {261, -22}}));
    connect(distillation1.port3, materialStream7.port1) annotation(Line(points = {{222, -72}, {260, -72}}));
    connect(materialStream5.port2, distillation1.port1) annotation(Line(points = {{170, -46}, {199, -46}}));
    connect(valve3.port2, materialStream5.port1) annotation(Line(points = {{130, -46}, {150, -46}}));
    connect(phFlashWithSizing1.port1, valve3.port1) annotation(Line(points = {{83, -27}, {83, -47}, {100, -47}, {100, -46.5}, {110, -46.5}, {110, -46}}));
    connect(phFlashWithSizing1.port2, valve2.port1) annotation(Line(points = {{84, 20}, {104, 20}, {104, 16}, {108, 16}}));
    connect(materialStream8.port2, phFlashWithSizing1.port3) annotation(Line(points = {{34, -2}, {61, -2}}));
    connect(materialStream2.port2, cSTR1.port2) annotation(Line(points = {{-275.5, 46}, {-247.25, 46}, {-247.25, 38}, {-225, 38}}));
    connect(materialStream1.port2, cSTR1.port1) annotation(Line(points = {{-274, -4}, {-235.5, -4}, {-235.5, 26}, {-225, 26}}));
    connect(cSTR1.port3, valve1.port1) annotation(Line(points = {{-181, 10}, {-181, 9}, {-125, 9}, {-125, -2}}));
    connect(valve2.port2, materialStream4.port1) annotation(Line(points = {{124, 16}, {144, 16}}));
    connect(heatExchanger1.port2, materialStream8.port1) annotation(Line(points = {{-8, -2}, {14, -2}}));
    connect(materialStream3.port2, heatExchanger1.port1) annotation(Line(points = {{-54, -2}, {-34, -2}}));
    connect(valve1.port2, materialStream3.port1) annotation(Line(points = {{-99, -2}, {-77, -2}}));
    annotation(Diagram(coordinateSystem(extent = {{-300, -100}, {300, 100}}, preserveAspectRatio = false)), Icon(coordinateSystem(extent = {{-300, -100}, {300, 100}}, preserveAspectRatio = false)), version = "", uses);end WithoutrecyclePhflash;

  model WithRecycle
  unitoperations.MaterialStream materialStream1(Flowrate = 22, Pressure = 24, Temperature = 600, molefraction = {0.9, 0, 0.1, 0}, step_value = 2, stepchange = true, unspecified = false) annotation(Placement(visible = true, transformation(origin = {-280, 18}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
  unitoperations.MaterialStream materialStream2(Flowrate = 66, Pressure = 24, Temperature = 600, molefraction = {0, 0.9, 0, 0.1}, unspecified = false) annotation(Placement(visible = true, transformation(origin = {-280, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  unitoperations.CSTR cSTR1(Ab = 0, Af = 5.1e11, Eab = 0, Eaf = 230e3, T_iso = 700, V_Total = 1, operation_mode = unitoperations.CSTR.operation_type.Isothermal, order_b = {0, 0, 0, 0}, order_f = {1, 0.5, 0, 0}) annotation(Placement(visible = true, transformation(origin = {-142.216, 18.9457}, extent = {{-19.7835, -21.582}, {19.7835, 17.985}}, rotation = 0)));
  unitoperations.MaterialStream materialStream3 annotation(Placement(visible = true, transformation(origin = {-58, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  unitoperations.valve valve1(OutletPfixed = false) annotation(Placement(visible = true, transformation(origin = {-94, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  unitoperations.FlashWithSizing flashWithSizing1 annotation(Placement(visible = true, transformation(origin = {1, 13.5154}, extent = {{-23, -42.0154}, {23, 30.011}}, rotation = 0)));
  unitoperations.valve valve2(OutletPfixed = true) annotation(Placement(visible = true, transformation(origin = {52, 28}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  unitoperations.valve valve3(OutletPfixed = true) annotation(Placement(visible = true, transformation(origin = {60, -26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  unitoperations.MaterialStream materialStream4 annotation(Placement(visible = true, transformation(origin = {92, 56}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  unitoperations.MaterialStream materialStream5 annotation(Placement(visible = true, transformation(origin = {92, -26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  unitoperations.Distillation distillation1(Override_Sizing_Calculations = false, h_weir = 0.01, specification1 = unitoperations.Distillation.spec1.RefluxRatio, specification1_value = 2, specification2 = unitoperations.Distillation.spec2.ProductMolarFlow, specification2_value = 15) annotation(Placement(visible = true, transformation(origin = {147, -22.8324}, extent = {{-31, -38.6205}, {31, 27.5861}}, rotation = 0)));
  unitoperations.MaterialStream materialStream6(flashCalculations = false) annotation(Placement(visible = true, transformation(origin = {225, -9}, extent = {{-17, -17}, {17, 17}}, rotation = 0)));
  unitoperations.MaterialStream materialStream7(flashCalculations = false) annotation(Placement(visible = true, transformation(origin = {207, -43}, extent = {{-17, -17}, {17, 17}}, rotation = 0)));
  unitoperations.Mixer mixer1 annotation(Placement(visible = true, transformation(origin = {-234, 4}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
  equation
    connect(materialStream7.port2, mixer1.port2) annotation(Line(points = {{222, -42}, {238, -42}, {238, -78}, {-280, -78}, {-280, -10}, {-248, -10}, {-248, -10}}));
    connect(distillation1.port3, materialStream7.port1) annotation(Line(points = {{169, -43}, {193, -43}}));
    connect(distillation1.port2, materialStream6.port1) annotation(Line(points = {{169, -9}, {211, -9}}));
    connect(materialStream5.port2, distillation1.port1) annotation(Line(points = {{100.5, -26}, {122, -26}}));
    connect(valve3.port2, materialStream5.port1) annotation(Line(points = {{68, -26}, {84, -26}}));
    connect(flashWithSizing1.port1, valve3.port1) annotation(Line(points = {{18, -11}, {18, -26}, {52, -26}}));
    connect(valve2.port2, materialStream4.port1) annotation(Line(points = {{60, 28}, {72, 28}, {72, 56}, {84, 56}}));
    connect(flashWithSizing1.port2, valve2.port1) annotation(Line(points = {{19, 29}, {44, 29}, {44, 28}}));
    connect(materialStream3.port2, flashWithSizing1.port3) annotation(Line(points = {{-49.5, 10}, {-17, 10}}));
    connect(valve1.port2, materialStream3.port1) annotation(Line(points = {{-86, 10}, {-66, 10}}));
    connect(cSTR1.port3, valve1.port1) annotation(Line(points = {{-125, 10}, {-102, 10}}));
    connect(materialStream2.port2, cSTR1.port2) annotation(Line(points = {{-271.5, 60}, {-203, 60}, {-203, 23.25}, {-159, 23.25}, {-159, 25}}));
    connect(mixer1.port3, cSTR1.port1) annotation(Line(points = {{-220, 4}, {-197.5, 4}, {-197.5, 18}, {-159, 18}}));
    connect(materialStream1.port2, mixer1.port1) annotation(Line(points = {{-270, 18}, {-248, 18}}));
    annotation(Diagram(coordinateSystem(extent = {{-300, -100}, {300, 100}}, preserveAspectRatio = false)), Icon(coordinateSystem(extent = {{-300, -100}, {300, 100}}, preserveAspectRatio = false)), version = "", uses);end WithRecycle;
end BenzeneProduction;