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

  model reactorFlashDistillation
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
  DistillationPackage.Distillationwithsizing distillationwithsizing1(Override_Sizing_Calculations = false, specification1 = DistillationPackage.Distillationwithsizing.spec1.RefluxRatio, specification1_value = 2, specification2 = DistillationPackage.Distillationwithsizing.spec2.ProductMolarFlow, specification2_value = 15)  annotation(Placement(visible = true, transformation(origin = {245, -67}, extent = {{-21, -21}, {21, 21}}, rotation = 0)));
  equation
    connect(materialStream5.port2, distillationwithsizing1.port1) annotation(Line(points = {{206, -68}, {227, -68}, {227, -67}}));
    connect(valve3.port2, materialStream5.port1) annotation(Line(points = {{152, -60}, {188, -60}, {188, -68}, {190, -68}}));
    connect(flashWithSizing1.port1, valve3.port1) annotation(Line(points = {{92, -44}, {94, -44}, {94, -60}, {136, -60}, {136, -60}}));
    connect(valve2.port2, materialStream4.port1) annotation(Line(points = {{154, -12}, {190, -12}, {190, -14}, {190, -14}}));
    connect(flashWithSizing1.port2, valve2.port1) annotation(Line(points = {{102, -26}, {138, -26}, {138, -12}, {138, -12}}));
    connect(materialStream3.port2, flashWithSizing1.port3) annotation(Line(points = {{54, -34}, {80, -34}, {80, -34}, {80, -34}}));
    connect(valve1.port2, materialStream3.port1) annotation(Line(points = {{22, -34}, {38, -34}}));
    connect(cSTR1.port3, valve1.port1) annotation(Line(points = {{-15, -34}, {6, -34}}));
    connect(materialStream2.port2, cSTR1.port2) annotation(Line(points = {{-71, 42}, {-29, 42}, {-29, -7}}));
    connect(materialStream1.port2, cSTR1.port1) annotation(Line(points = {{-72, -22}, {-45, -22}}));
  end reactorFlashDistillation;

  model WithRecycle
  end WithRecycle;
end BenzeneProduction;