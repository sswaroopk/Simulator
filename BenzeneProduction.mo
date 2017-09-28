package BenzeneProduction
  model Reactor
    unitoperations.MaterialStream materialStream1(Flowrate = 22, Pressure = 24, Temperature = 600, molefraction = {0.9, 0, 0.1, 0}, unspecified = false) annotation(Placement(visible = true, transformation(origin = {-82, -22}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
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
    unitoperations.PhFlashWithSizing phFlashWithSizing1(L(start = 22), OverrideSizeCalculations = false, Ti(start = 280), V(start = 66), connectedToInput = true, x(start = {0.7, 1e-12, 0.3, 0}))  annotation(Placement(visible = true, transformation(origin = {-56, 6}, extent = {{-26, -26}, {26, 26}}, rotation = 0)));
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
  
    unitoperations.MaterialStream materialStream1(Flowrate = 88, Pressure = 24, Temperature = 315, molefraction = {0.175, 0.625, 0.75, 0.125}, unspecified = false)  annotation(Placement(visible = true, transformation(origin = {-131, 5}, extent = {{-23, -23}, {23, 23}}, rotation = 0)));
  unitoperations.MaterialStream materialStream2 annotation(Placement(visible = true, transformation(origin = {73, 27}, extent = {{-19, -19}, {19, 19}}, rotation = 0)));
  unitoperations.valve valve1(OutletPfixed = true)  annotation(Placement(visible = true, transformation(origin = {11, 25}, extent = {{-19, -19}, {19, 19}}, rotation = 0)));
  unitoperations.valve valve2(OutletPfixed = true)  annotation(Placement(visible = true, transformation(origin = {15, -37}, extent = {{-17, -17}, {17, 17}}, rotation = 0)));
  unitoperations.MaterialStream materialStream3(Tbf(start = 350), Tdf(start = 353), zl(start = {0.2, 1e-15, 0.8, 0}))  annotation(Placement(visible = true, transformation(origin = {82, -38}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  unitoperations.FlashWithSizing flashWithSizing1(L(start = 22),OverrideSizeCalculations = false, V(start = 66), connectedToInput = true)  annotation(Placement(visible = true, transformation(origin = {-58, 6}, extent = {{-22, -22}, {22, 22}}, rotation = 0)));
  equation
    connect(flashWithSizing1.port1, valve2.port1) annotation(Line(points = {{-58, -10}, {-56, -10}, {-56, -38}, {2, -38}, {2, -38}}));
    connect(flashWithSizing1.port2, valve1.port1) annotation(Line(points = {{-42, 18}, {-6, 18}, {-6, 24}, {-4, 24}}));
    connect(materialStream1.port2, flashWithSizing1.port3) annotation(Line(points = {{-112, 6}, {-76, 6}, {-76, 6}, {-76, 6}}));
    connect(valve2.port2, materialStream3.port1) annotation(Line(points = {{28, -36}, {64, -36}, {64, -38}, {66, -38}}));
    connect(valve1.port2, materialStream2.port1) annotation(Line(points = {{26, 26}, {56, 26}, {56, 28}, {58, 28}}));
  end PTflashAfterReactor;
end BenzeneProduction;