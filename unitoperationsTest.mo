package unitoperationsTest
  model materialStreamTest
    unitoperations.MaterialStream materialStream1(Flowrate = 100, Pressure = 5, Temperature = 350, molefraction = {0.5, 0, 0.5, 0}, unspecified = false) annotation(Placement(visible = true, transformation(origin = {-24, 22}, extent = {{-24, -24}, {24, 24}}, rotation = 0)));
  end materialStreamTest;

  model CSTR_DynamicTest
    unitoperations.MaterialStream materialStream1(Flowrate = 100, Pressure = 24, Tbf(start = 350), Temperature = 600, molefraction = {0.9, 0, 0.1, 0}, step_value = 20, stepchange = false, stepchangetime = 0.01, unspecified = false) annotation(Placement(visible = true, transformation(origin = {-87, -1}, extent = {{-19, -19}, {19, 19}}, rotation = 0)));
    unitoperations.MaterialStream materialStream2(Flowrate = 100, Pressure = 24, Temperature = 350, molefraction = {0, 0.9, 0, 0.1}, unspecified = false) annotation(Placement(visible = true, transformation(origin = {-45, 65}, extent = {{-17, -17}, {17, 17}}, rotation = 0)));
    unitoperations.valve valve1(Control = false, OutletPfixed = true, OutletPressure = 5, valveCv = 0.4) annotation(Placement(visible = true, transformation(origin = {56, -20}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
    unitoperations.MaterialStream materialStream3 annotation(Placement(visible = true, transformation(origin = {110, -22}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    unitoperations.CSTR cSTR1(Ab = 0, Af = 5.1e11, Eab = 0, Eaf = 230e3, V_Total = 4, operation_mode = unitoperations.CSTR.operation_type.Isothermal, order_b = {0, 0, 0, 0}, order_f = {1, 0.5, 0, 0}) annotation(Placement(visible = true, transformation(origin = {-14, -10}, extent = {{-22, -22}, {22, 22}}, rotation = 0)));
  equation
    connect(cSTR1.port3, valve1.port1) annotation(Line(points = {{4, -26}, {38, -26}, {38, -20}, {40, -20}}));
    connect(materialStream2.port2, cSTR1.port2) annotation(Line(points = {{-30, 66}, {-14, 66}, {-14, 8}, {-14, 8}}));
    connect(materialStream1.port2, cSTR1.port1) annotation(Line(points = {{-70, 0}, {-34, 0}, {-34, -10}, {-32, -10}}));
    connect(valve1.port2, materialStream3.port1) annotation(Line(points = {{72, -20}, {102, -20}, {102, -22}, {102, -22}}));
  end CSTR_DynamicTest;

  model dynFlashtest
    unitoperations.MaterialStream materialStream1(Flowrate = 100, Pressure = 10, Temperature = 350, molefraction = {0.25, 0.25, 0.25, 0.25}, stepchange = true, unspecified = false) annotation(Placement(visible = true, transformation(origin = {-78, -2}, extent = {{-18, -18}, {18, 18}}, rotation = 0)));
    unitoperations.Flash flash1(connectedToInput = true) annotation(Placement(visible = true, transformation(origin = {-5, -1}, extent = {{-29, -29}, {29, 29}}, rotation = 0)));
    unitoperations.MaterialStream materialStream2 annotation(Placement(visible = true, transformation(origin = {101, 27}, extent = {{-19, -19}, {19, 19}}, rotation = 0)));
    unitoperations.MaterialStream materialStream3 annotation(Placement(visible = true, transformation(origin = {84, -38}, extent = {{-24, -24}, {24, 24}}, rotation = 0)));
    unitoperations.valve valve1(OutletPfixed = true, OutletPressure = 1) annotation(Placement(visible = true, transformation(origin = {51, 25}, extent = {{-17, -17}, {17, 17}}, rotation = 0)));
    unitoperations.valve valve2(OutletPfixed = true, OutletPressure = 1) annotation(Placement(visible = true, transformation(origin = {26, -48}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(materialStream1.port2, flash1.port3) annotation(Line(points = {{-62, -2}, {-28, -2}, {-28, 0}, {-28, 0}}));
    connect(valve2.port2, materialStream3.port1) annotation(Line(points = {{34, -48}, {64, -48}, {64, -38}, {64, -38}}));
    connect(flash1.port1, valve2.port1) annotation(Line(points = {{-4, -22}, {-4, -22}, {-4, -48}, {18, -48}, {18, -48}}));
    connect(valve1.port2, materialStream2.port1) annotation(Line(points = {{64, 26}, {84, 26}, {84, 28}, {86, 28}}));
    connect(flash1.port2, valve1.port1) annotation(Line(points = {{18, 14}, {36, 14}, {36, 24}, {38, 24}}));
  end dynFlashtest;

  model reactor_Flash_test
    unitoperations.MaterialStream materialStream1(Flowrate = 100, Pressure = 25e5, Temperature = 600, molefraction = {0.9, 0, 0.1, 0}, stepchange = true, unspecified = false) annotation(Placement(visible = true, transformation(origin = {-445, -30}, extent = {{-40, -40}, {40, 40}}, rotation = 0)));
    unitoperations.CSTR cSTR1(operation_mode = unitoperations.CSTR.operation_type.Isothermal)  annotation(Placement(visible = true, transformation(origin = {-310, -30}, extent = {{-40, -40}, {40, 40}}, rotation = 0)));
    unitoperations.MaterialStream materialStream2(Flowrate = 100, Pressure = 25e5, Temperature = 400, molefraction = {0, 0.9, 0, 0.1}, stepchange = false, unspecified = false) annotation(Placement(visible = true, transformation(origin = {-447.5, 77.5}, extent = {{-47.5, -47.5}, {47.5, 47.5}}, rotation = 0)));
    unitoperations.valve valve1(OutletPfixed = false) annotation(Placement(visible = true, transformation(origin = {-200, -60}, extent = {{-30, -30}, {30, 30}}, rotation = 0)));
    unitoperations.MaterialStream materialStream3 annotation(Placement(visible = true, transformation(origin = {-90, -60}, extent = {{-30, -30}, {30, 30}}, rotation = 0)));
    unitoperations.Flash flash1 annotation(Placement(visible = true, transformation(origin = {132.5, -62.5}, extent = {{-52.5, -52.5}, {52.5, 52.5}}, rotation = 0)));
    unitoperations.valve valve2(OutletPfixed = true, OutletPressure = 1e5) annotation(Placement(visible = true, transformation(origin = {262.5, 47.5}, extent = {{-37.5, -37.5}, {37.5, 37.5}}, rotation = 0)));
    unitoperations.valve valve3(OutletPfixed = true, OutletPressure = 1e5) annotation(Placement(visible = true, transformation(origin = {215, -155}, extent = {{-35, -35}, {35, 35}}, rotation = 0)));
    unitoperations.MaterialStream materialStream4 annotation(Placement(visible = true, transformation(origin = {367.5, 47.5}, extent = {{-47.5, -47.5}, {47.5, 47.5}}, rotation = 0)));
    unitoperations.MaterialStream materialStream5 annotation(Placement(visible = true, transformation(origin = {375, -155}, extent = {{-30, -30}, {30, 30}}, rotation = 0)));
    unitoperations.HeatExchanger heatExchanger1 annotation(Placement(visible = true, transformation(origin = {-22.5, -57.5}, extent = {{-22.5, -22.5}, {22.5, 22.5}}, rotation = 0)));
    unitoperations.MaterialStream materialStream6(zl(start = {0.442, 7e-17, 0.5547, 0.0033}), zv(start = {0.0062, 0.424, 0.0245, 0.5453})) annotation(Placement(visible = true, transformation(origin = {45, -60}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
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

  model reactor_Flash_Seperator
    unitoperations.MaterialStream materialStream1(Flowrate = 100, Pressure = 25e5, Temperature = 600, molefraction = {0.9, 0, 0.1, 0}, stepchange = true, unspecified = false) annotation(Placement(visible = true, transformation(origin = {-445, -30}, extent = {{-40, -40}, {40, 40}}, rotation = 0)));
    unitoperations.CSTR cSTR1 annotation(Placement(visible = true, transformation(origin = {-310, -30}, extent = {{-40, -40}, {40, 40}}, rotation = 0)));
    unitoperations.MaterialStream materialStream2(Flowrate = 140, Pressure = 25e5, Temperature = 400, molefraction = {0, 0.9, 0, 0.1}, stepchange = false, unspecified = false) annotation(Placement(visible = true, transformation(origin = {-447.5, 77.5}, extent = {{-47.5, -47.5}, {47.5, 47.5}}, rotation = 0)));
    unitoperations.valve valve1(Control = false, OutletPfixed = false, valveCv = 0.4) annotation(Placement(visible = true, transformation(origin = {-200, -60}, extent = {{-30, -30}, {30, 30}}, rotation = 0)));
    unitoperations.MaterialStream materialStream3 annotation(Placement(visible = true, transformation(origin = {-90, -60}, extent = {{-30, -30}, {30, 30}}, rotation = 0)));
    unitoperations.Flash flash1 annotation(Placement(visible = true, transformation(origin = {132.5, -62.5}, extent = {{-52.5, -52.5}, {52.5, 52.5}}, rotation = 0)));
    unitoperations.valve valve2(Control = false, OutletPfixed = true, OutletPressure = 1e5) annotation(Placement(visible = true, transformation(origin = {262.5, 47.5}, extent = {{-37.5, -37.5}, {37.5, 37.5}}, rotation = 0)));
    unitoperations.MaterialStream materialStream4 annotation(Placement(visible = true, transformation(origin = {367.5, 47.5}, extent = {{-47.5, -47.5}, {47.5, 47.5}}, rotation = 0)));
    unitoperations.MaterialStream materialStream5 annotation(Placement(visible = true, transformation(origin = {212.5, -157.5}, extent = {{-32.5, -32.5}, {32.5, 32.5}}, rotation = 0)));
    unitoperations.HeatExchanger heatExchanger1 annotation(Placement(visible = true, transformation(origin = {-22.5, -57.5}, extent = {{-22.5, -22.5}, {22.5, 22.5}}, rotation = 0)));
    unitoperations.MaterialStream materialStream6(zl(start = {0.7, 1e-16, 0.3, 0.001}), zv(start = {0.01, 0.75, 0.0145, 0.23})) annotation(Placement(visible = true, transformation(origin = {45, -60}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
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
    unitoperations.MaterialStream materialStream1(Flowrate = 100, Pressure = 5e5, Temperature = 310, molefraction = {0.5, 0.001, 0.4, 0.099}, stepchange = false, unspecified = false) annotation(Placement(visible = true, transformation(origin = {-68, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    unitoperations.CompoundSeperator compoundSeperator1 annotation(Placement(visible = true, transformation(origin = {-26, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    unitoperations.MaterialStream materialStream2 annotation(Placement(visible = true, transformation(origin = {54, 24}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    unitoperations.MaterialStream materialStream3 annotation(Placement(visible = true, transformation(origin = {16, -18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    unitoperations.valve valve1(Control = false, OutletPfixed = true, OutletPressure = 2e5) annotation(Placement(visible = true, transformation(origin = {10, 24}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(valve1.port2, materialStream2.port1) annotation(Line(points = {{18, 24}, {44, 24}, {44, 24}, {46, 24}}));
    connect(compoundSeperator1.port2, valve1.port1) annotation(Line(points = {{-16, 8}, {-14, 8}, {-14, 22}, {2, 22}, {2, 24}}));
    connect(materialStream1.port2, compoundSeperator1.port1) annotation(Line(points = {{-60, 0}, {-34, 0}, {-34, 0}, {-34, 0}}));
    connect(compoundSeperator1.port3, materialStream3.port1) annotation(Line(points = {{-16, -8}, {-16, -8}, {-16, -18}, {8, -18}, {8, -18}}));
  end separator_test;

  model Distillationtest
    unitoperations.MaterialStream materialStream1(Flowrate = 98, Pressure = 1, Temperature = 365, molefraction = {0.45, 0, 0.55, 0}, unspecified = false) annotation(Placement(visible = true, transformation(origin = {-70, 4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    unitoperations.MaterialStream materialStream2 annotation(Placement(visible = true, transformation(origin = {41, 27}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    unitoperations.Distillation distillation1(Override_Sizing_Calculations = false, specification1 = unitoperations.Distillation.spec1.RefluxRatio, specification1_value = 2, specification2 = unitoperations.Distillation.spec2.CompoundMolarFlow)  annotation(Placement(visible = true, transformation(origin = {-14, 6}, extent = {{-22, -22}, {22, 22}}, rotation = 0)));
  equation
    connect(materialStream1.port2, distillation1.port1) annotation(Line(points = {{-62, 4}, {-34, 4}, {-34, 6}, {-34, 6}}));
    connect(distillation1.port2, materialStream2.port1) annotation(Line(points = {{6, 24}, {30, 24}, {30, 28}, {32, 28}}));
  end Distillationtest;

  model flowsheet
    unitoperations.MaterialStream materialStream1(Flowrate = 100, Pressure = 25e5, Temperature = 600, molefraction = {0.9, 0, 0.1, 0}, step_value = 2, stepchange = true, stepchangetime = 0.01, unspecified = false) annotation(Placement(visible = true, transformation(origin = {-445, -30}, extent = {{-40, -40}, {40, 40}}, rotation = 0)));
    unitoperations.CSTR cSTR1 annotation(Placement(visible = true, transformation(origin = {-310, -30}, extent = {{-40, -40}, {40, 40}}, rotation = 0)));
    unitoperations.MaterialStream materialStream2(Flowrate = 140, Pressure = 25e5, Temperature = 400, molefraction = {0, 0.9, 0, 0.1}, stepchange = false, unspecified = false) annotation(Placement(visible = true, transformation(origin = {-447.5, 77.5}, extent = {{-47.5, -47.5}, {47.5, 47.5}}, rotation = 0)));
    unitoperations.valve valve1(Control = false, OutletPfixed = false, valveCv = 0.4) annotation(Placement(visible = true, transformation(origin = {-200, -60}, extent = {{-30, -30}, {30, 30}}, rotation = 0)));
    unitoperations.MaterialStream materialStream3 annotation(Placement(visible = true, transformation(origin = {-90, -60}, extent = {{-30, -30}, {30, 30}}, rotation = 0)));
    unitoperations.Flash flash1 annotation(Placement(visible = true, transformation(origin = {132.5, -62.5}, extent = {{-52.5, -52.5}, {52.5, 52.5}}, rotation = 0)));
    unitoperations.valve valve2(Control = false, OutletPfixed = true, OutletPressure = 1e5) annotation(Placement(visible = true, transformation(origin = {262.5, 47.5}, extent = {{-37.5, -37.5}, {37.5, 37.5}}, rotation = 0)));
    unitoperations.MaterialStream materialStream4 annotation(Placement(visible = true, transformation(origin = {367.5, 47.5}, extent = {{-47.5, -47.5}, {47.5, 47.5}}, rotation = 0)));
    unitoperations.MaterialStream materialStream5 annotation(Placement(visible = true, transformation(origin = {212.5, -157.5}, extent = {{-32.5, -32.5}, {32.5, 32.5}}, rotation = 0)));
    unitoperations.HeatExchanger heatExchanger1 annotation(Placement(visible = true, transformation(origin = {-22.5, -57.5}, extent = {{-22.5, -22.5}, {22.5, 22.5}}, rotation = 0)));
    unitoperations.MaterialStream materialStream6(zl(start = {0.7, 1e-16, 0.3, 0.001}), zv(start = {0.01, 0.75, 0.0145, 0.23})) annotation(Placement(visible = true, transformation(origin = {45, -60}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
    unitoperations.CompoundSeperator compoundSeperator1 annotation(Placement(visible = true, transformation(origin = {315, -160}, extent = {{-35, -35}, {35, 35}}, rotation = 0)));
    unitoperations.MaterialStream materialStream7 annotation(Placement(visible = true, transformation(origin = {502.5, -47.5}, extent = {{-27.5, -27.5}, {27.5, 27.5}}, rotation = 0)));
    unitoperations.MaterialStream materialStream8 annotation(Placement(visible = true, transformation(origin = {530, -180}, extent = {{-30, -30}, {30, 30}}, rotation = 0)));
    unitoperations.valve valve4(Control = false, OutletPfixed = false, OutletPressure = 1e5) annotation(Placement(visible = true, transformation(origin = {435, -135}, extent = {{-30, -30}, {30, 30}}, rotation = 0)));
    unitoperations.valve valve5(Control = false, OutletPfixed = true, OutletPressure = 1e5) annotation(Placement(visible = true, transformation(origin = {437.5, -182.5}, extent = {{-27.5, -27.5}, {27.5, 27.5}}, rotation = 0)));
    unitoperations.Distillation distillation1 annotation(Placement(visible = true, transformation(origin = {600, -40}, extent = {{-40, -40}, {40, 40}}, rotation = 0)));
    unitoperations.MaterialStream materialStream9 annotation(Placement(visible = true, transformation(origin = {712.5, 7.5}, extent = {{-27.5, -27.5}, {27.5, 27.5}}, rotation = 0)));
    unitoperations.MaterialStream materialStream10(Fl(start = 0), Fv(start = 50), Tbf(start = 392), Tdf(start = 395), zl(start = {0.9, 0, 0.1, 0}), zv(start = {0.86, 0, 0.14, 0})) annotation(Placement(visible = true, transformation(origin = {715, -85}, extent = {{-30, -30}, {30, 30}}, rotation = 0)));
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

  model flowsheet1
    unitoperations.MaterialStream materialStream1(Flowrate = 100, Pressure = 25e5, Temperature = 600, molefraction = {0.9, 0, 0.1, 0}, step_value = 2, stepchange = true, stepchangetime = 0.01, unspecified = false) annotation(Placement(visible = true, transformation(origin = {-445, -30}, extent = {{-40, -40}, {40, 40}}, rotation = 0)));
    unitoperations.CSTR cSTR1(Ab = 0, Af = 5.1e11, Eab = 0, Eaf = 230e3, V_Total = 4, operation_mode = modelicatest.CSTR.operation_type.Isothermal, order_b = {0, 0, 0, 0}, order_f = {1, 0.5, 0, 0}) annotation(Placement(visible = true, transformation(origin = {-310, -30}, extent = {{-40, -40}, {40, 40}}, rotation = 0)));
    unitoperations.MaterialStream materialStream2(Flowrate = 140, Pressure = 25e5, Temperature = 400, molefraction = {0, 0.9, 0, 0.1}, stepchange = false, unspecified = false) annotation(Placement(visible = true, transformation(origin = {-447.5, 77.5}, extent = {{-47.5, -47.5}, {47.5, 47.5}}, rotation = 0)));
    unitoperations.valve valve1(Control = false, OutletPfixed = false, valveCv = 0.4) annotation(Placement(visible = true, transformation(origin = {-200, -60}, extent = {{-30, -30}, {30, 30}}, rotation = 0)));
    unitoperations.MaterialStream materialStream3 annotation(Placement(visible = true, transformation(origin = {-90, -60}, extent = {{-30, -30}, {30, 30}}, rotation = 0)));
    unitoperations.Flash flash1 annotation(Placement(visible = true, transformation(origin = {132.5, -62.5}, extent = {{-52.5, -52.5}, {52.5, 52.5}}, rotation = 0)));
    unitoperations.valve valve2(Control = false, OutletPfixed = true, OutletPressure = 1e5) annotation(Placement(visible = true, transformation(origin = {262.5, 47.5}, extent = {{-37.5, -37.5}, {37.5, 37.5}}, rotation = 0)));
    unitoperations.MaterialStream materialStream4 annotation(Placement(visible = true, transformation(origin = {367.5, 47.5}, extent = {{-47.5, -47.5}, {47.5, 47.5}}, rotation = 0)));
    unitoperations.MaterialStream materialStream5 annotation(Placement(visible = true, transformation(origin = {212.5, -157.5}, extent = {{-32.5, -32.5}, {32.5, 32.5}}, rotation = 0)));
    unitoperations.HeatExchanger heatExchanger1 annotation(Placement(visible = true, transformation(origin = {-22.5, -57.5}, extent = {{-22.5, -22.5}, {22.5, 22.5}}, rotation = 0)));
    unitoperations.MaterialStream materialStream6(zl(start = {0.7, 1e-16, 0.3, 0.001}), zv(start = {0.01, 0.75, 0.0145, 0.23})) annotation(Placement(visible = true, transformation(origin = {45, -60}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
    unitoperations.CompoundSeperator compoundSeperator1 annotation(Placement(visible = true, transformation(origin = {315, -160}, extent = {{-35, -35}, {35, 35}}, rotation = 0)));
    unitoperations.MaterialStream materialStream7 annotation(Placement(visible = true, transformation(origin = {502.5, -47.5}, extent = {{-27.5, -27.5}, {27.5, 27.5}}, rotation = 0)));
    unitoperations.MaterialStream materialStream8 annotation(Placement(visible = true, transformation(origin = {530, -180}, extent = {{-30, -30}, {30, 30}}, rotation = 0)));
    unitoperations.valve valve4(Control = false, OutletPfixed = false, OutletPressure = 1e5) annotation(Placement(visible = true, transformation(origin = {435, -135}, extent = {{-30, -30}, {30, 30}}, rotation = 0)));
    unitoperations.valve valve5(Control = false, OutletPfixed = true, OutletPressure = 1e5) annotation(Placement(visible = true, transformation(origin = {437.5, -182.5}, extent = {{-27.5, -27.5}, {27.5, 27.5}}, rotation = 0)));
    unitoperations.Distillation1 distillation11 annotation(Placement(visible = true, transformation(origin = {630, -55}, extent = {{-35, -35}, {35, 35}}, rotation = 0)));
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

  model FlashWithSizingTest
    unitoperations.MaterialStream materialStream1(Flowrate = 200, Pressure = 10, Tdf(start = 400), Temperature = 300, molefraction = {0.25, 0.25, 0.25, 0.25}, step_value = 1, stepchange = true, unspecified = false) annotation(Placement(visible = true, transformation(origin = {-78, -2}, extent = {{-18, -18}, {18, 18}}, rotation = 0)));
    unitoperations.FlashWithSizing flash1(OverrideSizeCalculations = false, V(start = 100), connectedToInput = true) annotation(Placement(visible = true, transformation(origin = {-5, -1}, extent = {{-29, -29}, {29, 29}}, rotation = 0)));
    unitoperations.MaterialStream materialStream2(Tdf(start = 160)) annotation(Placement(visible = true, transformation(origin = {101, 27}, extent = {{-19, -19}, {19, 19}}, rotation = 0)));
    unitoperations.MaterialStream materialStream3(Tdf(start = 370)) annotation(Placement(visible = true, transformation(origin = {84, -38}, extent = {{-24, -24}, {24, 24}}, rotation = 0)));
    unitoperations.valve valve1(OutletPfixed = true, OutletPressure = 1) annotation(Placement(visible = true, transformation(origin = {51, 25}, extent = {{-17, -17}, {17, 17}}, rotation = 0)));
    unitoperations.valve valve2(OutletPfixed = true, OutletPressure = 1) annotation(Placement(visible = true, transformation(origin = {26, -48}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(materialStream1.port2, flash1.port3) annotation(Line(points = {{-62, -2}, {-28, -2}, {-28, 0}, {-28, 0}}));
    connect(valve2.port2, materialStream3.port1) annotation(Line(points = {{34, -48}, {64, -48}, {64, -38}, {64, -38}}));
    connect(flash1.port1, valve2.port1) annotation(Line(points = {{-4, -22}, {-4, -22}, {-4, -48}, {18, -48}, {18, -48}}));
    connect(valve1.port2, materialStream2.port1) annotation(Line(points = {{64, 26}, {84, 26}, {84, 28}, {86, 28}}));
    connect(flash1.port2, valve1.port1) annotation(Line(points = {{18, 14}, {36, 14}, {36, 24}, {38, 24}}));
  end FlashWithSizingTest;

  model PhFlashWithSizingTest
    unitoperations.MaterialStream materialStream1(Flowrate = 200, Pressure = 10, Temperature = 300, molefraction = {0.25, 0.25, 0.25, 0.25}, step_value = 2, stepchange = true, unspecified = false) annotation(Placement(visible = true, transformation(origin = {-78, -2}, extent = {{-18, -18}, {18, 18}}, rotation = 0)));
    unitoperations.PhFlashWithSizing flash1(OverrideSizeCalculations = false, connectedToInput = true) annotation(Placement(visible = true, transformation(origin = {-5, -1}, extent = {{-29, -29}, {29, 29}}, rotation = 0)));
    unitoperations.MaterialStream materialStream2(Tdf(start = 160)) annotation(Placement(visible = true, transformation(origin = {101, 27}, extent = {{-19, -19}, {19, 19}}, rotation = 0)));
    unitoperations.MaterialStream materialStream3(Tbf(start = 350)) annotation(Placement(visible = true, transformation(origin = {84, -38}, extent = {{-24, -24}, {24, 24}}, rotation = 0)));
    unitoperations.valve valve1(OutletPfixed = true, OutletPressure = 1) annotation(Placement(visible = true, transformation(origin = {51, 25}, extent = {{-17, -17}, {17, 17}}, rotation = 0)));
    unitoperations.valve valve2(OutletPfixed = true, OutletPressure = 1) annotation(Placement(visible = true, transformation(origin = {26, -48}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(materialStream1.port2, flash1.port3) annotation(Line(points = {{-62, -2}, {-28, -2}, {-28, 0}, {-28, 0}}));
    connect(valve2.port2, materialStream3.port1) annotation(Line(points = {{34, -48}, {64, -48}, {64, -38}, {64, -38}}));
    connect(flash1.port1, valve2.port1) annotation(Line(points = {{-4, -22}, {-4, -22}, {-4, -48}, {18, -48}, {18, -48}}));
    connect(valve1.port2, materialStream2.port1) annotation(Line(points = {{64, 26}, {84, 26}, {84, 28}, {86, 28}}));
    connect(flash1.port2, valve1.port1) annotation(Line(points = {{18, 14}, {36, 14}, {36, 24}, {38, 24}}));
  end PhFlashWithSizingTest;
end unitoperationsTest;
