package BenzeneProductionNewtest
  model ptflash_test
    unitoperationsModified.MaterialStream materialStream1(Flowrate = 88, molefraction = {0.25, 0.25, 0.25, 0.25}, pressure = 500000, specified_stream = true, step_value = 2, stepchange = true, temperature = 300) annotation(
      Placement(visible = true, transformation(origin = {-80, 8}, extent = {{-14, -14}, {14, 14}}, rotation = 0)));
    unitoperationsModified.PTFlash PTFlash1(Dynamic = false, L(start = 22), OverrideSizeCalculations = false, Pset = 500000, V(start = 66), connectedToInput = true, hset = 0.25) annotation(
      Placement(visible = true, transformation(origin = {-22, 10.1631}, extent = {{-30, -44.4964}, {30, 31.7832}}, rotation = 0)));
    unitoperationsModified.valve valve1(Dynamic = false,OutletPfixed = true, OutletPressure = 100000) annotation(
      Placement(visible = true, transformation(origin = {44, 32}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    unitoperationsModified.MaterialStream materialStream2 annotation(
      Placement(visible = true, transformation(origin = {80, 32}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  unitoperationsModified.valve valve2(Dynamic = false, OutletPfixed = true)  annotation(
      Placement(visible = true, transformation(origin = {34, -18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  unitoperationsModified.MaterialStream materialStream3 annotation(
      Placement(visible = true, transformation(origin = {88, -18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(valve2.port2, materialStream3.port1) annotation(
      Line(points = {{42, -18}, {78, -18}, {78, -18}, {80, -18}}));
    connect(PTFlash1.port1, valve2.port1) annotation(
      Line(points = {{0, -16}, {26, -16}, {26, -18}, {26, -18}}));
    connect(PTFlash1.port2, valve1.port1) annotation(
      Line(points = {{0, 28}, {36, 28}, {36, 32}, {36, 32}}));
    connect(valve1.port2, materialStream2.port1) annotation(
      Line(points = {{52, 32}, {70, 32}, {70, 32}, {72, 32}}));
    connect(materialStream1.port2, PTFlash1.port3) annotation(
      Line(points = {{-68, 8}, {-46, 8}}));
  end ptflash_test;


  model Reactor
    unitoperationsModified.MaterialStream toluene_feed(Flowrate = 22, molefraction = {1, 0, 0, 0}, pressure = 100000, specified_stream = true, step_value = 2, stepchange = true, stepchangetime = 0.5)  annotation(
      Placement(visible = true, transformation(origin = {-84, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    unitoperationsModified.MaterialStream hydrogenFeed(Flowrate = 66, Tdf(displayUnit = "K", start = 160), molefraction = {0, 1, 0, 0}, pressure = 2.5e+06, specified_stream = true)  annotation(
      Placement(visible = true, transformation(origin = {-76, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  unitoperationsModified.CSTR reactor(Ab = 0, Af = 5.1e11, Dynamic = false, Eab = 0, Eaf = 230e3,T_iso(displayUnit = "K") = 700, V_Total = 1, operation_mode = unitoperationsModified.types.operation_type.Isothermal, order_b = {0, 0, 0, 0}, order_f = {1, 0.5, 0, 0})  annotation(
      Placement(visible = true, transformation(origin = {-13.4704, 29.4704}, extent = {{-24.5296, -29.4355}, {24.5296, 24.5296}}, rotation = 0)));
  unitoperationsModified.valve valve1(Dynamic = false, OutletPfixed = true)  annotation(
      Placement(visible = true, transformation(origin = {44, 12}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  unitoperationsModified.MaterialStream reactorOutlet(Tdf(displayUnit = "K", start = 337))  annotation(
      Placement(visible = true, transformation(origin = {94, 12}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(valve1.port2, reactorOutlet.port1) annotation(
      Line(points = {{52, 12}, {84, 12}, {84, 12}, {86, 12}}));
    connect(reactor.port3, valve1.port1) annotation(
      Line(points = {{8, 16}, {34, 16}, {34, 12}, {36, 12}}));
    connect(toluene_feed.port2, reactor.port1) annotation(
      Line(points = {{-76, 0}, {-56, 0}, {-56, 28}, {-34, 28}, {-34, 28}}));
    connect(hydrogenFeed.port2, reactor.port2) annotation(
      Line(points = {{-68, 42}, {-34, 42}, {-34, 38}, {-34, 38}}));
  
  end Reactor;

  model Distillationtest
    unitoperationsModified.MaterialStream feed(Flowrate = 80, molefraction = {0.5, 0, 0.5, 0}, pressure = 100000, specified_stream = true, step_value = 2, stepchange = true, stepchangetime = 0.5)  annotation(
      Placement(visible = true, transformation(origin = {-78, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    unitoperationsModified.Distillation distillation1(Dynamic = true, Override_Sizing_Calculations = false, P_condenser = 100000, Pressure_drop = 700, specification1 = unitoperationsModified.types.Distillation_spec1.RefluxRatio, specification1_value = 2, specification2 = unitoperationsModified.types.Distillation_spec2.ProductMolarFlow, specification2_value = 30)  annotation(
      Placement(visible = true, transformation(origin = {-33, 0.787066}, extent = {{-16, -29.6981}, {16, 21.2129}}, rotation = 0)));
  unitoperationsModified.MaterialStream distillate annotation(
      Placement(visible = true, transformation(origin = {20, 18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  unitoperationsModified.MaterialStream bottoms annotation(
      Placement(visible = true, transformation(origin = {20, -28}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(distillation1.port3, bottoms.port1) annotation(
      Line(points = {{-22, -14}, {-4, -14}, {-4, -28}, {12, -28}, {12, -28}}));
    connect(distillation1.port2, distillate.port1) annotation(
      Line(points = {{-22, 12}, {10, 12}, {10, 18}, {12, 18}}));
    connect(feed.port2, distillation1.port1) annotation(
      Line(points = {{-70, 0}, {-46, 0}, {-46, -2}, {-46, -2}}));
  
  end Distillationtest;





end BenzeneProductionNewtest;
