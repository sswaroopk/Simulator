package NewTest
  model CSTR_Steady
    unitoperations.MaterialStream toluene_feed(Flowrate = 22, molefraction = {1, 0, 0, 0}, pressure = 100000, specified_stream = true, step_value = 2, stepchange = true, stepchangetime = 0.5, temperature = 573.15) annotation(
      Placement(visible = true, transformation(origin = {-84, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    unitoperations.MaterialStream hydrogenFeed(Flowrate = 66, Tdf(displayUnit = "K", start = 160), molefraction = {0, 1, 0, 0}, pressure = 2.5e+06, specified_stream = true) annotation(
      Placement(visible = true, transformation(origin = {-76, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    unitoperations.CSTR reactor(Ab = 0, Af = 5.1e11, Dynamic = false, Eab = 0, Eaf = 230e3, T_iso(displayUnit = "K") = 700, V_Total = 1, operation_mode = unitoperations.types.operation_type.Isothermal, order_b = {0, 0, 0, 0}, order_f = {1, 0.5, 0, 0}) annotation(
      Placement(visible = true, transformation(origin = {-13.4704, 29.4704}, extent = {{-24.5296, -29.4355}, {24.5296, 24.5296}}, rotation = 0)));
    unitoperations.valve valve1(OutletPfixed = true) annotation(
      Placement(visible = true, transformation(origin = {44, 12}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    unitoperations.MaterialStream reactorOutlet(Tdf(displayUnit = "K", start = 337)) annotation(
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
  
  end CSTR_Steady;



  model CSTR_Dynamic
    unitoperations.MaterialStream toluene_feed(Flowrate = 22, molefraction = {1, 0, 0, 0}, pressure = 100000, specified_stream = true, step_value = 2, stepchange = true, stepchangetime = 0.5, temperature = 573.15) annotation(
      Placement(visible = true, transformation(origin = {-84, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    unitoperations.MaterialStream hydrogenFeed(Flowrate = 66, Tdf(displayUnit = "K", start = 160), molefraction = {0, 1, 0, 0}, pressure = 2.5e+06, specified_stream = true) annotation(
      Placement(visible = true, transformation(origin = {-76, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    unitoperations.CSTR reactor(Ab = 0, Af = 5.1e11, Dynamic = true, Eab = 0, Eaf = 230e3, T_iso(displayUnit = "K") = 700, V_Total = 1, operation_mode = unitoperations.types.operation_type.Isothermal, order_b = {0, 0, 0, 0}, order_f = {1, 0.5, 0, 0}) annotation(
      Placement(visible = true, transformation(origin = {-13.4704, 29.4704}, extent = {{-24.5296, -29.4355}, {24.5296, 24.5296}}, rotation = 0)));
    unitoperations.valve valve1(OutletPfixed = true) annotation(
      Placement(visible = true, transformation(origin = {44, 12}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    unitoperations.MaterialStream reactorOutlet(Tdf(displayUnit = "K", start = 337)) annotation(
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
  
  end CSTR_Dynamic;


end NewTest;
