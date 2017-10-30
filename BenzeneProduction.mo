package BenzeneProduction
  model SteadyState
    unitoperations.MaterialStream HydrogenFeed(Flowrate = 66, Tdf(start = 423.15), molefraction = {0, 0.9, 0, 0.1}, pressure = 2.5e+06, specified_stream = true, stepchange = false, temperature(displayUnit = "K") = 500) annotation(
      Placement(visible = true, transformation(origin = {-182, 34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    unitoperations.MaterialStream TolueneFeed(Flowrate = 22, molefraction = {0.99, 0, 0.01, 0}, pressure = 2.5e+06, specified_stream = true, step_value = 2, stepchange = true, stepchangetime = 0.5, temperature = 623.15) annotation(
      Placement(visible = true, transformation(origin = {-183, 9}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    unitoperations.CSTR Reactor(Ab = 0, Af = 5e11, Dynamic = false, Eab = 0, Eaf = 230e3, T_iso(displayUnit = "K") = 700, V_Total = 1, operation_mode = unitoperations.types.operation_type.Isothermal, order_b = {0, 0, 0, 0}, order_f = {1, 0.5, 0, 0}) annotation(
      Placement(visible = true, transformation(origin = {-112.448, 8.66329}, extent = {{-19.5518, -23.4621}, {19.5518, 19.5518}}, rotation = 0)));
    unitoperations.valve valve1 annotation(
      Placement(visible = true, transformation(origin = {-73, -1}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    unitoperations.MaterialStream ReactorOutlet annotation(
      Placement(visible = true, transformation(origin = {-39, -1}, extent = {{-13, -13}, {13, 13}}, rotation = 0)));
    unitoperations.PTFlash FalshColumn(Dynamic = false) annotation(
      Placement(visible = true, transformation(origin = {5, 10.7335}, extent = {{-24, -37.0632}, {24, 26.4737}}, rotation = 0)));
    unitoperations.MaterialStream FlashVaporStream annotation(
      Placement(visible = true, transformation(origin = {84, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    unitoperations.MaterialStream FlashLiquidStream annotation(
      Placement(visible = true, transformation(origin = {76, -28}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    unitoperations.valve valve2(OutletPfixed = true) annotation(
      Placement(visible = true, transformation(origin = {55, 37}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    unitoperations.valve valve3(OutletPfixed = true) annotation(
      Placement(visible = true, transformation(origin = {48, -28}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    unitoperations.Distillation DistillationColumn(Dynamic = false, Override_Sizing_Calculations = false, specification1 = unitoperations.types.Distillation_spec1.RefluxRatio, specification1_value = 2, specification2 = unitoperations.types.Distillation_spec2.ProductMolarFlow, specification2_value = 15) annotation(
      Placement(visible = true, transformation(origin = {116, -26.0249}, extent = {{-23, -33.6349}, {23, 24.0249}}, rotation = 0)));
    unitoperations.MaterialStream Distillate annotation(
      Placement(visible = true, transformation(origin = {182, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    unitoperations.MaterialStream Bottoms annotation(
      Placement(visible = true, transformation(origin = {184, -58}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(ReactorOutlet.port2, FalshColumn.port3) annotation(
      Line(points = {{-28, -1}, {-18, -1}, {-18, 7}, {-12, 7}}));
    connect(FalshColumn.port2, valve2.port1) annotation(
      Line(points = {{24, 25}, {35, 25}, {35, 37}, {46, 37}}));
    connect(FalshColumn.port1, valve3.port1) annotation(
      Line(points = {{22, -10}, {34, -10}, {34, -28}, {40, -28}}));
    connect(DistillationColumn.port3, Bottoms.port1) annotation(
      Line(points = {{132, -44}, {160, -44}, {160, -58}, {176, -58}}));
    connect(DistillationColumn.port2, Distillate.port1) annotation(
      Line(points = {{132, -14}, {160, -14}, {160, 0}, {174, 0}, {174, 0}}));
    connect(FlashLiquidStream.port2, DistillationColumn.port1) annotation(
      Line(points = {{84.5, -28}, {98, -28}}));
    connect(valve3.port2, FlashLiquidStream.port1) annotation(
      Line(points = {{56, -28}, {68, -28}}));
    connect(valve2.port2, FlashVaporStream.port1) annotation(
      Line(points = {{64, 37}, {76, 37}, {76, 38}}));
    connect(valve1.port2, ReactorOutlet.port1) annotation(
      Line(points = {{-64, -1}, {-50, -1}}));
    connect(Reactor.port3, valve1.port1) annotation(
      Line(points = {{-95, -1}, {-82, -1}}));
    connect(HydrogenFeed.port2, Reactor.port2) annotation(
      Line(points = {{-173.5, 34}, {-158.25, 34}, {-158.25, 15}, {-129, 15}}));
    connect(TolueneFeed.port2, Reactor.port1) annotation(
      Line(points = {{-174, 9}, {-136.5, 9}, {-136.5, 8}, {-129, 8}}));
    annotation(
      Diagram(coordinateSystem(extent = {{-200, -100}, {200, 100}}, grid = {4, 4})),
      Icon(coordinateSystem(extent = {{-200, -100}, {200, 100}}, grid = {4, 4})),
      __OpenModelica_commandLineOptions = "",
  experiment(StartTime = 0, StopTime = 10000, Tolerance = 1e-06, Interval = 2));
  end SteadyState;

  model DynamicFlowSheet
    unitoperations.MaterialStream HydrogenFeed(Flowrate = 66, Tdf(start = 423.15), molefraction = {0, 0.9, 0, 0.1}, pressure = 2.5e+06, specified_stream = true, stepchange = false, temperature(displayUnit = "K") = 500) annotation(
      Placement(visible = true, transformation(origin = {-182, 34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    unitoperations.MaterialStream TolueneFeed(Flowrate = 22, molefraction = {0.99, 0, 0.01, 0}, pressure = 2.5e+06, specified_stream = true, step_value = 2, stepchange = true, stepchangetime = 0.5, temperature = 623.15) annotation(
      Placement(visible = true, transformation(origin = {-183, 9}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    unitoperations.CSTR Reactor(Ab = 0, Af = 5e11, Dynamic = true, Eab = 0, Eaf = 230e3, T_iso(displayUnit = "K") = 700, V_Total = 1, operation_mode = unitoperations.types.operation_type.Isothermal, order_b = {0, 0, 0, 0}, order_f = {1, 0.5, 0, 0}) annotation(
      Placement(visible = true, transformation(origin = {-112.448, 8.66329}, extent = {{-19.5518, -23.4621}, {19.5518, 19.5518}}, rotation = 0)));
    unitoperations.valve valve1 annotation(
      Placement(visible = true, transformation(origin = {-73, -1}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    unitoperations.MaterialStream ReactorOutlet annotation(
      Placement(visible = true, transformation(origin = {-39, -1}, extent = {{-13, -13}, {13, 13}}, rotation = 0)));
    unitoperations.PTFlash FlashColumn(Dynamic = true) annotation(
      Placement(visible = true, transformation(origin = {5, 10.7335}, extent = {{-24, -37.0632}, {24, 26.4737}}, rotation = 0)));
    unitoperations.MaterialStream FlashVaporStream annotation(
      Placement(visible = true, transformation(origin = {84, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    unitoperations.MaterialStream FlashLiquidStream annotation(
      Placement(visible = true, transformation(origin = {76, -28}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    unitoperations.valve valve2(OutletPfixed = true) annotation(
      Placement(visible = true, transformation(origin = {55, 37}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    unitoperations.valve valve3(OutletPfixed = true) annotation(
      Placement(visible = true, transformation(origin = {48, -28}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    unitoperations.Distillation DistillationColumn(Dynamic = true, Override_Sizing_Calculations = false, specification1 = unitoperations.types.Distillation_spec1.RefluxRatio, specification1_value = 2, specification2 = unitoperations.types.Distillation_spec2.ProductMolarFlow, specification2_value = 15) annotation(
      Placement(visible = true, transformation(origin = {116, -26.0249}, extent = {{-23, -33.6349}, {23, 24.0249}}, rotation = 0)));
    unitoperations.MaterialStream Distillate annotation(
      Placement(visible = true, transformation(origin = {182, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    unitoperations.MaterialStream Bottoms annotation(
      Placement(visible = true, transformation(origin = {184, -58}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(ReactorOutlet.port2, FlashColumn.port3) annotation(
      Line(points = {{-28, -1}, {-18, -1}, {-18, 7}, {-12, 7}}));
    connect(FlashColumn.port2, valve2.port1) annotation(
      Line(points = {{24, 25}, {35, 25}, {35, 37}, {46, 37}}));
    connect(FlashColumn.port1, valve3.port1) annotation(
      Line(points = {{22, -10}, {34, -10}, {34, -28}, {40, -28}}));
    connect(DistillationColumn.port3, Bottoms.port1) annotation(
      Line(points = {{132, -44}, {160, -44}, {160, -58}, {176, -58}}));
    connect(DistillationColumn.port2, Distillate.port1) annotation(
      Line(points = {{132, -14}, {160, -14}, {160, 0}, {174, 0}, {174, 0}}));
    connect(FlashLiquidStream.port2, DistillationColumn.port1) annotation(
      Line(points = {{84.5, -28}, {98, -28}}));
    connect(valve3.port2, FlashLiquidStream.port1) annotation(
      Line(points = {{56, -28}, {68, -28}}));
    connect(valve2.port2, FlashVaporStream.port1) annotation(
      Line(points = {{64, 37}, {76, 37}, {76, 38}}));
    connect(valve1.port2, ReactorOutlet.port1) annotation(
      Line(points = {{-64, -1}, {-50, -1}}));
    connect(Reactor.port3, valve1.port1) annotation(
      Line(points = {{-95, -1}, {-82, -1}}));
    connect(HydrogenFeed.port2, Reactor.port2) annotation(
      Line(points = {{-173.5, 34}, {-158.25, 34}, {-158.25, 15}, {-129, 15}}));
    connect(TolueneFeed.port2, Reactor.port1) annotation(
      Line(points = {{-174, 9}, {-136.5, 9}, {-136.5, 8}, {-129, 8}}));
    annotation(
      Diagram(coordinateSystem(extent = {{-200, -100}, {200, 100}}, grid = {4, 4})),
      Icon(coordinateSystem(extent = {{-200, -100}, {200, 100}}, grid = {4, 4})),
      __OpenModelica_commandLineOptions = "",
  experiment(StartTime = 0, StopTime = 10000, Tolerance = 1e-06, Interval = 2));
  end DynamicFlowSheet;

  model DynamicDistillationColumn "Dynamics are considered only for distillation column"
    unitoperations.MaterialStream HydrogenFeed(Flowrate = 66, Tdf(start = 423.15), molefraction = {0, 0.9, 0, 0.1}, pressure = 2.5e+06, specified_stream = true, stepchange = false, temperature(displayUnit = "K") = 500) annotation(
      Placement(visible = true, transformation(origin = {-182, 34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    unitoperations.MaterialStream TolueneFeed(Flowrate = 22, molefraction = {0.99, 0, 0.01, 0}, pressure = 2.5e+06, specified_stream = true, step_value = 2, stepchange = true, stepchangetime = 0.5, temperature = 623.15) annotation(
      Placement(visible = true, transformation(origin = {-183, 9}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    unitoperations.CSTR Reactor(Ab = 0, Af = 5e11, Dynamic = false, Eab = 0, Eaf = 230e3, T_iso(displayUnit = "K") = 700, V_Total = 1, operation_mode = unitoperations.types.operation_type.Isothermal, order_b = {0, 0, 0, 0}, order_f = {1, 0.5, 0, 0}) annotation(
      Placement(visible = true, transformation(origin = {-112.448, 8.66329}, extent = {{-19.5518, -23.4621}, {19.5518, 19.5518}}, rotation = 0)));
    unitoperations.valve valve1 annotation(
      Placement(visible = true, transformation(origin = {-73, -1}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    unitoperations.MaterialStream ReactorOutlet annotation(
      Placement(visible = true, transformation(origin = {-39, -1}, extent = {{-13, -13}, {13, 13}}, rotation = 0)));
    unitoperations.PTFlash FlashColumn(Dynamic = false) annotation(
      Placement(visible = true, transformation(origin = {5, 10.7335}, extent = {{-24, -37.0632}, {24, 26.4737}}, rotation = 0)));
    unitoperations.MaterialStream FlashVaporStream annotation(
      Placement(visible = true, transformation(origin = {84, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    unitoperations.MaterialStream FlashLiquidStream annotation(
      Placement(visible = true, transformation(origin = {76, -28}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    unitoperations.valve valve2(OutletPfixed = true) annotation(
      Placement(visible = true, transformation(origin = {55, 37}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    unitoperations.valve valve3(OutletPfixed = true) annotation(
      Placement(visible = true, transformation(origin = {48, -28}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    unitoperations.Distillation DistillationColumn(Dynamic = true, Override_Sizing_Calculations = false, specification1 = unitoperations.types.Distillation_spec1.RefluxRatio, specification1_value = 2, specification2 = unitoperations.types.Distillation_spec2.ProductMolarFlow, specification2_value = 15) annotation(
      Placement(visible = true, transformation(origin = {116, -26.0249}, extent = {{-23, -33.6349}, {23, 24.0249}}, rotation = 0)));
    unitoperations.MaterialStream Distillate annotation(
      Placement(visible = true, transformation(origin = {182, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    unitoperations.MaterialStream Bottoms annotation(
      Placement(visible = true, transformation(origin = {184, -58}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(ReactorOutlet.port2, FlashColumn.port3) annotation(
      Line(points = {{-28, -1}, {-18, -1}, {-18, 7}, {-12, 7}}));
    connect(FlashColumn.port2, valve2.port1) annotation(
      Line(points = {{24, 25}, {35, 25}, {35, 37}, {46, 37}}));
    connect(FlashColumn.port1, valve3.port1) annotation(
      Line(points = {{22, -10}, {34, -10}, {34, -28}, {40, -28}}));
    connect(DistillationColumn.port3, Bottoms.port1) annotation(
      Line(points = {{132, -44}, {160, -44}, {160, -58}, {176, -58}}));
    connect(DistillationColumn.port2, Distillate.port1) annotation(
      Line(points = {{132, -14}, {160, -14}, {160, 0}, {174, 0}, {174, 0}}));
    connect(FlashLiquidStream.port2, DistillationColumn.port1) annotation(
      Line(points = {{84.5, -28}, {98, -28}}));
    connect(valve3.port2, FlashLiquidStream.port1) annotation(
      Line(points = {{56, -28}, {68, -28}}));
    connect(valve2.port2, FlashVaporStream.port1) annotation(
      Line(points = {{64, 37}, {76, 37}, {76, 38}}));
    connect(valve1.port2, ReactorOutlet.port1) annotation(
      Line(points = {{-64, -1}, {-50, -1}}));
    connect(Reactor.port3, valve1.port1) annotation(
      Line(points = {{-95, -1}, {-82, -1}}));
    connect(HydrogenFeed.port2, Reactor.port2) annotation(
      Line(points = {{-173.5, 34}, {-158.25, 34}, {-158.25, 15}, {-129, 15}}));
    connect(TolueneFeed.port2, Reactor.port1) annotation(
      Line(points = {{-174, 9}, {-136.5, 9}, {-136.5, 8}, {-129, 8}}));
    annotation(
      Diagram(coordinateSystem(extent = {{-200, -100}, {200, 100}}, grid = {4, 4})),
      Icon(coordinateSystem(extent = {{-200, -100}, {200, 100}}, grid = {4, 4})),
      __OpenModelica_commandLineOptions = "",
  experiment(StartTime = 0, StopTime = 10000, Tolerance = 1e-06, Interval = 2));
  end DynamicDistillationColumn;

  model DynamicFlashAndDistillation
    unitoperations.MaterialStream HydrogenFeed(Flowrate = 66, Tdf(start = 423.15), molefraction = {0, 0.9, 0, 0.1}, pressure = 2.5e+06, specified_stream = true, stepchange = false, temperature(displayUnit = "K") = 500) annotation(
      Placement(visible = true, transformation(origin = {-182, 34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    unitoperations.MaterialStream TolueneFeed(Flowrate = 22, molefraction = {0.99, 0, 0.01, 0}, pressure = 2.5e+06, specified_stream = true, step_value = 2, stepchange = true, stepchangetime = 0.5, temperature = 623.15) annotation(
      Placement(visible = true, transformation(origin = {-183, 9}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    unitoperations.CSTR Reactor(Ab = 0, Af = 5e11, Dynamic = false, Eab = 0, Eaf = 230e3, T_iso(displayUnit = "K") = 700, V_Total = 1, operation_mode = unitoperations.types.operation_type.Isothermal, order_b = {0, 0, 0, 0}, order_f = {1, 0.5, 0, 0}) annotation(
      Placement(visible = true, transformation(origin = {-112.448, 8.66329}, extent = {{-19.5518, -23.4621}, {19.5518, 19.5518}}, rotation = 0)));
    unitoperations.valve valve1 annotation(
      Placement(visible = true, transformation(origin = {-73, -1}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    unitoperations.MaterialStream ReactorOutlet annotation(
      Placement(visible = true, transformation(origin = {-39, -1}, extent = {{-13, -13}, {13, 13}}, rotation = 0)));
    unitoperations.PTFlash FlashColumn(Dynamic = true) annotation(
      Placement(visible = true, transformation(origin = {5, 10.7335}, extent = {{-24, -37.0632}, {24, 26.4737}}, rotation = 0)));
    unitoperations.MaterialStream FlashVaporStream annotation(
      Placement(visible = true, transformation(origin = {84, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    unitoperations.MaterialStream FlashLiquidStream annotation(
      Placement(visible = true, transformation(origin = {76, -28}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    unitoperations.valve valve2(OutletPfixed = true) annotation(
      Placement(visible = true, transformation(origin = {55, 37}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    unitoperations.valve valve3(OutletPfixed = true) annotation(
      Placement(visible = true, transformation(origin = {48, -28}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    unitoperations.Distillation DistillationColumn(Dynamic = true, Override_Sizing_Calculations = false, specification1 = unitoperations.types.Distillation_spec1.RefluxRatio, specification1_value = 2, specification2 = unitoperations.types.Distillation_spec2.ProductMolarFlow, specification2_value = 15) annotation(
      Placement(visible = true, transformation(origin = {116, -26.0249}, extent = {{-23, -33.6349}, {23, 24.0249}}, rotation = 0)));
    unitoperations.MaterialStream Distillate annotation(
      Placement(visible = true, transformation(origin = {182, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    unitoperations.MaterialStream Bottoms annotation(
      Placement(visible = true, transformation(origin = {184, -58}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(ReactorOutlet.port2, FlashColumn.port3) annotation(
      Line(points = {{-28, -1}, {-18, -1}, {-18, 7}, {-12, 7}}));
    connect(FlashColumn.port2, valve2.port1) annotation(
      Line(points = {{24, 25}, {35, 25}, {35, 37}, {46, 37}}));
    connect(FlashColumn.port1, valve3.port1) annotation(
      Line(points = {{22, -10}, {34, -10}, {34, -28}, {40, -28}}));
    connect(DistillationColumn.port3, Bottoms.port1) annotation(
      Line(points = {{132, -44}, {160, -44}, {160, -58}, {176, -58}}));
    connect(DistillationColumn.port2, Distillate.port1) annotation(
      Line(points = {{132, -14}, {160, -14}, {160, 0}, {174, 0}, {174, 0}}));
    connect(FlashLiquidStream.port2, DistillationColumn.port1) annotation(
      Line(points = {{84.5, -28}, {98, -28}}));
    connect(valve3.port2, FlashLiquidStream.port1) annotation(
      Line(points = {{56, -28}, {68, -28}}));
    connect(valve2.port2, FlashVaporStream.port1) annotation(
      Line(points = {{64, 37}, {76, 37}, {76, 38}}));
    connect(valve1.port2, ReactorOutlet.port1) annotation(
      Line(points = {{-64, -1}, {-50, -1}}));
    connect(Reactor.port3, valve1.port1) annotation(
      Line(points = {{-95, -1}, {-82, -1}}));
    connect(HydrogenFeed.port2, Reactor.port2) annotation(
      Line(points = {{-173.5, 34}, {-158.25, 34}, {-158.25, 15}, {-129, 15}}));
    connect(TolueneFeed.port2, Reactor.port1) annotation(
      Line(points = {{-174, 9}, {-136.5, 9}, {-136.5, 8}, {-129, 8}}));
    annotation(
      Diagram(coordinateSystem(extent = {{-200, -100}, {200, 100}}, grid = {4, 4})),
      Icon(coordinateSystem(extent = {{-200, -100}, {200, 100}}, grid = {4, 4})),
      __OpenModelica_commandLineOptions = "",
  experiment(StartTime = 0, StopTime = 10000, Tolerance = 1e-06, Interval = 2));
  end DynamicFlashAndDistillation;

  model CSTR_Steady_test
    unitoperations.MaterialStream toluene_feed(Flowrate = 22, molefraction = {1, 0, 0, 0}, pressure = 100000, specified_stream = true, step_value = 2, stepchange = true, stepchangetime = 0.5, temperature = 573.15) annotation(
      Placement(visible = true, transformation(origin = {-84, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    unitoperations.MaterialStream hydrogenFeed(Flowrate = 66, Tdf(displayUnit = "K", start = 160), molefraction = {0, 1, 0, 0}, pressure = 2.5e+06, specified_stream = true) annotation(
      Placement(visible = true, transformation(origin = {-76, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    unitoperations.CSTR reactor(Ab = 0, Af = 5.1e11, ConstantTemperature = false, Dynamic = false, Eab = 0, Eaf = 230e3, T_iso(displayUnit = "K") = 700, V_Total = 1, operation_mode = unitoperations.types.operation_type.Isothermal, order_b = {0, 0, 0, 0}, order_f = {1, 0.5, 0, 0}) annotation(
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
    annotation(
      experiment(StartTime = 0, StopTime = 100, Tolerance = 1e-6, Interval = 0.2));
  end CSTR_Steady_test;

  model CSTR_Dynamic_test
    unitoperations.MaterialStream toluene_feed(Flowrate = 22, molefraction = {1, 0, 0, 0}, pressure = 100000, specified_stream = true, step_value = 2, stepchange = true, stepchangetime = 0.5, temperature = 573.15) annotation(
      Placement(visible = true, transformation(origin = {-84, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    unitoperations.MaterialStream hydrogenFeed(Flowrate = 66, Tdf(displayUnit = "K", start = 160), molefraction = {0, 1, 0, 0}, pressure = 2.5e+06, specified_stream = true) annotation(
      Placement(visible = true, transformation(origin = {-76, 42}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    unitoperations.CSTR reactor(Ab = 0, Af = 5.1e11, ConstantTemperature = false, Dynamic = true, Eab = 0, Eaf = 230e3, T_iso(displayUnit = "K") = 700, V_Total = 1, operation_mode = unitoperations.types.operation_type.Isothermal, order_b = {0, 0, 0, 0}, order_f = {1, 0.5, 0, 0}) annotation(
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
    annotation(
      experiment(StartTime = 0, StopTime = 100, Tolerance = 1e-06, Interval = 0.2));
  end CSTR_Dynamic_test;

  model Flash_Steady_test
    unitoperations.MaterialStream Feed(Flowrate = 88, molefraction = {0.25, 0.25, 0.25, 0.25}, pressure = 500000, specified_stream = true, step_value = 2, stepchange = true, temperature = 300) annotation(
      Placement(visible = true, transformation(origin = {-82, 12}, extent = {{-14, -14}, {14, 14}}, rotation = 0)));
    unitoperations.PTFlash FlashColumn(Dynamic = false, InputIsSpecifiedStream = true, L(start = 22), OverrideSizeCalculations = false, Pset = 500000, V(start = 66), hset = 1) annotation(
      Placement(visible = true, transformation(origin = {-20, 18.1631}, extent = {{-30, -44.4964}, {30, 31.7832}}, rotation = 0)));
    unitoperations.valve valve1(OutletPfixed = true, OutletPressure = 100000) annotation(
      Placement(visible = true, transformation(origin = {44, 32}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    unitoperations.MaterialStream VaporStream(Tdf(start = 277.15)) annotation(
      Placement(visible = true, transformation(origin = {80, 32}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    unitoperations.valve valve2(OutletPfixed = true) annotation(
      Placement(visible = true, transformation(origin = {34, -18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    unitoperations.MaterialStream LiquidStream(Tdf(start = 371.15)) annotation(
      Placement(visible = true, transformation(origin = {88, -18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(Feed.port2, FlashColumn.port3) annotation(
      Line(points = {{-70, 12}, {-57, 12}, {-57, 13}, {-41, 13}}));
    connect(FlashColumn.port2, valve1.port1) annotation(
      Line(points = {{3, 35}, {36, 35}, {36, 32}}));
    connect(FlashColumn.port1, valve2.port1) annotation(
      Line(points = {{2, -7}, {26, -7}, {26, -18}}));
    connect(valve2.port2, LiquidStream.port1) annotation(
      Line(points = {{42, -18}, {78, -18}, {78, -18}, {80, -18}}));
    connect(valve1.port2, VaporStream.port1) annotation(
      Line(points = {{52, 32}, {70, 32}, {70, 32}, {72, 32}}));
    annotation(
      experiment(StartTime = 0, StopTime = 10000, Tolerance = 1e-6, Interval = 2));
  end Flash_Steady_test;


  model Flash_Dynamic_test
    unitoperations.MaterialStream Feed(Flowrate = 88, molefraction = {0.25, 0.25, 0.25, 0.25}, pressure = 500000, specified_stream = true, step_value = 2, stepchange = true, temperature = 300) annotation(
      Placement(visible = true, transformation(origin = {-80, 8}, extent = {{-14, -14}, {14, 14}}, rotation = 0)));
    unitoperations.PTFlash FlashColumn(Dynamic = true, InputIsSpecifiedStream = true, L(start = 22), OverrideSizeCalculations = false, Pset = 500000, V(start = 66), hset = 1) annotation(
      Placement(visible = true, transformation(origin = {-22, 10.1631}, extent = {{-30, -44.4964}, {30, 31.7832}}, rotation = 0)));
    unitoperations.valve valve1(OutletPfixed = true, OutletPressure = 100000) annotation(
      Placement(visible = true, transformation(origin = {44, 32}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    unitoperations.MaterialStream VaporStream(Tdf(start = 277.15)) annotation(
      Placement(visible = true, transformation(origin = {80, 32}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    unitoperations.valve valve2(OutletPfixed = true) annotation(
      Placement(visible = true, transformation(origin = {34, -16}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    unitoperations.MaterialStream LiquidStream(Tdf(start = 371.15)) annotation(
      Placement(visible = true, transformation(origin = {88, -18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(FlashColumn.port1, valve2.port1) annotation(
      Line(points = {{0, -16}, {26, -16}}));
    connect(valve2.port2, LiquidStream.port1) annotation(
      Line(points = {{42, -16}, {70, -16}, {70, -18}, {80, -18}}));
    connect(FlashColumn.port2, valve1.port1) annotation(
      Line(points = {{0, 28}, {36, 28}, {36, 32}, {36, 32}}));
    connect(valve1.port2, VaporStream.port1) annotation(
      Line(points = {{52, 32}, {70, 32}, {70, 32}, {72, 32}}));
    connect(Feed.port2, FlashColumn.port3) annotation(
      Line(points = {{-68, 8}, {-46, 8}}));
    annotation(
      experiment(StartTime = 0, StopTime = 10000, Tolerance = 1e-06, Interval = 2));
  end Flash_Dynamic_test;


  model Distillation_Steady_test
    unitoperations.MaterialStream feed(Flowrate = 80, Tdf(start = 371.15), molefraction = {0.5, 0, 0.5, 0}, pressure = 100000, specified_stream = true, step_value = 2, stepchange = true, stepchangetime = 0.5) annotation(
      Placement(visible = true, transformation(origin = {-78, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    unitoperations.Distillation distillation1(Dynamic = false, M_eff = 1, Override_Sizing_Calculations = false, P_condenser = 100000, Pressure_drop = 700, specification1 = unitoperations.types.Distillation_spec1.RefluxRatio, specification1_value = 2, specification2 = unitoperations.types.Distillation_spec2.ProductMolarFlow, specification2_value = 30) annotation(
      Placement(visible = true, transformation(origin = {-33, 0.787066}, extent = {{-16, -29.6981}, {16, 21.2129}}, rotation = 0)));
    unitoperations.MaterialStream distillate annotation(
      Placement(visible = true, transformation(origin = {20, 18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    unitoperations.MaterialStream bottoms annotation(
      Placement(visible = true, transformation(origin = {20, -28}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(distillation1.port3, bottoms.port1) annotation(
      Line(points = {{-22, -14}, {-4, -14}, {-4, -28}, {12, -28}, {12, -28}}));
    connect(distillation1.port2, distillate.port1) annotation(
      Line(points = {{-22, 12}, {10, 12}, {10, 18}, {12, 18}}));
    connect(feed.port2, distillation1.port1) annotation(
      Line(points = {{-70, 0}, {-46, 0}, {-46, -2}, {-46, -2}}));
    annotation(
      experiment(StartTime = 0, StopTime = 1000, Tolerance = 1e-06, Interval = 2));
  end Distillation_Steady_test;

  model DistillationDynamic_test
    unitoperations.MaterialStream feed(Flowrate = 80, Tdf(start = 371.15), molefraction = {0.5, 0, 0.5, 0}, pressure = 100000, specified_stream = true, step_value = 2, stepchange = true, stepchangetime = 0.5) annotation(
      Placement(visible = true, transformation(origin = {-78, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    unitoperations.Distillation distillation1(Dynamic = true, M_eff = 1, Override_Sizing_Calculations = false, P_condenser = 100000, Pressure_drop = 700, specification1 = unitoperations.types.Distillation_spec1.RefluxRatio, specification1_value = 2, specification2 = unitoperations.types.Distillation_spec2.ProductMolarFlow, specification2_value = 30) annotation(
      Placement(visible = true, transformation(origin = {-33, 0.787066}, extent = {{-16, -29.6981}, {16, 21.2129}}, rotation = 0)));
    unitoperations.MaterialStream distillate annotation(
      Placement(visible = true, transformation(origin = {20, 18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    unitoperations.MaterialStream bottoms annotation(
      Placement(visible = true, transformation(origin = {20, -28}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(distillation1.port3, bottoms.port1) annotation(
      Line(points = {{-22, -14}, {-4, -14}, {-4, -28}, {12, -28}, {12, -28}}));
    connect(distillation1.port2, distillate.port1) annotation(
      Line(points = {{-22, 12}, {10, 12}, {10, 18}, {12, 18}}));
    connect(feed.port2, distillation1.port1) annotation(
      Line(points = {{-70, 0}, {-46, 0}, {-46, -2}, {-46, -2}}));
    annotation(
      experiment(StartTime = 0, StopTime = 1000, Tolerance = 1e-06, Interval = 2));
  end DistillationDynamic_test;
end BenzeneProduction;
