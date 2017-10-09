package Functions
  function Density
    input Real p[6], TC, T, Pressure;
    output Real density "units kmol/m3";
  protected
    Real Tr;
  protected
    parameter Real R = 8.314 "gas constant";
  algorithm
    Tr := T / TC;
    if T < TC then
      if p[1] == 105 then
        density := p[2] / p[3] ^ (1 + (1 - T / p[4]) ^ p[5]);
      elseif p[1] == 106 then
        density := p[2] * (1 - Tr) ^ (p[3] + p[4] * Tr + p[5] * Tr ^ 2 + p[6] * Tr ^ 3);
      end if;
    else
      density := Pressure / (R * T * 1000);
    end if;
  end Density;

  function VapCpId
    /*Calculates Vapor Specific Heat*/
    input Real VapCp[6] "from chemsep database";
    input Real T(unit = "K") "Temperature";
    output Real Cp(unit = "J/mol.K") "specific heat";
  algorithm
    Cp := (VapCp[2] + exp(VapCp[3] / T + VapCp[4] + VapCp[5] * T + VapCp[6] * T ^ 2)) / 1000;
  end VapCpId;

  function HV
    input Real HOV[6];
    input Real Tc;
    input Real T;
    output Real Hv;
  protected
    Real Tr = T / Tc;
  algorithm
    Hv := HOV[2] * (1 - Tr) ^ (HOV[3] + HOV[4] * Tr + HOV[5] * Tr ^ 2 + HOV[6] * Tr ^ 3) / 1000;
  end HV;

  function Psat
    input Real VP[6];
    input Real T;
    output Real Vp;
  algorithm
    Vp := exp(VP[2] + VP[3] / T + VP[4] * log(T) + VP[5] .* T .^ VP[6]);
  end Psat;

  function LiqCpId
    /*Calculates specific heat of liquid at given Temperature*/
    input Real LiqCp[6] "from chemsep database";
    input Real T(unit = "K") "Temperature";
    output Real Cp(unit = "J/mol") "Specific heat of liquid";
  algorithm
    Cp := (LiqCp[2] + exp(LiqCp[3] / T + LiqCp[4] + LiqCp[5] * T + LiqCp[6] * T ^ 2)) / 1000;
  end LiqCpId;

  function HVapId
    /* Calculates enthalpy of ideal vapor */
    input Real VapCp[6] "from chemsep database";
    input Real HOV[6] "from chemsep database";
    input Real Tc "critical temp, from chemsep database";
    input Real T(unit = "K") "Temperature";
    output Real Ent(unit = "J/mol") "Molar Enthalpy";
  protected
    Integer n = 20;
    Real Cp[n - 1];
  algorithm
    for i in 1:n - 1 loop
      Cp[i] := VapCpId(VapCp, 298.15 + i * (T - 298.15) / n);
    end for;
    if T > 298.15 then
      Ent := (T - 298.15) * (VapCpId(VapCp, T) / 2 + sum(Cp[:]) + VapCpId(VapCp, 298.15) / 2) / n;
    else
      Ent := -(T - 298.15) * (VapCpId(VapCp, T) / 2 + sum(Cp[:]) + VapCpId(VapCp, 298.15) / 2) / n;
    end if;
  end HVapId;

function HVapId1
  /* Calculates enthalpy of ideal vapor */
  input Real VapCp[6] "from chemsep database";
  input Real HOV[6] "from chemsep database";
  input Real Tc "critical temp, from chemsep database";
  input Real T(unit = "K") "Temperature";
  output Real Ent1(unit = "J/mol") "Molar Enthalpy";
protected
  Real Temp, Ent;// = 298.15;
algorithm
  Ent := 0.001;
  Temp := 298.15;
  if T > 298.15 then
    while Temp < T loop
      Ent := Ent + VapCpId(VapCp, Temp) * 1;
      Temp := Temp + 1;
    end while;
    
    Ent1 := Ent;
  else
    while Temp > T loop
      Ent := Ent + VapCpId(VapCp, Temp) * 1;
      Temp := Temp - 1;
    end while;
    Ent1 :=  -Ent;
  end if;
end HVapId1;

  function HLiqId1
    /* Calculates Enthalpy of Ideal Liquid*/
    //    input Real IGHF(unit = "J/kmol") "from chemsep database std. Heat of formation";
    input Real VapCp[6] "from chemsep database";
    input Real HOV[6] "from chemsep database";
    input Real Tc "critical temp, from chemsep database";
    input Real T(unit = "K") "Temperature";
    output Real Ent(unit = "J/mol") "Molar Enthalpy";
  algorithm
    Ent := 0.001;
    if T < Tc then
      Ent := Functions.HVapId1(VapCp, HOV, Tc, T) - HV(HOV, Tc, T);
    else
      Ent := Functions.HVapId1(VapCp, HOV, Tc, T);
    end if;
  end HLiqId1;

  function HLiqId
    /* Calculates Enthalpy of Ideal Liquid*/
    //    input Real IGHF(unit = "J/kmol") "from chemsep database std. Heat of formation";
    input Real VapCp[6] "from chemsep database";
    input Real HOV[6] "from chemsep database";
    input Real Tc "critical temp, from chemsep database";
    input Real T(unit = "K") "Temperature";
    output Real Ent(unit = "J/mol") "Molar Enthalpy";
  algorithm
    Ent := 0.001;
    if T < Tc then
      Ent := Functions.HVapId(VapCp, HOV, Tc, T) - HV(HOV, Tc, T);
    else
      Ent := Functions.HVapId(VapCp, HOV, Tc, T);
    end if;
  end HLiqId;
end Functions;
