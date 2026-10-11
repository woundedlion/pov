"""Routed rev 1.1 placement and footprint fixtures."""

QUILTER_FIXED = {
    "J1": (8.5, 18.8, 0),
    "J4": (8.5, 28.5, 90),
    "J2": (49.8, 2.77, 0),
    "J3A": (49.8, 11.7, 0),
    "J3B": (49.8, 20.6, 0),
    "JP_ID0": (56.0, 8.5, 90),
    "JP_ID1": (56.0, 11.9, 90),
    "JP_ID2": (56.0, 15.3, 90),
    "JP_SHLD": (56.0, 19.0, 90),
    "U_MCU": (29.0, 11.7, 0),
    "C_IN": (20.5, 27.0, 0),
    "C_DEC1": (12.5, 1.35, 0),
    "U1": (36.6, 26.5, 180),
    "C_DEC2": (31.3, 30.0, 180),
    "R_D1": (42.2, 27.77, 180),
    "R_D2": (42.2, 23.96, 180),
    "R_S": (31.7, 24.0, 90),
    "R1": (29.7, 26.2, 90),
    "R2": (29.7, 22.8, 90),
    "C_SYNC": (27.8, 22.8, 90),
    "D_BUS": (46.0, 26.5, 90),
    "R_PD": (46.0, 23.0, 90),
}

QUILTER_FIXED_FOOTPRINTS = {
    **{ref: "Connector_PinHeader_2.54mm:PinHeader_1x03_P2.54mm_Vertical"
       for ref in ("J2", "J3A", "J3B")},
    "J1": "Connector_PinHeader_2.54mm:PinHeader_1x02_P2.54mm_Vertical",
    "R1": "Resistor_SMD:R_0603_1608Metric",
    "R2": "Resistor_SMD:R_0603_1608Metric",
    "R_PD": "Resistor_SMD:R_0603_1608Metric",
    "R_S": "Resistor_SMD:R_0805_2012Metric",
}
