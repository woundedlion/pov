"""Manufacturer selections for the rev 1.3 differential interface."""


def part(manufacturer, mpn, datasheet):
    return {"Manufacturer": manufacturer, "MPN": mpn, "Datasheet": datasheet}


PULSE_RESISTOR_DATASHEET = "https://www.vishay.com/docs/20043/crcwhpe3.pdf"
PARTS = {
    "C_IN": part("Nichicon", "UPW1H101MPD",
                 "https://www.nichicon.com/en-us/part/upw1h101mpd/8471/"),
    "U_SYNC": part("Texas Instruments", "THVD2410DR",
                   "https://www.ti.com/lit/ds/symlink/thvd2410.pdf"),
    "D_SYNC": part("Bourns", "CDSOT23-SM712",
                   "https://www.bourns.com/docs/product-datasheets/cdsot23-sm712.pdf"),
    "R_A": part("Vishay", "CRCW060310R0FKEAHP", PULSE_RESISTOR_DATASHEET),
    "R_B": part("Vishay", "CRCW060310R0FKEAHP", PULSE_RESISTOR_DATASHEET),
    "R_TERM": part("Vishay", "CRCW0805120RFKEAHP", PULSE_RESISTOR_DATASHEET),
    "C_BULK3": part("YAGEO", "CC0603KRX7R8BB105",
                    "https://www.yageogroup.com/download/specsheet/CC0603KRX7R8BB105"),
    **{ref: part("YAGEO", "CC0603KRX7R9BB104",
                 "https://www.yageogroup.com/download/specsheet/CC0603KRX7R9BB104")
       for ref in ("C_DEC1", "C_DEC2", "C_DEC3")},
    **{ref: part("YAGEO", "RC0603FR-0710KL",
                 "https://www.yageogroup.com/component-documentation/download/specsheet/RC0603FR-0710KL")
       for ref in ("R_TX", "R_MEN", "R_DATA_PD", "R_CLK_PD")},
}
