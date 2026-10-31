"""Write the curated practicality table of the shortlist (data/curated/) and its README, hand-transcribed from sourced
research notes (each row keeps its source URL and access date). Usage: run this script with python from any directory;
the output paths are relative to the repository root."""
import csv
from pathlib import Path

D = "2026-10-25"
MTI = "https://www.mtixtl.com/products/"
SHALOM = "https://www.shalomeo.com/Shalomeo-lasercrystalsandcomponents-alexandritelasercrystal-product"
LUX = "https://luxiumsolutions.com/radiation-detection-scintillators/crystal-scintillators/bgo-bismuth-germanate"
GF = "https://www.goodfellow.com/uk/bismuth-disc-single-crystal-1000036376"
SAM = "https://www.samaterials.com/crystal-substrates/2430-yttria-stabilized-zirconia-ysz-crystal-substrates.html"
rows = [
    ("Al2O3 R-3c", "melting_point", "2040", "C", "https://www.crystran.com/optical-materials/sapphire-al2o3/", "MTI states the same"),
    ("Al2O3 R-3c", "growth_method", "Verneuil, Czochralski, Kyropoulos, ribbon growth", "", "https://www.crystran.com/optical-materials/sapphire-al2o3/", ""),
    ("Al2O3 R-3c", "largest_stated_crystal", "350 kg Kyropoulos boule", "", "https://semiconductor-today.com/news_items/2017/mar/monocrystal_130317.shtml", "trade-press report of a company claim (Monocrystal, 2017)"),
    ("Al2O3 R-3c", "water", "solubility 98e-6 g per 100 g water at 20 C", "", "https://www.crystran.com/optical-materials/sapphire-al2o3/", "no explicit hygroscopicity statement"),
    ("Al2O3 R-3c", "price", "24.99 USD for 10x10x0.5 mm", "10x10x0.5", MTI + "alc101005s1", "C-plane wafer ALC101005S1"),
    ("BeAl2O4 Pnma", "melting_point", "1870", "C", SHALOM, "alexandrite (Cr:BeAl2O4) data sheet"),
    ("BeAl2O4 Pnma", "growth_method", "improved Czochralski", "", SHALOM, "alexandrite"),
    ("BeAl2O4 Pnma", "largest_stated_crystal", "laser rods up to 9.53 mm dia x 117 mm", "", SHALOM, "no commercial boule size found; literature Cz boules about 110 x 25 mm (2008)"),
    ("BeAl2O4 Pnma", "water", "not found", "", "", ""),
    ("BeAl2O4 Pnma", "price", "not found", "", "", "quote-only vendors"),
    ("Bi R-3m", "melting_point", "271.3", "C", GF, ""),
    ("Bi R-3m", "growth_method", "not found", "", "", "Bridgman per the growable list (MaTecK)"),
    ("Bi R-3m", "largest_stated_crystal", "12 mm dia x 1.5 mm disc, (111)", "", GF, "only catalogue size found"),
    ("Bi R-3m", "water", "insoluble in water; surface oxidises in air", "", "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/5359367/JSON?heading=Solubility", "CRC Handbook via HSDB"),
    ("Bi R-3m", "price", "1970 GBP for 12 mm dia x 1.5 mm", "disc:12x1.5", GF, "BI002031, single piece"),
    ("Ge(Bi3O5)4 I23", "melting_point", "923", "C", "https://pmc.ncbi.nlm.nih.gov/articles/PMC12199007/", "1196 K reported melting temperature (ACS Omega 10, 25219, 2025)"),
    ("Ge(Bi3O5)4 I23", "growth_method", "modified Czochralski", "", "https://4lasers.com/components/bgo-crystals", ""),
    ("Ge(Bi3O5)4 I23", "largest_stated_crystal", "up to 3 in diameter aperture", "", "https://4lasers.com/components/bgo-crystals", ""),
    ("Ge(Bi3O5)4 I23", "water", "not found", "", "", ""),
    ("Ge(Bi3O5)4 I23", "price", "1495 USD for 2 in dia x 1 mm", "disc:50.8x1", MTI + "bgo20b50d10c1us", "BGO20b50D10C1US"),
    ("Ge3(BiO3)4 I-43d", "melting_point", "1050", "C", LUX, "1323 K"),
    ("Ge3(BiO3)4 I-43d", "growth_method", "low-thermal-gradient Czochralski", "", "https://www.niic.nsc.ru/science/razrabotki/619-materials/bgo/2853-bgo-en", ""),
    ("Ge3(BiO3)4 I-43d", "largest_stated_crystal", "ingots up to 130 mm dia x 450 mm", "", "https://www.niic.nsc.ru/science/razrabotki/619-materials/bgo/2853-bgo-en", "Luxium: detectors 1-5 in dia, up to 5 in thick"),
    ("Ge3(BiO3)4 I-43d", "water", "not hygroscopic", "", LUX, ""),
    ("Ge3(BiO3)4 I-43d", "price", "245 USD for 10x10x0.5 mm", "10x10x0.5", MTI + "bgo12b101005s2us", "BGO12b101005S2US"),
    ("MgO Fm-3m", "melting_point", "2800", "C", "https://www.crystran.com/optical-materials/magnesium-oxide-mgo/", "other sources: 2825 C (CRC), 2852 C (MTI)"),
    ("MgO Fm-3m", "growth_method", "arc melting", "", MTI + "mga101005s1", ""),
    ("MgO Fm-3m", "largest_stated_crystal", "2 in dia x 0.5 mm wafer; grown in relatively small sizes", "", MTI + "mga50d05c1us", "Crystran: grown in relatively small sizes"),
    ("MgO Fm-3m", "water", "moisture-sensitive (store in a vacuum box)", "", MTI + "mga101005s1", "solubility 0.00062 g/100 g water (Crystran)"),
    ("MgO Fm-3m", "price", "39 USD for 10x10x0.5 mm", "10x10x0.5", MTI + "mga101005s1", "MGa101005S1"),
    ("Si Fd-3m", "melting_point", "1420", "C", "https://www.crystran.com/optical-materials/silicon-si/", ""),
    ("Si Fd-3m", "growth_method", "Czochralski or float zone", "", "https://www.crystran.com/optical-materials/silicon-si/", ""),
    ("Si Fd-3m", "largest_stated_crystal", "300 mm wafers from CZ or FZ ingots", "", "https://www.siltronic.com/en/products.html", ""),
    ("Si Fd-3m", "water", "insoluble in water", "", "https://www.crystran.com/optical-materials/silicon-si/", ""),
    ("Si Fd-3m", "price", "79 USD for 6 in dia x 0.675 mm", "disc:152.4x0.675", MTI + "si-wafer-111-6-dia-x-0-65-0-7-mm-1sp-p-type-b-doped-with-resistivities-1-10-ohm-cm", "B-doped (111) wafer"),
    ("SiC P6_3mc", "melting_point", "2700 (sublimes)", "C", "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/9863/JSON?heading=Melting+Point", "NIOSH 4892 F sublimes; HSDB 2600 C; DOE PAC 2830 C; polytype not specified"),
    ("SiC P6_3mc", "growth_method", "not found on a vendor page", "", "", "physical vapour transport is the industry route but no opened vendor page states it; MTI says MOCVD (doubtful)"),
    ("SiC P6_3mc", "largest_stated_crystal", "6H: 100 and 150 mm wafers (reseller); 4H: 200 mm (Wolfspeed)", "", "https://universitywafer.com/sic-wafer-supplier.html", "Wolfspeed: https://www.wolfspeed.com/products/materials/"),
    ("SiC P6_3mc", "water", "insoluble in water", "", "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/9863/JSON?heading=Solubility", "CRC via HSDB"),
    ("SiC P6_3mc", "price", "199 USD for 10x10x0.26 mm (6H)", "10x10x0.26", MTI + "sic-6h-0001-10x10x0-26-mm-two-sides-polished", "SC6HZ1010026S2US"),
    ("SiO2 P3_221", "melting_point", "1710", "C", "https://www.crystran.com/optical-materials/crystal-quartz-sio2/", "MTI: 1610 C; alpha-beta transition 573.1 C"),
    ("SiO2 P3_221", "growth_method", "hydrothermal", "", MTI + "soz101005s1", ""),
    ("SiO2 P3_221", "largest_stated_crystal", "4 in dia wafers", "", MTI + "soz101d05c1", ""),
    ("SiO2 P3_221", "water", "insoluble in water", "", "https://www.crystran.com/optical-materials/crystal-quartz-sio2/", ""),
    ("SiO2 P3_221", "price", "25 USD for 10x10x0.5 mm", "10x10x0.5", MTI + "soz101005s1", "Z-cut SOZ101005S1"),
    ("Y3Al5O12 Ia-3d", "melting_point", "1940", "C", "https://www.crystran.com/optical-materials/yttrium-aluminium-garnet-yag/", "MTI: 1970 C"),
    ("Y3Al5O12 Ia-3d", "growth_method", "Czochralski", "", "https://www.crystran.com/optical-materials/yttrium-aluminium-garnet-yag/", ""),
    ("Y3Al5O12 Ia-3d", "largest_stated_crystal", "up to approximately 100 mm diameter", "", "https://www.crystran.com/optical-materials/yttrium-aluminium-garnet-yag/", ""),
    ("Y3Al5O12 Ia-3d", "water", "insoluble in water", "", "https://www.crystran.com/optical-materials/yttrium-aluminium-garnet-yag/", ""),
    ("Y3Al5O12 Ia-3d", "price", "214 USD for 10x10x0.5 mm", "10x10x0.5", MTI + "yaga101005s1", "YAGa101005S1"),
    ("ZrO2 P2_1/c", "melting_point", "2500 (YSZ, 8 mol% Y2O3)", "C", MTI + "ysza101005s1", "the grown crystal is cubic YSZ; pure ZrO2 2715 C (a-m.de FAQ, not peer reviewed)"),
    ("ZrO2 P2_1/c", "growth_method", "skull melting (YSZ)", "", SAM, "no monoclinic single-crystal offer found"),
    ("ZrO2 P2_1/c", "largest_stated_crystal", "2 in dia substrates (YSZ)", "", SAM, ""),
    ("ZrO2 P2_1/c", "water", "not found", "", "", ""),
    ("ZrO2 P2_1/c", "price", "49.95 USD for 10x10x0.5 mm (YSZ)", "10x10x0.5", MTI + "ysza101005s1", "YSZa101005S1"),
    # measurement search
    ("BeAl2O4 Pnma", "measurement_search", "none found in EXFOR (15 searches, validated on sapphire), OpenAlex (13 queries) or web searches (3); Google Scholar not searched", "", "neutron-filter-atlas-paper/literature/stage4/dossier_sources.md", "proposed measurement; the OSTI lead (Adib 2004) was opened: calculations for Si and sapphire only"),
    ("SiC P6_3mc", "measurement_search", "none found for any SiC polytype in EXFOR (15 searches, validated on sapphire), OpenAlex (13 queries) or web searches (3); Google Scholar not searched; one lead unopened", "", "neutron-filter-atlas-paper/literature/stage4/dossier_sources.md", "proposed measurement; unopened lead: Y. Zhu 2014 NC State thesis (repository HTTP 503; by its title most likely a thermal-scattering evaluation)"),
    # every single-piece USD catalogue price in the sourced notes, so the dossier can give a range
    ("Al2O3 R-3c", "price", "144 USD for 3 in dia x 0.5 mm", "disc:76.2x0.5", MTI + "al2o3-sapphire-wafer-c-plane-0001-3-dia-x0-5mm-2sp", "ALC76D05C2US"),
    ("Ge(Bi3O5)4 I23", "price", "925 USD for 10x10x1.0 mm", "10x10x1.0", MTI + "bgo20b101010s1us", "BGO20b101010S1US"),
    ("MgO Fm-3m", "price", "599 USD for 2 in dia x 0.5 mm", "disc:50.8x0.5", MTI + "mga50d05c1us", "MGa50D05C1US"),
    ("Si Fd-3m", "price", "155 USD for 4 in dia x 0.5 mm", "disc:101.6x0.5", "https://www.mtixtl.com/collections/crystals/products.json", "undoped FZ (111) SIUc101D05C1R10000"),
    ("SiC P6_3mc", "price", "1510 USD for 2 in dia x 0.33 mm (6H)", "disc:50.8x0.33", "https://www.pmoptics.com/images/silicon_carbide_semi-insulating_wafers.html", "PMOptics PWSC-31B23211, semi-insulating"),
    ("SiC P6_3mc", "price", "189 USD for 8 mm dia x 0.3 mm (6H)", "disc:8x0.3", MTI + "sic-6h-0001-8mm-dia-x-0-3-mm-th-one-side-polished", "SC6HZ08D03C1"),
    ("SiO2 P3_221", "price", "175 USD for 4 in dia x 0.5 mm", "disc:101.6x0.5", MTI + "soz101d05c1", "Z-cut SOZ101D05C1"),
    ("ZrO2 P2_1/c", "price", "179 USD for 25.4x25.4x0.5 mm (YSZ)", "25.4x25.4x0.5", MTI + "ysz-100-1x1x0-5mm-1sp-ysza252505s1", "YSZa252505S1"),
]
repo = Path(__file__).resolve().parents[1]
out = repo / "data" / "curated" / "stage4_dossier.csv"
with open(out, "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh, lineterminator="\n")
    w.writerow(["compound", "quantity", "value", "size_or_unit", "source_url", "access_date", "note"])
    for r in rows:
        w.writerow([r[0], r[1], r[2], r[3], r[4], D if r[4] else "", r[5]])
readme = """# stage4_dossier.csv

Practicality facts for the Tier-A shortlist (PREREG_STAGE4 §6), written by `scripts/make_stage4_dossier_csv.py`,
hand-transcribed from the sourced research notes in
`neutron-filter-atlas-paper/literature/stage4/dossier_sources.md` (verbatim quotes, URLs, access date 2026-10-25),
one primary value per quantity, except price, where every single-piece USD catalogue item in the notes is listed
(the dossier reports the range per cm3; referee m6). Disagreeing sources are named in `note`. Three rows were re-opened and checked by
hand (Shalom EO alexandrite page, Crystran sapphire page, Luxium BGO page). `size_or_unit` gives the product size of a
price in mm (`disc:<diameter>x<thickness>`), used for the price per cm3 in `stage4_dossier.json`. These are thin
substrates, so that figure says little about the price of a filter-sized block.
"""
(repo / "data" / "curated" / "stage4_dossier_README.md").write_text(readme, encoding="utf-8", newline="\n")
print(len(rows), "rows ->", out)
