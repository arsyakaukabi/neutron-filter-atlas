# growable_crystals.csv: Tier A "commercially growable" single crystals

Built on 2026-10-20. The file has 374 evidence rows covering 166 distinct compounds. Elements count as compounds. The two polymorphs of BaB2O4 (alpha and beta) count separately, as do diamond and graphite.

## How the list was built

1. Vendor catalogue pages were found by browsing vendor home pages, with a few web searches. Each page was opened with WebFetch on 2026-10-20 and read for: the material listed, whether it is called a single crystal, the largest size listed, the growth method and the price.
2. Each row records one piece of evidence for one compound from one vendor. `catalogue_url` is always a page that was actually opened and loaded in this session. Where the evidence came from a list or category page, that page is cited, not a product page that was not opened.
3. `max_size_stated`, `growth_method_stated` and `price_stated` are filled only with what the page said, paraphrased closely. They are empty when the page said nothing. Nothing was estimated.
   - MTI: `max_size_stated` is the largest item visible on the collection page. The price is one example item, with its size or part number in brackets. Prices are in USD as listed on 2026-10-20.
   - MaTecK metals: sizes are the "standard size" printed on each element page. Custom sizes are offered on request.
4. `formula` is the host formula. Dopants (Ce, Tl, Nd, Yb, Cr, Eu and so on) go in `common_name`/`notes`. Solid solutions keep a composite formula, for example `(LaAlO3)0.3(Sr2AlTaO6)0.7` for LSAT, `(Lu,Y)2SiO5` for LYSO, `Tl(Br,I)` for KRS-5, and `Cd1-xZnxTe` and `SrxBa1-xNb2O6`. YSZ is written as `ZrO2`, with the stabilizer noted.
5. `crystal_system_or_phase` comes from standard crystallography, not necessarily from the vendor page. Crystran pages do state the structure, and it matched.
6. Neutron-absorber flags were added automatically from tabulated natural-element thermal (2200 m/s) absorption cross sections (Sears 1992):
   - "STRONG" means sigma_a >= 100 b.
   - "moderate" means 10 b <= sigma_a < 100 b.
   - These flags are for filtering later. Compounds were not excluded because of them.
   - Hydrogenous crystals (KDP) are noted separately because of incoherent scattering.

## Vendors used (14)

| Vendor | Pages used |
|---|---|
| Crystran (UK) | optical-materials index and per-material pages |
| MTI Corporation (USA) | `mtixtl.com/collections/...` crystal-substrate and metal single-crystal collections |
| MaTecK (Germany) | metal, oxide, semiconductor and optical single-crystal pages, plus per-element pages |
| Princeton Scientific (USA) | metal, oxide and semiconductor single-crystal catalogues |
| CASTECH (China) | NLO, laser and birefringent crystal pages |
| Shalom EO (China) | scintillator and laser/NLO crystal pages |
| Epic Crystal (China) | scintillation crystal pages |
| Luxium Solutions (formerly Saint-Gobain Crystals) | crystal scintillator pages |
| Crytur (Czech Republic) | materials pages |
| SurfaceNet (Germany) | CaWO4 page and crystal data page |
| Korth Kristalle (Germany) | home-page material list ("Eigenzucht" = own-grown) and SrF2 page |
| Wolfspeed | SiC substrates |
| HexaTech | AlN substrates |
| Isomet | PbMoO4 data sheet |

Tried but not cited:
- Crystal GmbH: the product page returned 404, and the substrate page names no materials.
- MolTech / mt-berlin.com: connection refused.
- korth.de detail pages other than SrF2: the fetch was blocked.
- Hellma Materials: redirect only; no product page reached.
- MTI collections with no single-crystal products: LaF3, ScMgAlO4, Co, Pb, V, Pd, Zr, W, Mo, Ti and Fe. They were empty or polycrystalline only.
- MTI AlN: ceramic or template on sapphire only, so excluded.

Not attempted: Alineason, Hilger, EKSMA, Rubicon, Monocrystal, Kyocera, Coherent/II-VI, Del Mar, PI-KEM, Biotain, Gamdan and Dayoptics.

## Seed compounds

All seed compounds except Be have at least one row. Weak cases:

- **AlN**: confirmed only by HexaTech ("single crystalline AlN substrates", c- and m-plane). No size or price is stated on the page. MTI's AlN is ceramic or a template on sapphire, so it was excluded.
- **Diamond (C)**: confirmed only as small single crystals. Crystran says "economically up to about 4 mm"; MTI sells 3x3 mm CVD single crystal. Both fail the ~10 mm criterion. Large diamond optics are polycrystalline CVD.
- **GaN**: only MTI, with free-standing wafers up to 50.8 mm diameter but only 0.25 mm thick.
- **Bi and Pb**: only metal single crystals of about 10-12 mm diameter and 1-2 mm thickness (MaTecK standard sizes). Princeton Scientific lists both with no size stated.
- **Spinel MgAl2O4**: Crystran's "spinel" is a polycrystalline ceramic, so it was excluded. Single-crystal spinel is confirmed by MTI (2 in), MaTecK and Princeton.
- **ZnSe/ZnS**: most optical stock (Crystran, MaTecK "Cleartran"/CVD, MTI windows) is polycrystalline CVD and was excluded. Single-crystal ZnSe/ZnS is confirmed by MaTecK (Bridgman/Markov, 30 mm diameter), Princeton (CVT/PVT) and MTI (oriented 10x10 mm).

Not confirmed:
- **Be metal**: not listed by MaTecK or Princeton Scientific. No catalogue page listing Be single crystals was opened.
- **CaMoO4**: no vendor page found.
- **PbF2**: Crystran states "We cannot supply this material".

These three are not in the CSV.

## Caveats

- **HOPG is not a true single crystal.** It is a c-axis-oriented mosaic polycrystal, randomly oriented in-plane. No vendor page stated a mosaic spread. Its rows (MTI, Princeton) are flagged in `crystal_system_or_phase`.
- **Many rows are below the ~10 mm criterion.** They are kept, with the stated size, so they can be filtered:
  - Metal single crystals: Zr 3-4 mm; In 7 mm; Re, Ru, Hf about 6-8 mm. Princeton says Ti and Zr are "limited to a few millimeters".
  - Diamond, YIG (5 mm) and the MTI Mg substrate (5 mm).
- **Size limits.** Metal crystals at MaTecK are typically 6-20 mm discs, 1-2 mm thick (Cu up to 60 mm). MTI substrate sizes (often 10x10x0.5 mm) are catalogue wafer formats. They are not boule limits; boules are often far larger. The listed sizes show only that the material is grown, not the largest crystal possible.
- **"Single crystal" is often implied, not stated.**
  - Many MTI collection pages and Crystran pages never use the words "single crystal". Single-crystal status is implied by oriented cuts, cleavage or optic-axis statements, or growth from a seed.
  - The MaTecK optical-crystals page and the Korth home page also do not state single crystal per item.
  - Each such case is noted in the row's `notes`.
- **Natural rather than grown.** Calcite (CaCO3) from Crystran is natural mined material; vendors still sell it as single-crystal pieces.
- **Composition not stated.** SBN, PMN-PT, CZT, KRS-6 and DKDP deuteration have no composition on the pages.
- **Absorber-bearing compounds are included on purpose.** Examples: Gd garnets/scandates/GdVO4/KGW, Li/B compounds (LiF, LBO, BBO, CLBO, LiNbO3, LiYF4), Cd compounds, In compounds, SmScO3, DyScO3, Hf/Ir/Rh/Re metals and CLLB.
- **Multiple rows per vendor.** A compound can appear more than once for one vendor when it was found on different pages, for example MaTecK per-material and list pages, or CASTECH's beta- and alpha-BBO.
- **Intermetallic and alloy single crystals were not included.** MaTecK lists many, such as NiAl, Ni3Al, FeAl, Cu3Au, NiTi, LaB6, TiC and TaC, all offered "grown to spec". They could be added from `https://mateck.com/material/metal-single-crystals` if needed.
