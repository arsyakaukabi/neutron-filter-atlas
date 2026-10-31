// Builds growable_crystals.csv from hand-entered evidence rows (all from vendor pages opened on 2026-10-20).
const fs = require('fs');
const OUT = process.argv[2];
const DATE = '2026-10-20';

// Thermal (2200 m/s) absorption cross sections, natural abundance, barns (Sears 1992). Used only to flag absorbers.
const SA = {H:0.3326,D:0.000519,Li:70.5,Be:0.0076,B:767,C:0.0035,N:1.9,O:0.00019,F:0.0096,Na:0.53,Mg:0.063,Al:0.231,Si:0.171,P:0.172,S:0.53,Cl:33.5,K:2.1,Ca:0.43,Sc:27.5,Ti:6.09,V:5.08,Cr:3.05,Mn:13.3,Fe:2.56,Co:37.18,Ni:4.49,Cu:3.78,Zn:1.11,Ga:2.75,Ge:2.2,As:4.5,Se:11.7,Br:6.9,Rb:0.38,Sr:1.28,Y:1.28,Zr:0.185,Nb:1.15,Mo:2.48,Ru:2.56,Rh:144.8,Pd:6.9,Ag:63.3,Cd:2520,In:193.8,Sn:0.626,Sb:4.91,Te:4.7,I:6.15,Cs:29,Ba:1.1,La:8.97,Ce:0.63,Nd:50.5,Sm:5922,Eu:4530,Gd:49700,Tb:23.4,Dy:994,Ho:64.7,Er:159,Tm:100,Yb:34.8,Lu:74,Hf:104.1,Ta:20.6,W:18.3,Re:89.7,Ir:425,Pt:10.3,Au:98.65,Tl:3.43,Pb:0.171,Bi:0.0338};

function absorberNote(formula){
  const els = [...new Set((formula.match(/[A-Z][a-z]?/g)||[]))];
  const strong = [], mod = [];
  for (const e of els){
    const s = SA[e];
    if (s === undefined) continue;
    if (s >= 100) strong.push(`${e}(${s} b)`); else if (s >= 10) mod.push(`${e}(${s} b)`);
  }
  const parts = [];
  if (strong.length) parts.push('STRONG neutron absorber: ' + strong.join(', '));
  if (mod.length) parts.push('moderate neutron absorber: ' + mod.join(', '));
  return parts.join('; ');
}

const CR = 'https://www.crystran.com/optical-materials/';
const MTI = 'https://mtixtl.com/collections/';
const MK = 'https://mateck.com/material/';
const PS = 'https://princetonscientific.com/materials/';
const CT = 'https://www.castech.com/product/';
const SH = 'https://www.shalomeo.com/';
const EP = 'https://www.epic-crystal.com/scintillation-crystals/';
const LX = 'https://luxiumsolutions.com/radiation-detection-scintillators/crystal-scintillators';
const CY = 'https://www.crytur.com/materials/';

const OS = 'oriented substrate listed under MTI Crystal Substrate A-Z; listing page does not use the words "single crystal" (single-crystal nature implied by crystallographic orientation)';
const OSsc = 'listing explicitly says "single crystal"';

// [formula, common_name, phase, vendor, url, max_size, method, price, notes]
const R = [
// ---------------- Crystran ----------------
['BaF2','barium fluoride','cubic (fluorite, Fm-3m)','Crystran',CR+'barium-fluoride-baf2/','','vacuum Stockbarger','','page does not use the words "single crystal"; cleaves on (111)'],
['CaF2','calcium fluoride','cubic (fluorite, Fm-3m)','Crystran',CR+'calcium-fluoride-caf2/','diameters up to about 250 mm','vacuum Stockbarger','',''],
['LiF','lithium fluoride','cubic (rock salt)','Crystran',CR+'lithium-fluoride-lif/','ingots approximately 100 mm diameter','vacuum Stockbarger','',''],
['MgF2','magnesium fluoride','tetragonal (rutile type)','Crystran',CR+'magnesium-fluoride-mgf2/','ingots of various diameters (no max stated)','vacuum Stockbarger','','birefringent; optic-axis discussion implies single crystal'],
['SrF2','strontium fluoride','cubic (fluorite, Fm-3m)','Crystran',CR+'strontium-fluoride-srf2/','','vacuum Stockbarger','',''],
['Al2O3','sapphire','trigonal (corundum)','Crystran',CR+'sapphire-al2o3/','','Verneuil and Czochralski (standard grade); Kyropoulos (higher quality); ribbon growth (thin sheets)','',''],
['Ge','germanium','cubic (diamond)','Crystran',CR+'germanium-ge/','','Czochralski','',''],
['Si','silicon','cubic (diamond)','Crystran',CR+'silicon-si/','','Czochralski (CZ)','',''],
['SiO2','alpha-quartz (crystal quartz)','trigonal (alpha-quartz, P3121/P3221)','Crystran',CR+'crystal-quartz-sio2/','','','','page: synthetic quartz "produced ... in large, long-faceted crystals"; fused silica (amorphous) is a separate, excluded product'],
['MgO','magnesium oxide (periclase)','cubic (rock salt)','Crystran',CR+'magnesium-oxide-mgo/','','','','page: "grown in relatively small sizes, mainly in Japan and China" and "quite difficult to source"'],
['TiO2','rutile','tetragonal (rutile)','Crystran',CR+'rutile-tio2/','','','','listed as optical material with data sheet; page states no size/method'],
['Y3Al5O12','YAG (yttrium aluminium garnet)','cubic (garnet)','Crystran',CR+'yttrium-aluminium-garnet-yag/','up to approximately 100 mm diameter','Czochralski','',''],
['CsI','caesium iodide','cubic (CsCl type)','Crystran',CR+'caesium-iodide-csi/','ingots of approximately 70 mm diameter','sealed ampoule Stockbarger','',''],
['CsBr','caesium bromide','cubic (CsCl type)','Crystran',CR+'caesium-bromide-csbr/','','sealed ampoule Stockbarger','','not in Crystran stock-optics menu; on enquiry'],
['KBr','potassium bromide','cubic (rock salt)','Crystran',CR+'potassium-bromide-kbr/','"large ingots" (no dimension stated)','Kyropoulos','',''],
['KCl','potassium chloride','cubic (rock salt)','Crystran',CR+'potassium-chloride-kcl/','"large ingots" (no dimension stated)','Kyropoulos','',''],
['NaCl','sodium chloride (rock salt)','cubic (rock salt)','Crystran',CR+'sodium-chloride-nacl/','"large ingots" (no dimension stated)','Kyropoulos','',''],
['KI','potassium iodide','cubic (rock salt)','Crystran',CR+'potassium-iodide-ki/','"large ingots" (no dimension stated)','Kyropoulos','','page: only useful in controlled lab conditions (very soft, very water soluble)'],
['NaF','sodium fluoride','cubic (rock salt)','Crystran',CR+'sodium-fluoride-naf/','useful size of pieces limited to about 80 mm','vacuum Stockbarger','',''],
['LaF3','lanthanum fluoride','trigonal/hexagonal (tysonite)','Crystran',CR+'lanthanum-fluoride-laf3eu3/','small ingots of about 10 mm diameter','','','small; borderline for 10 mm criterion'],
['CaCO3','calcite','trigonal (calcite, R-3c)','Crystran',CR+'calcite-caco3/','','','','page: calcite "is often mined naturally"; Crystran stocks small rhombs (natural, not grown)'],
['GaAs','gallium arsenide','cubic (zinc blende)','Crystran',CR+'gallium-arsenide-gaas/','','Czochralski or horizontal Bridgman','',''],
['AgCl','silver chloride','cubic (rock salt)','Crystran',CR+'silver-chloride-agcl/','"small ingots"','sealed-ampoule Stockbarger','','soft, cold-flows'],
['AgBr','silver bromide','cubic (rock salt)','Crystran',CR+'silver-bromide-agbr/','"small ingots"','sealed ampoule Stockbarger','','not in stock menus; data sheet only'],
['Tl(Br,I)','KRS-5 (thallium bromo-iodide)','cubic (CsCl type)','Crystran',CR+'krs5-thallium-bromo-iodide/','','sealed-ampoule Stockbarger','','composition 42 mol% TlBr / 58 mol% TlI (mixed crystal)'],
['RbCl','rubidium chloride','cubic (rock salt)','Crystran',CR+'rubidium-chloride-rbcl/','','sealed-ampoule Stockbarger','','page: "only specialist applications"'],
['C','diamond','cubic (diamond)','Crystran',CR+'diamond-c-cubic-carbon/','single crystal "economically up to about 4 mm diameter" (CVD polycrystalline to 75 mm)','','','single-crystal diamond is far below the 10 mm criterion; large pieces are polycrystalline CVD'],
// ---------------- MTI Corporation ----------------
['Al2O3','sapphire','trigonal (corundum)','MTI Corporation',MTI+'al2o3','4 in dia x 0.5 mm wafer (ALC100D05C1)','','$12.75 (5x5x0.5 mm C-plane wafer ALC050505S1)','listed as "Single Crystal Substrate" / "Sapphire Wafer"'],
['BaF2','barium fluoride','cubic (fluorite)','MTI Corporation',MTI+'baf2','1 in dia x 1.0 mm (BFc25D10C1)','','$235.00 (BFc25D10C1)',OS],
['BaTiO3','barium titanate','tetragonal perovskite (RT)','MTI Corporation',MTI+'batio3','10x10x1.0 mm','','$1,049.00 (BTOb101010S1US)','"substrate grade (with domains)"; ferroelectric domains'],
['Bi4Ge3O12','BGO (bismuth germanate, eulytite)','cubic (eulytite)','MTI Corporation',MTI+'bi4ge3o12','10x10x0.5 mm','','$245.00 (BGO12b101005S2US)',OS],
['Bi12GeO20','BGO20 (sillenite bismuth germanate)','cubic (sillenite)','MTI Corporation',MTI+'bi12geo20','2 in dia x 1 mm (BGO20b50D10C1US)','','$475.00 (10x5x1.0 mm)',OS],
['CaF2','calcium fluoride','cubic (fluorite)','MTI Corporation',MTI+'caf2','10x10x1.0 mm','','$125.00 (CFa101010S1)',OS],
['CaCO3','calcite','trigonal (calcite)','MTI Corporation',MTI+'caco3','10x10x1.0 mm (CCOa101010S1)','','$79.00 (5x5x0.5 mm)','listed as "Single Crystal CaCO3"'],
['CdWO4','cadmium tungstate','monoclinic (wolframite)','MTI Corporation',MTI+'cdwo4','10x10x0.5 mm','','$269.00 (CdWO4b101005S2)',OS],
['CdS','cadmium sulfide','hexagonal (wurtzite)','MTI Corporation',MTI+'cds','25.4 mm dia x 1.0 mm','','$499.00 (CdSc101005S1, 10x10x0.5 mm)',OS],
['CdSe','cadmium selenide','hexagonal (wurtzite)','MTI Corporation',MTI+'cdse','10x10x1.0 mm','','$1,499.00',OSsc],
['CdTe','cadmium telluride','cubic (zinc blende)','MTI Corporation',MTI+'cdte-czt','10x10x1.0 mm','','$415.00 (5x5x1.0 mm)',OS],
['Cd1-xZnxTe','CZT (cadmium zinc telluride)','cubic (zinc blende) solid solution','MTI Corporation',MTI+'cdte-czt','10x10x1.0 mm','','$425.00 (5x5x1.0 mm)','alloy (solid solution), Zn fraction not stated; '+OS],
['CsI','caesium iodide (incl. CsI:Tl)','cubic (CsCl type)','MTI Corporation',MTI+'csitl','10x10x5.0 mm block (CSI101050SN)','','$249.00 (CSI101050SN)','listed as "Single Crystal Block"/"Single Crystal Substrate"; Tl-doped and undoped'],
['DyScO3','dysprosium scandate','orthorhombic perovskite','MTI Corporation',MTI+'dysco3','10x10x0.5 mm','','$259.00 (5x5x0.5 mm)',OS],
['C','diamond','cubic (diamond)','MTI Corporation',MTI+'diamond','3x3x0.25 mm','CVD','$429.00','"CVD Single Crystal Diamond"; far below 10 mm criterion'],
['Ga2O3','beta-gallium oxide','monoclinic (beta-Ga2O3)','MTI Corporation',MTI+'ga2o3','10x10x0.6 mm','','$1,089.00 (5x5x0.6 mm)',OSsc],
['GaAs','gallium arsenide','cubic (zinc blende)','MTI Corporation',MTI+'gaas','4 in dia wafer','VGF (stated for several products)','$129.00 (5x5x0.5 mm SI)',OS],
['GaN','gallium nitride','hexagonal (wurtzite)','MTI Corporation',MTI+'gan','50.8 mm dia x 0.25 mm','','$1,095.00 (10x10.5x0.27 mm)','"GaN-Single Crystal Substrate" (GaN-on-sapphire templates on same page excluded); thin wafers only'],
['GaP','gallium phosphide','cubic (zinc blende)','MTI Corporation',MTI+'gap','3 in dia x 0.45 mm','','$169.00 (10x10x0.5 mm)',OS],
['GaSb','gallium antimonide','cubic (zinc blende)','MTI Corporation',MTI+'gasb','2 in dia','','$135.00 (5x5x0.5 mm)',OS],
['Ge','germanium','cubic (diamond)','MTI Corporation',MTI+'ge','4 in dia','','$769.00 (3 in dia undoped wafer)','"Ge Single Crystal wafer"; some items prefixed "VGF-Ge"'],
['Gd3Ga5O12','GGG (gadolinium gallium garnet)','cubic (garnet)','MTI Corporation',MTI+'ggg','3 in dia x 0.5 mm','','$189.00 (5x5x0.5 mm)',OS],
['Gd2.6Ca0.4Ga4.1Mg0.25Zr0.65O12','SGGG (substituted GGG)','cubic (garnet)','MTI Corporation',MTI+'ggg','10x10x0.5 mm','','',OS],
['GdScO3','gadolinium scandate','orthorhombic perovskite','MTI Corporation',MTI+'gdsco3','10x10x0.5 mm','','$335.00 (5x5x0.5 mm)',OS],
['C','graphite (HOPG / pyrolytic graphite)','hexagonal graphite, mosaic (not a true single crystal)','MTI Corporation',MTI+'pyrolytic-graphite','','','$49.95 (pyrolytic graphite, c-axis textured, 10x(4-5)x0.5 mm)','HOPG SPI-2 grade listed; mosaic spread not stated; HOPG is a highly oriented polycrystal (c-axis aligned, random in-plane)'],
['InAs','indium arsenide','cubic (zinc blende)','MTI Corporation',MTI+'inas','2 in dia x 0.5 mm','','$1,150.00 (2 in wafer)',OS],
['InP','indium phosphide','cubic (zinc blende)','MTI Corporation',MTI+'inp-1','3 in dia x 0.625 mm','VGF','$129.00 (5x5x0.35 mm)',OS],
['InSb','indium antimonide','cubic (zinc blende)','MTI Corporation',MTI+'insb-1','2 in dia x 0.5 mm','','$279.00 (10x10x0.3 mm)',OS],
['KCl','potassium chloride','cubic (rock salt)','MTI Corporation',MTI+'kcl','10x10x(2.2-2.5) mm','','$119.00',OS],
['KH2PO4','KDP (potassium dihydrogen phosphate)','tetragonal (KDP)','MTI Corporation',MTI+'kh2po4','10x10x1.0 mm','','$125.00 (10x5x1.0 mm)','hydrogenous (incoherent H scattering); '+OS],
['KTaO3','potassium tantalate','cubic perovskite','MTI Corporation',MTI+'ktao3','10x10 mm','','$629.00 (10x10x0.5 mm)',OSsc],
['LaAlO3','lanthanum aluminate','rhombohedral perovskite','MTI Corporation',MTI+'laalo3','3 in dia x 0.5 mm','','$25.00 (5x5x0.5 mm)','twinned perovskite; '+OS],
['LiAlO2','gamma-lithium aluminate','tetragonal (gamma-LiAlO2)','MTI Corporation',MTI+'lialo2','2 in dia x 0.5 mm','','$75.00 (5x5x0.5 mm)',OS],
['LiGaO2','lithium gallate','orthorhombic','MTI Corporation',MTI+'ligao2','10x10x0.5 mm','','$379.00 (10x10x0.5 mm)','"with domain"; '+OS],
['La3Ga5SiO14','langasite (LGS)','trigonal (langasite)','MTI Corporation',MTI+'la3ga5sio14','3 in dia x 0.5 mm','','$155.00 (5x5x0.5 mm)',OS],
['(LaAlO3)0.3(Sr2AlTaO6)0.7','LSAT','cubic perovskite (solid solution)','MTI Corporation',MTI+'lsat','2 in dia','','$49.00 (5x5x0.5 mm)',OSsc],
['LiNbO3','lithium niobate','trigonal (LiNbO3, R3c)','MTI Corporation',MTI+'linbo3','4 in dia wafer','','$39.00 (10x10x0.5 mm Z-cut)',OS],
['LiTaO3','lithium tantalate','trigonal (LiNbO3 type)','MTI Corporation',MTI+'litao3','4 in dia x 0.35 mm','','$69.00 (10x10x0.5 mm)',OS],
['Lu2SiO5','LSO:Ce (lutetium oxyorthosilicate)','monoclinic','MTI Corporation',MTI+'lu2sio5-ce','10x10x0.5 mm','','$525.00','Ce-doped; '+OS],
['LiF','lithium fluoride','cubic (rock salt)','MTI Corporation',MTI+'lif','10x10x1.0 mm','','$129.00',OS],
['MgAl2O4','spinel','cubic (spinel)','MTI Corporation',MTI+'mgal2o4','2 in dia x 0.5 mm','','$89.00 (10x10x0.5 mm)',OSsc],
['MgF2','magnesium fluoride','tetragonal (rutile type)','MTI Corporation',MTI+'mgf2','10x10x1.0 mm','','$119.00',OS],
['MgO','magnesium oxide','cubic (rock salt)','MTI Corporation',MTI+'mgo','12.7x12.7x0.5 mm (visible items); page also shows a "2 inch dia" size filter','','$19.95 (5x5x0.5 mm)',OS],
['NaCl','sodium chloride','cubic (rock salt)','MTI Corporation',MTI+'nacl','10x10x2.5 mm','','$199.00 (10x10x2.0 mm (111))',OS],
['NdAlO3','neodymium aluminate','rhombohedral perovskite','MTI Corporation',MTI+'ndalo3','10x10x0.5 mm','','$359.00 (5x5x0.5 mm)',OS],
['NdCaAlO4','neodymium calcium aluminate','tetragonal (K2NiF4 type)','MTI Corporation',MTI+'ndcaalo4','10x10x0.5 mm','','$185.00',OS],
['NdGaO3','neodymium gallate','orthorhombic perovskite','MTI Corporation',MTI+'ndgao3','2 in dia x 0.5 mm','','$145.00 (5x5x0.5 mm)',OS],
['Nd3Ga5O12','NGG (neodymium gallium garnet)','cubic (garnet)','MTI Corporation',MTI+'ngg-nd3ga5o12','10x10x0.5 mm','','$279.00 (5x5x0.5 mm)',OS],
['PbWO4','PWO (lead tungstate)','tetragonal (scheelite)','MTI Corporation',MTI+'pbwo4','10x10x0.5 mm','','$102.00 (10x3x0.45 mm random orientation)',OSsc],
['PbS','lead sulfide (galena)','cubic (rock salt)','MTI Corporation',MTI+'pbs','10x10x0.5 mm','','$799.00 (5x5x0.5 mm)',OSsc],
['PbSe','lead selenide','cubic (rock salt)','MTI Corporation',MTI+'pbse','10x10x0.5 mm','','$999.00 (5x5x0.5 mm)',OS],
['PbTe','lead telluride','cubic (rock salt)','MTI Corporation',MTI+'pbte-1','10x10x1.0 mm','','$2,695.00',OS],
['SrTiO3','strontium titanate','cubic perovskite','MTI Corporation',MTI+'srtio3','1 in dia x 0.5 mm (first of 4 listing pages)','','$125.00 (10x10x0.5 mm)',OS],
['SrLaAlO4','strontium lanthanum aluminate','tetragonal (K2NiF4 type)','MTI Corporation',MTI+'srlaalo4','10x20x1.0 mm','','$129.00 (5x5x0.5 mm)',OS],
['SrLaGaO4','strontium lanthanum gallate','tetragonal (K2NiF4 type)','MTI Corporation',MTI+'srlagao4','10x10x0.5 mm','','$525.00 (10x5x0.5 mm)',OSsc],
['SiO2','alpha-quartz','trigonal (alpha-quartz)','MTI Corporation',MTI+'sio2-quartz','100 mm dia x 0.5 mm (SOZ100D05C2)','','$25.00 (10x10x0.5 mm Z-cut)','"Quartz single crystal"'],
['SiC','silicon carbide (4H/6H)','hexagonal (4H, 6H polytypes)','MTI Corporation',MTI+'sic-4h-6h-3c','1 in dia x 0.26 mm','','$475.00 (1 in 4H wafer)',OS],
['Si','silicon','cubic (diamond)','MTI Corporation',MTI+'si','6 in dia wafer','FZ (stated for some products)','$79.00 (6 in (111) B-doped wafer)',OS],
['SmScO3','samarium scandate','orthorhombic perovskite','MTI Corporation',MTI+'smsco3','10x10x0.5 mm','','$799.00',OS],
['TbScO3','terbium scandate','orthorhombic perovskite','MTI Corporation',MTI+'tbsco3','10x10x0.5 mm','','$345.00 (5x5x0.5 mm)',OS],
['TiO2','rutile','tetragonal (rutile)','MTI Corporation',MTI+'tio2-rutile','10x10x1.0 mm (first listing page)','','$159.00 (10x5x0.5 mm)',OS],
['Tb3Ga5O12','TGG (terbium gallium garnet)','cubic (garnet)','MTI Corporation',MTI+'tgg','10x10x0.5 mm','','$299.00',OSsc],
['TeO2','paratellurite','tetragonal (paratellurite)','MTI Corporation',MTI+'teo2','10x10x0.5 mm','','$115.00',OS],
['YAlO3','YAP (yttrium aluminium perovskite)','orthorhombic perovskite','MTI Corporation',MTI+'yalo3','10x10x1.0 mm','','$309.00 (10x10x0.5 mm)',OS],
['Y3Al5O12','YAG','cubic (garnet)','MTI Corporation',MTI+'yag','10x10x1.0 mm','','$214.00 (10x10x0.5 mm (111))',OS],
['Y3Fe5O12','YIG (yttrium iron garnet)','cubic (garnet)','MTI Corporation',MTI+'yig','5 mm dia x 0.5 mm','','$1,019.00 (5 mm dia x 0.35 mm)','"YIG single crystal substrate"; below 10 mm criterion'],
['YVO4','yttrium orthovanadate','tetragonal (zircon)','MTI Corporation',MTI+'yvo4','10x10x0.5 mm','','$105.00 (5x5x0.5 mm)',OS],
['ZrO2','YSZ (yttria-stabilized zirconia)','cubic fluorite (Y2O3-stabilized)','MTI Corporation',MTI+'ysz','1 in x 1 in x 0.5 mm','','$25.95 (5x5x0.5 mm)','"YSZ Single Crystal Substrate"; contains Y2O3 stabilizer (exact mol% not stated)'],
['Y3Sc2Ga3O12','YSGG','cubic (garnet)','MTI Corporation',MTI+'ysgg-y3sc2ga3o12','10x10x0.5 mm','','$459.00',OS],
['Gd3Sc2Ga3O12','GSGG','cubic (garnet)','MTI Corporation',MTI+'gsgggd3sc2ga3o12','2 in dia x 0.5 mm','','$489.00 (10x10x0.5 mm)',OS],
['ZnO','zinc oxide','hexagonal (wurtzite)','MTI Corporation',MTI+'zno','10x10x1.0 mm','','$179.00 (5x5x0.5 mm)',OS],
['ZnS','zinc sulfide','cubic (zinc blende)','MTI Corporation',MTI+'zns-1','10x10x0.5 mm ((111) oriented)','','$1,495.00 (ZSc101005S1US)','(111)-oriented item implies single crystal; poly-crystal multispectral window on same page excluded'],
['ZnSe','zinc selenide','cubic (zinc blende)','MTI Corporation',MTI+'znse','10x10x1.0 mm ((110) oriented)','','','oriented (100)/(110)/(111) items imply single crystal; ZnSe poly-crystal windows on same page excluded; price of oriented item not captured'],
['ZnTe','zinc telluride','cubic (zinc blende)','MTI Corporation',MTI+'znte','10x10x0.5 mm','','$699.00',OSsc],
['Bi2Te3','bismuth telluride','trigonal (tetradymite)','MTI Corporation',MTI+'bi2te3','about 10x10x0.1 mm (irregular, as cleaved)','','$1,395.00','"Highly Oriented Crystal Substrate" (not stated as single crystal)'],
['SrxBa1-xNb2O6','SBN (strontium barium niobate)','tetragonal tungsten bronze','MTI Corporation',MTI+'sbn','10x10x0.5 mm','','$1,195.00','composition x not stated; '+OS],
['Fe2O3','hematite (alpha-Fe2O3)','trigonal (corundum)','MTI Corporation',MTI+'fe2o3','10x10x0.5 mm','','$1,998.00','"Alpha-Fe2O3 Crystal (0001)"'],
['Cu','copper','cubic fcc','MTI Corporation',MTI+'cu','2 in dia x 1.0 mm','','$99.00 (5x5x1.0 mm)',OSsc],
['Al','aluminium','cubic fcc','MTI Corporation',MTI+'aluminum','2 in dia x 0.95 mm','','$195.00 (10x10x0.5 mm)',OSsc],
['Ni','nickel','cubic fcc','MTI Corporation',MTI+'nickelsinglextl','30 mm dia x 1.0 mm','','$1,230.25 (30 mm dia <110>)',OSsc],
['Ag','silver','cubic fcc','MTI Corporation',MTI+'ag-silver','10x10x0.5 mm','','$695.00 (5x5x0.5 mm <111>)',OSsc],
['Au','gold','cubic fcc','MTI Corporation',MTI+'au','10x10x0.5 mm','','$1,496.00 (10x5x0.5 mm <110>)',OSsc],
['Mg','magnesium','hexagonal hcp','MTI Corporation',MTI+'mg-1','5x5x(0.9-1.0) mm','','$99.96','"Magnesium (Mg) Single Crystal Substrate"; below 10 mm criterion'],
['Zn','zinc','hexagonal hcp','MTI Corporation',MTI+'zn','10x10x0.5 mm','','RFQ (price on request)','"Zn Single Crystal Substrate (0001)"'],
// ---------------- MaTecK: metals (element pages) ----------------
['Bi','bismuth','trigonal (A7, rhombohedral)','MaTecK',MK+'elements/bismuth','standard size 12 mm dia x 1-2 mm (custom on request)','Bridgman','','purity 99.999%'],
['Pb','lead','cubic fcc','MaTecK',MK+'elements/lead','standard size 10-12 mm dia x 1-2 mm','Bridgman','','purity 99.99%'],
['Cu','copper','cubic fcc','MaTecK',MK+'elements/copper','standard size 10-60 mm dia x 1-2 mm','Czochralski','','purity 99.999%'],
['Nb','niobium','cubic bcc','MaTecK',MK+'elements/niobium','standard size 8-12 mm dia x 1-2 mm','floating zone','','purity 99.99%'],
['W','tungsten','cubic bcc','MaTecK',MK+'elements/tungsten','standard size 8-12 mm dia x 1-2 mm','floating zone','','purity 99.999%'],
['Al','aluminium','cubic fcc','MaTecK',MK+'elements/aluminum','standard size 12-20 mm dia x 1-2 mm','Bridgman or Czochralski','','purity 99.9999%'],
['Ni','nickel','cubic fcc','MaTecK',MK+'elements/nickel','standard size 10-20 mm dia x 1-2 mm','Czochralski','','purity 99.99%'],
['Mo','molybdenum','cubic bcc','MaTecK',MK+'elements/molybdenum','standard size 8-12 mm dia x 1-2 mm','floating zone','','purity 99.99%'],
['Ta','tantalum','cubic bcc','MaTecK',MK+'elements/tantalum','standard size 6-12 mm dia x 1-2 mm','floating zone','','purity 99.99%'],
['Fe','iron (alpha)','cubic bcc','MaTecK',MK+'elements/iron','standard size 6-10 mm dia x 1-2 mm','strain annealing','','purity 99.98%'],
['Ti','titanium (alpha)','hexagonal hcp','MaTecK',MK+'elements/titanium','standard size 8-10 mm dia x 1-2 mm','floating zone','','purity 99.97%'],
['Ag','silver','cubic fcc','MaTecK',MK+'elements/silver','standard size about 8-20 mm dia x 1-2 mm','Czochralski','','purity 99.999%'],
['Au','gold','cubic fcc','MaTecK',MK+'elements/gold','standard size 10-20 mm dia x 1-2 mm','Czochralski','','purity 99.999%'],
['Zn','zinc','hexagonal hcp','MaTecK',MK+'elements/zinc','standard size 10-20 mm dia x 1-2 mm','Bridgman','','purity 99.999%'],
['Sn','tin (beta)','tetragonal (beta-Sn)','MaTecK',MK+'elements/tin','standard size 9-10 mm dia x 1-2 mm','Bridgman','','purity 99.999%'],
['Mg','magnesium','hexagonal hcp','MaTecK',MK+'elements/magnesium','standard size 10-20 mm dia x 1-2 (thickness unit not stated)','Bridgman','','purity 99.999%'],
['Zr','zirconium (alpha)','hexagonal hcp','MaTecK',MK+'elements/zirconium','standard size 3-4 mm dia x 1 mm','floating zone','','small (below 10 mm criterion); 99.99% "ex Hf"'],
['V','vanadium','cubic bcc','MaTecK',MK+'elements/vanadium','standard size 12 mm dia x 1-2 mm','floating zone','','purity 99.99%'],
['Pt','platinum','cubic fcc','MaTecK',MK+'elements/platinum','standard size 10-12 mm dia x 1-2 mm','Czochralski','','purity 99.99%'],
['Pd','palladium','cubic fcc','MaTecK',MK+'elements/palladium','standard size 6-12 mm dia x 1-2 mm','Czochralski','','purity 99.999%'],
['Co','cobalt','hexagonal hcp (RT)','MaTecK',MK+'elements/cobalt','standard size 8-10 mm dia x 1-2 mm','Czochralski','','purity 99.99%'],
['Cr','chromium','cubic bcc','MaTecK',MK+'elements/chromium','standard size 6-10 mm dia x 1-2 mm','floating zone','','purity 99.999%'],
['Sb','antimony','trigonal (A7)','MaTecK',MK+'elements/antimony','standard size "diameter 12, thickness 1" (units not stated)','Bridgman','','purity 99.999%'],
['Te','tellurium','trigonal','MaTecK',MK+'elements/tellurium','standard size 10 mm dia x 1-2 mm','Czochralski','','purity 99.999%'],
['In','indium','tetragonal','MaTecK',MK+'elements/indium','standard size 7 mm dia x 3 mm','Bridgman','','purity 99.9999%; below 10 mm'],
['Cd','cadmium','hexagonal hcp','MaTecK',MK+'elements/cadmium','standard size 12 mm dia x 1-2 mm','Bridgman','','purity 99.999%'],
['Re','rhenium','hexagonal hcp','MaTecK',MK+'elements/rhenium','standard size about 6 mm dia x 1-2 mm','floating zone','','below 10 mm'],
['Ir','iridium','cubic fcc','MaTecK',MK+'elements/iridium','standard size 8-10 mm dia x 1-2 (unit not stated)','floating zone','','purity 99.99%'],
['Hf','hafnium','hexagonal hcp','MaTecK',MK+'elements/hafnium','standard size 7-8 mm dia x 1 mm','floating zone','','below 10 mm'],
['Ru','ruthenium','hexagonal hcp','MaTecK',MK+'elements/ruthenium','standard size about 6-8 mm dia x 1-2 mm','floating zone','','below 10 mm'],
['Rh','rhodium','cubic fcc','MaTecK',MK+'metal-single-crystals','','','','listed in pure-metals single crystal table; element page not opened'],
['Y','yttrium','hexagonal hcp','MaTecK',MK+'metal-single-crystals','','','','listed in pure-metals single crystal table; element page not opened'],
['Li','lithium','cubic bcc (RT)','MaTecK',MK+'metal-single-crystals','','','','listed in pure-metals single crystal table; element page not opened'],
['Gd','gadolinium','hexagonal hcp','MaTecK',MK+'metal-single-crystals','','','','listed in pure-metals single crystal table; element page not opened'],
['Dy','dysprosium','hexagonal hcp','MaTecK',MK+'metal-single-crystals','','','','listed in pure-metals single crystal table; element page not opened'],
['Ho','holmium','hexagonal hcp','MaTecK',MK+'metal-single-crystals','','','','listed in pure-metals single crystal table; element page not opened'],
// ---------------- MaTecK: semiconductors ----------------
['Si','silicon','cubic (diamond)','MaTecK',MK+'semiconductor-single-crystals','dia 2 in - 6 in','CZ / FZ','',''],
['Ge','germanium','cubic (diamond)','MaTecK',MK+'semiconductor-single-crystals','dia 1 in - 4 in','CZ / Br (as written)','',''],
['GaAs','gallium arsenide','cubic (zinc blende)','MaTecK',MK+'semiconductor-single-crystals','dia 3 in','CZ, LEC & HB','',''],
['GaP','gallium phosphide','cubic (zinc blende)','MaTecK',MK+'semiconductor-single-crystals','dia 50 mm','CZ & LEC','',''],
['GaSb','gallium antimonide','cubic (zinc blende)','MaTecK',MK+'semiconductor-single-crystals','dia 3 in','CZ & LEC','',''],
['InAs','indium arsenide','cubic (zinc blende)','MaTecK',MK+'semiconductor-single-crystals','dia 30 mm','LEC','',''],
['InP','indium phosphide','cubic (zinc blende)','MaTecK',MK+'semiconductor-single-crystals','dia 2 in','LEC','',''],
['CdS','cadmium sulfide','hexagonal (wurtzite)','MaTecK',MK+'semiconductor-single-crystals','dia 2 in','Bridgman & Markov','',''],
['CdSe','cadmium selenide','hexagonal (wurtzite)','MaTecK',MK+'semiconductor-single-crystals','dia 2 in','Bridgman & Markov','',''],
['CdTe','cadmium telluride','cubic (zinc blende)','MaTecK',MK+'semiconductor-single-crystals','','Bridgman','',''],
['ZnS','zinc sulfide','cubic (zinc blende)','MaTecK',MK+'semiconductor-single-crystals','dia 30 mm','Bridgman & Markov','','single crystal (separate from CVD Cleartran optical grade on optical-crystals page)'],
['ZnSe','zinc selenide','cubic (zinc blende)','MaTecK',MK+'semiconductor-single-crystals','dia 30 mm','Bridgman & Markov','','single crystal (separate from CVD optical grade)'],
['ZnTe','zinc telluride','cubic (zinc blende)','MaTecK',MK+'semiconductor-single-crystals','dia 30 mm','Bridgman & Markov','',''],
// ---------------- MaTecK: oxides / compounds ----------------
['Al2O3','sapphire','trigonal (corundum)','MaTecK',MK+'compounds/aluminum-oxide','standard 10x10 mm (0.5 or 1.0 mm thick); custom on request','Czochralski, HEM, Verneuil (possible methods listed)','',''],
['MgO','magnesium oxide','cubic (rock salt)','MaTecK',MK+'compounds/magnesium-oxide','wafers up to 3 or 4 in dia; standard 10x10 mm','arc fusion (special arc melting)','',''],
['SrTiO3','strontium titanate','cubic perovskite','MaTecK',MK+'compounds/strontium-titanate','standard 10x10 mm, 0.5/1.0 mm thick','Verneuil (twin-free)','',''],
['ZnO','zinc oxide','hexagonal (wurtzite)','MaTecK',MK+'compounds/zinc-oxide','10x10 mm, 10x5 mm, 0.4 mm thick; max dia < 1 in','hydrothermal','',''],
['SiO2','alpha-quartz','trigonal (alpha-quartz)','MaTecK',MK+'compounds/quartz','','','','"Naturally birefringent hexagonal single crystal"'],
['ZrO2','YSZ (yttria-stabilized zirconia)','cubic fluorite (Y2O3-stabilized)','MaTecK',MK+'compounds/yttria-stabilized-zirconia','standard 10x10x0.5 mm','skull melting','','"supplied twin-free"; Y2O3 stabilizer'],
['LiF','lithium fluoride','cubic (rock salt)','MaTecK',MK+'compounds/lithium-fluoride','','','','listed on optical crystals page; size on request'],
['CdTe','cadmium telluride','cubic (zinc blende)','MaTecK',MK+'compounds/cadmium-telluride','','','','offered as optical window / IR lens / CdZnTe detector substrate basis'],
['BaB2O4','BBO (barium borate)','trigonal (beta-BBO, R3c)','MaTecK',MK+'oxide-single-crystals','','','','listed in oxide single crystal table (phase not specified there)'],
['Bi12SiO20','BSO (bismuth silicon oxide)','cubic (sillenite)','MaTecK',MK+'oxide-single-crystals','','','',''],
['Bi4Ge3O12','BGO','cubic (eulytite)','MaTecK',MK+'oxide-single-crystals','','','',''],
['Bi12GeO20','BGO20','cubic (sillenite)','MaTecK',MK+'oxide-single-crystals','','','',''],
['BaTiO3','barium titanate','tetragonal perovskite (RT)','MaTecK',MK+'oxide-single-crystals','','','',''],
['CdWO4','cadmium tungstate','monoclinic (wolframite)','MaTecK',MK+'oxide-single-crystals','','','',''],
['CaCO3','calcite','trigonal (calcite)','MaTecK',MK+'oxide-single-crystals','','','',''],
['Cr2O3','chromium(III) oxide (eskolaite)','trigonal (corundum)','MaTecK',MK+'oxide-single-crystals','','','',''],
['CoO','cobalt(II) oxide','cubic (rock salt)','MaTecK',MK+'oxide-single-crystals','','','',''],
['Gd3Ga5O12','GGG','cubic (garnet)','MaTecK',MK+'oxide-single-crystals','','','',''],
['FeO','wustite','cubic (rock salt, non-stoichiometric Fe1-xO)','MaTecK',MK+'oxide-single-crystals','','','',''],
['Fe2O3','hematite','trigonal (corundum)','MaTecK',MK+'oxide-single-crystals','','','',''],
['Fe3O4','magnetite','cubic (inverse spinel)','MaTecK',MK+'oxide-single-crystals','','','',''],
['LaAlO3','lanthanum aluminate','rhombohedral perovskite','MaTecK',MK+'oxide-single-crystals','','','',''],
['(LaAlO3)0.3(Sr2AlTaO6)0.7','LSAT','cubic perovskite (solid solution)','MaTecK',MK+'oxide-single-crystals','','','',''],
['PbWO4','lead tungstate','tetragonal (scheelite)','MaTecK',MK+'oxide-single-crystals','','','',''],
['LiAlO2','lithium aluminate','tetragonal (gamma-LiAlO2)','MaTecK',MK+'oxide-single-crystals','','','',''],
['LiGaO2','lithium gallate','orthorhombic','MaTecK',MK+'oxide-single-crystals','','','',''],
['LiNbO3','lithium niobate','trigonal (R3c)','MaTecK',MK+'oxide-single-crystals','','','',''],
['LiTaO3','lithium tantalate','trigonal (R3c)','MaTecK',MK+'oxide-single-crystals','','','',''],
['LiB3O5','LBO (lithium triborate)','orthorhombic (Pna21)','MaTecK',MK+'oxide-single-crystals','','','',''],
['Lu2SiO5','LSO:Ce','monoclinic','MaTecK',MK+'oxide-single-crystals','','','','Ce-doped scintillator'],
['MgAl2O4','spinel','cubic (spinel)','MaTecK',MK+'oxide-single-crystals','','','',''],
['MnO','manganese(II) oxide','cubic (rock salt)','MaTecK',MK+'oxide-single-crystals','','','',''],
['NdCaAlO4','neodymium calcium aluminate','tetragonal (K2NiF4 type)','MaTecK',MK+'oxide-single-crystals','','','',''],
['NdGaO3','neodymium gallate','orthorhombic perovskite','MaTecK',MK+'oxide-single-crystals','','','',''],
['NiO','nickel(II) oxide','cubic (rock salt; rhombohedral distortion below T_N)','MaTecK',MK+'oxide-single-crystals','','','',''],
['Pb(Mg1/3Nb2/3)O3-PbTiO3','PMN-PT','perovskite relaxor solid solution','MaTecK',MK+'oxide-single-crystals','','','','composition (PT fraction) not stated'],
['KNbO3','potassium niobate','orthorhombic perovskite (RT)','MaTecK',MK+'oxide-single-crystals','','','',''],
['KTaO3','potassium tantalate','cubic perovskite','MaTecK',MK+'oxide-single-crystals','','','',''],
['KTiOPO4','KTP (potassium titanyl phosphate)','orthorhombic (Pna21)','MaTecK',MK+'oxide-single-crystals','','','',''],
['SrxBa1-xNb2O6','SBN','tetragonal tungsten bronze','MaTecK',MK+'oxide-single-crystals','','','',''],
['SrLaAlO4','strontium lanthanum aluminate','tetragonal (K2NiF4 type)','MaTecK',MK+'oxide-single-crystals','','','',''],
['SrLaGaO4','strontium lanthanum gallate','tetragonal (K2NiF4 type)','MaTecK',MK+'oxide-single-crystals','','','',''],
['TeO2','paratellurite','tetragonal','MaTecK',MK+'oxide-single-crystals','','','',''],
['Tb3Ga5O12','TGG','cubic (garnet)','MaTecK',MK+'oxide-single-crystals','','','',''],
['TbScO3','terbium scandate','orthorhombic perovskite','MaTecK',MK+'oxide-single-crystals','','','',''],
['TiO2','rutile','tetragonal (rutile)','MaTecK',MK+'oxide-single-crystals','','','',''],
['Y3Al5O12','YAG','cubic (garnet)','MaTecK',MK+'oxide-single-crystals','','','',''],
['YAlO3','YAP','orthorhombic perovskite','MaTecK',MK+'oxide-single-crystals','','','',''],
['YVO4','yttrium orthovanadate','tetragonal (zircon)','MaTecK',MK+'oxide-single-crystals','','','',''],
// ---------------- MaTecK: optical crystals ----------------
['BaF2','barium fluoride','cubic (fluorite)','MaTecK',MK+'optical-crystals','','','','optical-crystals page does not state each item is single crystal'],
['CaF2','calcium fluoride','cubic (fluorite)','MaTecK',MK+'optical-crystals','','','','optical-crystals page does not state each item is single crystal'],
['MgF2','magnesium fluoride','tetragonal (rutile type)','MaTecK',MK+'optical-crystals','','','','optical-crystals page does not state each item is single crystal'],
['KBr','potassium bromide','cubic (rock salt)','MaTecK',MK+'optical-crystals','','','','optical-crystals page does not state each item is single crystal'],
['KCl','potassium chloride','cubic (rock salt)','MaTecK',MK+'optical-crystals','','','','optical-crystals page does not state each item is single crystal'],
['NaCl','sodium chloride','cubic (rock salt)','MaTecK',MK+'optical-crystals','','','','optical-crystals page does not state each item is single crystal'],
['CsI','caesium iodide','cubic (CsCl type)','MaTecK',MK+'optical-crystals','','','','optical-crystals page does not state each item is single crystal'],
['Tl(Br,I)','KRS-5','cubic (CsCl type)','MaTecK',MK+'optical-crystals','','','','TlBr-TlI mixed crystal'],
// ---------------- Princeton Scientific: metals ----------------
...['Al','Bi','Cd','Cr','Co','Cu','Dy','Gd','Au','Hf','Ho','In','Ir','Fe','Pb','Li','Mg','Mo','Ni','Nb','Pd','Pt','Re','Rh','Ru','Ag','Ta','Sn','Ti','W','V','Y','Zn','Zr'].map(el => {
  const names = {Al:'aluminium',Bi:'bismuth',Cd:'cadmium',Cr:'chromium',Co:'cobalt',Cu:'copper',Dy:'dysprosium',Gd:'gadolinium',Au:'gold',Hf:'hafnium',Ho:'holmium',In:'indium',Ir:'iridium',Fe:'iron',Pb:'lead',Li:'lithium',Mg:'magnesium',Mo:'molybdenum',Ni:'nickel',Nb:'niobium',Pd:'palladium',Pt:'platinum',Re:'rhenium',Rh:'rhodium',Ru:'ruthenium',Ag:'silver',Ta:'tantalum',Sn:'tin',Ti:'titanium',W:'tungsten',V:'vanadium',Y:'yttrium',Zn:'zinc',Zr:'zirconium'};
  const phase = {Al:'cubic fcc',Bi:'trigonal (A7)',Cd:'hexagonal hcp',Cr:'cubic bcc',Co:'hexagonal hcp (RT)',Cu:'cubic fcc',Dy:'hexagonal hcp',Gd:'hexagonal hcp',Au:'cubic fcc',Hf:'hexagonal hcp',Ho:'hexagonal hcp',In:'tetragonal',Ir:'cubic fcc',Fe:'cubic bcc',Pb:'cubic fcc',Li:'cubic bcc (RT)',Mg:'hexagonal hcp',Mo:'cubic bcc',Ni:'cubic fcc',Nb:'cubic bcc',Pd:'cubic fcc',Pt:'cubic fcc',Re:'hexagonal hcp',Rh:'cubic fcc',Ru:'hexagonal hcp',Ag:'cubic fcc',Ta:'cubic bcc',Sn:'tetragonal (beta-Sn)',Ti:'hexagonal hcp',W:'cubic bcc',V:'cubic bcc',Y:'hexagonal hcp',Zn:'hexagonal hcp',Zr:'hexagonal hcp'};
  let note = 'page-level statement: "Depending on the metal the rod diameter can be set to values between several millimeters to 100 mm"; per-element size not given';
  let max = '';
  if (el==='Ti'||el==='Zr') { max = 'crystal size limited to a few millimeters (page footnote)'; note = 'flagged on page: "Crystal size limited to a few millimeters"'; }
  return [el, names[el], phase[el], 'Princeton Scientific', PS+'metal-single-crystal-catalog', max, 'Bridgman, Czochralski and zone melting (page-level, not per element)', '', note];
}),
// ---------------- Princeton Scientific: semiconductors ----------------
['GaAs','gallium arsenide','cubic (zinc blende)','Princeton Scientific',PS+'semiconductor-single-crystals','4 in dia x 0.6 mm','CZ, LEC, HB','',''],
['GaP','gallium phosphide','cubic (zinc blende)','Princeton Scientific',PS+'semiconductor-single-crystals','2 in dia x 0.5 mm','CZ & LEC','',''],
['GaSb','gallium antimonide','cubic (zinc blende)','Princeton Scientific',PS+'semiconductor-single-crystals','3 in dia x 0.45 mm','CZ & LEC','',''],
['Ge','germanium','cubic (diamond)','Princeton Scientific',PS+'semiconductor-single-crystals','4 in dia x 0.5 mm','CZ/Br (as written)','',''],
['InAs','indium arsenide','cubic (zinc blende)','Princeton Scientific',PS+'semiconductor-single-crystals','2 in dia x 0.5 mm','LEC','',''],
['InP','indium phosphide','cubic (zinc blende)','Princeton Scientific',PS+'semiconductor-single-crystals','2 in dia x 0.5 mm','LEC','',''],
['InSb','indium antimonide','cubic (zinc blende)','Princeton Scientific',PS+'semiconductor-single-crystals','2 in dia x 0.5 mm','','',''],
['Si','silicon','cubic (diamond)','Princeton Scientific',PS+'semiconductor-single-crystals','6 in dia x 0.6 mm','CZ / FZ','',''],
['CdS','cadmium sulfide','hexagonal (wurtzite)','Princeton Scientific',PS+'semiconductor-single-crystals','10x10x1 mm','PVT or Bridgman','',''],
['CdSe','cadmium selenide','hexagonal (wurtzite)','Princeton Scientific',PS+'semiconductor-single-crystals','10x10x1 mm','Bridgman or PVT','',''],
['CdTe','cadmium telluride','cubic (zinc blende)','Princeton Scientific',PS+'semiconductor-single-crystals','15x15x0.5 mm','Bridgman (vertical or horizontal)','',''],
['ZnS','zinc sulfide','cubic (zinc blende)','Princeton Scientific',PS+'semiconductor-single-crystals','10x10x1 mm','CVT or PVT','',''],
['ZnSe','zinc selenide','cubic (zinc blende)','Princeton Scientific',PS+'semiconductor-single-crystals','4 in dia x 1 mm (as listed)','CVT or PVT','',''],
['ZnTe','zinc telluride','cubic (zinc blende)','Princeton Scientific',PS+'semiconductor-single-crystals','10x10x1 mm','Bridgman or PVT','',''],
// ---------------- Princeton Scientific: oxides ----------------
...[
 ['BaTiO3','barium titanate','tetragonal perovskite (RT)'],['CoO','cobalt(II) oxide','cubic (rock salt)'],['Cr2O3','chromium(III) oxide','trigonal (corundum)'],
 ['Gd3Ga5O12','GGG','cubic (garnet)'],['FeO','wustite','cubic (rock salt, Fe1-xO)'],['LaAlO3','lanthanum aluminate','rhombohedral perovskite'],
 ['(LaAlO3)0.3(Sr2AlTaO6)0.7','LSAT','cubic perovskite (solid solution)'],['LiAlO2','lithium aluminate','tetragonal (gamma)'],['LiNbO3','lithium niobate','trigonal (R3c)'],
 ['LiTaO3','lithium tantalate','trigonal (R3c)'],['MgO','magnesium oxide','cubic (rock salt)'],['MnO','manganese(II) oxide','cubic (rock salt)'],
 ['NdCaAlO4','neodymium calcium aluminate','tetragonal (K2NiF4 type)'],['NdGaO3','neodymium gallate','orthorhombic perovskite'],['NiO','nickel(II) oxide','cubic (rock salt)'],
 ['KTaO3','potassium tantalate','cubic perovskite'],['SiO2','alpha-quartz','trigonal (alpha-quartz)'],['MgAl2O4','spinel','cubic (spinel)'],
 ['SrLaAlO4','strontium lanthanum aluminate','tetragonal (K2NiF4 type)'],['SrLaGaO4','strontium lanthanum gallate','tetragonal (K2NiF4 type)'],['SrTiO3','strontium titanate','cubic perovskite'],
 ['TiO2','rutile (page: "Titanium Dioxide")','tetragonal (rutile assumed; polymorph not stated on page)'],['Y3Al5O12','YAG','cubic (garnet)'],['YAlO3','YAP','orthorhombic perovskite'],
 ['YVO4','yttrium orthovanadate','tetragonal (zircon)'],['ZnO','zinc oxide','hexagonal (wurtzite)'],['ZrO2','YSZ','cubic fluorite (Y-stabilized)']
].map(([f,n,p]) => [f,n,p,'Princeton Scientific',PS+'oxide-single-crystals','','','','"available according to your specifications"; no sizes on page']),
['C','graphite (HOPG)','hexagonal graphite, mosaic (not a true single crystal)','Princeton Scientific',PS+'oxide-single-crystals','','','','HOPG listed on oxide single crystal page; mosaic spread not stated'],
// ---------------- CASTECH ----------------
['KTiOPO4','KTP','orthorhombic (Pna21)','CASTECH',CT+'KTP---Potassium-Titanyl-Phosphate-186.html','up to 20x20x40 mm3, max length 60 mm','','',''],
['LiB3O5','LBO (lithium triborate)','orthorhombic (Pna21)','CASTECH',CT+'LBO---Lithium-Triborate-121.html','aperture up to 55x55 mm2, length up to 60 mm','flux growth','',''],
['BaB2O4','beta-BBO (beta-barium borate)','trigonal (R3c)','CASTECH',CT+'--BBO---Beta-Barium-Borate-120.html','up to 15x15x15 mm3; up to 25 mm long','','',''],
['BaB2O4','alpha-BBO (alpha-barium borate)','trigonal (R-3c, high-T phase)','CASTECH',CT+'--BBO---Alpha-Barium-Borate-90.html','aperture up to dia 50 mm, length up to 40 mm; ~1 kg blanks reported','Czochralski (stated as typical method for alpha-BBO)','','different polymorph from beta-BBO (same formula)'],
['BiB3O6','BIBO (bismuth triborate)','monoclinic (C2)','CASTECH',CT+'BIBO---Bismuth-Trborate-184.html','up to 10x10x15 mm3','top-seeded solution growth (TSSG)','',''],
['CsLiB6O10','CLBO (caesium lithium borate)','tetragonal (I-42d)','CASTECH',CT+'CLBO---Cesium-Lithium-Borate-185.html','','','','hygroscopic'],
['RbTiOPO4','RTP (rubidium titanyl phosphate)','orthorhombic (Pna21)','CASTECH',CT+'RTP---Rubidium-Titanyl-Phosphate-117.html','max length 25 mm at 5x5 mm2 aperture','','',''],
['KTiOAsO4','KTA (potassium titanyl arsenate)','orthorhombic (Pna21)','CASTECH',CT+'KTA---Potassium-Titanyl-Arsenate--KTiOAsO4--189.html','up to 10x10x30 mm3','','',''],
['KH2PO4','KDP','tetragonal (I-42d)','CASTECH',CT+'KDP-DKDP---Potassium-Dihydrogen-Dideuterium-Phosphate-112.html','','','','"supplies high quality KDP and DKDP crystals in large quantities"; hydrogenous (incoherent H scattering)'],
['KD2PO4','DKDP (deuterated KDP)','tetragonal (I-42d)','CASTECH',CT+'KDP-DKDP---Potassium-Dihydrogen-Dideuterium-Phosphate-112.html','','','','deuterated analogue of KDP (deuteration level not stated)'],
['LiIO3','lithium iodate','hexagonal (P63)','CASTECH',CT+'NLO-Crystals-33','','','','listed on NLO crystals list page; product page not opened'],
['LiNbO3','lithium niobate (incl. MgO:LiNbO3)','trigonal (R3c)','CASTECH',CT+'NLO-Crystals-33','','','','listed on NLO crystals list page'],
['YVO4','yttrium orthovanadate','tetragonal (zircon)','CASTECH',CT+'YVO4---Yttrium-Orthovanadate-92.html','bulk up to dia 35 x 50 mm; finished up to dia 20 x 20 mm','Czochralski','',''],
['BeAl2O4','alexandrite (Cr:BeAl2O4, chrysoberyl host)','orthorhombic (olivine type, Pnma)','CASTECH',CT+'Cr3--3ABeAl2O4-Alexandrite--238.html','rods up to 9.53 mm dia (x117 mm); longest 8x130 mm','','','Cr3+-doped; beryllium-containing'],
['LiYF4','YLF (Nd:YLF host)','tetragonal (scheelite)','CASTECH',CT+'Nd-3AYLF---Neodymium-Doped-Yttrium-Lithium-Fluoride-100.html','rods 2-10 mm dia, 1-150 mm long; larger on request','Czochralski','','Nd-doped'],
['GdVO4','gadolinium orthovanadate (Nd:GdVO4 host)','tetragonal (zircon)','CASTECH',CT+'Nd-3AGdVO4---Neodymium-Doped-Gadolinium-Orthovanadate-107.html','','','','Nd-doped'],
['KGd(WO4)2','KGW (Yb:KGW host)','monoclinic (C2/c)','CASTECH',CT+'Yb-3AKGW---Ytterbium-Doped-Potassium-Gadolinum-Tungstate-96.html','maximum length 50 mm','','','Yb-doped'],
['LiSrAlF6','LiSAF (Cr:LiSAF host)','trigonal (colquiriite)','CASTECH',CT+'Cr-3ALiSAF---Cr-doped-Colquiriite-98.html','rods up to 16 mm dia, 60 mm length','Czochralski','','Cr3+-doped'],
['CaGdAlO4','CALGO (Yb:CALGO host)','tetragonal (K2NiF4 type)','CASTECH',CT+'Yb-3ACALGO---Ytterbium-Doped-Calcium-Gadolinium-Aluminate-110.html','boules up to dia 30 x 50 mm3','','','Yb-doped'],
['CaF2','calcium fluoride (Yb:CaF2)','cubic (fluorite)','CASTECH',CT+'Laser-Crystals-34','','','','Yb-doped laser crystal on laser crystal list page'],
// ---------------- Shalom EO ----------------
['LaCl3','lanthanum chloride (LaCl3:Ce)','hexagonal (UCl3 type)','Shalom EO',SH+'Scintillators/LaCl3-Ce/product-396.html','dia 75 mm x 75 mm','','','"Single crystals with Aluminum housing"; Ce-doped; hygroscopic'],
['Y2SiO5','YSO:Ce (yttrium oxyorthosilicate)','monoclinic','Shalom EO',SH+'Scintillators/YSO-Ce/product-400.html','dia 80 mm x 160 mm','','','Ce-doped; formula not printed on page'],
['Bi4Ge3O12','BGO','cubic (eulytite)','Shalom EO',SH+'Scintillators/BGO/product-388.html','maximum diameter 200 mm','Bridgman','',''],
['ZnGeP2','ZGP (zinc germanium phosphide)','tetragonal (chalcopyrite)','Shalom EO',SH+'Laser-Crystals-and-Components/ZGP-Crystals/product-334.html','7x5x25 mm (largest listed element)','','',''],
['YCa4O(BO3)3','YCOB (yttrium calcium oxyborate)','monoclinic (Cm)','Shalom EO',SH+'Laser-Crystals-and-Components/YCOB-Crystals/product-333.html','maximum diameter 100 mm','Czochralski','',''],
['YAlO3','YAP:Ce','orthorhombic perovskite','Shalom EO',SH+'Scintillators/YAP-Ce/product-399.html','dia 50 mm x 160 mm','Czochralski','','Ce 0.2-1.2 at%'],
['Lu3Al5O12','LuAG:Ce','cubic (garnet)','Shalom EO',SH+'Scintillators/LuAG-Ce/product-397.html','screens up to 4 in dia','','','Ce-doped'],
['LiIO3','lithium iodate','hexagonal (P63)','Shalom EO',SH+'Laser-Crystals-and-Components/LiIO3-Crystals/product-329.html','','','',''],
['Gd3Al2Ga3O12','GAGG:Ce','cubic (garnet)','Shalom EO',SH+'Scintillators','','','','listed on scintillator list page'],
['(Lu,Y)2SiO5','LYSO:Ce','monoclinic','Shalom EO',SH+'Scintillators','','','','listed on scintillator list page; Lu/Y ratio not stated'],
['LaBr3','lanthanum bromide (LaBr3:Ce)','hexagonal (UCl3 type)','Shalom EO',SH+'Scintillators','','','','listed on scintillator list page'],
['CeBr3','cerium bromide','hexagonal (UCl3 type)','Shalom EO',SH+'Scintillators','','','','listed on scintillator list page'],
['NaI','sodium iodide (NaI:Tl)','cubic (rock salt)','Shalom EO',SH+'Scintillators','encapsulated NaI(Tl) series maximum size 950x950x300 mm','','','size refers to encapsulated detector series'],
['CdWO4','cadmium tungstate','monoclinic (wolframite)','Shalom EO',SH+'Scintillators','','','','listed on scintillator list page'],
['BaF2','barium fluoride','cubic (fluorite)','Shalom EO',SH+'Scintillators','','','','listed on scintillator list page'],
// ---------------- Epic Crystal ----------------
['NaI','sodium iodide (NaI:Tl)','cubic (rock salt)','Epic Crystal',EP+'nai-tl-crystal.html','largest ingot 150 mm dia x 400 mm','Bridgman','','Tl-doped'],
['(Lu,Y)2SiO5','LYSO:Ce','monoclinic','Epic Crystal',EP+'lysoce-crystal.html','ingot 90 mm dia x 200 mm (commercial volume)','Czochralski (CZ)','','Ce-doped'],
['LaBr3','lanthanum bromide (LaBr3:Ce)','hexagonal (UCl3 type)','Epic Crystal',EP+'labr3ce-crystal.html','max ingot 80 mm dia x 200 mm','Bridgman','','Ce-doped; hygroscopic'],
['CeBr3','cerium bromide','hexagonal (UCl3 type)','Epic Crystal',EP+'cebr3-crystal.html','max ingot about 100 mm dia x 150 mm','Bridgman','','hygroscopic'],
['Gd3Al2Ga3O12','GAGG:Ce','cubic (garnet)','Epic Crystal',EP+'gaggce-crystal.html','ingots around 60 mm dia x 200 mm','','','Ce-doped'],
['CdWO4','cadmium tungstate','monoclinic (wolframite)','Epic Crystal',EP+'cdwo4-crystal.html','largest boule 80 mm dia x 200 mm','Bridgman','',''],
['CsI','caesium iodide (CsI:Tl, CsI:Na, pure)','cubic (CsCl type)','Epic Crystal',EP+'csi-crystal.html','up to 150 mm dia x 400 mm (standard 90 x 320 mm)','Bridgman','',''],
// ---------------- Luxium Solutions (ex Saint-Gobain Crystals) ----------------
['Cs2LiLaBr6','CLLB:Ce (elpasolite)','cubic (elpasolite)','Luxium Solutions',LX+'/cllb','2 in dia x 2 in long (largest typical size)','','','Ce-doped'],
['NaI','sodium iodide (NaI:Tl)','cubic (rock salt)','Luxium Solutions',LX+'/naitl-scintillation-crystals','typical: cylinders up to 5 in dia x 5 in; blocks 4x4x16 in; plates 24x18x0.375 in','','','Polyscin variant is polycrystalline (recrystallized single-crystal ingot) - not single crystal'],
['LaBr3','lanthanum bromide (LaBr3:Ce)','hexagonal (UCl3 type)','Luxium Solutions',LX+'/lanthanum-bromide-labr3','standard to 3 in dia x 3 in; achieved 97 mm dia x 244 mm','','',''],
['Bi4Ge3O12','BGO','cubic (eulytite)','Luxium Solutions',LX+'/bgo-bismuth-germanate','typical 1-5 in dia, thickness up to 5 in','','',''],
['CsI','caesium iodide','cubic (CsCl type)','Luxium Solutions',LX,'','','','listed on crystal scintillator page'],
['(Lu,Y)2SiO5','LYSO:Ce','monoclinic','Luxium Solutions',LX,'','','','listed on crystal scintillator page'],
// ---------------- Crytur ----------------
['Y3Al5O12','YAG (undoped and doped)','cubic (garnet)','Crytur',CY,'','Czochralski (company-level statement for garnet and perovskite crystals)','',''],
['YAlO3','YAP (undoped and Ce-doped)','orthorhombic perovskite','Crytur',CY+'yap-undoped/','','','','formula YAlO3 given on page'],
['Lu3Al5O12','LuAG (Ce, Pr doped; Yb:LuAG)','cubic (garnet)','Crytur',CY,'','Czochralski (company-level statement for garnet and perovskite crystals)','',''],
['LuAlO3','LuAP:Ce (lutetium aluminium perovskite)','orthorhombic perovskite','Crytur',CY+'luap-ce/','','','','Ce-doped'],
['Gd3Al2Ga3O12','GAGG+ (Ce)','cubic (garnet)','Crytur',CY+'gagg-1/','','','','"Gadolinium aluminium gallium garnet doped with cerium"; Crytur page does not state the Al/Ga ratio (nominal GAGG formula used)'],
['Bi4Ge3O12','BGO','cubic (eulytite)','Crytur',CY+'bgo/','','','',''],
['PbWO4','PWO','tetragonal (scheelite)','Crytur',CY+'pwo/','','','',''],
['Al2O3','sapphire (Sapphit)','trigonal (corundum)','Crytur',CY+'sapphire/','','','','"grown as single crystal in a strictly controlled environment"'],
// ---------------- SurfaceNet ----------------
['CaWO4','calcium tungstate (scheelite)','tetragonal (scheelite)','SurfaceNet',"https://surfacenet.de/files/kr_Calcium_Tungstate.php",'','','','standard orientations (010),(100),(001),(110),(111); doped variants (Nd, Ce, Pr, Sm, Er, Eu) also listed; page does not use the words "single crystal"'],
['TiO2','rutile','tetragonal (rutile)','SurfaceNet','https://www.surfacenet.de/kristall_daten.php','','floating zone','','"Floating Zone grown TiO2 Rutile (substrate)"'],
['Cr2O3','chromium(III) oxide','trigonal (corundum)','SurfaceNet','https://www.surfacenet.de/kristall_daten.php','10x75 mm','','',''],
['KNbO3','potassium niobate','orthorhombic perovskite (RT)','SurfaceNet','https://www.surfacenet.de/kristall_daten.php','','','',''],
['GdScO3','gadolinium scandate','orthorhombic perovskite','SurfaceNet','https://www.surfacenet.de/kristall_daten.php','','grown in-house (method not named)','','formula taken from image caption/file name on page'],
['DyScO3','dysprosium scandate','orthorhombic perovskite','SurfaceNet','https://www.surfacenet.de/kristall_daten.php','','grown in-house (method not named)','','formula taken from image caption/file name on page'],
['Er2SiO5','erbium oxyorthosilicate','monoclinic','SurfaceNet','https://www.surfacenet.de/kristall_daten.php','','','','listed in current-crystals data page; no specs'],
// ---------------- Korth Kristalle ----------------
['SrF2','strontium fluoride','cubic (fluorite)','Korth Kristalle','https://korth.de/material/detail/Strontiumfluorid','','','','Korth grows crystals itself ("Wir zuechten Kristalle"); page mentions seeded growth'],
...[
 ['BaF2','barium fluoride','cubic (fluorite)'],['CsI','caesium iodide','cubic (CsCl type)'],['KBr','potassium bromide','cubic (rock salt)'],['KCl','potassium chloride','cubic (rock salt)'],
 ['CaF2','calcium fluoride','cubic (fluorite)'],['Tl(Br,I)','KRS-5','cubic (CsCl type)'],['LiF','lithium fluoride','cubic (rock salt)'],['MgF2','magnesium fluoride','tetragonal (rutile type)'],
 ['NaCl','sodium chloride','cubic (rock salt)'],['NaF','sodium fluoride','cubic (rock salt)'],['RbI','rubidium iodide','cubic (rock salt)'],['AgBr','silver bromide','cubic (rock salt)'],['AgCl','silver chloride','cubic (rock salt)']
].map(([f,n,p]) => [f,n,p,'Korth Kristalle','https://www.korth.de/','','grown in-house ("Eigenzucht"; method not named on home page)','','marked as own-grown material on Korth home page material list; single-crystal wording not on page']),
['Tl(Br,Cl)','KRS-6 (thallium bromo-chloride)','cubic (CsCl type)','Korth Kristalle','https://www.korth.de/','','','','listed on material list (not marked own-grown); composition not stated on page'],
// ---------------- Others ----------------
['SiC','silicon carbide (4H)','hexagonal (4H polytype)','Wolfspeed','https://www.wolfspeed.com/materials/products/sic-substrates','150 mm and 200 mm dia wafers','','','4H n-type, 4 deg off-axis substrates'],
['AlN','aluminium nitride','hexagonal (wurtzite)','HexaTech','https://www.hexatechinc.com/aln-wafer-sales.html','','','','"single crystalline AlN substrates in c-plane and m-plane orientation"; sizes/prices on request; home page cites proprietary growth of single crystalline AlN boules'],
['PbMoO4','lead molybdate (wulfenite)','tetragonal (scheelite, I41/a)','Isomet','https://mail.isomet.com/App-Manual_pdf/Lead%20Molyb_DS902.pdf','typical boule 1 in dia x 4 in long','Czochralski','','Isomet data sheet (acousto-optic device maker); not a price catalogue'],
];

const header = ['formula','common_name','crystal_system_or_phase','vendor','catalogue_url','date_accessed','max_size_stated','growth_method_stated','price_stated','notes'];
const q = s => { s = String(s ?? ''); return /[",\n]/.test(s) ? '"' + s.replace(/"/g,'""') + '"' : s; };
const lines = [header.join(',')];
const compounds = new Set();
for (const r of R) {
  if (r.length !== 9) throw new Error('bad row ' + JSON.stringify(r));
  const [f,n,p,v,u,mx,gm,pr,nt] = r;
  const ab = absorberNote(f);
  const notes = [nt, ab].filter(Boolean).join('; ');
  lines.push([f,n,p,v,u,DATE,mx,gm,pr,notes].map(q).join(','));
  let key = f;
  if (f==='C') key = /graphite/i.test(n) ? 'C-graphite' : 'C-diamond';
  if (f==='BaB2O4') key = /alpha/i.test(n) ? 'BaB2O4-alpha' : 'BaB2O4-beta';
  compounds.add(key);
}
fs.writeFileSync(OUT, lines.join('\r\n') + '\r\n', 'utf8');
const vendors = [...new Set(R.map(r=>r[3]))];
console.log('rows', R.length, 'compounds', compounds.size, 'vendors', vendors.length, vendors.join(' | '));
// per-compound vendor counts
const per = {};
for (const r of R){ let key=r[0]; if(key==='C') key=/graphite/i.test(r[1])?'C-graphite':'C-diamond'; (per[key]=per[key]||new Set()).add(r[3]); }
const single = Object.entries(per).filter(([k,s])=>s.size===1).map(([k,s])=>k+':'+[...s][0]);
console.log('single-vendor compounds', single.length); console.log(single.join('; '));
