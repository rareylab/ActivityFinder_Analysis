from PIL import ImageDraw, ImageFont
from rdkit import Chem
from rdkit.Chem import Draw

font = ImageFont.truetype("/usr/share/fonts/dejavu/DejaVuSans.ttf", 25)


def plot_two_inchis(
    inchi_chembl, inchi_pdb, output_path, title="this is a test", size=400
):
    chembl_inchi = Chem.MolFromInchi(inchi_chembl)
    chembl_inchi.SetProp("_Name", f"ChEMBL")
    pdb_inchi = Chem.MolFromInchi(inchi_pdb)
    pdb_inchi.SetProp("_Name", f"PDB")
    mols = [chembl_inchi, pdb_inchi]
    img = Draw.MolsToGridImage(
        mols,
        molsPerRow=2,
        subImgSize=(size, size),
        legends=[x.GetProp("_Name") for x in mols],
    )
    print(type(img))
    draw = ImageDraw.Draw(img)
    draw.text((size / 1.5, 0), title, (0, 0, 0), font=font)
    img.save(output_path)


# Beispiele von Molekülen die identisch sind aber nicht perfekten inchi match haben
# ZSB ist ein isotop von wasserstoffen unterschied
chembl_inchi = "InChI=1S/C23H24N8O3/c1-24-20-12-18(27-16-6-5-11-30(23(16)33)19-7-3-4-10-25-19)29-21-14(13-26-31(20)21)22(32)28-15-8-9-17(15)34-2/h3-7,10-13,15,17,24H,8-9H2,1-2H3,(H,27,29)(H,28,32)/t15-,17-/m1/s1/i2D3"
pdb_inchi = "InChI=1S/C23H24N8O3/c1-24-20-12-18(27-16-6-5-11-30(23(16)33)19-7-3-4-10-25-19)29-21-14(13-26-31(20)21)22(32)28-15-8-9-17(15)34-2/h3-7,10-13,15,17,24H,8-9H2,1-2H3,(H,27,29)(H,28,32)/t15-,17-/m1/s1"
plot_two_inchis(
    chembl_inchi,
    pdb_inchi,
    "examples_ehmki_isotope_difference2.png",
    "ABCD ZSB CHEMBL5401575",
)


# 2357 ist interessant, da ist chembl ein stabiles radikal und wir sagen es ist ein ungeladenes hydroxid ... naomi problem incoming ich will nicht mehr
chembl_inchi = "InChI=1S/C28H31N6O3/c1-27(2)16-21(17-28(3,4)34(27)37)31-26(36)19-8-5-7-18(13-19)23-15-25-29-12-11-24(33(25)32-23)30-20-9-6-10-22(35)14-20/h5-15,21,30,35H,16-17H2,1-4H3,(H,31,36)"
pdb_inchi =    "InChI=1S/C28H32N6O3/c1-27(2)16-21(17-28(3,4)34(27)37)31-26(36)19-8-5-7-18(13-19)23-15-25-29-12-11-24(33(25)32-23)30-20-9-6-10-22(35)14-20/h5-15,21,30,35,37H,16-17H2,1-4H3,(H,31,36)"
plot_two_inchis(
    chembl_inchi,
    pdb_inchi,
    "examples_ehmki_radikal_difference.png",
    "3KXZ 925 CHEMBL1163016",
)


# 2422 ist ein beispiel wo wir ein adamantan haben und wir schreiben das mit chiralen wasserstoffen und sie schreiben es ohne, wenn wir das dann beide mal mit unseren USMILES machen ist es natürlich identisch
chembl_inchi = "InChI=1S/C13H18N4O3S2/c14-22(19,20)12-17-16-11(21-12)15-10(18)13-4-7-1-8(5-13)3-9(2-7)6-13/h7-9H,1-6H2,(H2,14,19,20)(H,15,16,18)"
pdb_inchi = "InChI=1S/C13H18N4O3S2/c14-22(19,20)12-17-16-11(21-12)15-10(18)13-4-7-1-8(5-13)3-9(2-7)6-13/h7-9H,1-6H2,(H2,14,19,20)(H,15,16,18)/t7-,8+,9-,13-"
plot_two_inchis(
    chembl_inchi,
    pdb_inchi,
    "examples_ehmki_rs_adamantan_difference.png",
    "3MHC ARZ CHEMBL1209038",
)

# 2354 PDB 3kvw IRB ist ein schönes beispiel, unterschiedlich kekulisierte strukturen zwischen PDB und chembl

chembl_inchi = "InChI=1S/C17H10BrN3O4/c18-10-3-1-2-8-12(16(22)20-13(8)10)15-14(21-25)9-6-7(17(23)24)4-5-11(9)19-15/h1-6,19,25H,(H,20,22)(H,23,24)/b15-12-,21-14+"
pdb_inchi = "InChI=1S/C17H10BrN3O4/c18-10-3-1-2-8-12(16(22)20-13(8)10)15-14(21-25)9-6-7(17(23)24)4-5-11(9)19-15/h1-6,19,21,25H,(H,23,24)"
plot_two_inchis(
    chembl_inchi,
    pdb_inchi,
    "examples_ehmki_kekulisiert_difference.png",
    "3KVW IRB CHEMBL1233659",
)

# 8767 Es gibt viele beispiele von doppelbindungsstereo bei einer R1R2C=NNHR3 (Oder R1R2C=NR3) doppelindung. Ich glaube die kann über konjugation invertieren, haben wir da vlt. einen fall für in naomi?

chembl_inchi = "InChI=1S/C18H13Cl2N5O4S/c19-11-3-2-10(7-12(11)20)15-8-30-18(22-15)24-23-14(17(26)27)5-9-1-4-13(21)16(6-9)25(28)29/h1-4,6-8H,5,21H2,(H,22,24)(H,26,27)/b23-14-"
pdb_inchi = "InChI=1S/C18H13Cl2N5O4S/c19-11-3-2-10(7-12(11)20)15-8-30-18(22-15)24-23-14(17(26)27)5-9-1-4-13(21)16(6-9)25(28)29/h1-4,6-8H,5,21H2,(H,22,24)(H,26,27)/b23-14+"
plot_two_inchis(
    chembl_inchi,
    pdb_inchi,
    "examples_ehmki_ez_difference.png",
    "4TQC 34J CHEMBL3310215",
)

# 5656 Chembl hat ein isotop, wir nicht. Ist ein unterschied im inchi und wird im SMILES dann wahrscheinlich rausgeschmissen

chembl_inchi = "InChI=1S/C22H22F3N5O2/c1-28-5-7-29(8-6-28)18-12-19-20(10-16(22(23,24)25)11-21(19)27-14-18)26-13-15-3-2-4-17(9-15)30(31)32/h2-4,9-12,14,26H,5-8,13H2,1H3/i13D"
pdb_inchi = "InChI=1S/C22H22F3N5O2/c1-28-5-7-29(8-6-28)18-12-19-20(10-16(22(23,24)25)11-21(19)27-14-18)26-13-15-3-2-4-17(9-15)30(31)32/h2-4,9-12,14,26H,5-8,13H2,1H3"
plot_two_inchis(
    chembl_inchi,
    pdb_inchi,
    "examples_ehmki_isotope_difference.png",
    "4GG5 OJ3 CHEMBL3221131",
)

# 3705, CHEMBL2069945 mit hcl gemessen in CHEMBL2072943 wo alles mit hcl gemessen wurde. Es gibt natürlich die alternate form ohne gegenion, aber wenn man einfach nur blind INCHI nimmt verliert man diesen Datenpunkt
chembl_inchi = "InChI=1S/C11H16N6.ClH/c12-10-9-5-16-17(11(9)15-7-14-10)6-8-1-3-13-4-2-8;/h5,7-8,13H,1-4,6H2,(H2,12,14,15);1H"
pdb_inchi = "InChI=1S/C11H16N6/c12-10-9-5-16-17(11(9)15-7-14-10)6-8-1-3-13-4-2-8/h5,7-8,13H,1-4,6H2,(H2,12,14,15)"
plot_two_inchis(
    chembl_inchi,
    pdb_inchi,
    "examples_ehmki_ions_difference.png",
    "3UQG B5A CHEMBL2069937",
)

# Hier beginnen nun beispiele von molekülen die NICHT identisch sind

# Example of hydroxyl->carbonyl where i checked that the pdb is a hydroxyl
chembl_inchi = "InChI=1S/C10H11N5O4/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(18)6(17)4(1-16)19-10/h1-4,6-7,10,17-18H,(H2,11,12,13)/t4-,6+,7+,10+/m0/s1"
pdb_inchi = "InChI=1S/C10H13N5O4/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(18)6(17)4(1-16)19-10/h2-4,6-7,10,16-18H,1H2,(H2,11,12,13)/t4-,6-,7-,10-/m1/s1"
plot_two_inchis(
    chembl_inchi, pdb_inchi, "examples_hyroxyl_carbonyl.png", "1D4F ADN CHEMBL77518"
)

# Example of double_bond -> single bond
chembl_inchi = "InChI=1S/C16H18O/c1-12-10-15(11-13(2)16(12)17)9-8-14-6-4-3-5-7-14/h3-7,10-11,17H,8-9H2,1-2H3"
pdb_inchi = "InChI=1S/C16H16O/c1-12-10-15(11-13(2)16(12)17)9-8-14-6-4-3-5-7-14/h3-11,17H,1-2H3/b9-8+"
plot_two_inchis(
    chembl_inchi,
    pdb_inchi,
    "examples_double_bond_difference.png",
    "3CN0 LJ1 CHEMBL493267",
)

# Example of nitrile -> amine bond
chembl_inchi = "InChI=1S/C20H21F3N4/c21-20(22,23)16-7-4-6-15(12-16)18-13-17(25-19(14-24)26-18)8-5-11-27-9-2-1-3-10-27/h4,6-7,12-13H,1-3,5,8-11H2"
pdb_inchi = "InChI=1S/C20H25F3N4/c21-20(22,23)16-7-4-6-15(12-16)18-13-17(25-19(14-24)26-18)8-5-11-27-9-2-1-3-10-27/h4,6-7,12-13H,1-3,5,8-11,14,24H2"
plot_two_inchis(
    chembl_inchi, pdb_inchi, "examples_nitril_amine.png", "3KWZ KWZ CHEMBL1084400"
)

# Example of halogen difference
chembl_inchi = "InChI=1S/C23H16ClNO2/c1-14-21(23(26)27)19-13-18(24)11-12-20(19)25-22(14)17-9-7-16(8-10-17)15-5-3-2-4-6-15/h2-13H,1H3,(H,26,27)"
pdb_inchi = "InChI=1S/C23H16FNO2/c1-14-21(23(26)27)19-13-18(24)11-12-20(19)25-22(14)17-9-7-16(8-10-17)15-5-3-2-4-6-15/h2-13H,1H3,(H,26,27)"
plot_two_inchis(
    chembl_inchi, pdb_inchi, "examples_halogen_difference.png", "1D3G BRE CHEMBL398986"
)

# Example of oxygen-> sulfur exchange
chembl_inchi = "InChI=1S/C16H18N4OS/c1-11(2)13-7-3-4-8-14(13)17-15-18-16(20-19-15)22-10-12-6-5-9-21-12/h3-9,11H,10H2,1-2H3,(H2,17,18,19,20)"
pdb_inchi = "InChI=1S/C16H18N4S2/c1-11(2)13-7-3-4-8-14(13)17-15-18-16(20-19-15)22-10-12-6-5-9-21-12/h3-9,11H,10H2,1-2H3,(H2,17,18,19,20)"
plot_two_inchis(
    chembl_inchi, pdb_inchi, "examples_hetatom_exchange.png", "2OAZ I96 CHEMBL444369"
)

# Example of tautomerism apparently not accounted for
# chembl_inchi = "InChI=1S/C10H8O4/c11-8(6-9(12)10(13)14)7-4-2-1-3-5-7/h1-6,11H,(H,13,14)/b8-6-"
# pdb_inchi = "InChI=1S/C10H8O4/c11-8(6-9(12)10(13)14)7-4-2-1-3-5-7/h1-6,12H,(H,13,14)/b9-6-"
# plot_two_inchis(chembl_inchi, pdb_inchi, "examples_identical_molecules_badmatch.png", "4AWF XI7 CHEMBL16326")
# InChI=1S/C10H8O4/c11-8(6-9(12)10(13)14)7-4-2-1-3-5-7/h1-5H,6H2,(H,13,14)
# InChI=1S/C10H8O4/c11-8(6-9(12)10(13)14)7-4-2-1-3-5-7/h1-6,11H,(H,13,14)/b8-6-

# Nun kommen beispiele von 1/2 matches die aber gut sind

chembl_inchi = "InChI=1S/C15H10O8/c16-6-3-7(17)11-10(4-6)23-15(14(22)13(11)21)5-1-8(18)12(20)9(19)2-5/h1-4,15-20H"
pdb_inchi = "InChI=1S/C15H10O8/c16-6-3-7(17)11-10(4-6)23-15(14(22)13(11)21)5-1-8(18)12(20)9(19)2-5/h1-4,16-20,22H"
plot_two_inchis(
    chembl_inchi,
    pdb_inchi,
    "examples_identical_molecules_badmatch.png",
    "2O63 MYC CHEMBL198468",
)


chembl_inchi = (
    "InChI=1S/C11H14O2S/c1-5-7(2)6-11(4)9(12)8(3)10(13)14-11/h5-6,13H,1H2,2-4H3/b7-6+"
)
pdb_inchi = "InChI=1S/C11H14O2S/c1-5-7(2)6-11(4)9(12)8(3)10(13)14-11/h5-6,12H,1H2,2-4H3/b7-6+/t11-/m1/s1"
plot_two_inchis(
    chembl_inchi,
    pdb_inchi,
    "examples_identical_molecules_badmatch_2.png",
    "2VB8 TLM CHEMBL322740",
)


# InChI=1S/C30H31N3O3/c1-35-27-18-23-17-26-29(22-7-5-20(6-8-22)21-9-11-24(34)12-10-21)31-32-30(26)25(23)19-28(27)36-16-15-33-13-3-2-4-14-33/h5-12,18-19,34H,2-4,13-17H2,1H3,(H,31,32) # chembl
# InChI=1S/C30H31N3O3/c1-35-27-18-23-17-26-29(22-7-5-20(6-8-22)21-9-11-24(34)12-10-21)31-32-30(26)25(23)19-28(27)36-16-15-33-13-3-2-4-14-33/h5-12,17-19,31-32,34H,2-4,13-16H2,1H3
chembl_inchi = "InChI=1S/C30H31N3O3/c1-35-27-18-23-17-26-29(22-7-5-20(6-8-22)21-9-11-24(34)12-10-21)31-32-30(26)25(23)19-28(27)36-16-15-33-13-3-2-4-14-33/h5-12,18-19,34H,2-4,13-17H2,1H3,(H,31,32)"
pdb_inchi = "InChI=1S/C30H31N3O3/c1-35-27-18-23-17-26-29(22-7-5-20(6-8-22)21-9-11-24(34)12-10-21)31-32-30(26)25(23)19-28(27)36-16-15-33-13-3-2-4-14-33/h5-12,17-19,31-32,34H,2-4,13-16H2,1H3"
plot_two_inchis(
    chembl_inchi,
    pdb_inchi,
    "examples_identical_molecules_badmatch_3.png",
    "4FTN 2HK CHEMBL550479",
)
