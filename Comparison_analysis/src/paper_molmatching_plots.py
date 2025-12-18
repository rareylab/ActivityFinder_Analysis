from src.utils.visualisation import create_and_plot_mol_array_align, create_and_plot_mol_array
from rdkit.Chem import Draw
from rdkit import Chem
from src.utils import getter
from pathlib import Path

chembl_entries = ["CHEMBL4861071", "CHEMBL5206501", "CHEMBL5183847", "CHEMBL5183847"]
chembl_to_usmiles = getter.get_chembl_to_usmiles(recalc=False)

names = ["SDA" , "CHEMBL306195"]
mols = ["Clc1c(N)c(Cl)cc(S(=O)(=O)N)c1", "Nc1c(Br)cc(S(N)(=O)=O)cc1Br"]

create_and_plot_mol_array_align(mols, names, output_name="molmatching_level2_example.png", align=True, size_molecule=500)

names = ["UUJ" , "CHEMBL33719"]
mols = ["Brc1cc(C#N)c(O)cc1", "NCc1cc(Cl)ccc1O"]

create_and_plot_mol_array_align(mols, names, output_name="molmatching_level1_example.png", align=True, size_molecule=500)



names = ["2WE hetcode", "2WE in 4PTC"]
mols = [" COc1c(sc(n1)c2ccnc(c2)NC(=O)C3CC3)C(=O)N", "S1C(=NC(O)=C1C(=O)N)c2ccnc(NC(=O)C3CC3)c2"]
create_and_plot_mol_array_align(mols, names, output_name="2we_incomplete.png", align=True, size_molecule=400)

names = ["MYC hetcode", "MYC_A_414", "MYC_A_408"]
mols = [" c1c(cc(c(c1O)O)O)C2=C(C(=O)c3c(cc(cc3O2)O)O)O", "CCCCCCCCCCCC", "CCCCCCCCCCC"]

create_and_plot_mol_array_align(mols, names, output_name="myc_alkanes.png", align=False, size_molecule=500)


names = ["BRH from 1M5C", "CHEMBL73698"]
mols = ["BrC=1C(=O)NOC1C[C@H](N)C(=O)O", "CHEMBL73698"]
#create_and_plot_mol_array_align(mols, names, output_name="brh.svg", align=True)
create_and_plot_mol_array_align(mols, names, output_name="brh.png", align=True, size_molecule=300)

names = ["OJ3", "CHEMBL3221131"]
mols = ["C(Nc1cc(C(F)(F)F)cc2ncc(N3CCN(C)CC3)cc12)c1cccc([N+](=O)[O-])c1", "[2H]C(Nc1cc(C(F)(F)F)cc2ncc(N3CCN(C)CC3)cc12)c1cccc([N+](=O)[O-])c1"]
create_and_plot_mol_array_align(
    mols, names, output_name="oj3_solo.png", align=True
)
molecules = [Chem.MolFromSmiles(chembl_to_usmiles[x]) for x in chembl_entries]
for i,mol in enumerate(molecules):
    mol.SetProp("_Name", chembl_entries[i])
create_and_plot_mol_array_align(chembl_entries, chembl_entries, align=False, output_name="molekuele_mpro.svg")
create_and_plot_mol_array_align(chembl_entries, chembl_entries, align=False, output_name="molekuele_mpro.png", size_molecule=600)


names = ["H4B and CHEMBL1201774", "CHEMBL1201775"]
mols = ["C[C@@H]([C@@H]([C@H]1CNC2=C(N1)C(=O)NC(=N2)N)O)O", "CHEMBL1201775"]
create_and_plot_mol_array_align(mols, names, output_name="h4b.svg", align=True)
create_and_plot_mol_array_align(mols, names, output_name="h4b.png", align=True, size_molecule=800)

names = ["Rebastinib"]
mol = Chem.MolFromSmiles("CNC(=O)c1cc(Oc2ccc(NC(=O)Nc3cc(C(C)(C)C)nn3-c3ccc4ncccc4c3)c(F)c2)ccn1")

img = Draw.MolToImage(mol, size=[400,400])
img.save(Path("/work/gutermuth/activityfinderanalysis/all_output_files/molecule_difference_data") / "pictures" / "rebastinib.png")
names = ["P8T built from 6YQN and CHEMBL4303293", "Hetcode P8T and CHEMBL5197137"]
mols = [
    "Clc1ccc(cc1)C2=N[C@H](C3=NN=C(N3C=4SC(=C(C42)C)C)C)CC(=O)N",
    "Cc1c(sc-2c1C(=N[C@H](c3n2c(nn3)C)CC(=O)Nc4ccc(cc4)C(=O)Nc5ccccc5N)c6ccc(cc6)Cl)C",
]
create_and_plot_mol_array_align(
    mols, names, output_name="molecule_not_fully_modelled_p8t.png", size_molecule=400
)
names = ["826 built from 1I30 and CHEMBL123617", "Hetcode 826 and CHEMBL72232"]
mols = ["O=C(N)c1ccc(O)cc1", "CHEMBL72232"]
create_and_plot_mol_array_align(
    mols, names, output_name="molecule_not_fully_modelled_826.svg"
)
names = ["5LO built from 5EC9", "CHEMBL3781132"]
mols = [
    "O=C(O)c1ccc(cc1)[C@]23c4c(O[C@@H]3CC=C2)cc5c(c4)C(CCC5(C)C)(C)C",
    "CHEMBL3781132",
]
create_and_plot_mol_array_align(
    mols, names, output_name="molecule_differently_built_5lo.svg", align=False
)
names = ["Hetcode CQD and CHEMBL337173", "CQD built from 4LMT"]
mols = ["CHEMBL337173", "Cl[C@@H]1[C@H](O)[C@H]2[C@H](NCCC2)[C@H](O)[C@H]1NCCN3CCOCC3"]
create_and_plot_mol_array_align(
    mols, names, output_name="molecule_differently_built_cqd.svg", align=False
)
names = [
    "CHEMBL198468 found only by ActivityFinder",
    "CHEMBL164 found by both approaches",
]
mols = ["CHEMBL198468", "CHEMBL164"]
create_and_plot_mol_array_align(
    mols, names, output_name="MYC_found_tautomer.svg", align=True
)
names = ["CHEMBL374731", "LLT in 3HP1 or 3QEO"]
mols = ["CHEMBL374731", "O=C1N(C=C(C(=O)N1)C)[C@H]2O[C@H]([C@H](O)C2)CO"]
create_and_plot_mol_array_align(
    mols, names, output_name="LLT_only_we_found_this.svg", align=True, size_molecule=800
)
create_and_plot_mol_array_align(
    mols, names, output_name="LLT_only_we_found_this.png", align=True
)
names = ["3BT", "CHEMBL1527769"]
mols = ['[H]/N=C(\\N)/SCCc1cccc(c1)CCSC(=N)N', "CHEMBL1527769"]
create_and_plot_mol_array_align(mols, names, output_name="3bt.svg", align=True)

names = ["MC", "CHEMBL3350473"]
mols = ["CC1=C(C(=O)c2c(c3n(c2C1=O)C[C@@H]([C@@H]3O)N)COC(=O)N)N", "CHEMBL3350473"]
create_and_plot_mol_array_align(mols, names, output_name="mc.svg", align=True)


# Fälle die wir kategorisch zu einem anderen liganden mappen, da sie falsch gebaut wurden
# 826 1
create_and_plot_mol_array(
    chembl_id_us="CHEMBL123617",
    chembl_id_them="CHEMBL72232",
    usmiles_chembl_built="O=C(N)c1ccc(O)cc1",
    usmiles_pdb_built="O=C(N)c1ccc(O)cc1",
    output_name="826_1.png",
)
names = ["CHEMBL123617 us", "CHEMBL72232"]
molecules = ["CHEMBL123617", "CHEMBL72232"]
# create_and_plot_mol_array_lists()
# 826 2
create_and_plot_mol_array(
    chembl_id_us="CHEMBL2441358",
    chembl_id_them="CHEMBL72232",
    usmiles_chembl_built="S=C(N)c1ccc(O)cc1",
    usmiles_pdb_built="O=C(N)c1ccc(O)cc1",
    output_name="826_2.png",
)
# P8T
create_and_plot_mol_array(
    chembl_id_us="CHEMBL4303293",
    chembl_id_them="CHEMBL5197137",
    usmiles_chembl_built="Clc1ccc(cc1)C2=NC(C3=NN=C(N3C=4SC(=C(C42)C)C)C)CC(=O)N",
    usmiles_pdb_built="Clc1ccc(cc1)C2=N[C@H](C3=NN=C(N3C=4SC(=C(C42)C)C)C)CC(=O)N",
    output_name="P8T_1.png",
)
# Fälle in denen wir das richtige matchen, aber nicht mit nem inchi match und ich nicht verstehe warum
# 5SS ringsymetrie 1
# create_and_plot_mol_array(chembl_id_us="CHEMBL3780189", chembl_id_them="CHEMBL3780324", usmiles_chembl_built="O=C1NNc2c1cc(c3cc(ncc3)C(NC(=O)C4CCC(CN)CC4)Cc5ccccc5)cc2",usmiles_pdb_built="O=C1NNc2c1cc(c3cc(ncc3)[C@@H](NC(=O)C4CCC(CN)CC4)Cc5ccccc5)cc2", output_name="5SS_1_nochiral.png")
# 5SS ringsymetrie 2
create_and_plot_mol_array(
    chembl_id_us="CHEMBL3780324",
    chembl_id_them="CHEMBL3780324",
    usmiles_chembl_built="O=C1NNc2c1cc(c3cc(ncc3)C(NC(=O)C4CCC(CN)CC4)Cc5ccccc5)cc2",
    usmiles_pdb_built="O=C1NNc2c1cc(c3cc(ncc3)[C@@H](NC(=O)C4CCC(CN)CC4)Cc5ccccc5)cc2",
    output_name="5SS_2_withchiral.png",
)
# QMW 1 nochiral
create_and_plot_mol_array(
    chembl_id_us="CHEMBL4793170",
    chembl_id_them="CHEMBL4793170",
    usmiles_chembl_built="O=C1C(=C2NNC(=N)C(N3CCC(O)CC3)=C2)C=CC=C1",
    usmiles_pdb_built="Oc1c(c2nnc(N)c(N3CCC(O)CC3)c2)cccc1",
    output_name="QMW_difference.png",
)
# QMW 2 clipped connection
# create_and_plot_mol_array(chembl_id_us="CHEMBL4762019", chembl_id_them="CHEMBL4793170", usmiles_chembl_built="O=C1C(=C2NNC(=N)C(N3CCC(=O)CC3)=C2)C=CC=C1",usmiles_pdb_built="Oc1c(c2nnc(N)c(N3CCC(O)CC3)c2)cccc1", output_name="QMW_2_clippedconnection.png")
# J2Y 1
# create_and_plot_mol_array(chembl_id_us="CHEMBL4280731", chembl_id_them="CHEMBL4288634", usmiles_chembl_built="[O+](=C1C=NC(=C2c3ncccc3NC2=C4C=C(NC(=O)C)NC=C4)C=C1)C",usmiles_pdb_built="O=C(Nc1nccc(c1)C=2Nc3c(nccc3)C2c4ncc(OC)cc4)C", output_name="J2Y_1.png")
# J2Y 1
create_and_plot_mol_array(
    chembl_id_us="CHEMBL4288634",
    chembl_id_them="CHEMBL4288634",
    usmiles_chembl_built="[O+](=C1C=NC(=C2c3ncccc3NC2=C4C=C(NC(=O)C)NC=C4)C=C1)C",
    usmiles_pdb_built="O=C(Nc1nccc(c1)C=2Nc3c(nccc3)C2c4ncc(OC)cc4)C",
    output_name="J2Y_difference.png",
)

# Fälle wo wir auch andere level finden da wir sie anders bauen als im hetcode
# 5LO (Doppelbindung mehr da sie da ist)
create_and_plot_mol_array(
    chembl_id_us="CHEMBL3781132",
    usmiles_chembl_built="O=C(O)c1ccc(cc1)C23c4c(OC2CCC3)cc5c(c4)C(CCC5(C)C)(C)C",
    usmiles_pdb_built="O=C(O)c1ccc(cc1)[C@]23c4c(O[C@@H]3CC=C2)cc5c(c4)C(CCC5(C)C)(C)C",
    output_name="5LO_difference_doublebond.png",
)
# 71C doppelbindung mehr
create_and_plot_mol_array(
    chembl_id_us="CHEMBL596234",
    usmiles_chembl_built="O=C(O)c1c(cc(c(c1)C(=O)O)C(=O)N(C2c3c(cccc3)CCC2)Cc4cc(Oc5ccccc5)ccc4)C(=O)O",
    usmiles_pdb_built="O=C(O)c1c(cc(c(c1)C(=O)O)C(=O)N([C@@H]2c3c(cccc3)C=CC2)Cc4cc(Oc5ccccc5)ccc4)C(=O)O",
    output_name="71C_difference_doublebond.png",
)
# G5R ringe nicht aromatisch
create_and_plot_mol_array(
    chembl_id_us="CHEMBL2441358",
    chembl_id_them="CHEMBL72232",
    usmiles_chembl_built="S=C(N)c1ccc(O)cc1",
    usmiles_pdb_built="O=C(N)c1ccc(O)cc1",
    output_name="g5r_2.png",
)
# Z60 doppelbindung
create_and_plot_mol_array(
    chembl_id_us="CHEMBL2377859",
    usmiles_chembl_built="ClC=1SC(C(=O)C=2SC(=NCC=C)NC2N)=CC1",
    usmiles_pdb_built="ClC=1SC(C(=O)C=2SC(=NC2N)NCCC)=CC1",
    output_name="z60_difference_doublebond.png",
)
# 4JG ringsystem komisch
create_and_plot_mol_array(
    chembl_id_us="CHEMBL3410576",
    usmiles_chembl_built="[S+]1=C(Nc2cc3c(cc2)cccc3)NC(=C1C(=O)C)C",
    usmiles_pdb_built="S1C(=Nc2cc3c(cc2)cccc3)N=C([C@@H]1C(=O)C)C",
    output_name="4jg_difference_weird_ringsystem.png",
)
# BRH (Das ist ein falscher perfekter inchi match den wir richtig einordnen)
create_and_plot_mol_array(
    chembl_id_us="CHEMBL71073",
    chembl_id_them="CHEMBL73698",
    usmiles_chembl_built="BrC=1C(=O)NOC1CC(N)C(=O)O",
    usmiles_pdb_built="BrC=1C(=O)NOC1C[C@H](N)C(=O)O",
    output_name="BRH_difference_wefindstereo.png",
)

# Fall in denen wir einen weiteren richtigen liganden finden, da wir nicht nur inchi machen
# MYC finden wir noch ein tautomer und 2 andere liganden
create_and_plot_mol_array(
    chembl_id_us="CHEMBL198468",
    chembl_id_them="CHEMBL164",
    usmiles_chembl_built="O=C1c2c(OC(c3cc(O)c(O)c(O)c3)C1=O)cc(O)cc2O",
    usmiles_pdb_built="O=C1c2c(OC(c3cc(O)c(O)c(O)c3)=C1O)cc(O)cc2O",
    output_name="MYC_difference_tautomer.png",
)

# Fälle in denen nur wir die richtigen liganden finden
# LLT
create_and_plot_mol_array(
    chembl_id_us="CHEMBL374731",
    usmiles_chembl_built="O=C1N(C=C(C(=O)N1)C)C2OC(C(O)C2)CO",
    usmiles_pdb_built="O=C1N(C=C(C(=O)N1)C)[C@H]2O[C@H]([C@H](O)C2)CO",
    output_name="LLT_onlywefindthis.png",
)


# JDALSSGOBMTZEP-NWJWHWDBSA-N
# JDALSSGOBMTZEP-NWJWHWDBSA-N
# PQVTLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF
# PQVTLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF
