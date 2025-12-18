from collections import Counter

from PIL import ImageDraw
from matplotlib import pyplot as plt
from matplotlib_venn import venn2
from rdkit import Chem
from rdkit.Chem import Draw, AllChem, rdFMCS
from rdkit.Chem.Draw import rdMolDraw2D
from supervenn import supervenn
import plotnine as p9
import pandas as pd
from src import definitions
from src.utils import getter

def prepare_own_venn_diagramm(sets_of_sets, labels, output_name, title, overall_names):
    for i, _ in enumerate(sets_of_sets):
        assert len(sets_of_sets[i]) == len(labels)
        assert len(sets_of_sets[i]) == 2
    new_data_structure = []
    for i, _ in enumerate(sets_of_sets):
        a = sets_of_sets[i][0]
        b = sets_of_sets[i][1]
        c = a.intersection(b)
        total = len(a) + len(b) - len(c)
        new_data_structure.append([100*(len(a)-len(c))/total,len(a)-len(c), labels[0], overall_names[i],0])
        new_data_structure.append([100*(len(b)-len(c))/total,len(b)-len(c), labels[1], overall_names[i],2])
        new_data_structure.append([100*len(c)/total,len(c), "both", overall_names[i], 1])
    print(new_data_structure)
    new_data_structure = pd.DataFrame(new_data_structure, columns=["Percent","Total","Label", "Type", "Order"])
    new_data_structure.to_csv("/work/gutermuth/activityfinderanalysis/src/ownvenn_data.csv")
    plot_own_venn_diagramm(new_data_structure, output_name, title, order=overall_names)

def plot_own_venn_diagramm(dataframe, output_filename, title="", order=None,thousandspace=True, color=False):
    plot = p9.ggplot(dataframe, p9.aes(x="Type", y="Percent", fill="reorder(Label, Order)"))
    plot += p9.geom_col(alpha=1, width=0.8)
    plot += p9.coord_flip()
    column_label = "Total"
    if thousandspace:
        dataframe['label_formatted'] = dataframe[column_label].apply(lambda x: f"{x:,}".replace(",", " "))
        column_label = "label_formatted"
    if color:
        plot += p9.geom_text(p9.aes(label=column_label),color=dataframe["color"], show_legend=False, position=p9.positions.position_stack(vjust=0.5), size=7)
        #plot += p9.scale_color_manual({"high": "black", "low": "white"})  # new
    else:
        plot += p9.geom_text(p9.aes(label=column_label), position=p9.positions.position_stack(vjust=0.5), size=7)
    #plot += p9.guide_legend(title=title)
    plot += p9.labs(fill=" ", title=title, x=" ",y="Percentage")
    plot += p9.scale_fill_cmap_d()
    plot += p9.scale_x_discrete(limits=order)
    #plot += p9.theme(text=p9.element_text(size=8))
    plot.save(output_filename, dpi=800, width=6, height=3)



def plot_supervenn_diagram(set1, set2, name1, name2, output_filename, title=""):
    supervenn([set1, set2],[name1, name2], side_plots=False)
    plt.title(title)
    plt.savefig(output_filename)
    plt.close()


def plot_venn_diagram(
    set1, set2, name1, name2, output_filename, title_addition="", shift=False
):
    """
    Plots a Venn diagram for two sets and saves it to a file.

    Parameters:
    set1 (set): The first set.
    set2 (set): The second set.
    name1 (str): The name for the first set.
    name2 (str): The name for the second set.
    output_filename (str): The name of the output file (including extension, e.g., 'venn_diagram.png').
    """
    # Create Venn diagram
    venn_diagram = venn2(
        [set1, set2], (name1, name2), set_colors=("#fde725", "#3e4989"), alpha=0.7
    )
    if shift:
        venn_diagram.get_label_by_id("10").set_x(-0.65)

    # Add title
    #plt.title(f"Venn Diagram of {name1} and {name2} {title_addition}")

    # Save the figure
    plt.savefig(output_filename)

    plt.close()


def generate_piechart_hetcodes(hetcode_list, output_name):
    # Zähle die Frequenzen der Werte
    frequenzen = Counter(hetcode_list)

    # Extrahiere die Daten für das Tortendiagramm
    labels, sizes = zip(
        *sorted(frequenzen.items(), key=lambda item: item[1], reverse=True)
    )

    # Erstelle das Tortendiagramm
    plt.figure(figsize=(8, 6))  # Größe des Diagramms
    plt.pie(sizes, labels=labels, autopct="%1.1f%%", startangle=140)
    plt.axis("equal")  # Damit das Tortendiagramm rund aussieht
    plt.title(f"Frequenz der Werte total {len(hetcode_list)}")
    plt.savefig(
        definitions.plot_path / output_name, bbox_inches="tight"
    )  # Zeige das Diagramm an


def visualize_mols_in_df(
    dataframe,
    picture_name,
    number_of_entries_per_side=20,
    subimgsize=300,
    sort_by="confidence_level",
    skip_their_chembl=False
):
    molecules_to_showcase = []
    dataframe.sort_values(by=sort_by, inplace=True)
    molecule_difference_path = definitions.general_output_files / "molecule_difference_data"
    chembl_to_usmiles = getter.get_chembl_to_usmiles()
    picture_path = molecule_difference_path / "pictures"
    if not picture_path.exists():
        picture_path.mkdir()
    for index, row in dataframe.iterrows():
        chembl_id_we = row["chembl_id"]
        chembl_usmiles_us = row["usmiles_from_chembl"]
        chembl_id_them = row["CHEMBL_ID_MOL"]
        our_smiles = row["usmiles_from_pdb"]
        hetcode = row["ligand_id"]
        matching_level = row["confidence_level"]
        if (
            not chembl_id_them in chembl_to_usmiles
            or not chembl_id_we in chembl_to_usmiles
        ):
            if chembl_id_we in chembl_to_usmiles and skip_their_chembl:
                pass
            else:
                print("Skipping one entry here.....")
                continue
        our_chembl_mol = Chem.MolFromSmiles(chembl_to_usmiles[chembl_id_we])
        our_chembl_mol.SetProp("_Name", f"our_{chembl_id_we}")
        our_chembl_mol_built = Chem.MolFromSmiles(chembl_usmiles_us)
        our_chembl_mol_built.SetProp("_Name", f"built_{chembl_id_we}")
        our_pdb_mol = Chem.MolFromSmiles(our_smiles)
        our_pdb_mol.SetProp("_Name", f"our_smiles_{hetcode} {matching_level}")
        if not skip_their_chembl:
            their_chembl_mol = Chem.MolFromSmiles(chembl_to_usmiles[chembl_id_them])
            their_chembl_mol.SetProp("_Name", f"their_{chembl_id_them}")
            molecules_to_showcase.extend(
                [our_pdb_mol, our_chembl_mol, our_chembl_mol_built, their_chembl_mol]
            )
        else:
            molecules_to_showcase.extend([our_pdb_mol, our_chembl_mol, our_chembl_mol_built])
    dopts = rdMolDraw2D.MolDrawOptions()
    dopts.addStereoAnnotation = True
    if skip_their_chembl:
        molsperrow = 3
    else:
        molsperrow = 4

    n_mols = molsperrow * number_of_entries_per_side
    for i in range(0, int(len(molecules_to_showcase) / n_mols)):
        img = Draw.MolsToGridImage(
            molecules_to_showcase[i * n_mols : (i + 1) * n_mols],
            molsPerRow=molsperrow,
            subImgSize=(subimgsize, subimgsize),
            legends=[
                x.GetProp("_Name")
                for x in molecules_to_showcase[i * n_mols : (i + 1) * n_mols]
            ],
            drawOptions=dopts,
        )
        img.save(picture_path / f"{picture_name}_{i}.png")


def align_and_visualise_mols_mcs(
    molecules_to_showcase, subimgsize=400, molperrow=3, output_name="test", align=True
):
    molecule_difference_path = definitions.molecule_difference_path
    dopts = rdMolDraw2D.MolDrawOptions()
    dopts.addStereoAnnotation = True
    dopts.minFontSize = 20
    for mol in molecules_to_showcase:
        AllChem.Compute2DCoords(mol)
    if align:
        res = rdFMCS.FindMCS(molecules_to_showcase)
        mcs_mol = Chem.MolFromSmarts(res.smartsString)
        AllChem.Compute2DCoords(mcs_mol)
        for m in molecules_to_showcase:
            _ = AllChem.GenerateDepictionMatching2DStructure(m, mcs_mol)
    if output_name.endswith(".png"):
        img = Draw.MolsToGridImage(
            molecules_to_showcase,
            molsPerRow=molperrow,
            subImgSize=(subimgsize, subimgsize),
            legends=[x.GetProp("_Name") for x in molecules_to_showcase],
            drawOptions=dopts
        )
        draw = ImageDraw.Draw(img)
        img.save(molecule_difference_path / "pictures" / f"paper_example_{output_name}")
    elif output_name.endswith(".svg"):
        img = Draw.MolsToGridImage(
            molecules_to_showcase,
            molsPerRow=molperrow,
            legends=[x.GetProp("_Name") for x in molecules_to_showcase],
            drawOptions=dopts,
            useSVG=True
        )
        with open(molecule_difference_path / "pictures" / f"paper_example_{output_name}", "w") as f:
            f.write(img)


def create_and_plot_mol_array_align(molecule_names, names, output_name, align=True, size_molecule=800):
    molecules = []
    chembl_to_usmiles = getter.get_chembl_to_usmiles()
    if not len(molecule_names) == len(names):
        print(
            f"Number of molecules {len(molecule_names)} and number of names {len(names)}"
        )
        raise RuntimeError("Number of molecules uneuqal to number of names")
    for index, name in enumerate(molecule_names):
        if name.startswith("CHEMBL"):
            smiles = chembl_to_usmiles[name]
        else:
            smiles = name
        mol = Chem.MolFromSmiles(smiles)
        mol.SetProp("_Name", names[index])
        molecules.append(mol)
    align_and_visualise_mols_mcs(
        molecules, molperrow=len(molecules), output_name=output_name, align=align, subimgsize=size_molecule
    )


def create_and_plot_mol_array(
    chembl_id_us,
    usmiles_chembl_built,
    usmiles_pdb_built,
    output_name,
    chembl_id_them=None,
):
    chembl_to_usmiles = getter.get_chembl_to_usmiles()
    our_chembl = Chem.MolFromSmiles(chembl_to_usmiles[chembl_id_us])
    our_chembl.SetProp("_Name", f"our_{chembl_id_us}")
    built_chembl = Chem.MolFromSmiles(usmiles_chembl_built)
    built_chembl.SetProp("_Name", f"our_chembl_built")
    our_pdb = Chem.MolFromSmiles(usmiles_pdb_built)
    our_pdb.SetProp("_Name", f"our_pdb_built")
    if not chembl_id_them == None:
        their_chembl = Chem.MolFromSmiles(chembl_to_usmiles[chembl_id_them])
        their_chembl.SetProp("_Name", f"their_{chembl_id_them}")
        mols = [our_chembl, built_chembl, our_pdb, their_chembl]
        molsperrow = 4
    else:
        mols = [our_chembl, built_chembl, our_pdb]
        molsperrow = 3
    align_and_visualise_mols_mcs(mols, molperrow=molsperrow, output_name=output_name)