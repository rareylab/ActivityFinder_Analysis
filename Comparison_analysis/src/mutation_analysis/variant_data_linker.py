from src import definitions
import pandas as pd
import json
import numpy as np

from src.utils import db


def get_target_to_assay():
    file_path = definitions.automated_path / "assay_to_target.csv"
    if file_path.exists():
        data = pd.read_csv(file_path)
    else:
        sql_statament = """SELECT 
                a.chembl_id AS assay_id, 
                t.chembl_id AS target_id 
                FROM assay AS a
                JOIN target as t on t.tid = a.tid
        """
        data = db.query_activitydb(sql_statament)
        data.to_csv(file_path)
    data_dict = {}
    for index, row in data.iterrows():
        assay_id = row["assay_id"]
        target_id = row["target_id"]
        data_dict[assay_id] = target_id
    return data_dict


def get_hetcode_to_chembl():
    file_path = definitions.automated_path / "hetcode_to_chembl.csv"
    if file_path.exists():
        data = pd.read_csv(file_path)
    else:
        sql_statement = """SELECT 
            pl.ligand_id, 
            cl.chembl_id
        FROM 
            pdb_ligand AS pl
        JOIN 
            pdb_ligand_chembl_ligand AS pcl ON pl.pdblid = pcl.pdblid
        JOIN 
        chembl_ligand AS cl ON pcl.clid = cl.clid;"""
        data = db.query_activitydb(sql_statement)
        data.to_csv(file_path)
    data_dict = {}
    for index, row in data.iterrows():
        hetcode = row["ligand_id"]
        chembl_id = row["chembl_id"]
        data_dict[hetcode] = chembl_id
    return data_dict


def get_component_sequence_of_assay_chembl(assay_id):
    sql_statement = f"""select cs.sequence
                        from assays as a
                        JOIN target_dictionary td on a.tid = td.tid
                        JOIN target_components tc on td.tid = tc.tid
                        JOIN component_sequences cs on tc.component_id = cs.component_id
                        WHERE a.chembl_id = '{assay_id}';"""
    data = db.query_chembl(sql_statement)
    all_sequences = [x for x in data["sequence"]]
    if len(all_sequences) == 0:
        raise RuntimeError(
            f"No component sequence here for assay {assay_id} in the chembl"
        )
    return all_sequences


def get_component_sequence_of_assay_activitydb(assay_id):
    """
    WARNING: This function is not working correctly currently
    :param chembl_id:
    :return:
    """
    sql_statement = f"""SELECT
        a.chembl_id AS assay_id, t.chembl_id AS target_id, s.sequence AS sequence, s.sid AS sequence_id
        FROM assay AS a
        JOIN target AS t on a.tid = t.tid
        JOIN sequence_target as seq on seq.chembl_tid = t.chembl_tid
        JOIN sequence as s on s.sid = seq.sid
        WHERE a.chembl_id = '{assay_id}'
    """
    component_data = db.query_activitydb(sql_statement)
    component_data.to_csv("test.csv")
    component_sequences = [x for x in component_data["sequence"]]
    return component_sequences


def get_variant_sequence_activitydb(assay_id):
    sql_statament = f"""SELECT 
            s.sequence, a.chembl_id, v.variant_id, v.mutation
            FROM assay AS a
            JOIN chembl_variant_sequence as v on v.variant_id = a.variant_id
            JOIN sequence as s on s.sid = v.sid
            WHERE a.chembl_id = '{assay_id}'
    """

    data = db.query_activitydb(sql_statament)
    variant_sequence = [x for x in data["sequence"]]
    if len(variant_sequence) != 1:
        raise RuntimeError(
            f"Multiple or no variant sequences found {len(variant_sequence)} for assay {assay_id}"
        )
    variant_sequence = variant_sequence[0]
    mutation = [x for x in data["mutation"]][0]
    return variant_sequence, mutation


def get_variant_sequence_chembl(assay_id):
    sql_statament = f"""select vs.sequence, vs.mutation
                        FROM assays as a
                        JOIN variant_sequences vs on a.variant_id = vs.variant_id
                        WHERE a.chembl_id = '{assay_id}';
                            """
    data = db.query_chembl(sql_statament)
    variant_sequence = [x for x in data["sequence"]]
    if len(variant_sequence) != 1:
        raise RuntimeError(
            f"Multiple or no variant sequences found {len(variant_sequence)} for assay {assay_id}"
        )
    variant_sequence = variant_sequence[0]
    mutation = [x for x in data["mutation"]][0]
    return variant_sequence, mutation


def get_sequences_assay_chembl(assay_id):
    # first variant sequence again
    variant_sequence, mutation = get_variant_sequence_chembl(assay_id)
    component_sequence = get_component_sequence_of_assay_chembl(assay_id)
    return variant_sequence, mutation, component_sequence


def get_sequences_assay(assay_id):
    # First we get the variant sequence of the assay
    variant_sequence, mutation = get_variant_sequence_activitydb(assay_id)
    # Now we get the component sequence(s) of the target
    component_sequences = get_component_sequence_of_assay_activitydb(assay_id)
    return variant_sequence, mutation, component_sequences


def is_ligand_in_assay(assay_id, ligand_id):
    sql_statement = f"""select distinct md.chembl_id
                        FROM molecule_dictionary AS md
                        JOIN public.activities a on md.molregno = a.molregno
                        JOIN public.assays AS ass on ass.assay_id = a.assay_id
                        WHERE ass.chembl_id = '{assay_id}'"""
    data = db.query_chembl(sql_statement)
    all_molecules = [x for x in data["chembl_id"]]
    return ligand_id in all_molecules


# Variant
# PISPIXTVPVKLKPGMDGPKVKQWPLTEEKIKALTEICTEMEKEGKIEKIGPENPYNTPVFAIKKKDSTKWRKVVDFRELNKRTQDFWEVQLGIPHPAGLKKKKSVTVLDVGDAYFSVPLDKDFRKYTAFTIPSINNETPGIRYQYNVLPQGWKGSPAIFQSSMTKILEPFRKQNPDIVIIQYMDDLYVGSDLEIEQHRAKIEELRQHLLRWGFTTPDKKHQKEPPFLWMGYELHPDKWTVQPIVLPEKDSWTVNXXXX
# Component
# PISPIETVPVKLKPGMDGPKVKQWPLTEEKIKALVEICTEMEKEGKISKIGPENPYNTPVFAIKKKDSTKWRKLVDFRELNKRTQDFWEVQLGIPHPAGLKKRKSVTVLDVGDAYFSVPLDEDFRKYTAFTIPSINNETPGIRYQYNVLPQGWKGSPAIFQSSMTKILEPFRKQNPDIVIYQYMDDLYVGSDLEIGQHRTKIEELRQHLLRWGLTTPDKKHQKEPPFLWMGYELHPDKWTVQPIVLPEKDSWTVNDIQKLVGKLNWASQIYPGIRVRQLCKLLRGTKALTEVIPLTEEAELELAENREILKEPVHGVYYDPSKDLIAEIQKQGQGQWTYQIYQEPFKNLRTGKYARMRGAHTNDVKQLTEAVQKITTESIVIWGKTPKFKLPIQKETWETWWTEYWQATWIPEWEFVNTPPLVKLWYQLEKEPIVGAETFYVDGAANRETKLGKAGYVTNRGRQKVVTLTDTTNQKTELQAIYLALQDSGLEVNIVTDSQYALGIIQAQPDQSESELVNQIIEQLIKKEKVYLAWVPAHKGIGGNEQVDKLVSAGIRKVLFLDGID


recalc = True
hetcode_to_chembl = get_hetcode_to_chembl()
target_data = get_target_to_assay()

with open(
        definitions.output_files_path / "activityfinder_full.json", "r", encoding="utf-8"
) as file:
    activityfinder_dict = json.load(file)

with open(
        definitions.output_files_path / "activityfinder_exact.json", "r", encoding="utf-8"
) as file:
    activityfinder_dict_exact = json.load(file)

with open(
        definitions.output_files_path / "chembiograph.json", "r", encoding="utf-8"
) as file:
    chembiodata_dict = json.load(file)

found_combination_path = (
        definitions.output_files_path / "variant_advantage_found_combinations.json"
)
target_to_assay_dict_path = (
        definitions.output_files_path / "variant_advantage_target_to_assay.json"
)
relevant_mutation_data_path = (
        definitions.automated_path / "difference_variant_mutation_data_relevant.csv"
)

if (
    not found_combination_path.exists()
    or not target_to_assay_dict_path.exists()
    or not relevant_mutation_data_path.exists()
    or recalc
):
    data_with_activities = []
    found_combinations = set()
    target_to_assay_used = {}
    for chunk in pd.read_csv(
            definitions.automated_path / "difference_variant_mutation_data_nogaps.csv",
        chunksize=10000,
    ):
        for index, row in chunk.iterrows():
            pdb = row["PDB"]
            ligand = row["name"]
            assay = row["assay"]
            if ligand[0:3] not in hetcode_to_chembl:
                continue
            chembl_id = hetcode_to_chembl[ligand[0:3]]
            target_id = target_data[assay]
            if pdb in activityfinder_dict:
                if target_id in activityfinder_dict[pdb]:
                    if chembl_id in activityfinder_dict[pdb][target_id]:
                        data_with_activities.append([x for x in row])
                        current_combination = f"{pdb}:{target_id}:{chembl_id}"
                        found_combinations.add(current_combination)
                        key = f"{target_id}:{chembl_id}"
                        if key in target_to_assay_used:
                            target_to_assay_used[key].add(assay)
                        else:
                            target_to_assay_used[key] = {assay}

    with open(found_combination_path, "w") as f:
        json.dump([x for x in found_combinations], f)

    with open(target_to_assay_dict_path, "w") as f:
        json.dump(
            {x: [y for y in target_to_assay_used[x]] for x in target_to_assay_used}, f
        )
    data_with_activities = pd.DataFrame(
        data_with_activities, columns=[x for x in chunk.columns.values]
    )
    data_with_activities.to_csv(relevant_mutation_data_path)
else:
    with open(target_to_assay_dict_path, "r", encoding="utf-8") as file:
        target_to_assay_used = json.load(file)
    with open(found_combination_path, "r", encoding="utf-8") as file:
        found_combinations = json.load(file)
    data_with_activities = pd.read_csv(relevant_mutation_data_path)


print(len(found_combinations))
n_found = 0
n_not_found = 0
coupled_combinations = []
# final_data = []
final_data_columns = [
    "PDB",
    "Target",
    "Ligand",
    "BioChemGraph",
    "Exact",
    "Assays",
    "N_assays",
]
mutation_data_columns = [
    "PDB",
    "Target",
    "Ligand",
    "BioChemGraph",
    "Assay",
    "Mutations",
    "n_component_sequences",
    "variant",
    "component",
]
# detailed_mutation_data = []

for index, combination in enumerate(found_combinations):
    final_data = []
    detailed_mutation_data = []
    pdb, target, ligand = combination.split(":")
    recovered = False
    exact = False
    pdb = pdb.lower()
    key = f"{target}:{ligand}"
    assays = target_to_assay_used[key]
    kept_assays = []
    for assay in assays:
        if is_ligand_in_assay(assay, ligand):
            kept_assays.append(assay)
    if len(kept_assays) == 0:
        continue
    print(pdb, target, ligand, len(kept_assays))
    if (
        pdb in chembiodata_dict.keys()
        and target in chembiodata_dict[pdb].keys()
        and ligand in chembiodata_dict[pdb][target]
    ):
        n_found += 1
        coupled_combinations.append(combination)
        recovered = True
    else:
        n_not_found += 1
    pdb_upper = pdb.upper()
    if (
        pdb_upper in activityfinder_dict_exact.keys()
        and target in activityfinder_dict_exact[pdb_upper].keys()
        and ligand in activityfinder_dict_exact[pdb_upper][target]
    ):
        exact = True
    final_data.append(
        [pdb, target, ligand, recovered, exact, kept_assays, len(kept_assays)]
    )
    for assay in kept_assays:
        variant_sequence, mutation, component_sequences = get_sequences_assay_chembl(
            assay
        )
        detailed_mutation_data.append(
            [
                pdb,
                target,
                ligand,
                recovered,
                assay,
                mutation,
                len(component_sequences),
                variant_sequence,
                component_sequences,
            ]
        )
    final_data = pd.DataFrame(final_data, columns=final_data_columns)
    detailed_mutation_data = pd.DataFrame(
        detailed_mutation_data, columns=mutation_data_columns
    )
    print(final_data.shape, detailed_mutation_data.shape)
    if index == 0:
        final_data.to_csv(
            definitions.automated_path / "combinations_with_better_variant_matches.csv"
        )
        detailed_mutation_data.to_csv(
            definitions.automated_path / "detailed_mutation_data.csv"
        )
    else:
        final_data.to_csv(
            definitions.automated_path / "combinations_with_better_variant_matches.csv",
            mode="a",
            header=False,
        )
        detailed_mutation_data.to_csv(
            definitions.automated_path / "detailed_mutation_data.csv",
            mode="a",
            header=False,
        )


print(coupled_combinations)
print(n_found, n_not_found, len(found_combinations))
"""final_data = pd.DataFrame(final_data, columns=["PDB", "Target", "Ligand", "BioChemGraph", "Assays", "N_assays"])
final_data.to_csv(definitions.automated_path / "combinations_with_better_variant_matches.csv")
detailed_mutation_data = pd.DataFrame(detailed_mutation_data, columns=["PDB", "Target", "Ligand", "BioChemGraph", "Assay", "Mutations", "n_component_sequences", "variant", "component"])
detailed_mutation_data.to_csv(definitions.automated_path / "detailed_mutation_data.csv")"""
final_data = pd.read_csv(
    definitions.automated_path / "combinations_with_better_variant_matches.csv"
)
detailed_mutation_data = pd.read_csv(
    definitions.automated_path / "detailed_mutation_data.csv"
)
all_pdbs = set(x for x in final_data["PDB"])
ligands = set(x for x in final_data["Ligand"])
targets = set(x for x in final_data["Target"])
assays = [x for x in final_data["N_assays"]]
n_pdbs = len(all_pdbs)
n_ligands = len(ligands)
n_targets = len(targets)
print(
    n_pdbs,
    n_ligands,
    n_targets,
    np.mean(assays),
    np.median(assays),
    final_data.shape,
    data_with_activities.shape,
)


# Variant
# PQVTLWKRPFVTVKVGGQLKEALLDTGADVTIFEDINLPGRWKPKMVGGIGGFLKVREYDQVPIEICGHKVIGTVLIGPTPVNVIGRDTMTQIGCTLNF
# Component
# PQVTLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF
# Mutation: I84V
# PDB : 2qci

# Assay : CHEMBL5160886
# Variant
# MRPSGTAGAALLALLAALCPASRALEEKKVCQGTSNKLTQLGTFEDHFLSLQRMFNNCEVVLGNLEITYVQRNYDLSFLKTIQEVAGYVLIALNTVERIPLENLQIIRGNMYYENSYALAVLSNYDANKTGLKELPMRNLQEILHGAVRFSNNPALCNVESIQWRDIVSSDFLSNMSMDFQNHLGSCQKCDPSCPNGSCWGAGEENCQKLTKIICAQQCSGRCRGKSPSDCCHNQCAAGCTGPRESDCLVCRKFRDEATCKDTCPPLMLYNPTTYQMDVNPEGKYSFGATCVKKCPRNYVVTDHGSCVRACGADSYEMEEDGVRKCKKCEGPCRKVCNGIGIGEFKDSLSINATNIKHFKNCTSISGDLHILPVAFRGDSFTHTPPLDPQELDILKTVKEITGFLLIQAWPENRTDLHAFENLEIIRGRTKQHGQFSLAVVSLNITSLGLRSLKEISDGDVIISGNKNLCYANTINWKKLFGTSGQKTKIISNRGENSCKATGQVCHALCSPEGCWGPEPRDCVSCRNVSRGRECVDKCNLLEGEPREFVENSECIQCHPECLPQAMNITCTGRGPDNCIQCAHYIDGPHCVKTCPAGVMGENNTLVWKYADAGHVCHLCHPNCTYGCTGPGLEGCPTNGPKIPSIATGMVGALLLLLVVALGIGLFMRRRHIVRKRTLRRLLQERELVEPLTPSGEAPNQALLRILKETEFKKIKVLGSGAFGTVYKGLWIPEGEKVKIPVAIKELREATSPKANKEILDEAYVMASVDNPHVCRLLGICLTSTVQLITQLMPFGCLLDYVREHKDNIGSQYLLNWCVQIAKGMNYLEDRRLVHRDLAARNVLVKTPQHVKITDFGRAKLLGAEEKEYHAEGGKVPIKWMALESILHRIYTHQSDVWSYGVTVWELMTFGSKPYDGIPASEISSILEKGERLPQPPICTIDVYMIMVKCWMIDADSRPKFRELIIEFSKMARDPQRYLVIQGDERMHLPSPTDSNFYRALMDEEDMDDVVDADEYLIPQQGFFSSPSTSRTPLLSSLSATSNNSTVACIDRNGLQSCPIKEDSFLQRYSSDPTGALTEDSIDDTFLPVPEYINQSVPKRPAGSVQNPVYHNQPLNPAPSRDPHYQDPHSTAVGNPEYLNTVQPTCVNSTFDSPAHWAQKGSHQISLDNPDYQQDFFPKEAKPNGIFKGSTAENAEYLRVAPQSSEFIGA
# Component
# MPRRAENWDEAEVGAEEAGVEEYGPEEDGGEESGAEESGPEESGPEELGAEEEMEAGRPRPVLRSVNSREPSQVIFCNRSPRVVLPVWLNFDGEPQPYPTLPPGTGRRIHSYRGHLWLFRDAGTHDGLLVNQTELFVPSLNVDGQPIFANITLPVYTLKERCLQVVRSLVKPENYRRLDIVRSLYEDLEDHPNVQKDLERLTQERIAHQRMGD
# MRPSGTAGAALLALLAALCPASRALEEKKVCQGTSNKLTQLGTFEDHFLSLQRMFNNCEVVLGNLEITYVQRNYDLSFLKTIQEVAGYVLIALNTVERIPLENLQIIRGNMYYENSYALAVLSNYDANKTGLKELPMRNLQEILHGAVRFSNNPALCNVESIQWRDIVSSDFLSNMSMDFQNHLGSCQKCDPSCPNGSCWGAGEENCQKLTKIICAQQCSGRCRGKSPSDCCHNQCAAGCTGPRESDCLVCRKFRDEATCKDTCPPLMLYNPTTYQMDVNPEGKYSFGATCVKKCPRNYVVTDHGSCVRACGADSYEMEEDGVRKCKKCEGPCRKVCNGIGIGEFKDSLSINATNIKHFKNCTSISGDLHILPVAFRGDSFTHTPPLDPQELDILKTVKEITGFLLIQAWPENRTDLHAFENLEIIRGRTKQHGQFSLAVVSLNITSLGLRSLKEISDGDVIISGNKNLCYANTINWKKLFGTSGQKTKIISNRGENSCKATGQVCHALCSPEGCWGPEPRDCVSCRNVSRGRECVDKCNLLEGEPREFVENSECIQCHPECLPQAMNITCTGRGPDNCIQCAHYIDGPHCVKTCPAGVMGENNTLVWKYADAGHVCHLCHPNCTYGCTGPGLEGCPTNGPKIPSIATGMVGALLLLLVVALGIGLFMRRRHIVRKRTLRRLLQERELVEPLTPSGEAPNQALLRILKETEFKKIKVLGSGAFGTVYKGLWIPEGEKVKIPVAIKELREATSPKANKEILDEAYVMASVDNPHVCRLLGICLTSTVQLITQLMPFGCLLDYVREHKDNIGSQYLLNWCVQIAKGMNYLEDRRLVHRDLAARNVLVKTPQHVKITDFGLAKLLGAEEKEYHAEGGKVPIKWMALESILHRIYTHQSDVWSYGVTVWELMTFGSKPYDGIPASEISSILEKGERLPQPPICTIDVYMIMVKCWMIDADSRPKFRELIIEFSKMARDPQRYLVIQGDERMHLPSPTDSNFYRALMDEEDMDDVVDADEYLIPQQGFFSSPSTSRTPLLSSLSATSNNSTVACIDRNGLQSCPIKEDSFLQRYSSDPTGALTEDSIDDTFLPVPEYINQSVPKRPAGSVQNPVYHNQPLNPAPSRDPHYQDPHSTAVGNPEYLNTVQPTCVNSTFDSPAHWAQKGSHQISLDNPDYQQDFFPKEAKPNGIFKGSTAENAEYLRVAPQSSEFIGA
# Mutation : L858R
# PDB :6lud

# Other variant
# Assay : CHEMBL5160887
# L858R,T790M
# MRPSGTAGAALLALLAALCPASRALEEKKVCQGTSNKLTQLGTFEDHFLSLQRMFNNCEVVLGNLEITYVQRNYDLSFLKTIQEVAGYVLIALNTVERIPLENLQIIRGNMYYENSYALAVLSNYDANKTGLKELPMRNLQEILHGAVRFSNNPALCNVESIQWRDIVSSDFLSNMSMDFQNHLGSCQKCDPSCPNGSCWGAGEENCQKLTKIICAQQCSGRCRGKSPSDCCHNQCAAGCTGPRESDCLVCRKFRDEATCKDTCPPLMLYNPTTYQMDVNPEGKYSFGATCVKKCPRNYVVTDHGSCVRACGADSYEMEEDGVRKCKKCEGPCRKVCNGIGIGEFKDSLSINATNIKHFKNCTSISGDLHILPVAFRGDSFTHTPPLDPQELDILKTVKEITGFLLIQAWPENRTDLHAFENLEIIRGRTKQHGQFSLAVVSLNITSLGLRSLKEISDGDVIISGNKNLCYANTINWKKLFGTSGQKTKIISNRGENSCKATGQVCHALCSPEGCWGPEPRDCVSCRNVSRGRECVDKCNLLEGEPREFVENSECIQCHPECLPQAMNITCTGRGPDNCIQCAHYIDGPHCVKTCPAGVMGENNTLVWKYADAGHVCHLCHPNCTYGCTGPGLEGCPTNGPKIPSIATGMVGALLLLLVVALGIGLFMRRRHIVRKRTLRRLLQERELVEPLTPSGEAPNQALLRILKETEFKKIKVLGSGAFGTVYKGLWIPEGEKVKIPVAIKELREATSPKANKEILDEAYVMASVDNPHVCRLLGICLTSTVQLIMQLMPFGCLLDYVREHKDNIGSQYLLNWCVQIAKGMNYLEDRRLVHRDLAARNVLVKTPQHVKITDFGRAKLLGAEEKEYHAEGGKVPIKWMALESILHRIYTHQSDVWSYGVTVWELMTFGSKPYDGIPASEISSILEKGERLPQPPICTIDVYMIMVKCWMIDADSRPKFRELIIEFSKMARDPQRYLVIQGDERMHLPSPTDSNFYRALMDEEDMDDVVDADEYLIPQQGFFSSPSTSRTPLLSSLSATSNNSTVACIDRNGLQSCPIKEDSFLQRYSSDPTGALTEDSIDDTFLPVPEYINQSVPKRPAGSVQNPVYHNQPLNPAPSRDPHYQDPHSTAVGNPEYLNTVQPTCVNSTFDSPAHWAQKGSHQISLDNPDYQQDFFPKEAKPNGIFKGSTAENAEYLRVAPQSSEFIGA


# PDB 2qlq
# Engineered mutation S->C
# Assay CHEMBL978650
# S345C
# Variant
# MGSSKSKPKDPSQRRRSLEPPDSTHHGGFPASQTPNKTAAPDTHRTPSRSFGTVATEPKLFGGFNTSDTVTSPQRAGALAGGVTTFVALYDYESRTETDLSFKKGERLQIVNNTEGDWWLAHSLTTGQTGYIPSNYVAPSDSIQAEEWYFGKITRRESERLLLNPENPRGTFLVRESETTKGAYCLSVSDFDNAKGLNVKHYKIRKLDSGGFYITSRTQFSSLQQLVAYYSKHADGLCHRLTNVCPTSKPQTQGLAKDAWEIPRESLRLEVKLGQGCFGEVWMGTWNGTTRVAIKTLKPGTMSPEAFLQEAQVMKKLRHEKLVQLYAVVSEEPIYIVTEYMSKGSLLDFLKGEMGKYLRLPQLVDMAAQIASGMAYVERMNYVHRDLRAANILVGENLVCKVADFGLARLIEDNEYTARQGAKFPIKWTAPEAALYGRFTIKSDVWSFGILLTELTTKGRVPYPGMVNREVLDQVERGYRMPCPPECPESLHDLMCQCWRKDPEERPTFEYLQAFLEDYFTSTEPQYQPGENL

# Component
# MGSSKSKPKDPSQRRRSLEPPDSTHHGGFPASQTPNKTAAPDTHRTPSRSFGTVATEPKLFGGFNTSDTVTSPQRAGALAGGVTTFVALYDYESRTETDLSFKKGERLQIVNNTEGDWWLAHSLTTGQTGYIPSNYVAPSDSIQAEEWYFGKITRRESERLLLNPENPRGTFLVRESETTKGAYCLSVSDFDNAKGLNVKHYKIRKLDSGGFYITSRTQFSSLQQLVAYYSKHADGLCHRLTNVCPTSKPQTQGLAKDAWEIPRESLRLEVKLGQGCFGEVWMGTWNGTTRVAIKTLKPGTMSPEAFLQEAQVMKKLRHEKLVQLYAVVSEEPIYIVMEYMSKGCLLDFLKGEMGKYLRLPQLVDMAAQIASGMAYVERMNYVHRDLRAANILVGENLVCKVADFGLARLIEDNEYTARQGAKFPIKWTAPEAALYGRFTIKSDVWSFGILLTELTTKGRVPYPGMVNREVLDQVERGYRMPCPPECPESLHDLMCQCWRKDPEERPTFEYLQAFLEDYFTSTEPQYQPGENL

# Other assay
# Assay CHEMBL978651
# S345C,T338M
# Variant
# MGSSKSKPKDPSQRRRSLEPPDSTHHGGFPASQTPNKTAAPDTHRTPSRSFGTVATEPKLFGGFNTSDTVTSPQRAGALAGGVTTFVALYDYESRTETDLSFKKGERLQIVNNTEGDWWLAHSLTTGQTGYIPSNYVAPSDSIQAEEWYFGKITRRESERLLLNPENPRGTFLVRESETTKGAYCLSVSDFDNAKGLNVKHYKIRKLDSGGFYITSRTQFSSLQQLVAYYSKHADGLCHRLTNVCPTSKPQTQGLAKDAWEIPRESLRLEVKLGQGCFGEVWMGTWNGTTRVAIKTLKPGTMSPEAFLQEAQVMKKLRHEKLVQLYAVVSEEPIYIVTEYMSKGSLLDFLKGEMGKYLRLPQLVDMAAQIASGMAYVERMNYVHRDLRAANILVGENLVCKVADFGLARLIEDNEYTARQGAKFPIKWTAPEAALYGRFTIKSDVWSFGILLTELTTKGRVPYPGMVNREVLDQVERGYRMPCPPECPESLHDLMCQCWRKDPEERPTFEYLQAFLEDYFTSTEPQYQPGENL


# 8 assays, 5 haben eine variant id, 2 haben die richtige mutation
# CHEMBL1108061
# CHEMBL1108057
# CHEMBL2217258
# CHEMBL2217259

# CHEMBL978651 S345C,T338M
# (CHEMBL997160) T338M (ligand nicht drin)

# CHEMBL978650 100% S345C
# CHEMBL4024872 100% S345C
# CHEMBL4024870 100% S345C


# 3f3w
# Assays mit variant {'CHEMBL4024870', 'CHEMBL4024872', 'CHEMBL997160', 'CHEMBL978650', 'CHEMBL978651'}
# Assays mit 100%  CHEMBL997159 CHEMBL997160 CHEMBL997163

# 2vuk
# assays mit mutation ['CHEMBL5107906', 'CHEMBL2184458', 'CHEMBL2184459', 'CHEMBL2183185', 'CHEMBL4154033', 'CHEMBL3796661']

# Alle assays :
# CHEMBL2183185
# CHEMBL2184458
# CHEMBL2184459
# CHEMBL3796661
# CHEMBL4154033
# CHEMBL5107906

# Leider alle assays haben ne variant sequence


# 3qrj
# Assay CHEMBL5108948
# Mutation T315I
# Mutation im PDB ist dieselbe
# Andere assays sind CHEMBL3991603  CHEMBL3991604 CHEMBL5060457 (CHEMBL5108948) CHEMBL5108949 CHEMBL5133794
# seq id nur mit variant seq 100 % oder mit seqres/seqadv
# Biochemgraph annotiert zwar nicht die assays selbst, aber deren description und er hat sogar sowohl assays mit als auch ohne

# Component
# MLEICLKLVGCKSKKGLSSSSSCYLEEALQRPVASDFEPQGLSEAARWNSKENLLAGPSENDPNLFVALYDFVASGDNTLSITKGEKLRVLGYNHNGEWCEAQTKNGQGWVPSNYITPVNSLEKHSWYHGPVSRNAAEYLLSSGINGSFLVRESESSPGQRSISLRYEGRVYHYRINTASDGKLYVSSESRFNTLAELVHHHSTVADGLITTLHYPAPKRNKPTVYGVSPNYDKWEMERTDITMKHKLGGGQYGEVYEGVWKKYSLTVAVKTLKEDTMEVEEFLKEAAVMKEIKHPNLVQLLGVCTREPPFYIITEFMTYGNLLDYLRECNRQEVNAVVLLYMATQISSAMEYLEKKNFIHRDLAARNCLVGENHLVKVADFGLSRLMTGDTYTAHAGAKFPIKWTAPESLAYNKFSIKSDVWAFGVLLWEIATYGMSPYPGIDLSQVYELLEKDYRMERPEGCPEKVYELMRACWQWNPSDRPSFAEIHQAFETMFQESSISDEVEKELGKQGVRGAVSTLLQAPELPTKTRTSRRAAEHRDTTDVPEMPHSKGQGESDPLDHEPAVSPLLPRKERGPPEGGLNEDERLLPKDKKTNLFSALIKKKKKTAPTPPKRSSSFREMDGQPERRGAGEEEGRDISNGALAFTPLDTADPAKSPKPSNGAGVPNGALRESGGSGFRSPHLWKKSSTLTSSRLATGEEEGGGSSSKRFLRSCSASCVPHGAKDTEWRSVTLPRDLQSTGRQFDSSTFGGHKSEKPALPRKRAGENRSDQVTRGTVTPPPRLVKKNEEAADEVFKDIMESSPGSSPPNLTPKPLRRQVTVAPASGLPHKEEAGKGSALGTPAAAEPVTPTSKAGSGAPGGTSKGPAEESRVRRHKHSSESPGRDKGKLSRLKPAPPPPPAASAGKAGGKPSQSPSQEAAGEAVLGAKTKATSLVDAVNSDAAKPSQPGEGLKKPVLPATPKPQSAKPSGTPISPAPVPSTLPSASSALAGDQPSSTAFIPLISTRVSLRKTRQPPERIASGAITKGVVLDSTEALCLAISRNSEQMASHSAVLEAGKNLYTFCVSYVDSIQQMRNKFAFREAINKLENNLRELQICPATAGSGPAATQDFSKLLSSVKEISDIVQR
# Variant
# MLEICLKLVGCKSKKGLSSSSSCYLEEALQRPVASDFEPQGLSEAARWNSKENLLAGPSENDPNLFVALYDFVASGDNTLSITKGEKLRVLGYNHNGEWCEAQTKNGQGWVPSNYITPVNSLEKHSWYHGPVSRNAAEYLLSSGINGSFLVRESESSPGQRSISLRYEGRVYHYRINTASDGKLYVSSESRFNTLAELVHHHSTVADGLITTLHYPAPKRNKPTVYGVSPNYDKWEMERTDITMKHKLGGGQYGEVYEGVWKKYSLTVAVKTLKEDTMEVEEFLKEAAVMKEIKHPNLVQLLGVCTREPPFYIIIEFMTYGNLLDYLRECNRQEVNAVVLLYMATQISSAMEYLEKKNFIHRDLAARNCLVGENHLVKVADFGLSRLMTGDTYTAHAGAKFPIKWTAPESLAYNKFSIKSDVWAFGVLLWEIATYGMSPYPGIDLSQVYELLEKDYRMERPEGCPEKVYELMRACWQWNPSDRPSFAEIHQAFETMFQESSISDEVEKELGKQGVRGAVSTLLQAPELPTKTRTSRRAAEHRDTTDVPEMPHSKGQGESDPLDHEPAVSPLLPRKERGPPEGGLNEDERLLPKDKKTNLFSALIKKKKKTAPTPPKRSSSFREMDGQPERRGAGEEEGRDISNGALAFTPLDTADPAKSPKPSNGAGVPNGALRESGGSGFRSPHLWKKSSTLTSSRLATGEEEGGGSSSKRFLRSCSASCVPHGAKDTEWRSVTLPRDLQSTGRQFDSSTFGGHKSEKPALPRKRAGENRSDQVTRGTVTPPPRLVKKNEEAADEVFKDIMESSPGSSPPNLTPKPLRRQVTVAPASGLPHKEEAGKGSALGTPAAAEPVTPTSKAGSGAPGGTSKGPAEESRVRRHKHSSESPGRDKGKLSRLKPAPPPPPAASAGKAGGKPSQSPSQEAAGEAVLGAKTKATSLVDAVNSDAAKPSQPGEGLKKPVLPATPKPQSAKPSGTPISPAPVPSTLPSASSALAGDQPSSTAFIPLISTRVSLRKTRQPPERIASGAITKGVVLDSTEALCLAISRNSEQMASHSAVLEAGKNLYTFCVSYVDSIQQMRNKFAFREAINKLENNLRELQICPATAGSGPAATQDFSKLLSSVKEISDIVQR
