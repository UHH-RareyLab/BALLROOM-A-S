from PIL import ImageDraw
from rdkit import Chem
from rdkit.Chem import Draw, AllChem, rdFMCS
from rdkit.Chem.Draw import rdMolDraw2D

from src.utils import definitions


def align_and_visualise_mols_mcs(
    molecules_to_showcase, subimgsize=400, molperrow=3, output_name="test", align=True
):
    molecule_difference_path = definitions.output_path
    dopts = rdMolDraw2D.MolDrawOptions()
    dopts.addStereoAnnotation = True
    dopts.minFontSize = 45  # Ändert molekül aber nicht caption
    dopts.legendFontSize = 45
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
            drawOptions=dopts,
        )
        draw = ImageDraw.Draw(img)
        img.save(molecule_difference_path / f"paper_example_{output_name}")
    elif output_name.endswith(".svg"):
        img = Draw.MolsToGridImage(
            molecules_to_showcase,
            molsPerRow=molperrow,
            legends=[x.GetProp("_Name") for x in molecules_to_showcase],
            drawOptions=dopts,
            useSVG=True,
        )
        with open(molecule_difference_path / f"paper_example_{output_name}", "w") as f:
            f.write(img)


def create_and_plot_mol_array_align(
    molecule_names, names, output_name, align=True, size_molecule=800
):
    molecules = []
    if not len(molecule_names) == len(names):
        print(
            f"Number of molecules {len(molecule_names)} and number of names {len(names)}"
        )
        raise RuntimeError("Number of molecules uneuqal to number of names")
    for index, name in enumerate(molecule_names):
        smiles = name
        mol = Chem.MolFromSmiles(smiles)
        mol.SetProp("_Name", names[index])
        molecules.append(mol)
    align_and_visualise_mols_mcs(
        molecules,
        molperrow=len(molecules),
        output_name=output_name,
        align=align,
        subimgsize=size_molecule,
    )


size_molecule = 800

names = [
    "CHEMBL442450 / 6F5L_CQB_A_823 \t 0.65nM",
    "CHEMBL4513258 \t 8530nM",
]  # ASSAY CHEMBL4311581
mols = [
    "O=C(O[C@H](C(=O)O)CCC(=O)O)N[C@H](C(=O)O)CC(C)C",
    "CC(C)C[C@@H](NC(=O)O[C@H](CCC(=O)O)C(=O)O)C(=O)O",
]

create_and_plot_mol_array_align(
    mols,
    names,
    output_name="ballrooms_highdiff_pair.png",
    align=True,
    size_molecule=size_molecule,
)

names = ["CHEMBL4096813", "CHEMBL4087054"]  # ASSAY CHEMBL4037917
mols = ["O=C(NCC)C=1Nc2nccc(N3[C@H](COCC3)C)c2C1", "O=C(OCC)c1c2c(ncc1)NC=C2"]

create_and_plot_mol_array_align(
    mols,
    names,
    output_name="ballrooma_highdiff_pair_1.png",
    align=True,
    size_molecule=size_molecule,
)

names = [
    "CHEMBL5205374 / 6DY7_HH7_A_401 \t 66000nM",
    "CHEMBL5172311 / 6E22_HLS_A_401 \t 1.3nM",
]  # ASSAY CHEMBL5232467
mols = [
    "Clc1c(OCCCNC2=NCCN2)ccc(Cl)c1",
    "Fc1c(cc(cc1)C(=O)NCc2cc(OC)cc(OC)c2)CNC3=NCCN3",
]
create_and_plot_mol_array_align(
    mols,
    names,
    output_name="ballrooma_highdiff_pair_2.png",
    align=True,
    size_molecule=size_molecule,
)

# Ballroom A

# diff 3
# 3B28_B2X_A_237,3B27_B2T_A_1,Clc1c(c2nc(SC)nc(n2)N)c3c4c(c1)COCc4ccc3,Clc1c(c2nc(SC)nc(n2)N)cccc1,3.4,6900.0,3,Kd,nM,CHEMBL1838880
names = [
    "CHEMBL1834096 / 3B28_B2X_A_237 \t 3.4nM",
    "CHEMBL1834095 / 3B27_B2T_A_1 \t 6900nM",
]
mols = ["Clc1c(c2nc(SC)nc(n2)N)c3c4c(c1)COCc4ccc3", "Clc1c(c2nc(SC)nc(n2)N)cccc1"]
create_and_plot_mol_array_align(
    mols,
    names,
    output_name="ballrooma_diff_3_pair.png",
    align=True,
    size_molecule=size_molecule,
)

# diff 2
# 255,254,4LY1_20Y_B_408,5IX0_6EZ_A_404,S1C(c2cc(NC(=O)c3ccc(NC(=O)C)cc3)c(N)cc2)=CC=C1,Fc1ccc(c2cc(NC(=O)[C@@H]3C[C@@H]4O[C@H](C3)CC4)c(N)cc2)cc1,1.5,519.0,2,Ki,nM,CHEMBL3830979,CHEMBL1937,4LY1_20Y_B_408,5IX0_6EZ_A_404,16646511.0,16646509.0
names = [
    "CHEMBL235842 / 4LY1_20Y_B_408 \t 1.5nM",
    "CHEMBL3827611 / 5IX0_6EZ_A_404 \t 519nM",
]
mols = [
    "S1C(c2cc(NC(=O)c3ccc(NC(=O)C)cc3)c(N)cc2)=CC=C1",
    "Fc1ccc(c2cc(NC(=O)[C@@H]3C[C@@H]4O[C@H](C3)CC4)c(N)cc2)cc1",
]
create_and_plot_mol_array_align(
    mols,
    names,
    output_name="ballrooma_diff_2_pair.png",
    align=False,
    size_molecule=size_molecule,
)

# diff 1
# 164,163,8D6C_QGR_B_405,8D6E_QGI_B_402,Brc1cc2nc3c(nc2cc1)N(c4c(c(O)ccc4C)C)C(N)=C3C(=O)N,O=C(N)C=1c2c(nc(c(c2)C)C)N(c3c(c(O)ccc3C)C)C1N,360.0,14.0,1,IC50,nM,CHEMBL5141607,CHEMBL3984,8D6C_QGR_B_405,8D6E_QGI_B_402,24826080.0,24826105.0

names = [
    "CHEMBL5185146 / 8D6C_QGR_B_405 \t 360nM",
    "CHEMBL5199076 / 8D6E_QGI_B_402 \t 14nM",
]
mols = [
    "Brc1cc2nc3c(nc2cc1)N(c4c(c(O)ccc4C)C)C(N)=C3C(=O)N",
    "O=C(N)C=1c2c(nc(c(c2)C)C)N(c3c(c(O)ccc3C)C)C1N",
]
create_and_plot_mol_array_align(
    mols,
    names,
    output_name="ballrooma_diff_1_pair.png",
    align=False,
    size_molecule=size_molecule,
)

# Ballroom S

# diff 3
# 649,648,CHEMBL4461007,CHEMBL4573978,CS(=O)(=O)N1CCC(NC(=O)c2cc(C3CC3)on2)CC1,N[C@H]1CC[C@H](S(=O)(=O)N2CCC(NC(=O)c3cc(C4CC4)on3)CC2)CC1,5010.0,1.6,3,IC50,nM,CHEMBL4338950,CHEMBL2321643,6P6K_L0J_A_508;6PAF_O6A_A_509;,6P6K_L0J_A_508;6PAF_O6A_A_509;,18979782.0,18979791.0

names = ["CHEMBL4461007 \t 5010nM", "CHEMBL4573978 \t 1.6nM"]
mols = [
    "CS(=O)(=O)N1CCC(NC(=O)c2cc(C3CC3)on2)CC1",
    "N[C@H]1CC[C@H](S(=O)(=O)N2CCC(NC(=O)c3cc(C4CC4)on3)CC2)CC1",
]
create_and_plot_mol_array_align(
    mols,
    names,
    output_name="ballrooms_diff_3_pair.png",
    align=True,
    size_molecule=size_molecule,
)

# diff 2
# 1124,1123,CHEMBL2420638,CHEMBL2420654,O=C(NCc1ccc(S(=O)(=O)c2ccccc2)cc1)c1cnc2nccn2c1,Cn1c(C(=O)NCc2ccc(S(=O)(=O)c3ccccc3)cc2)cc2ccncc21,2.0,410.0,2,IC50,nM,CHEMBL2423680,CHEMBL1744525,4LWW_LWW_A_601;,4LWW_LWW_A_601;,13428905.0,13428924.0
names = ["CHEMBL2420638 \t 2nM", "CHEMBL2420654 \t 410nM"]
mols = [
    "O=C(NCc1ccc(S(=O)(=O)c2ccccc2)cc1)c1cnc2nccn2c1",
    "Cn1c(C(=O)NCc2ccc(S(=O)(=O)c3ccccc3)cc2)cc2ccncc21",
]
create_and_plot_mol_array_align(
    mols,
    names,
    output_name="ballrooms_diff_2_pair.png",
    align=True,
    size_molecule=size_molecule,
)

# diff 1
# 253,252,CHEMBL4285444,CHEMBL4283434,Cn1c(=O)n(CC2CC2)c(=O)c2cc(S(=O)(=O)NC3(C)CC3)ccc21,Cn1c(=O)n(Cc2cnn(CC#N)c2)c(=O)c2cc(S(=O)(=O)NC3(C)CC3)ccc21,4600.0,370.0,1,EC50,nM,CHEMBL4265322,CHEMBL1795143,6HMK_7JC_A_1008;,6HMK_7JC_A_1008;,18756760.0,18756796.0
names = ["CHEMBL4285444 \t 4600nM", "CHEMBL4283434 \t 370nM"]
mols = [
    "Cn1c(=O)n(CC2CC2)c(=O)c2cc(S(=O)(=O)NC3(C)CC3)ccc21",
    "Cn1c(=O)n(Cc2cnn(CC#N)c2)c(=O)c2cc(S(=O)(=O)NC3(C)CC3)ccc21",
]
create_and_plot_mol_array_align(
    mols,
    names,
    output_name="ballrooms_diff_1_pair.png",
    align=True,
    size_molecule=size_molecule,
)
