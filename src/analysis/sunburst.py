import json
import shutil
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go
from matplotlib import cm

from src.utils import definitions

my_viridis = [
    "#453581",
    "#3d4d8a",
    "#34618d",
    "#2b748e",
    "#24878e",
    "#1f998a",
    "#25ac82",
    "#40bd72",
    "#67cc5c",
    "#98d83e",
    "#cde11d",
    "#fde725",
]


def generate_lookup_numbers(lookup_path):
    """
    Generates a lookup class id -> parent id
    @param lookup_path: path to the lookup data
    @return:
    """
    lookup_table = {}
    data = pd.read_csv(lookup_path, dtype=str)
    for index, row in data.iterrows():
        class_id = int(float(row["protein_class_id"]))
        raw_parent_class_id = row["parent_id"]
        if pd.isna(raw_parent_class_id):
            parent_class_id = ""
        else:
            parent_class_id = int(float(raw_parent_class_id))
        lookup_table[class_id] = parent_class_id
    return lookup_table


def generate_lookup_names(lookup_path, short_name):
    """
    Generates lookup id -> name
    @param lookup_path: path to the lookup data
    @param short_name:
    @return:
    """
    lookup_table = {}
    data = pd.read_csv(lookup_path, dtype=str)
    for index, row in data.iterrows():
        class_id = int(float(row["protein_class_id"]))
        if short_name:
            name = row["short_name"]
        else:
            name = row["pref_name"]
        lookup_table[class_id] = name
    lookup_table[""] = ""
    # I have to override some names because they are identical to their parents without the subgroup
    # I was able to drop this, but had to adjust it in the raw data which is messy...
    # Todo: automatically handle this case
    # lookup_table[571] = "Protease-activated receptor subgroup"
    # lookup_table[407] = "Parathyroid hormone receptor subgroup"
    # lookup_table[503] = "Corticotropin releasing factor receptor subgroup"
    return lookup_table


def generate_only_name_lookup_table(lookup_path, short_name):
    """
    Generates lookup name -> parent name
    @param lookup_path: path to the lookup data
    @param short_name: Flag what name to use
    @return:
    """
    id_name_lookup = generate_lookup_names(lookup_path, short_name)
    lookup_table = {}
    data = pd.read_csv(lookup_path, dtype=str)
    for index, row in data.iterrows():
        class_id = int(float(row["protein_class_id"]))
        if pd.isna(row["parent_id"]):
            parent_name = ""
        else:
            parent_name = id_name_lookup[int(float(row["parent_id"]))]
        if short_name:
            name = row["short_name"]
        else:
            name = row["pref_name"]
        lookup_table[name] = parent_name
    lookup_table[""] = ""
    return lookup_table


def extract_data_pdb_ids_nonredundant(
    raw_data_path, lookup_path, short_name, relevant_pdbs
):
    """
    This function extracts the raw data from the csv file
    This function is way too long, I should rewrite it and split it up....
    @param raw_data_path: Path to the csv file
    @param lookup_path: Path to the chembl lookup
    @param short_name: Flag if we want to extract the preferred or short name
    @return:
    """
    pdb_ids = []
    new_chembl_target_dict = {}
    name_lookup_table = generate_lookup_names(lookup_path, short_name)
    used_target_names = []
    target_id_lookups = {}
    labels = {}
    # First run through the data, we extract the pdb, chembl targets and all parent class ids
    for chunk in pd.read_csv(raw_data_path, dtype=str, chunksize=10000, skiprows=1):
        for index, row in chunk.iterrows():
            pdb = row["pdb.pdb_id"]
            if relevant_pdbs and pdb not in relevant_pdbs:
                continue
            target_id = row["target.chembl_id"]
            normal_name = row["target.preferred_name"] + "ctarget"
            if target_id in target_id_lookups:
                target_name = target_id_lookups[target_id]
            else:
                while normal_name in used_target_names:
                    normal_name = normal_name + " "
                target_name = normal_name
                target_id_lookups[target_id] = normal_name
                used_target_names.append(normal_name)
            pdb_target = tuple([pdb, target_name])
            class_level = row["Class_Level"]
            class_id = row["Protein_Class_ID"]
            if pdb_target in pdb_ids or pd.isna(class_id):
                continue
            split_class_ids = list(set([int(float(x)) for x in class_id.split(";")]))
            if not target_name in new_chembl_target_dict.keys():
                new_chembl_target_dict[target_name] = [
                    class_level,
                    [pdb],
                    [name_lookup_table[x] for x in split_class_ids],
                ]
            else:
                if new_chembl_target_dict[target_name][0] != class_level:
                    # print(class_level,pdb, class_id, new_chembl_target_dict[target_name], target_name, target_id)
                    raise RuntimeError(
                        "Class level stimmt nicht überein, irgendwas ist hier komisch..."
                    )
                if not pdb in new_chembl_target_dict[target_name][1]:
                    new_chembl_target_dict[target_name][1].append(pdb)
                to_extend = [
                    name_lookup_table[x]
                    for x in split_class_ids
                    if not name_lookup_table[x]
                    in new_chembl_target_dict[target_name][2]
                ]
                new_chembl_target_dict[target_name][2].extend(to_extend)
                pdb_ids.append(pdb_target)

    # Second run through the data, now we go through the extracted data and rewrite it using all pdbs
    # In addition, here the magic with the whitespaces happens to be able to plot it in the end
    # If two names are identical, the sunburst chart dies. But as one PDB can have multiple targets and
    # one target can have multiple parents we have to enforce this so every single PDB and target is unique
    # We do this using whitespaces before and after the name
    # An improvement on this would be using specific and unique ids and then use the names as labels
    all_parents = []
    used_pdbs = {}
    final_data = []
    counter = 1000
    for target in new_chembl_target_dict:
        parents = new_chembl_target_dict[target][2]
        pdbs = new_chembl_target_dict[target][1]
        for pdb in pdbs:
            if not pdb in used_pdbs:
                used_pdbs[pdb] = 0
            else:
                used_pdbs[pdb] += 1
        for i, parent in enumerate(parents):
            if not parent in all_parents:
                all_parents.append(parent)
            current_label_target = "torben" + str(counter)
            labels[current_label_target] = target.strip()
            counter += 1
            final_data.append([current_label_target, parent, 0])
            for pdb in pdbs:
                current_label_pdb = "torben" + str(counter)
                labels[current_label_pdb] = pdb.strip()
                counter += 1
                final_data.append([current_label_pdb, current_label_target, 1])

    # Now we have all PDB and chembl targets in the data, we need to add the complete parent tree to the data
    grandparent_dict = {}
    only_name_lookup = generate_only_name_lookup_table(lookup_path, short_name)
    for parent in all_parents:
        current_parent = parent
        while current_parent != "":
            grandparent = only_name_lookup[current_parent]
            key = tuple([current_parent, grandparent])
            if key in grandparent_dict.keys():
                break
            else:
                grandparent_dict[key] = 0
                current_parent = grandparent

    # Last thing is to slightly rewrite the data and we are finished
    for key in grandparent_dict.keys():
        final_data.append([key[0], key[1], 0])
        if not key[0] in labels.keys():
            labels[key[0]] = key[0]
    final_data = pd.DataFrame(final_data, columns=["Name", "Parent", "Number"])
    return final_data, labels


def finish_up_sunburst_data(
    sunburst_data, labels, remove_targets=False, remove_pdbs=True
):
    """
    This function alters the sunburst data so it is ready to be plotted and removes unwanted data
    If you remove targets but not PDBS the plot will not work anymore
    @param sunburst_data: Raw sunburst data
    @param remove_targets: Flag to remove targets
    @param remove_pdbs:  Flag to remove PDBS
    @return:
    """
    if not remove_pdbs and remove_targets:
        raise RuntimeError(
            "Removing targets without removing PDBS not possible, exiting..."
        )
    final_data = {}
    for index, row in sunburst_data.iterrows():
        name = row["Name"]
        number = row["Number"]
        parent = row["Parent"]
        final_data[name] = [parent, number]  # Generate lookup data

    for index, row in sunburst_data.iterrows():
        name = row["Name"]
        parent = row["Parent"]
        number = int(float(row["Number"]))
        current_parent = parent
        if number > 0:
            while not pd.isna(current_parent):
                # Now we go the tree until we reach the stem
                if current_parent == "":
                    break
                final_data[current_parent][
                    1
                ] += number  # Increase the parent by its child
                current_parent = final_data[current_parent][0]

    plotting_data = []
    for key in final_data:
        if len(str(labels[key]).strip()) == 4 and remove_pdbs:  # Remove PDB if wanted
            continue
        elif (
            str(labels[key]).strip().endswith(" ctarget") and remove_targets
        ):  # Remove target if wanted
            continue
        plotting_data.append(
            [key, final_data[key][0], final_data[key][1]]
        )  # Save data in necessary format
    plotting_data = pd.DataFrame(plotting_data, columns=["Name", "Parent", "Number"])
    return plotting_data


def get_bottom_level_labels(plotting_data):
    bottom_level_labels = []
    bottom_level_innerworkings = []
    for index, row in plotting_data.iterrows():
        parent = row["Parent"]
        name = row["Name"]
        if parent == "Protein Class":
            bottom_level_innerworkings.append(name)
    return bottom_level_innerworkings


def generate_sunburst_chart(
    raw_data_path,
    lookup_path,
    output_directory,
    short_name=True,
    remove_pdbs=True,
    remove_targets=False,
    restart=False,
    relevant_pdbs=[],
):
    """
    General function to generate a sunburst plot including all data
    @param raw_data_path: Path to the data to be analysed. Both parsed activity and complete master table works, also filtered and unfiltered both work
    @param lookup_path: Path to the chembl lookup table (is written during mt creation)
    @param output_directory: Directory to output the data and plots to
    @param short_name: If you want to use abbreviated or preferred chembl names
    @param remove_pdbs: Remove pdb file layer (necessary for big plots)
    @param remove_targets: Remove chembl target layer
    @param restart: Do not calculate the raw data again but only redo the plot
    @return:
    """
    if restart:
        if output_directory.exists():
            shutil.rmtree(output_directory)
        output_directory.mkdir()
        new_data, label_data = extract_data_pdb_ids_nonredundant(
            raw_data_path, lookup_path, short_name, relevant_pdbs
        )
        new_data.to_csv(output_directory / "sunburst_data_raw.csv")
        with open(output_directory / "labels.json", "w") as f:
            json.dump(label_data, f)
    else:
        new_data = pd.read_csv(output_directory / "sunburst_data_raw.csv", index_col=0)
        with open((output_directory / "labels.json"), "r") as f:
            label_data = json.load(f)
    keepers = []
    used_names = []
    for index, row in new_data.iterrows():
        name = row["Name"]
        parent = row["Parent"]
        if name == parent:
            raise RuntimeError(f"name {name} and parent {parent} identical")
        if name in used_names:
            raise RuntimeError(f"Name ,{name}, already used")
        used_names.append(name)
        if len(row["Name"].strip()) == 4:
            keepers.append(index)
    plotting_data = finish_up_sunburst_data(
        new_data, label_data, remove_pdbs=remove_pdbs, remove_targets=remove_targets
    )
    plotting_data.to_csv(output_directory / "sunburst_data_plot_rmsd.csv")
    names = plotting_data["Name"]
    parents = plotting_data["Parent"]
    numbers = plotting_data["Number"]
    labels_corr = [label_data[x].replace("ctarget", "") for x in names]
    fig = go.Figure(
        go.Sunburst(
            ids=names,
            labels=labels_corr,
            parents=parents,
            values=numbers,
            branchvalues="total",
            maxdepth=3,
            # domain=dict(column=0),
            # marker=dict(colorscale=px.colors.sequential.Viridis_r)
            # textfont=dict(size=18)  # << HIER Schriftgröße einstellen
        )
    )
    # fig.update_layout(uniformtext=dict(minsize=16, mode="hide"))
    offset = 1
    color_map = cm.get_cmap("viridis", 15 + offset)

    fig.update_layout(sunburstcolorway=my_viridis, font_family="Sans Serif")
    fig.update_layout(margin=dict(t=50, l=25, r=25, b=25))
    fig.write_html(output_directory / "sunburst_plot_pdb.html")
    fig.write_image(output_directory / "sunburst_plot_pdb.svg")


def execute_sunburst(restart=True):
    new_data_path = definitions.base_stractable_path / "filtered_activities.csv"
    lookup_path = "/scratch/gutermuth/sge_trialerror/tool-wrappers/output/protein_classification_lookup.csv"
    plot_output_path = definitions.output_path / Path("sunburst")
    plot_output_path.mkdir(exist_ok=True)
    ballrooma_data = pd.read_csv(definitions.ballrooma / "complete_data.csv")
    ballrooms_data = pd.read_csv(definitions.ballrooms / "complete_data.csv")
    ballrooma_pdbs = [x.split("_")[0] for x in ballrooma_data["Dockinto1"]]
    ballrooma_pdbs.extend([x.split("_")[0] for x in ballrooma_data["Dockinto1"]])
    ballrooma_pdbs = list(set(ballrooma_pdbs))
    ballrooms_pdbs = [x.split("_")[0] for x in ballrooms_data["Dockinto1"]]
    ballrooms_pdbs.extend([x.split("_")[0] for x in ballrooms_data["Dockinto1"]])
    ballrooms_pdbs = list(set(ballrooms_pdbs))
    ballrooma_output_path = definitions.output_path / Path("sunburst") / "ballrooma"
    ballrooma_output_path.mkdir(exist_ok=True)
    ballrooms_output_path = definitions.output_path / Path("sunburst") / "ballrooms"
    ballrooms_output_path.mkdir(exist_ok=True)

    generate_sunburst_chart(
        new_data_path,
        lookup_path,
        ballrooma_output_path,
        False,
        remove_pdbs=False,
        remove_targets=False,
        restart=restart,
        relevant_pdbs=ballrooma_pdbs,
    )
    generate_sunburst_chart(
        new_data_path,
        lookup_path,
        ballrooms_output_path,
        False,
        remove_pdbs=False,
        remove_targets=False,
        restart=restart,
        relevant_pdbs=ballrooms_pdbs,
    )


def main():
    execute_sunburst(restart=True)


if __name__ == "__main__":
    main()
