from pathlib import Path

import pandas as pd

from calc_activity_pairs import calc_difference_between_two_activities, ActivityPair
from check_pair_benchmark_contradictions import get_noncontradicting_set
from src.utils import definitions, db_interaction

min_difference = 1


def calc_pairs_in_assay(assay_data, assay, min_difference=1):
    new_pairs = []
    for index1, row1 in assay_data.iterrows():
        for index2, row2 in assay_data.iterrows():
            if index1 >= index2:
                continue
            affinity1 = row1["standard_value"]
            affinity2 = row2["standard_value"]
            if (
                pd.isna(affinity1)
                or pd.isna(affinity2)
                or affinity1 <= 0
                or affinity2 <= 0
            ):
                continue
            identifier1 = row1["mol_chembl"]
            identifier2 = row2["mol_chembl"]
            smiles1 = row1["canonical_smiles"]
            smiles2 = row2["canonical_smiles"]
            activity_type = row1["standard_type"]
            activity_unit = row1["standard_units"]
            activity_id1 = row1["activity_id"]
            activity_id2 = row2["activity_id"]
            difference = calc_difference_between_two_activities(affinity1, affinity2)
            if difference >= min_difference:
                new_pairs.append(
                    ActivityPair(
                        identifier1,
                        identifier2,
                        smiles1,
                        smiles2,
                        affinity1,
                        affinity2,
                        difference,
                        assay,
                        target,
                        activity_type,
                        activity_unit,
                        identifier1,
                        identifier2,
                        activity_id1,
                        activity_id2,
                    )
                )
    return new_pairs


def get_activities_for_assay(assay):
    sql_query = f"""SELECT
    activities.standard_value,
    activities.activity_id,
    activities.standard_type,
    activities.standard_units,
    molecule_dictionary.chembl_id AS mol_chembl,
    compound_structures.canonical_smiles
    FROM activities
    JOIN assays ON activities.assay_id = assays.assay_id
    JOIN molecule_dictionary ON activities.molregno = molecule_dictionary.molregno
    JOIN compound_structures ON activities.molregno = compound_structures.molregno
    WHERE assays.chembl_id = '{assay}'
    AND activities.standard_relation = '='
    AND activities.standard_type IN ('IC50', 'EC50', 'Ki', 'Kd')
    AND (activities.data_validity_comment IS NULL OR activities.data_validity_comment IN ('Manually validated','Potential missing data','Potential author error','Potential transcription error'))"""
    return db_interaction.query_chembl(sql_query)


output_path = Path("/gutermuth/output/chembl_contradiction_analysis")
# First step get all targets/assays combinations from chembl
assay_target_data = db_interaction.query_chembl(
    "select assays.chembl_id as assay_id, target_dictionary.chembl_id as target_id from assays join target_dictionary on assays.tid = target_dictionary.tid"
)
target_assay_dict = {
    target: list(group["assay_id"])
    for target, group in assay_target_data.groupby("target_id")
}
# Second step, iterate over all targets and gather assay data
target_contradiction_data = []
for target in target_assay_dict:
    all_assays = target_assay_dict[target]
    print(all_assays)
    if len(all_assays) < 2:
        continue
    all_pairs = []
    for assay in all_assays:
        current_assay_data = get_activities_for_assay(assay)
        # Third step, create pairs for all assays and initiatie matrix
        assay_pairs = calc_pairs_in_assay(
            current_assay_data, assay, min_difference=min_difference
        )
        all_pairs.extend(assay_pairs)
    data = [x.return_pair_as_list() for x in all_pairs]
    df = pd.DataFrame(
        data,
        columns=definitions.pair_data_columns,
    )
    assays_in_pairs = set(x for x in df["Assay"])
    print(assays_in_pairs)
    if len(all_pairs) < 2 or len(assays_in_pairs) < 2:
        print("target skipped")
        continue
        # Fourth step, get contradiction measurements
    (
        contradicting_assays,
        _,
        _,
        _,
        assay_labels,
        contradiction_matrix,
    ) = get_noncontradicting_set(df, target)
    contradiction_matrix = pd.DataFrame(
        contradiction_matrix, index=assay_labels, columns=assay_labels
    )
    target_path = output_path / "specific_target_data" / f"{target}"
    target_path.mkdir(parents=True, exist_ok=True)
    target_contradiction_matrix_path = (
        target_path / f"contradiction_matrix_difference_{min_difference}.csv"
    )
    assays_to_discard_path = (
        target_path / f"assays_to_discard_{target}_difference_{min_difference}.csv"
    )
    df.to_csv(target_path / "pair_data.csv")
    contradiction_matrix.to_csv(target_contradiction_matrix_path)
    if len(contradicting_assays) > 0:
        with open(assays_to_discard_path, "w") as f:
            for assay in contradicting_assays:
                f.write(f"{assay}\n")
    target_contradiction_data.append(
        [
            target,
            len(contradicting_assays),
            len(assay_labels),
            contradicting_assays,
            assay_labels,
        ]
    )
    print(target, len(target_contradiction_data))

target_contradiction_data = pd.DataFrame(
    target_contradiction_data,
    columns=[
        "Target",
        "N contradictions",
        "N total",
        "Assays to discard",
        "All assays",
    ],
)
target_contradiction_data.to_csv(
    output_path / f"target_contradiction_data_difference_{min_difference}.csv"
)
