import math
import random

import numpy as np
import pandas as pd

from src.generation.retrieve_bioactivity_info_from_chembl import (
    check_validity_comment,
)
from src.utils.definitions import (
    allowed_activities,
    matchings_proteins,
    matchings_molecules,
)


def orderOfMagnitude(number):
    print(number)
    return math.floor(math.log(number, 10))


class DatasetFilter:
    """
    This DatasetFilter shall solve multiple problems that can be present in an assay after filtering
    Problems to fix :
    1. Multiple activity types per assay (WHY THE FUCK IS THIS A THING)
    2. Multiple identical molecules per pdb
    3. Multiple activities per molecule in a pdb (unique molecule)
    4. Multiple molecules per activity
    5. Check that now no molecules are doubled and no activities are doubled
    After checking for that we should only have unique pdbs, unique molecules and unique activities per assay
    """

    def __init__(self):
        random.seed(1994)

    def run_prefiltering(self, dataframe, min_datapoints, filter_uncertain=True):
        if filter_uncertain:
            first_step = self.filter_uncertain_datapoints(dataframe)
            print(
                f"Filtered {dataframe.shape[0] - first_step.shape[0]} datapoints from uncertainty, "
                f"old {dataframe.shape[0]}, new {first_step.shape[0]}, "
                f"assay filtered {first_step.shape[0] < min_datapoints and dataframe.shape[0] > min_datapoints} "
                f"assay {dataframe.iloc[0]['assay.chembl_id']}"
            )
        else:
            first_step = dataframe
        step2 = self.make_activity_types_unique(first_step)
        step3 = self.decide_multiple_activity_values_per_ligand(step2)
        if step3.shape[0] >= min_datapoints:
            return step3
        else:
            raise RuntimeError("Dataframe too small even in prefiltering")

    def filter_uncertain_datapoints(self, dataframe):
        keepers = []
        for index, row in dataframe.iterrows():
            comment = row["chembl_activity.data_validity_comment"]
            potential_duplicate = int(float(row["chembl_activity.potential_duplicate"]))
            if check_validity_comment(comment) and potential_duplicate == 0:
                keepers.append(index)
            else:
                print("Filtered one for uncertainty")
        return dataframe.loc[keepers]

    def run(
        self, dataframe: pd.DataFrame, min_datapoints, min_difference
    ) -> pd.DataFrame:
        """
        This function executes the filter cascade on the dataframe
        @param dataframe: The dataframe that shall be filtered
        @return: filtered dataframe
        """
        # Remove data not part of the most prevalent activity type
        step1 = self.make_activity_types_unique(
            dataframe
        )  # Das sollte ich auf jeden Fall machen
        # Remove multiple identical ligands in a pdb file
        step2 = self.decide_multiple_ligands_per_pdb(
            step1
        )  # Das kann ich mit Siena auflösen und ansonsten machen
        # Remove multiple Activity values per ligand (e.g. R/S enantiomer)
        step3 = self.decide_multiple_activity_values_per_ligand(
            step2
        )  # Das sollte ich auf jeden Fall machen
        # Remove multiple Ligands per Activity (e.g. two identical ligands in multiple PDBs)
        step4 = self.decide_multiple_rows(
            step3, "activity.chembl_activity_id"
        )  # Das könnte ich unter umständen auch mit SIENA auflösen
        # Remove multiple ligands with the same USMILES ? Shouldnt this already be resolved with activities? different protonations?
        step5 = self.decide_multiple_rows(step4, "usmiles")  # Same
        worked = self.check_if_dataframe_is_ok(step5, min_datapoints, min_difference)
        if worked:
            return step5
        else:
            raise RuntimeError("Happy little error")

    def make_activity_types_unique(self, dataframe: pd.DataFrame) -> pd.DataFrame:
        """
        Filters the dataframe for the activity type with the most unique ligands
        In most cases, only one activity type should be in an assay...
        If not, the one with most unique ligands is used, if it is a draw the allowed activities are
        looped through in the order of importance and the first one with most unique ligands is chosen
        @param dataframe: dataframe to be filtered
        @return: filtered dataframe
        """
        activity_type_map = {}
        for index, row in dataframe.iterrows():
            activity_type_row = row["activity.standard_type"]
            usmiles = row["usmiles"]
            if activity_type_row in activity_type_map.keys():
                activity_type_map[activity_type_row].add(usmiles)
            else:
                activity_type_map[activity_type_row] = {usmiles}
        if len(activity_type_map.keys()) == 1:
            # Everything is as expected
            return dataframe
        else:
            # More than one activity type present
            found = False
            longest_key = ""
            length = 0
            for key in activity_type_map.keys():
                if len(activity_type_map[key]) > length:
                    # if activity type with more ligands is found use it
                    length = len(activity_type_map[key])
                    longest_key = [key]
                    found = True
                elif len(activity_type_map[key]) == length:
                    # if more than one activity type with most ligands is found
                    found = False
                    longest_key.append(key)
            if found:
                # Good solution found
                dataframe = dataframe[
                    dataframe["activity.standard_type"] == longest_key[0]
                ]
                return dataframe
            else:
                # No good solution found, use the first one in allowed activities... avoids random
                for gen_type in allowed_activities:
                    if gen_type in longest_key:
                        return dataframe[
                            dataframe["activity.standard_type"] == gen_type
                        ]
                raise RuntimeError

    def decide_multiple_ligands_per_pdb(self, dataframe: pd.DataFrame) -> pd.DataFrame:
        """
        Filters the dataframe for multiple ligands in the same PDB File
        @param dataframe: dataframe to be filtered
        @return: filtered dataframe
        """
        pdb_usmiles_map = {}
        all_indexes = []
        indexes_to_drop = []
        # Calc all pdb usmiles combos
        for index, row in dataframe.iterrows():
            all_indexes.append(index)
            usmiles = row["usmiles"]
            pdb = row["pdb"]
            name = row["name"]
            key = (pdb, usmiles)
            if not key in pdb_usmiles_map.keys():
                pdb_usmiles_map[key] = {name: [index]}
            else:
                if not name in pdb_usmiles_map[key].keys():
                    pdb_usmiles_map[key][name] = [index]
                else:
                    pdb_usmiles_map[key][name].append(index)
        for key in pdb_usmiles_map.keys():
            all_names = pdb_usmiles_map[key].keys()
            # If there is more than one name then more than one row with the same pdb/usmiles is present
            if len(all_names) > 1:
                new_indexes = []
                for name in all_names:
                    new_indexes.extend([index for index in pdb_usmiles_map[key][name]])
                # Create dataset with all non unique pdb/usmiles
                new_data = dataframe.loc[new_indexes]
                best_edia = 0
                best_edia_name = ""
                found = False
                for index, row in new_data.iterrows():
                    edia = float(row["EDIAm"])
                    name = row["name"]
                    if name in best_edia_name:
                        continue
                    if edia > best_edia:
                        # Found a solution using edia
                        best_edia = edia
                        best_edia_name = [name]
                        found = True
                    elif edia == best_edia:
                        # Mhh, didint work in the end
                        found = False
                        best_edia_name.append(name)
                if found:
                    # if found use all molecules with the best edia
                    blub = list(set(pdb_usmiles_map[key].keys()) - set(best_edia_name))
                else:
                    # else decide between the ones with the same edia using random
                    choice = random.choice(best_edia_name)
                    blub = list(set(pdb_usmiles_map[key].keys()) - {choice})
                discard = [pdb_usmiles_map[key][test] for test in blub]
                indexes_to_drop.extend([x[0] for x in discard])
        return dataframe.loc[list(set(all_indexes) - set(indexes_to_drop))]

    def decide_multiple_activity_values_per_ligand(
        self, dataframe: pd.DataFrame
    ) -> pd.DataFrame:
        """
        Filters the dataframe for multiple activity values per ligand
        @param dataframe: dataframe to be filtered
        @return: filtered dataframe
        """
        pdb_usmiles_map = {}
        all_indexes = []
        indexes_to_drop = []
        for index, row in dataframe.iterrows():
            all_indexes.append(index)
            activity_id = row["activity.chembl_activity_id"]
            pdb = row["pdb"]
            name = row["name"]
            key = (pdb, name)
            if not key in pdb_usmiles_map.keys():
                pdb_usmiles_map[key] = {activity_id: [index]}
            else:
                if not activity_id in pdb_usmiles_map[key].keys():
                    pdb_usmiles_map[key][activity_id] = [index]
                else:
                    pdb_usmiles_map[key][activity_id].append(index)
        for key in pdb_usmiles_map.keys():
            all_names = pdb_usmiles_map[key].keys()
            if len(all_names) > 1:
                new_indexes = []
                for name in all_names:
                    new_indexes.extend([index for index in pdb_usmiles_map[key][name]])
                new_data = dataframe.loc[new_indexes]
                best_matching = 0
                best_matching_name = ""
                found = False
                for index, row in new_data.iterrows():
                    molecule_matching = row[
                        "small_molecule_matching_confidence_level.comment"
                    ]
                    activity_id = row["activity.chembl_activity_id"]
                    if best_matching == 0:
                        # Found an initial matching
                        best_matching = molecule_matching
                        best_matching_name = [activity_id]
                        found = True
                        continue
                    for matching_level in matchings_molecules:
                        if (
                            matching_level == best_matching
                            and matching_level == molecule_matching
                        ):
                            found = False
                            best_matching_name.append(activity_id)
                        elif matching_level == best_matching:
                            continue
                        elif matching_level == molecule_matching:
                            found = True
                            best_matching = matching_level
                            best_matching_name = [activity_id]
                if found:
                    blub = list(
                        set(pdb_usmiles_map[key].keys()) - set(best_matching_name)
                    )
                else:
                    choice = random.choice(best_matching_name)
                    blub = list(set(pdb_usmiles_map[key].keys()) - set([choice]))
                discard = [pdb_usmiles_map[key][test] for test in blub]
                indexes_to_drop.extend([x[0] for x in discard])
        return dataframe.loc[list(set(all_indexes) - set(indexes_to_drop))]

    def get_single_best_small_molecule(self, dataframe: pd.DataFrame) -> pd.DataFrame:
        """
        Filters a dataframe for the best matching molecules
        @param dataframe:
        @return:
        """
        best_match = ""
        best_match_indexes = 0
        for index, row in dataframe.iterrows():
            current_match = row["small_molecule_matching_confidence_level.comment"]
            if best_match == "":
                best_match = current_match
                best_match_indexes = [index]
                continue
            for matching in matchings_molecules:
                if current_match == best_match:
                    best_match_indexes.append(index)
                    break
                elif current_match == matching:
                    best_match = current_match
                    best_match_indexes = [index]
                    break
                elif best_match == matching:
                    break
        return dataframe.loc[best_match_indexes]

    def get_single_best_protein_matching(self, dataframe: pd.DataFrame) -> pd.DataFrame:
        """
        Filters a dataframe for the best matching protein
        @param dataframe: dataframe to be filtered
        @return: Filtered dataframe
        """
        best_match = ""
        best_match_indexes = 0
        for index, row in dataframe.iterrows():
            current_match = row["Level_pm"]
            if best_match == "":
                best_match = current_match
                best_match_indexes = [index]
                continue
            for matching in matchings_proteins:
                if current_match == best_match:
                    best_match_indexes.append(index)
                    break
                elif current_match == matching:
                    best_match = current_match
                    best_match_indexes = [index]
                    break
                elif best_match == matching:
                    break
        return dataframe.loc[best_match_indexes]

    def get_single_best_edia_match(self, dataframe: pd.DataFrame) -> pd.DataFrame:
        """
        Gets the rows with the best edia in a dataframe
        @param dataframe: dataframe to be filtered
        @return: filtered dataframe
        """
        best_match = 0
        best_match_indexes = 0
        for index, row in dataframe.iterrows():
            edia = float(row["EDIAm"])
            if edia > best_match:
                best_match = edia
                best_match_indexes = [index]
            elif edia == best_match:
                best_match_indexes.append(index)
        return dataframe.loc[best_match_indexes]

    def check_if_dataframe_is_ok(
        self, dataframe: pd.DataFrame, min_molecules: int, min_difference: int
    ) -> bool:
        """
        Checks a dataframe if the correct format is present (unique usmiles and activity ids)
        Also checks if the correct number of molecules are present and the difference in values is big enough
        @param dataframe: dataframe to be filtered
        @param min_molecules: Minimal number of molecules present
        @param min_difference: Minimal difference in order of magnitude of activities
        @return: True if dataframe is ok, False if not
        """
        activities = []
        all_usmiles = []
        if dataframe.shape[0] < min_molecules:
            return False
        assay_min_level = {x: np.finfo(float).max for x in matchings_proteins}
        assay_max_level = {x: np.finfo(float).min for x in matchings_proteins}
        all_usmiles_level = {x: set() for x in matchings_proteins}
        found = False
        found_levels = set()
        for index, row in dataframe.iterrows():
            activity = row["activity.chembl_activity_id"]
            usmiles = row["usmiles"]
            score = row["Level_pm"]
            value = float(row["activity.standard_value"])
            if math.isclose(value, 0):
                # Somehow there are some examples of 0 as value which we obviously want to omit
                continue
            found_levels.add(
                score
            )  # We have to know which levels to not use the max/min floats at a later stage
            unit = row["activity.standard_units"]
            for sublevel in matchings_proteins:
                if sublevel == score or found:
                    found = True
                    all_usmiles_level[sublevel].add(usmiles)
                    if value < assay_min_level[sublevel]:
                        assay_min_level[sublevel] = value
                    if value > assay_max_level[sublevel]:
                        assay_max_level[sublevel] = value
            if activity in activities:
                print(f"Activity {activity} is at least twice present....")
                return False
            if usmiles in all_usmiles:
                print(f"Usmiles {usmiles} at least twice present...")
                return False
            activities.append(activity)
            all_usmiles.append(usmiles)
        works = False
        for sublevel in matchings_proteins:
            if sublevel not in found_levels:
                # Filter for found levels to not use max/min floats
                continue
            diff_order_of_magnitude = orderOfMagnitude(
                assay_max_level[sublevel]
            ) - orderOfMagnitude(assay_min_level[sublevel])
            if (
                len(all_usmiles_level[sublevel]) >= min_molecules
                and diff_order_of_magnitude >= min_difference
            ):
                works = True
        return works

    def decide_multiple_rows(
        self, dataframe: pd.DataFrame, columnname: str
    ) -> pd.DataFrame:
        """
        Decides which row to use when there are multiple rows with the same values in columname
        Decides this using the best molecule matching, protein matching, edia in this order
        If there are rows that are identical in all three, a random one is used (rare)
        @param dataframe: dataframe to be filtered
        @param columnname: row that shall be used as discriminator
        @return: filtered dataframe
        """
        activities_index_map = {}
        all_indexes = []
        indexes_drop_final = []
        for index, row in dataframe.iterrows():
            all_indexes.append(index)
            activity = row[columnname]
            if activity in activities_index_map:
                activities_index_map[activity].append(index)
            else:
                activities_index_map[activity] = [index]
        for key in activities_index_map.keys():
            if len(activities_index_map[key]) > 1:
                new_dataframe = dataframe.loc[activities_index_map[key]]
                newer_dataframe = self.get_single_best_small_molecule(new_dataframe)
                even_newer_dataframe = self.get_single_best_protein_matching(
                    newer_dataframe
                )
                newest_dataframe = self.get_single_best_edia_match(even_newer_dataframe)
                if newest_dataframe.shape[0] > 1:
                    final_dataframe = newest_dataframe.iloc[
                        [random.randint(0, newest_dataframe.shape[0] - 1)]
                    ]
                else:
                    final_dataframe = newest_dataframe
                indexes_drop_final.extend(
                    list(set(new_dataframe.index) - set(final_dataframe.index))
                )
        indexes_to_keep = list(set(all_indexes) - set(indexes_drop_final))
        return dataframe.loc[indexes_to_keep]
