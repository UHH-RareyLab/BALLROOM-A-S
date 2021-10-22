from pathlib import Path

import pandas as pd

from src.utils import definitions, db_interaction as db


def get_timecode_pdb_dict():
    timecode_pdb_dict = {}
    timecode_file_path = definitions.mirror_base_path / "derived_data/index/entries.idx"
    with open(timecode_file_path, "r") as timecode_file:
        for line in timecode_file.readlines():
            line = line.strip()
            if line.startswith("-") or line.startswith("IDCODE"):
                continue
            (
                pdb,
                header,
                accession_date,
                compound,
                source,
                author_list,
                resolution,
                experiment_type,
            ) = line.split("\t")
            timecode_pdb_dict[pdb] = pd.to_datetime(
                accession_date.strip(), format="%m/%d/%y"
            )
    return timecode_pdb_dict


def get_chembl_release_and_date_for_activity(activity_id):
    activity_id = int(float(activity_id))
    query = (
        f"SELECT docs.chembl_release_id, docs.year from assays "
        f"JOIN docs on assays.doc_id = docs.doc_id "
        f"JOIN activities on assays.assay_id = activities.assay_id "
        f"WHERE activities.activity_id = '{activity_id}'"
    )
    data = db.query_chembl(query)
    return data.loc[0, "chembl_release_id"], data.loc[0, "year"]


def get_chembl_release_and_date_for_assay(assay_id):
    query = (
        f"SELECT docs.chembl_release_id, docs.year from assays "
        f"JOIN docs on assays.doc_id = docs.doc_id "
        f"WHERE assays.chembl_id = '{assay_id}'"
    )
    data = db.query_chembl(query)
    assert data.shape[0] == 1
    return data.loc[0, "chembl_release_id"], data.loc[0, "year"]


def get_chembl_release_dict():
    query = f"SELECT * FROM chembl_release"
    data = db.query_chembl(query)
    release_dict = {}
    for index, row in data.iterrows():
        release_dict[row["chembl_release_id"]] = pd.to_datetime(row["creation_date"])
    return release_dict


def annotate_timecodes_on_benchmark(benchmark_path):
    ballrooma_path = benchmark_path / "ballrooma" / "complete_data.csv"
    ballrooms_path = benchmark_path / "ballrooms" / "complete_data.csv"
    benchmarks = {"ballrooma": ballrooma_path, "ballrooms": ballrooms_path}
    chembl_release_dict = get_chembl_release_dict()
    timecode_pdb_dict = get_timecode_pdb_dict()
    for benchmark in benchmarks:
        earliest_times = []
        latest_times = []
        assay_times = []
        data = pd.read_csv(benchmarks[benchmark], index_col=0)
        for index, row in data.iterrows():
            all_dates = []
            assay = row["Assay"]
            dockinto1 = row["Dockinto1"]
            dockinto2 = row["Dockinto2"]
            raw_pdbs = dockinto1.split(";") + dockinto2.split(";")
            pdbs = set(x[0:4] for x in raw_pdbs if not x == "")
            dates = [timecode_pdb_dict[x] for x in pdbs]
            chembl_realease, date = get_chembl_release_and_date_for_assay(assay)
            assay_date = chembl_release_dict[chembl_realease]
            dates.append(assay_date)
            dates = pd.Series(dates)
            earliest_times.append(dates.min())
            latest_times.append(dates.max())
            assay_times.append(assay_date)
        data["earliest_available_time"] = earliest_times
        data["latest_available_time"] = latest_times
        data["assay_time"] = assay_times
        data.to_csv(benchmarks[benchmark])
        print(earliest_times, latest_times)


annotate_timecodes_on_benchmark(
    Path(
        "/work/gutermuth/bob_analyse_skripte/gutermuth/benchmark_directory/presubmission_ballroom"
    )
)
