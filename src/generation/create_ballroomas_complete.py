import os
import shutil
from pathlib import Path

from calc_activity_pairs import create_pair_benchmark_crystal, create_pair_benchmark_mmp
from create_mmp_directory import make_mmp_benchmark
from download_benchmark import download_pair_benchmark, download_mmp_pair_benchmark
from make_benchmark_final import ActivityBenchmarkMaker
from make_siena_clusters import make_siena_clusters
from src.utils.definitions import siena_path, sienatools_path, mm_strict
from src.utils.helper_functions import cleanup_pair_benchmark
from src.utils.helper_functions import create_stats_for_benchmark


def main():
    general_output_path = (
        Path(__file__).parent.parent.parent
        / "gutermuth"
        / "benchmark_directory"
        / "presubmission_ballroom"
    )
    ballrooma_path = general_output_path / "ballrooma"
    ballrooms_path = general_output_path / "ballrooms"
    raw_benchmark_path = general_output_path / "raw_benchmark_start"
    mmp_input_path = general_output_path / "ballrooms_startingpoint"
    if general_output_path.exists():
        shutil.rmtree(general_output_path)
    general_output_path.mkdir()
    os.chdir(general_output_path)
    benchmark = ActivityBenchmarkMaker()
    benchmark.run(
        min_datapoints=1,
        min_difference=1,
        ediam=0.4,
        allowed_protein_matches=["Gold", "Silver", "Bronze"],
        allowed_molecule_matches=mm_strict,
        benchmark_name=raw_benchmark_path.name,
        filter_uncertain=True,
    )
    create_pair_benchmark_crystal(raw_benchmark_path, ballrooma_path, 1)
    make_mmp_benchmark(
        start_path=raw_benchmark_path, end_path=mmp_input_path, rdkit=True, mmpdb=True
    )
    create_pair_benchmark_mmp(mmp_input_path, ballrooms_path, 1)

    cleanup_pair_benchmark(ballrooma_path)
    download_mmp_pair_benchmark(ballrooms_path)
    download_pair_benchmark(ballrooma_path)

    make_siena_clusters(ballrooms_path, siena_path, sienatools_path, "pairs")
    make_siena_clusters(ballrooma_path, siena_path, sienatools_path, "pairs")
    # old_ballrooma = Path("/home/gutermuth/master_table/benchmarks/BALLROOMA")
    # old_ballrooms = Path("/home/gutermuth/master_table/benchmarks/BALLROOMS")
    create_stats_for_benchmark(ballrooma_path)
    create_stats_for_benchmark(ballrooms_path)


if __name__ == "__main__":
    main()
