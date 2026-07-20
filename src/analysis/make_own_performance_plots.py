import numpy as np
import pandas as pd
import plotnine as p9
from statsmodels.stats.proportion import proportion_confint
from statsmodels.stats.proportion import proportions_ztest

from src.utils.definitions import output_path


def annotate_significance_by_benchmark(df, alpha=0.05):
    df = df.copy()
    df["Rate"] = df["Positives"] / df["Total"]
    df["significance"] = None

    # Process each benchmark independently
    for benchmark, sub in df.groupby("Benchmark"):
        # Identify the best method within this benchmark
        best_idx = sub["Rate"].idxmax()
        best_row = df.loc[best_idx]

        # Assign 'best'
        df.at[best_idx, "significance"] = "Best method"

        # Compare all other methods in the same benchmark to the best one
        for i, row in sub.iterrows():
            if i == best_idx:
                continue

            # Z-test for two proportions
            count = np.array([row["Positives"], best_row["Positives"]])
            nobs = np.array([row["Total"], best_row["Total"]])

            z, p = proportions_ztest(count, nobs)

            if p >= alpha:
                df.at[i, "significance"] = "Similar method"
            else:
                df.at[i, "significance"] = "Worse method"

    return df


def create_performance_plot(dataframe, output_name, method_list):
    plot = p9.ggplot(dataframe, p9.aes(x="Method", y="Ratio", color="significance"))
    plot += p9.geom_point()
    plot += p9.coord_flip()
    plot += p9.theme(
    text=p9.element_text(size=14),
    axis_title=p9.element_text(size=16),
    axis_text=p9.element_text(size=12),
    legend_title=p9.element_text(size=13),
    legend_text=p9.element_text(size=11),
    plot_title=p9.element_text(size=18))
    plot += p9.scale_color_manual(values=colors)
    plot += p9.geom_hline(yintercept=50, color="black", linetype="dashed", size=1)
    plot += p9.ylab("Percentage correct activity pairs")
    plot += p9.scale_x_discrete(limits=method_list)
    plot += p9.theme(
        subplots_adjust={"wspace": 10, "hspace": 10}, legend_position="none"
    )
    plot += p9.facet_wrap("Benchmark", ncol=2)
    # plot += p9.ylim([40, 100])

    # Fehlerbalken inklusive der Ästhetik für Fehlergrößen
    plot += p9.geom_errorbar(p9.aes(ymin="lower_bound", ymax="upper_bound", width=0.2))
    plot.save(base_output_path / output_name, units="mm")


base_output_path = output_path

mmp_data = pd.read_csv(base_output_path / "ballrooms_rmsd.csv")
mmp_data = mmp_data[
    ~mmp_data["Method"].isin(["jamda_rescoring", "jamda_workingrmsd"])
]  # We do not want these entries for BALLLROOM S
mmp_data.loc[mmp_data["Method"].astype(str) == "Vina", "Method"] = f"Autodock Vina"
mmp_data.loc[mmp_data["Method"].astype(str) == "tpsa", "Method"] = "TPSA"
mmp_data.loc[mmp_data["Method"].astype(str) == "masse", "Method"] = "Weight"
mmp_data.loc[mmp_data["Method"].astype(str) == "logp", "Method"] = "cLogP"
mmp_data.loc[mmp_data["Method"].astype(str) == "jamda", "Method"] = "JAMDA"
mmp_data.loc[mmp_data["Method"].astype(str) == "DOCK", "Method"] = "DOCK"
mmp_data.loc[
    mmp_data["Method"].astype(str) == "boltz_pred", "Method"
] = f"Boltz2 activity"
mmp_data.loc[
    mmp_data["Method"].astype(str) == "boltz_affinity", "Method"
] = f"Boltz2 affinity"
mmp_data.loc[
    mmp_data["Method"].astype(str) == "nn_tobi", "Method"
] = f"Nearest Neighbor"
mmp_data.loc[
    mmp_data["Method"].astype(str) == "jamda_rescoring", "Method"
] = f"JAMDA Rescoring"
mmp_data.loc[
    mmp_data["Method"].astype(str) == "jamda_workingrmsd", "Method"
] = f"JAMDA good poses"

redocking_data = pd.read_csv(base_output_path / "ballrooma_rmsd.csv")
redocking_data.loc[
    redocking_data["Method"].astype(str) == "Vina", "Method"
] = "Autodock Vina"
redocking_data.loc[redocking_data["Method"].astype(str) == "tpsa", "Method"] = "TPSA"
redocking_data.loc[redocking_data["Method"].astype(str) == "masse", "Method"] = "Weight"
redocking_data.loc[redocking_data["Method"].astype(str) == "logp", "Method"] = "cLogP"
redocking_data.loc[redocking_data["Method"].astype(str) == "jamda", "Method"] = "JAMDA"
redocking_data.loc[redocking_data["Method"].astype(str) == "DOCK", "Method"] = "DOCK"
redocking_data.loc[
    redocking_data["Method"].astype(str) == "boltz_pred", "Method"
] = f"Boltz2 activity"
redocking_data.loc[
    redocking_data["Method"].astype(str) == "boltz_affinity", "Method"
] = f"Boltz2 affinity"
redocking_data.loc[
    redocking_data["Method"].astype(str) == "nn_tobi", "Method"
] = f"Nearest Neighbor"
redocking_data.loc[
    redocking_data["Method"].astype(str) == "jamda_rescoring", "Method"
] = f"JAMDA Rescoring"
redocking_data.loc[
    redocking_data["Method"].astype(str) == "jamda_workingrmsd", "Method"
] = f"JAMDA good poses"

all_relevant_methods = [
    "TPSA",
    "cLogP",
    "Weight",
    "Boltz2 affinity",
    "Nearest Neighbor",
]
traditional_relevant_methods = [
    "TPSA",
    "cLogP",
    "Weight",
    "DOCK",
    "Autodock Vina",
    "JAMDA",
]
rmsd_relevant_methods = ["Weight", "JAMDA", "JAMDA good poses", "JAMDA Rescoring"]
significane = 0.05

datasets = [redocking_data, mmp_data]
for dataset in datasets:
    upper_bound_list = []
    lower_bound_list = []
    best_method_trad_list = []
    best_method_ml_list = []
    for index, row in dataset.iterrows():
        lower_bound, upper_bound = proportion_confint(
            count=row["Positives"], nobs=row["Total"], alpha=significane
        )
        upper_bound_list.append(upper_bound * 100)
        lower_bound_list.append(lower_bound * 100)
    dataset["upper_bound"] = upper_bound_list
    dataset["lower_bound"] = lower_bound_list

redocking_data["Benchmark"] = "BALLROOM A"
redocking_data_ges = redocking_data[redocking_data["Target"] == "ges"]

mmp_data["Benchmark"] = "BALLROOM S"
mmp_data_ges = mmp_data[mmp_data["Target"] == "ges"]

print(mmp_data.head())
print(redocking_data.head())
colors = {"BALLROOM A": "#666a90", "BALLROOM S": "#00461e"}
colors = {"Best method": "blue", "Similar method": "grey", "Worse method": "red"}
complete_data_ges = pd.concat([mmp_data_ges, redocking_data_ges])
complete_data_ml = complete_data_ges[
    complete_data_ges["Method"].isin(all_relevant_methods)
]
complete_data_trad = complete_data_ges[
    complete_data_ges["Method"].isin(traditional_relevant_methods)
]
complete_data_rmsd = complete_data_ges[
    complete_data_ges["Method"].isin(rmsd_relevant_methods)
]
complete_data_rmsd = complete_data_rmsd[
    ~complete_data_rmsd["Benchmark"].isin(["BALLROOM S"])
]
complete_data_trad = annotate_significance_by_benchmark(complete_data_trad, significane)
complete_data_ml = annotate_significance_by_benchmark(complete_data_ml, significane)
complete_data_rmsd = annotate_significance_by_benchmark(complete_data_rmsd, significane)


create_performance_plot(
    complete_data_trad,
    "traditional_methods_performance_rmsd.png",
    traditional_relevant_methods,
)
create_performance_plot(
    complete_data_ml, "ml_methods_performance_rmsd.png", all_relevant_methods
)
create_performance_plot(
    complete_data_rmsd, "rmsd_performance.png", rmsd_relevant_methods
)
