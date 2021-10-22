import pandas as pd
import plotnine as p9

from src.utils.definitions import output_path

base_output_path = output_path
colors = {"BALLROOM A": "#666a90", "BALLROOM S": "#00461e"}


pair_data_normal = pd.read_csv(
    "/work/gutermuth/bob_analyse_skripte/gutermuth/benchmark_directory/presubmission_ballroom/ballrooma/complete_data.csv"
)
pair_data_mmp = pd.read_csv(
    "/work/gutermuth/bob_analyse_skripte/gutermuth/benchmark_directory/presubmission_ballroom/ballrooms/complete_data.csv"
)
pair_data_normal["Benchmark"] = "BALLROOM A"
pair_data_mmp["Benchmark"] = "BALLROOM S"
full_data = pd.concat([pair_data_mmp, pair_data_normal])
print(pair_data_normal.head())

plot = p9.ggplot(pair_data_mmp, p9.aes(x="Difference"))
plot += p9.geom_histogram()
plot += p9.scale_y_log10()
plot += p9.ggtitle("Histogram of difference [Orders of magnitude] for mmp benchmark")
plot.save(base_output_path / "Difference_histogram_mmp.png")
plot = p9.ggplot(pair_data_normal, p9.aes(x="Difference"))
plot += p9.geom_histogram()
plot += p9.scale_y_log10()
plot += p9.ggtitle(
    "Histogram of difference [Orders of magnitude] for crystal benchmark"
)
plot.save(base_output_path / "Difference_histogram_cryst.png")


plot = p9.ggplot(pair_data_mmp, p9.aes(x="Activity_Type"))
plot += p9.geom_histogram()
plot += p9.scale_y_log10()
plot += p9.ggtitle("Histogram of difference [Orders of magnitude] for mmp benchmark")
plot.save(base_output_path / "activity_types_histogram_mmp.png")
plot = p9.ggplot(pair_data_normal, p9.aes(x="Activity_Type"))
plot += p9.geom_histogram()
plot += p9.scale_y_log10()
plot += p9.ggtitle(
    "Histogram of difference [Orders of magnitude] for crystal benchmark"
)
plot.save(base_output_path / "activity_types_histogram_cryst.png")

plot = p9.ggplot(
    full_data, p9.aes(x="Difference", fill="Benchmark", y=p9.after_stat("count"))
)
plot += p9.geom_histogram(binwidth=1)
plot += p9.facet_wrap("Benchmark", ncol=1)
# plot += p9.geom_density(alpha=0.5)
plot += p9.scale_y_log10()
plot += p9.ylab("Count")
# plot += p9.xlim([0,5])
plot += p9.xlab("Difference [10^x]")
plot += p9.theme(text=p9.element_text(size=16), legend_position="none")
plot += p9.ggtitle("Difference between the orders of magnitude of activities in pairs")
plot += p9.scale_fill_manual(values=colors)
plot.save(
    base_output_path / "Difference_histogram_both.png",
    units="mm",
    height=120,
    width=175,
    dpi=1200,
)

n_ballrooma = 0
n_ballroomb = 0
ic50_a = 0
ic50_b = 0
ec50_a = 0
ec50_b = 0
kd_a = 0
kd_b = 0
ki_a = 0
ki_b = 0
for index, row in full_data.iterrows():
    if row["Benchmark"] == "BALLROOM A":
        if row["Activity_Type"] == "IC50":
            ic50_a += 1
        elif row["Activity_Type"] == "EC50":
            ec50_a += 1
        elif row["Activity_Type"] == "Kd":
            kd_a += 1
        elif row["Activity_Type"] == "Ki":
            ki_a += 1
        else:
            raise RuntimeError
        n_ballrooma += 1
    elif row["Benchmark"] == "BALLROOM S":
        if row["Activity_Type"] == "IC50":
            ic50_b += 1
        elif row["Activity_Type"] == "EC50":
            ec50_b += 1
        elif row["Activity_Type"] == "Kd":
            kd_b += 1
        elif row["Activity_Type"] == "Ki":
            ki_b += 1
        else:
            raise RuntimeError
        n_ballroomb += 1

benchmark_data = []
benchmark_data.append(["IC50", "BALLROOM A", ic50_a / n_ballrooma])
benchmark_data.append(["EC50", "BALLROOM A", ec50_a / n_ballrooma])
benchmark_data.append(["Ki", "BALLROOM A", ki_a / n_ballrooma])
benchmark_data.append(["Kd", "BALLROOM A", kd_a / n_ballrooma])
benchmark_data.append(["IC50", "BALLROOM S", ic50_b / n_ballroomb])
benchmark_data.append(["EC50", "BALLROOM S", ec50_b / n_ballroomb])
benchmark_data.append(["Ki", "BALLROOM S", ki_b / n_ballroomb])
benchmark_data.append(["Kd", "BALLROOM S", kd_b / n_ballroomb])
benchmark_data = pd.DataFrame(
    benchmark_data, columns=["Activity Type", "Benchmark", "Percentage"]
)
plot = p9.ggplot(
    benchmark_data, p9.aes(x="Activity Type", fill="Benchmark", y="Percentage")
)
plot += p9.geom_col()
plot += p9.facet_wrap("Benchmark", ncol=1)
# plot += p9.geom_density(alpha=0.5)
# plot += p9.scale_y_log10()
plot += p9.ylab("Ratio")
plot += p9.xlab("Activity Type")
plot += p9.ggtitle("Activity Type Distribution")
plot += p9.theme(text=p9.element_text(size=16), legend_position="none")
plot += p9.scale_fill_manual(values=colors)
plot.save(
    base_output_path / "Difference_activity_types_both.png",
    units="mm",
    height=120,
    width=175,
    dpi=1400,
)
