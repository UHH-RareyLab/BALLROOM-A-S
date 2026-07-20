import json

import pandas as pd
import plotnine as p9

from src.utils.definitions import score_path, output_path, ballrooma

data_ballroom = pd.read_csv(ballrooma / "complete_data.csv")
relevant_complexes = set()
for index, row in data_ballroom.iterrows():
    complex1 = row["ID1"]
    complex2 = row["ID2"]
    relevant_complexes.add(complex1)
    relevant_complexes.add(complex2)


rmsd_normal_path = score_path / "jamda_rmsds.json"
rmsd_rescoring_path = score_path / "jamda_rmsds_rescoring.json"

with open(rmsd_normal_path) as f:
    rmsd_normal = json.load(f)
with open(rmsd_rescoring_path) as f:
    rmsd_rescoring = json.load(f)

# This is necessary, as some complexes are only present in BS
rmsd_rescoring = {
    x: rmsd_rescoring[x] for x in rmsd_rescoring if x in relevant_complexes
}
rmsd_normal = {x: rmsd_normal[x] for x in rmsd_normal if x in relevant_complexes}


all_data = []
found_keys = set()

for key in rmsd_rescoring:
    if key in found_keys:
        continue
    found_keys.add(key)
    rescoring_rmsd = rmsd_rescoring[key]
    if key in rmsd_normal:
        normal_rmsd = rmsd_normal[key]
    else:
        normal_rmsd = +10000
        print(f"Key only present in rescoring dict {key}")
    if not key.startswith("CHEMBL"):
        all_data.append([key, normal_rmsd, rescoring_rmsd])

for key in rmsd_normal:
    if key in found_keys:
        continue
    found_keys.add(key)
    rescoring_rmsd = rmsd_normal[key]
    if key in rmsd_rescoring:
        normal_rmsd = rmsd_rescoring[key]
    else:
        normal_rmsd = +10000
        print(f"Key only present in normal dict {key}")
    if not key.startswith("CHEMBL"):
        all_data.append([key, normal_rmsd, rescoring_rmsd])

all_data = pd.DataFrame(all_data, columns=["key", "normal_rmsd", "rescoring_rmsd"])
print(all_data[all_data["normal_rmsd"] < 2].shape)
cutoff = 2
print(f"There are {all_data[all_data["normal_rmsd"] < cutoff].shape[0]} entries under {cutoff} Angstroem and {all_data[all_data["normal_rmsd"] >= cutoff].shape[0]} above it")
print(f"This results in a correct docking of {100*all_data[all_data["normal_rmsd"] < cutoff].shape[0]/(all_data.shape[0])} percent")
print("rescoring")
print(f"There are {all_data[all_data["rescoring_rmsd"] < cutoff].shape[0]} entries under {cutoff} Angstroem and {all_data[all_data["rescoring_rmsd"] >= cutoff].shape[0]} above it")
print(f"This results in a correct docking of {100*all_data[all_data["rescoring_rmsd"] < cutoff].shape[0]/(all_data.shape[0])} percent")
all_data.to_csv("test.csv")

df = pd.DataFrame(
    {"Docking": all_data["normal_rmsd"], "Rescoring": all_data["rescoring_rmsd"]}
)

long = df.melt(
    value_vars=["Docking", "Rescoring"], var_name="series", value_name="rmsd"
).dropna()
binwidth = 1  # anpassen nach Bedarf
plot = (
    p9.ggplot(long)
    + p9.geom_histogram(
        p9.aes(x="rmsd", fill="series"),
        binwidth=binwidth,
        position=p9.position_dodge(width=0.9),
        color="black",
    )
    + p9.theme(
    text=p9.element_text(size=14),
axis_title = p9.element_text(size=16),
axis_text = p9.element_text(size=12),
legend_title = p9.element_text(size=13),
legend_text = p9.element_text(size=11),
plot_title = p9.element_text(size=18))
    + p9.xlim(-1, 20)
    + p9.scale_fill_manual(values=["#fde725", "#440154"])
    + p9.labs(fill="")
    + p9.xlab("RMSD")
)
plot.save(output_path / "rmsd_comparison_jamda.png")
plot = (
    p9.ggplot(long)
    + p9.geom_density(p9.aes(x="rmsd", fill="series"), color="black")
    + p9.xlim(-1, 20)
    + p9.scale_fill_manual(values=["#fde725", "#440154"])
    + p9.labs(fill="")
)
plot.save(output_path / "rmsd_comparison_jamda_density.png")
plot = (
    p9.ggplot(long)
    + p9.geom_histogram(
        p9.aes(x="rmsd", fill="series"),
        binwidth=binwidth,
        position=p9.position_dodge(width=0.9),
        color="black",
    )
    + p9.xlim(-1, 20)
    + p9.scale_fill_manual(values=["#fde725", "#440154"])
    + p9.labs(fill="")
    + p9.scale_x_log10()
)
plot.save(output_path / "rmsd_comparison_jamda_log.png")
