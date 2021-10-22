import ast

import pandas as pd
import plotnine as p9

from src.utils.db_interaction import query_chembl
from src.utils.definitions import base_path

output_data_path = (
    base_path
    / "gutermuth"
    / "output"
    / "chembl_contradiction_analysis"
    / "target_contradiction_data_difference_1_annotated.csv"
)
data = pd.read_csv(output_data_path)
output_path = base_path / "gutermuth" / "output"

print(data.head())
counts = data["target_type"].value_counts()
print(counts)
n = 50
values_to_keep = counts[counts > n].index
filtered_df = data[data["target_type"].isin(values_to_keep)]
plot = p9.ggplot(filtered_df, p9.aes(x="N total", y="N contradictions"))
plot += p9.geom_bin_2d()
plot += p9.scale_y_log10()
plot += p9.scale_x_log10()
plot += p9.labs(x="Number of total assays", y="Number of assays to delete for harmony")
plot.save(output_path / "chembl_contradictions_2d.png")

plot = p9.ggplot(filtered_df, p9.aes(x="N total", y="ratio"))
plot += p9.geom_bin_2d()
plot += p9.scale_x_log10()
plot += p9.labs(x="Number of total assays", y="Ratio of assays to delete for harmony")
plot.save(output_path / "chembl_contradictions_2d_ratio.png")


plot = p9.ggplot(filtered_df, p9.aes(x="N total", y="N contradictions"))
plot += p9.geom_bin_2d()
plot += p9.scale_y_log10()
plot += p9.scale_x_log10()
plot += p9.facet_wrap("target_type")
plot += p9.labs(x="Number of total assays", y="Number of assays to delete for harmony")
plot.save(output_path / "chembl_contradictions_2d_target_type.png")

plot = p9.ggplot(filtered_df, p9.aes(x="target_type", y="ratio"))
plot += p9.geom_violin()
plot += p9.theme(axis_text_x=p9.element_text(angle=45, hjust=1))
plot += p9.scale_y_log10()
plot.save(output_path / "chembl_contradictions_violin_target_type.png")

plot = p9.ggplot(filtered_df, p9.aes(x="target_type", y="ratio"))
plot += p9.geom_boxplot()
plot += p9.theme(axis_text_x=p9.element_text(angle=45, hjust=1))
plot += p9.scale_y_log10()
plot += p9.ylab("Ratio of assays to delete to total assays")
plot.save(output_path / "chembl_contradictions_boxplot_target_type.png")

grouped_stats = (
    filtered_df.groupby("target_type")["ratio"].agg(["mean", "std"]).reset_index()
)
print(grouped_stats)

assays_to_discard = []
all_assays = []
for index, row in data.iterrows():
    current_discard = ast.literal_eval(row["Assays to discard"])
    current_all = ast.literal_eval(row["All assays"])
    all_assays.extend(current_all)
    assays_to_discard.extend(current_discard)
print(len(assays_to_discard))
print(len(all_assays))
new_data = []
for assay in all_assays:
    if assay in assays_to_discard:
        discard = True
    else:
        discard = False
    new_data.append([assay, discard])
new_data = pd.DataFrame(new_data, columns=["Assay", "Discard"])
query = f"""SELECT * FROM assays WHERE assays.chembl_id IN ({','.join([f"'{id}'" for id in all_assays])})"""
assay_data = query_chembl(query)
print(assay_data.head())

new_data = new_data.merge(assay_data, how="left", left_on="Assay", right_on="chembl_id")
new_data["has_variant"] = new_data["variant_id"].notna()
print(new_data.head())

current_column = "assay_type"
counts = new_data[current_column].value_counts()
print(counts)
n = 50
values_to_keep = counts[counts > n].index
filtered_df = new_data[new_data[current_column].isin(values_to_keep)]

plot = p9.ggplot(
    filtered_df, p9.aes(x=f"factor({current_column})", fill="factor(Discard)")
)
plot += p9.geom_bar(position="fill")
plot += p9.geom_label(
    p9.aes(label=p9.after_stat("count")),
    stat="count",
    position="fill",
    size=9,
)
plot += p9.guides(fill=p9.guide_legend(title="Discarded"))
plot += p9.xlab("Assay Type")
plot += p9.ylab("Percentage")
plot += p9.scale_y_continuous(labels=lambda x: ["{:.0f}%".format(v * 100) for v in x])
plot.save(output_path / f"chembl_contradictions_{current_column}.png")

current_column = "confidence_score"
counts = new_data[current_column].value_counts()
print(counts)
n = 10
values_to_keep = counts[counts > n].index
filtered_df = new_data[new_data[current_column].isin(values_to_keep)]

plot = p9.ggplot(
    filtered_df, p9.aes(x=f"factor({current_column})", fill="factor(Discard)")
)
plot += p9.geom_bar(position="fill")
plot += p9.geom_label(
    p9.aes(label=p9.after_stat("count")),
    stat="count",
    position="fill",
    size=9,
)
plot += p9.scale_y_continuous(labels=lambda x: ["{:.0f}%".format(v * 100) for v in x])
plot += p9.guides(fill=p9.guide_legend(title="Discarded"))
plot += p9.xlab("Confidence Score")
plot += p9.ylab("Percentage")
plot.save(output_path / f"chembl_contradictions_{current_column}.png")

current_column = "curated_by"
counts = new_data[current_column].value_counts()
print(counts)
n = 50
values_to_keep = counts[counts > n].index
filtered_df = new_data[new_data[current_column].isin(values_to_keep)]

plot = p9.ggplot(
    filtered_df, p9.aes(x=f"factor({current_column})", fill="factor(Discard)")
)
plot += p9.geom_bar(position="fill")
plot += p9.geom_label(
    p9.aes(label=p9.after_stat("count")),
    stat="count",
    position="fill",
    size=9,
)
plot += p9.scale_y_continuous(labels=lambda x: ["{:.0f}%".format(v * 100) for v in x])
plot += p9.guides(fill=p9.guide_legend(title="Discarded"))
plot += p9.xlab("Curated by")
plot += p9.ylab("Percentage")
# plot += p9.geom_text(p9.aes(label=p9.after_stat("count")), stat="count", va="bottom", nudge_y=0.125, position=p9.position_dodge2(width=0.9))
plot.save(output_path / f"chembl_contradictions_{current_column}.png")

current_column = "has_variant"
counts = new_data[current_column].value_counts()
n = 50
values_to_keep = counts[counts > n].index
filtered_df = new_data[new_data[current_column].isin(values_to_keep)]

plot = p9.ggplot(
    filtered_df, p9.aes(x=f"factor({current_column})", fill="factor(Discard)")
)
plot += p9.geom_bar(position="fill")
plot += p9.geom_label(
    p9.aes(label=p9.after_stat("count")),
    stat="count",
    position="fill",
    size=9,
)
plot += p9.scale_y_continuous(labels=lambda x: ["{:.0f}%".format(v * 100) for v in x])
plot += p9.guides(fill=p9.guide_legend(title="Discarded"))
plot += p9.xlab("Has variant")
plot += p9.ylab("Percentage")
# plot += p9.geom_text(p9.aes(label=p9.after_stat("count")), stat="count", va="bottom", nudge_y=0.125, position=p9.position_dodge2(width=0.9))
plot.save(output_path / f"chembl_contradictions_{current_column}.png")
