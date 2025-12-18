import sys

import pandas as pd

import definitions
import plotnine as p9


def prepare_dataframe_for_2d_plot(
    df_path,
    column_x,
    column_size,
    column_y,
    suffix,
    column_wrap=None,
    order_y=None,
    order_x=None,
    nrow=None,
    ncol=None,
    color_cutoff=0.5,
    rowwise=True,
    totals=False,
    order_facet=None,
    label_y=None,
    rotation=30,
    percentage=False,
    textsize=9
):
    data = pd.read_csv(df_path)
    if rowwise:
        if column_wrap is None:
            data["Category Sum"] = data.groupby([column_x])[column_size].transform("sum")
        else:
            data["Category Sum"] = data.groupby([column_x, column_wrap])[
            column_size
            ].transform("sum")
        if totals:
            new_rows = []
            for index, row in data.iterrows():
                row_to_copy = data.loc[index].copy()
                row_to_copy[column_y] = "Total"
                row_to_copy[column_size] = row_to_copy["Category Sum"]
                row_to_copy = [x for x in row_to_copy]
                row_to_copy.extend([0, "high"])
                new_rows.append(row_to_copy)
            unique_rows = [list(t) for t in set(tuple(inner) for inner in new_rows)]
    else:
        if column_wrap is None:
            data["Category Sum"] = data.groupby([column_y])[
                column_size
            ].transform("sum")
        else:
            data["Category Sum"] = data.groupby([column_y, column_wrap])[
                column_size
            ].transform("sum")
        if totals:
            new_rows = []
            for index, row in data.iterrows():
                row_to_copy = data.loc[index].copy()
                row_to_copy[column_x] = "Total"
                row_to_copy[column_size] = row_to_copy["Category Sum"]
                row_to_copy = [x for x in row_to_copy]
                row_to_copy.extend([0, "high"])
                new_rows.append(row_to_copy)
            unique_rows = [list(t) for t in set(tuple(inner) for inner in new_rows)]
    df_with_totals = data
    df_with_totals["Percentage"] = (
        df_with_totals[column_size] / df_with_totals["Category Sum"]
    )
    df_with_totals["p_group"] = pd.cut(
        df_with_totals["Percentage"],
        (
            0,
            df_with_totals["Percentage"].max() * color_cutoff,
            df_with_totals["Percentage"].max(),
        ),
        labels=("low", "high"),
    )
    if totals:
        print(df_with_totals.head())
        print(unique_rows)
        df_with_totals = pd.concat(
            [
                df_with_totals,
                pd.DataFrame(unique_rows, columns=[x for x in data.columns.values]),
            ],
            ignore_index=True,
        )
    if (not totals or not rowwise) and order_y :
        order_y = [x for x in order_y if x != "Total"]
    new_df_path = df_path.parent / (df_path.stem + f"_{suffix}.csv")
    if order_facet is not None:
        df_with_totals[column_wrap] = df_with_totals[column_wrap].astype("category").cat.reorder_categories(order_facet)
    df_with_totals.to_csv(new_df_path)
    make_2d_quality_plot_of_results(
        df_with_totals,
        order_y=order_y,
        column_y=column_y,
        column_x=column_x,
        column_size="Percentage",
        column_wrap=column_wrap,
        suffix=suffix,
        order_x=order_x,
        nrow=nrow,
        ncol=ncol,
        rowwise=rowwise,
        rotation=rotation,
        textsize=textsize,
        percentage=percentage,
    )


def make_2d_quality_plot_of_results(
    df,
    column_y,
    column_x,
    column_size,
    column_wrap,
    suffix,
    order_y=None,
    order_x=None,
    nrow=None,
    ncol=None,
    rowwise=True,
    label_y=None,
    rotation=30,
    textsize=11,
    percentage=False,
    thousandspace=True
):
    if rowwise:
        df_small = df[df[column_y] != "Total"]
    else:
        df_small = df[df[column_x] != "Total"]
    plot = p9.ggplot()
    plot += p9.geom_tile(
        df_small,
        p9.aes(width=0.95, height=0.85, x=column_x, y=column_y, fill=column_size)
    )
    if percentage:
        df["label"] = df["Size"].astype(str) + "(" + (df[column_size] * 100).round(2).astype(str) + "%)"
        column_label = "label"
    else:
        column_label = "Size"
    if thousandspace:
        df['label_formatted'] = df[column_label].apply(lambda x: f"{x:,}".replace(",", " "))
        column_label = "label_formatted"
    plot += p9.geom_text(
        df,
        p9.aes(x=column_x, y=column_y, label=column_label, color="p_group"),
        show_legend=False,size=textsize
    )
    if column_wrap is not None:
        plot += p9.facet_wrap(column_wrap, nrow=nrow, ncol=ncol)
    plot += p9.theme(panel_grid=p9.element_blank(), panel_background=p9.element_blank())
    if order_y:
        plot += p9.scale_y_discrete(limits=order_y)
    if order_x:
        plot += p9.scale_x_discrete(limits=order_x[::1])
    if label_y:
        plot += p9.ylab(label_y)
    if rotation > 0:
        plot += p9.theme(axis_text_x=p9.element_text(rotation=rotation, hjust=1))
    plot += p9.scale_color_manual({"high": "black", "low": "white"})  # new
    #plot += p9.scale_fill_continuous("summer",na_value="#00FFFFFF")
    #plot += p9.scale_fill_manual(values=["#21918c", "#5ec962"])
    plot += p9.theme(text=p9.element_text(family="sans-serif"))
    plot += p9.labs(x=" ", y=" ")
    plot.save(definitions.plot_path / f"2d_quality_plot_{suffix}.png", dpi=800)

def main():
    """input_path = definitions.input_files_path / "paula_file_new.csv"
    order_y = ["belowtwo", "twotothree", "overthree"]
    order_y = ["RMSD <= 2Å", "2Å < RMSD <= 3Å", "RMSD > 3Å"]
    order_y.reverse()
    prepare_dataframe_for_2d_plot(
        df_path=input_path,
        column_x="Type",
        column_size="Size",
        suffix="paula",
        column_y="Cutoff",
        totals=True,
        rowwise=True,
        order_y=order_y,
        rotation=0,
        textsize=12,
        percentage=True
    )"""


    order_y = ["[100]", "[95,100)", "[90,95)", "[85,90)", "[80,85)", "Total"]
    order_y.reverse()
    order_x = [
        "InChi-Str Connection Layer",
        "InChi-Str Connection & Hydrogen Layer",
        "SMI w/o Stereo",
        "SMI w Stereo",
        "InChI-Key",
    ]
    order_x_total = [
        "InChi-Str Connection Layer",
        "InChi-Str Connection & Hydrogen Layer",
        "SMI w/o Stereo",
        "SMI w Stereo",
        "InChI-Key",
        "Total",
    ]
    order_sequences = ["ATOMSEQ", "SEQRES", "SEQADV"]
    input_path = definitions.input_files_path / "example_data_result_plot.csv"
    prepare_dataframe_for_2d_plot(
        df_path=input_path,
        order_y=order_y,
        column_x="PDB",
        column_wrap="ChEMBL",
        column_size="Size",
        suffix="automated_test",
        column_y="Percent Identity Interval",
        totals=True,
        rowwise=True,
        order_x=order_sequences
    )

    input_path = definitions.input_files_path / "component_sequence_ligand_matches.csv"
    prepare_dataframe_for_2d_plot(
        df_path=input_path,
        order_y=order_y,
        column_x="Mol Matching",
        column_wrap="PDB",
        column_size="Size",
        suffix="compenent_ligand_matches",
        column_y="Percent identity interval",
        ncol=1,
        order_x=order_x,
        color_cutoff=0.65,
    )
    prepare_dataframe_for_2d_plot(
        df_path=input_path,
        order_y=order_y,
        column_x="Mol Matching",
        column_wrap="PDB",
        column_size="Size",
        suffix="compenent_ligand_matches_columnwise",
        column_y="Percent identity interval",
        ncol=1,
        order_x=order_x,
        color_cutoff=0.65,
        rowwise=False,
        order_facet=order_sequences
    )
    prepare_dataframe_for_2d_plot(
        df_path=input_path,
        order_y=order_y,
        column_x="Mol Matching",
        column_wrap="PDB",
        column_size="Size",
        suffix="compenent_ligand_matches_columnwise_totals",
        column_y="Percent identity interval",
        ncol=1,
        order_x=order_x_total,
        color_cutoff=0.65,
        rowwise=False,
        totals=True,
        order_facet=order_sequences
    )
    prepare_dataframe_for_2d_plot(
        df_path=input_path,
        order_y=order_y,
        column_x="PDB",
        column_wrap="Mol Matching",
        column_size="Size",
        suffix="compenent_ligand_matches_consistent",
        column_y="Percent identity interval",
        order_x=order_x,
    )
    prepare_dataframe_for_2d_plot(
        df_path=input_path,
        order_y=order_x,
        column_wrap="PDB",
        column_y="Mol Matching",
        column_size="Size",
        suffix="compenent_ligand_matches_novel",
        column_x="Percent identity interval",
        order_x=order_y,
        ncol=1,
    )
    input_path = definitions.input_files_path / "variant_sequence_ligand_matches.csv"
    prepare_dataframe_for_2d_plot(
        df_path=input_path,
        order_y=order_y,
        column_x="Mol Matching",
        column_wrap="PDB",
        column_size="Size",
        suffix="variant_ligand_matches",
        column_y="Percent identity interval",
        ncol=1,
        order_x=order_x,
    )
    prepare_dataframe_for_2d_plot(
        df_path=input_path,
        order_y=order_x,
        column_wrap="PDB",
        column_y="Mol Matching",
        column_size="Size",
        suffix="variant_ligand_matches_novel",
        column_x="Percent identity interval",
        order_x=order_y,
        ncol=1,
    )

if __name__ == "__main__":
    main()
