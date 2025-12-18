import plotnine as p9
import pandas as pd

from src import definitions

def make_histogram(data, columnname, plot_output_path, title, binwidth=1, logy=False):
    plot = p9.ggplot(data, p9.aes(x=columnname))
    plot += p9.geom_histogram(binwidth=binwidth, fill="#2a788e")
    plot += p9.ggtitle(title=title)
    #plot += p9.scale_fill_manual("#2a788e")
    plot += p9.xlab("Maximum sequence identity")
    plot += p9.ylab("count (log10-scale)")
    if logy:
        plot += p9.scale_y_log10()
    plot.save(plot_output_path, dpi=2000)


def make_density(data, columnname, plot_output_path, title):
    plot = p9.ggplot(data, p9.aes(x=columnname))
    plot += p9.geom_density()
    plot += p9.ggtitle(title=title)
    plot.save(plot_output_path)


sequence_data = pd.read_csv(
     definitions.output_files_path / "sequence_identity_data_biochem_full.csv",
    index_col=0,
)
sequence_data.drop_duplicates(inplace=True)
make_histogram(
    sequence_data,
    "Best sequence identity",
     definitions.plot_path / "best_sequence_identity_full.png",
    " ",
)
make_histogram(
    sequence_data,
    "Best sequence identity",
     definitions.plot_path / "best_sequence_identity_full_log.png",
    " ",
    logy=True
)
print(sequence_data.shape)
print(sequence_data[sequence_data["Best sequence identity"] < 90].shape)
