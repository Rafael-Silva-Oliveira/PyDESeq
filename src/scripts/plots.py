import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn import decomposition
from sklearn.cluster import KMeans
import json
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from gseapy import dotplot
import matplotlib.pyplot as plt
from datetime import datetime
from loguru import logger
date = datetime.now().strftime("%Y_%m_%d-%I_%M_%S_%p")


class Plots:

    def __init__(self, plot_volcano: bool, plot_pca: bool, plot_heatmap: bool, plot_gsea: bool, plot_dotplot: bool, config: dict):
        self.plot_volcano = plot_volcano
        self.plot_pca = plot_pca
        self.plot_heatmap = plot_heatmap
        self.plot_gsea = plot_gsea
        self.plot_dotplot = plot_dotplot
        self.config = config

    def volcano_fn(self, dds_stats_df: pd.DataFrame):
        from sanbomics.plots import volcano
        volcano(dds_stats_df, symbol="Symbols", to_label=20, save=True)

    def pca_fn(self, dds):
        # pca_data = pca_fn(count_matrix_processed,
        #                  filter_by_quantile=False, quant=0.95)
        sc.tl.pca(dds)
        sc.pl.pca(dds, color="Condition", size=300)
        # sc.pp.neighbors(dds, n_neighbors=5)
        # sc.tl.leiden(dds, resolution=0.5)
        # sc.tl.umap(dds)
        # sc.pl.umap(dds, color=["leiden", "Condition"])

    def dotplot_fn(self, go_ontologies):
        ax = dotplot(go_ontologies.res2d, column="FDR q-val", title="DisGeNET Ontologies",
                     cmap=plt.cm.viridis, size=6, figsize=(10, 10), cutoff=0.50, show_ring=True, ofname="dotplot.png")

    def gseaplot_fn(self, go_ontologies):
        # Make a GSEA plot based on a given term
        terms = go_ontologies.res2d.Term
        plt.figure(figsize=(20, 20))  # Set the desired figure size
        axs = go_ontologies.plot(terms=terms[1:10],
                                 show_ranking=True,
                                 figsize=(20, 20)
                                 )
        axs.figure.set_size_inches(25, 25)  # Set a bigger figure size
        axs.savefig("gsea.png")

        tt = 2

    def heatmap_fn(self, dds, significant_genes, ensg_to_gene_mapper):

        dds.layers["log1p"] = np.log1p(dds.layers["normed_counts"])

        # Filter only by the significant genes by a given column order showing only the top n genes.
        if self.config["plotting"]["plot_heatmap"]["filter"]["use_filter"]:
            show_top_n = self.config["plotting"]["plot_heatmap"]["filter"]["show_top_n"]
            order_significant_genes_by = self.config["plotting"]["plot_heatmap"][
                "filter"]["order_significant_genes_by"]
            ascending = self.config["plotting"]["plot_heatmap"]["filter"]["ascending"]

            dds_sigs = dds[:, significant_genes.sort_values(
                order_significant_genes_by, ascending=ascending).index[:show_top_n]]
            heatmap_name = "filtered_heatmap"
        else:
            dds_sigs = dds[:, significant_genes.index]
            heatmap_name = "unfiltered_heatmap"

        heatmap_df = pd.DataFrame(
            dds_sigs.layers["log1p"].T, columns=dds_sigs.obs_names, index=dds_sigs.var_names)
        heatmap_df.index = heatmap_df.index.map(
            ensg_to_gene_mapper.mapper)

        heatmap = sns.clustermap(heatmap_df, z_score=0, cmap="RdYlBu_r")
        heatmap.savefig(f"{heatmap_name}.png")

    def run_plots(self, dds, dds_stats_df,
                  significant_genes,  ensg_to_gene_mapper,  go_ontologies=None):
        logger.info("Creating plots.")
        if self.plot_volcano:
            self.volcano_fn(dds_stats_df=dds_stats_df)
        if self.plot_pca:
            self.pca_fn(dds)
        if self.plot_heatmap:
            self.heatmap_fn(dds, significant_genes, ensg_to_gene_mapper)
        if self.plot_gsea:
            self.gseaplot_fn(go_ontologies)
        if self.plot_dotplot:
            self.dotplot_fn(go_ontologies)


def pca_fn(df, filter_by_quantile=True, quant=0.95):
    data = df.copy()
    pca = decomposition.PCA(n_components=2)

    if filter_by_quantile == True:

        data = data[data.Sum > data.Sum.quantile(quant)]

    comps = pca.fit(data).transform(data)
    # labels = {
    #     str(i): f"PC {i+1} ({var:.1f}%)"
    #     for i, var in enumerate(pca.explained_variance_ratio_ * 100)
    # }
    per_var = np.round(pca.explained_variance_ratio_*100, decimals=1)
    labels = ['PC' + str(x) for x in range(1, len(per_var)+1)]

    # plt.bar(x=range(1, len(per_var)+1), height=per_var, tick_label=labels)
    # plt.ylabel('Percentage of explained variability')
    # plt.xlabel('Principal component')
    # plt.title('Scree plot')
    # plt.show()

    pca_df = pd.DataFrame(comps, index=data.index, columns=labels)
    pca_df.drop("Sum", inplace=True)

    plt.scatter(pca_df.PC1, pca_df.PC2, alpha=0.5, c=[
        "red", "red", "red", "red", "red", "cyan", "cyan", "cyan", "cyan", "cyan"])
    plt.title("PCA graph")
    plt.xlabel('PC1- {0}%'.format(per_var[0]))
    plt.ylabel('PC2- {0}%'.format(per_var[1]))
    # plt.xlim(0, 10000)
    # plt.ylim(-2500, 10000)
    plt.show()

    # Apply KMeans
    # Coloring the components
    sample_color = {}

    for index, sample_name in enumerate(data.columns):
        if sample_name not in sample_color:
            sample_color[sample_name] = 0
    vect_palette = sns.color_palette(None, len(sample_color))

    # como estamos a olhar para 2 PC, o n_clusters vai ser 2. O metodo
    kmeans = KMeans(n_clusters=2, random_state=0)

    # Compute cluster centers and predict cluster indices
    cluster_labels = kmeans.fit_predict(pca_df)

    pca_df["cluster_labels"] = cluster_labels

    # filter rows of original data
    filtered_label0 = pca_df[cluster_labels == 0]
    filtered_label1 = pca_df[cluster_labels == 1]

    # Plotting the results
    plt.scatter(filtered_label0.iloc[:, 0],
                filtered_label0.iloc[:, 1], color='red')
    plt.scatter(filtered_label1.iloc[:, 0],
                filtered_label1.iloc[:, 1], color='black')

    plt.show()
    return pca_df
