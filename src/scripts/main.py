from deseq import DESeq
from plots import Plots
import json
import numpy as np
import pandas as pd
from utils import replace_colnames, txt_2_dataframe
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from sanbomics.tools import id_map
import scanpy as sc
import gseapy as gp
from gseapy.plot import gseaplot
import seaborn as sns
from gseapy import dotplot
import matplotlib.pyplot as plt
from datetime import datetime
from loguru import logger
date = datetime.now().strftime("%Y_%m_%d-%I_%M_%S_%p")


def run_pipeline(CONFIG_PATH: str, diff_exp_matrix_path: str):
    """_summary_

    Args:
        CONFIG_PATH (str): _description_
        diff_exp_matrix_path (str): _description_
    """

    config = json.load(open(CONFIG_PATH, "r"))
    logger.info("Loading configuration parameters.")
    control = config["differential_expression"]["Condition"]["control"]
    experiment = config["differential_expression"]["Condition"]["experiment"]
    species = config["differential_expression"]["species"]
    padj = config["differential_expression"]["thresholds_for_significant_genes"]["padj"]
    log2FoldChange = config["differential_expression"]["thresholds_for_significant_genes"]["log2FoldChange"]
    baseMean = config["differential_expression"]["thresholds_for_significant_genes"]["baseMean"]
    count_matrix_path = config["differential_expression"]["count_matrix"]["path"]
    count_matrix_gene_column = config["differential_expression"]["count_matrix"]["gene_id_column_name"]

    # Gene Set enrichment analysis parameters
    perform_GSEA = config["gene_set_enrichment_analysis"]["perform_GSEA"]
    gene_sets = config["gene_set_enrichment_analysis"]["gene_sets"]
    permutations_num = config["gene_set_enrichment_analysis"]["permutations"]
    filter_by_organ = config["gene_set_enrichment_analysis"]["filter_by_organ"]

    # Plotting parameters
    plot_volcano = config["plotting"]["plot_volcano"]["plot"]
    plot_pca = config["plotting"]["plot_pca"]["plot"]
    plot_heatmap = config["plotting"]["plot_heatmap"]["plot"]
    plot_dotplot = config["plotting"]["plot_dotplot"]["plot"]
    plot_gsea = config["plotting"]["plot_gsea"]["plot"]

    ##########################################################
    # Read input data from bash processing
    count_matrix = pd.read_csv(r"{}".format(count_matrix_path), sep="\t")
    # diff_exp_matrix = txt_2_dataframe(diff_exp_matrix_path)

    # Convert ENSG to gene names
    ensg_to_gene_mapper = id_map(species=species)

    # Instatiate DESeq class
    deseq = DESeq(gene_sets, permutations_num,
                  config, ensg_to_gene_mapper)

    count_matrix_processed, metadata = deseq.pre_processing(
        count_matrix=count_matrix, gene_id_col_name=count_matrix_gene_column)

    dds, stat_res, dds_stats_df, significant_genes, ranked_by_stat = deseq.DEG(
        count_matrix_processed=count_matrix_processed, metadata=metadata, experiment=experiment, control=control, padj=padj, log2FoldChange=log2FoldChange, baseMean=baseMean)

    if perform_GSEA:
        go_ontologies, go_ontologies_df = deseq.GSEA(
            ranked_by_stat, filter_by_organ)

        vis = Plots(plot_volcano, plot_pca,
                    plot_heatmap, plot_gsea, plot_dotplot, config)
        vis.run_plots(dds, dds_stats_df,
                      significant_genes, ensg_to_gene_mapper, go_ontologies)
    else:
        vis = Plots(plot_volcano, plot_pca,
                    plot_heatmap, plot_gsea=False, plot_dotplot=False, config=config)
        vis.run_plots(dds, dds_stats_df,
                      significant_genes, ensg_to_gene_mapper)

    logger.info("Pipeline completed.")


if __name__ == "__main__":
    run_pipeline(CONFIG_PATH=r"C:\Users\rafaelo\OneDrive - NTNU\Documents\Projects\PyDESeq\src\scripts\config.json",
                 diff_exp_matrix_path=r"C:\Users\rafaelo\OneDrive - NTNU\Documents\Projects\BIOMOL8008\Voom_diffexp_gnms.txt")
