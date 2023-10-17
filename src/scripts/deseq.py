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


class DESeq:

    __slots__ = ["gene_sets", "gsea_permutations",
                 "config", "ensg_to_gene_mapper"]

    def __init__(self, gene_sets: dict, gsea_permutations: int, config: dict, ensg_to_gene_mapper):

        self.gene_sets = gene_sets
        self.gsea_permutations = gsea_permutations
        self.config = config
        self.ensg_to_gene_mapper = ensg_to_gene_mapper

    def pre_processing(self, count_matrix, gene_id_col_name: str):
        """_summary_

        Args:
            count_matrix (_type_): _description_
            gene_id_col_name (str): _description_
            diff_exp_matrix (_type_): _description_

        Returns:
            _type_: _description_
        """
        logger.info("Pre-processing count matrix text file.")

        count_matrix_copy = count_matrix.copy()

        count_matrix_processed = pd.DataFrame(
            replace_colnames(count_matrix_copy))
        count_matrix_processed.set_index(gene_id_col_name, inplace=True)

        # Filter go_ontologies_list all rows that only contain 0 counts
        count_matrix_processed = count_matrix_processed[(
            count_matrix_processed.iloc[:, :] != 0).all(axis=1)]

        # Transpose so our genes become columns and samples rows
        count_matrix_processed = count_matrix_processed.transpose()

        # Add group or condition. In this case, IDs are ordered and so we can do this. Ideally, a sample would have a dictinary mapping their condition.
        labels = ["Smoker"]*5 + ["Non-Smoker"]*5

        # diff_exp_matrix_processed = diff_exp_matrix_copy[(diff_exp_matrix_copy["adj.P.Val"] < 0.05) & (
        #     abs(diff_exp_matrix_copy["logFC"]) > 0.5)]

        metadata = pd.DataFrame(zip(count_matrix_processed.index, labels), columns=[
            "Sample", "Condition"]).set_index("Sample")

        return count_matrix_processed,  metadata

    def DEG(self, count_matrix_processed, metadata, experiment: str, control: str, padj: float, log2FoldChange: float, baseMean: int):

        logger.info("Creating Deseq Dataset for PyDESeq2.")
        dds = DeseqDataSet(counts=count_matrix_processed,
                           metadata=metadata, design_factors="Condition")
        dds.deseq2()

        logger.info(
            f"Generating differential expression statistics with Experimental : {experiment} vs Control : {control}.")
        stat_res = DeseqStats(dds, n_cpus=4, contrast=(
            "Condition", experiment, control))
        stat_res.summary()
        # Create dataframe with DESeq statistics
        dds_stats_df = stat_res.results_df

        # Do DEG
        # sc.tl.rank_genes_groups(dds, groupby="Condition",
        #                         method="t-test_overestim_var", n_genes=10000)
        # diff_exp_genes = sc.get.rank_genes_groups_df(dds, group=None)

        # sc.pl.rank_genes_groups(dds, n_genes=20, sharey=False)

        dds_stats_df["Symbols"] = dds_stats_df.index.map(
            self.ensg_to_gene_mapper.mapper)
        # Drop any symols that weren't recognized as genes (most likely are pseudo genes)
        dds_stats_df.dropna(subset=["Symbols"], inplace=True)

        logger.info("Extracting significantly expressed genes.")
        # Remove genes that are lowly expressed and rank the significant genes
        significant_genes = dds_stats_df[(dds_stats_df["padj"] < padj) & (
            abs(dds_stats_df["log2FoldChange"]) > log2FoldChange) & (dds_stats_df.baseMean >= baseMean)]
        ranked_by_stat = dds_stats_df[["Symbols", "stat"]].dropna().sort_values(
            "stat", ascending=False).drop_duplicates("Symbols")

        logger.info(
            f"\nShape of dataframe with every gene before ranking and filtering: {dds_stats_df.shape}.\nShape of dataframe with every gene after ranking by Stat: {ranked_by_stat.shape}.\nShape of significant genes (after filtering using the following threshold parameters -> p-value: {padj} || log2FC: {log2FoldChange} || baseMean: {baseMean}): {significant_genes.shape} ")

        # Save ranked significant genes to excel
        logger.info(
            "Saving significant genes and raneked genes by stat to Excel.")
        significant_genes.to_excel("Significant_genes_{}.xlsx".format(date))
        ranked_by_stat.to_excel("Ranked_genes_by_stat_{}.xlsx".format(date))

        return dds, stat_res, dds_stats_df, significant_genes, ranked_by_stat

    def GSEA(self, ranked_by_stat: pd.DataFrame, filter_by_organ: str):
        """_summary_

        Args:
            ranked_by_stat (pd.DataFrame): _description_

        Returns:
            _type_: _description_
        """
        # Do GSEA
        # Add GO ontologies
        logger.info(
            "Performing Gene Set Enrichment Analysis (GSEA)")
        gene_sets_list = [key for key,
                          value in self.gene_sets.items() if value]

        go_ontologies = gp.prerank(
            rnk=ranked_by_stat, gene_sets=gene_sets_list, seed=6, permutation_num=self.gsea_permutations)

        go_ontologies_list = []

        # If we want to filter for those that contain lung
        if filter_by_organ != "":
            logger.warning(
                f"You're currently filtering for the {filter_by_organ} organ. ")
            go_ontologies_dict = {
                k: v for k, v in go_ontologies.results.items() if filter_by_organ in k}
            updated_list_of_terms = [el for el in go_ontologies_dict.keys()]
            go_ontologies.results.clear()  # Clear the existing dictionary
            # Update with the filtered dictionary
            go_ontologies.results.update(go_ontologies_dict)

            # Update the res2d attribute to make the gsea plot with just lung
            go_ontologies.res2d = go_ontologies.res2d.loc[go_ontologies.res2d.Term.isin(
                updated_list_of_terms)]

        for term in list(go_ontologies.results):
            go_ontologies_list.append([term, go_ontologies.results[term]["fdr"],
                                       go_ontologies.results[term]["es"], go_ontologies.results[term]["nes"], go_ontologies.results[term]["tag %"], go_ontologies.results[term]["gene %"], go_ontologies.results[term]["lead_genes"], go_ontologies.results[term]["matched_genes"]])

        # Create a version of the data in the list as a dataframe
        go_ontologies_df = pd.DataFrame(go_ontologies_list, columns=["Term", "ffr", "es", "nes", "tag%", "gene%", "lead_genes", "matched_genes"]).sort_values(
            "ffr").reset_index(drop=True)  # ffr = p-value, es = enrichment score, nes = normalized enrichment score
        # We can now sort via nes to get the most negative or positively enriched
        go_ontologies_df.sort_values("nes")
        logger.info(
            f"Saving GEAS file with the following database(s): {[k for k,v in  self.gene_sets.items() if v]}.")
        go_ontologies_df.to_excel("GEAS_{}.xlsx".format(date))
        return go_ontologies,  go_ontologies_df
