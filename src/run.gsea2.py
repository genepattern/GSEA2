import os
import sys
import subprocess
from optparse import OptionParser
from datetime import datetime
from zipfile import ZipFile
from os.path import basename
import argparse
import shutil
import json
import random
import pandas
import numpy
import dominate
from dominate.tags import *
from dominate.util import raw
from scipy.integrate import simps


# Better boolean command line parsing
def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def main():
    usage = "%prog [options]" + "\n"
    ap = argparse.ArgumentParser()
    ap.add_argument("--libdir", action="store",
                    dest="libdir", help="Working directory to load support library from.")
    ap.add_argument("--dataset", action="store", dest="dataset",
                    help="Input Expression Dataset.")
    ap.add_argument("--gsdb", action="store", dest="gsdb",
                    help="Gene Set Database File.")
    ap.add_argument("--nperm", action="store", dest="nperm",
                    default=1000, type=int, help="Number of permutations.")
    ap.add_argument("--cls", action="store", dest="cls", help="CLS file.")
    ap.add_argument("--reverse", action="store", type=str2bool, nargs='?', const=True, dest="reverse",
                    default=False, help="Reverse the phenotype comparison defined in the CLS file.")
    ap.add_argument("--perm", action="store", dest="perm", default="sample",
                    help="Permutation mode. Options are 'sample' (phenotype) and 'set' (gene set).")
    ap.add_argument("--collapse", action="store", dest="collapse", default="none",
                    help="Method for computing mathematical collapse. Supports 'none' ('no collapse'), 'sum', 'mean', 'median', 'max', 'absmax'")
    ap.add_argument("--chip", action="store", dest="chip",
                    default="none", help="Chip file used for performing collapse.")
    ap.add_argument("--metric", action="store",
                    dest="rank_metric", help="Metric for ranking genes.")
    ap.add_argument("--alg", action="store", dest="method",
                    help="Enrichment Method. 'ks' (Kolmogorov Smirnov, old GSEA), 'ksa' (Kolmogorov Smirnov area), and 'cidac' (cumulative information divergence with antisymmetricity and complementation, next generation GSEA) supported.")
    ap.add_argument("--exponent", action="store", dest="exponent", default=1.0,
                    type=float, help="Weight for ks or ksa enrichment method.")
    ap.add_argument("--max", action="store", dest="max",
                                    default=500, type=int, help="Max gene set size.")
    ap.add_argument("--min", action="store", dest="min",
                                    default=15, type=int, help="Min gene set size.")
    ap.add_argument("--seed", action="store", dest="seed",
                    default="timestamp", help="Random seed used for permutations.")
    ap.add_argument("--ogllv", action="store", type=str2bool, nargs='?', const=True, dest="override",
                    default=False, help="Override reasonableness check for input dataset gene list size.")
    ap.add_argument("--nplot", action="store", dest="nplot", default=25,
                    type=int, help="Number of enrichment results to plot.")
    ap.add_argument("--zip", action="store", type=str2bool, nargs='?', const=True,
                                    dest="zip", default=True, help="Create ZIP bundle of results.")
    ap.add_argument("--cpu", action="store", dest="cpu",
                                    default=1, type=int, help="Job CPU Count.")
    options = ap.parse_args()

    sys.path.insert(1, options.libdir)
    import GSEAlib

    # Make a directory to store processed input files
    os.mkdir("input")

    # Generate and set the random seed at the Python level and save it to pass to GSEA
    if options.seed == "timestamp":
        options.seed = int(round(datetime.now().timestamp()))
        random.seed(options.seed)
    else:
        options.seed = int(round(float(options.seed)))
        random.seed(options.seed)

    # Parse GCT file
    if options.dataset.split(".")[-1] == "gct":
        if options.collapse != "none":
            chip_file = GSEAlib.read_chip(options.chip)
            input_ds = GSEAlib.collapse_dataset(
                options.dataset, chip_file, method=options.collapse, drop=True)
            input_ds['mappings'].to_csv(
                'input/collapse_dataset_mapping_details.tsv', sep="\t", na_rep="No Symbol Mapping")
            collapse_length = input_ds['collapse_length']
        else:
            input_ds = GSEAlib.read_gct(options.dataset)
        input_length = input_ds['input_length']
        input_ds = input_ds['data']
    # elif options.dataset.split(".")[-1] == "rnk":
    #     input_ds = pandas.read_csv(
    #         options.dataset, sep='\t', index_col=0, skip_blank_lines=True, header=None)
    #     if any(input_ds.index.str.startswith('#')):
    #         input_ds = input_ds.rename(columns=input_ds.iloc[int(
    #             numpy.where(input_ds.index.str.startswith('#'))[0])])
    #         input_ds = input_ds[input_ds.index.str.startswith('#') != True]
    #     else:
    #         input_ds = input_ds.rename(columns={1: "Preranked Metric"})
    #     input_ds.index.name = "Name"
    #     input_length = len(input_ds.index)
    #     if options.collapse != "none":
    #         chip_file = GSEAlib.read_chip(options.chip)
    #         input_ds = GSEAlib.collapse_dataset(
    #             input_ds, chip_file, method=options.collapse, drop=True)
    #         input_ds['mappings'].to_csv(
    #             'input/collapse_dataset_mapping_details.tsv', sep="\t", na_rep="No Symbol Mapping")
    #         input_length = input_ds['input_length']
    #         collapse_length = input_ds['collapse_length']
    #         input_ds = input_ds['data']
    else:
        input_ds = pandas.read_csv(
            options.dataset, sep='\t', index_col=0, skip_blank_lines=True)
        if "description" in input_ds.columns.str.lower():
            description_loc = input_ds.columns.str.lower().to_list().index('description')
            input_ds.drop(
                input_ds.columns[[description_loc]], axis=1, inplace=True)
        input_ds.index.name = "Name"
        input_length = len(input_ds.index)
        if options.collapse != "none":
            chip_file = GSEAlib.read_chip(options.chip)
            input_ds = GSEAlib.collapse_dataset(
                input_ds, chip_file, method=options.collapse, drop=True)
            input_ds['mappings'].to_csv(
                'input/collapse_dataset_mapping_details.tsv', sep="\t", na_rep="No Symbol Mapping")
            input_length = input_ds['input_length']
            collapse_length = input_ds['collapse_length']
            input_ds = input_ds['data']

    if len(input_ds) < 10000 and options.override == False and options.collapse == "none":
        sys.exit(print("Only ", len(input_ds), "genes were identified in the dataset.\nEither the dataset did not contain all expressed genes, or collapse dataset may need to be run with an appropriate chip file.\n\nIf this was intentional, to bypass this check you can set 'override gene list length validation' (--ogllv) to 'True' but this is not recommended."))
    if len(input_ds) < 10000 and options.override == False and options.collapse != "none":
        sys.exit(print("Only ", len(input_ds), "genes were identified in the dataset.\nEither the dataset did not contain all expressed genes, or there was possibly a problem with the chip selected for collapse dataset.\n\nIf this was intentional, to bypass this check you can set 'override gene list length validation' (--ogllv) to 'True' but this is not recommended."))
    if len(input_ds) < 10000 and options.override == True:
        print("Only", len(input_ds), "genes were identified in the dataset, but the user specified overriding this check. Continuing analysis, as-is however this is not recommended. The input dataset should include all expressed genes.")

    # Parse CLS file
    labels, phenotypes = GSEAlib.read_cls(options.cls)
    phenotypes = GSEAlib.match_phenotypes(input_ds, phenotypes)
    if options.reverse == True and phenotypes.columns[0] == "Labels":
        phenotypes["Phenotypes"] = numpy.where((phenotypes["Phenotypes"] == 0) | (
            phenotypes["Phenotypes"] == 1), phenotypes["Phenotypes"] ^ 1, phenotypes["Phenotypes"])
        labels = {0: labels[1], 1: labels[0]}
    phenotypes = phenotypes.sort_values('Phenotypes')

    # Order the dataset using the phenotypes and write out both files
    input_ds = input_ds.reindex(columns=phenotypes.index)
    input_ds.to_csv('input/gene_by_sample.tsv', sep="\t")
    tbs_df = pandas.DataFrame(phenotypes['Phenotypes']).transpose()
    tbs_df.index.name = "Target"
    tbs_df = tbs_df.rename(index={tbs_df.index[0]: labels[0]})
    tbs_df.to_csv('input/target_by_sample.tsv', sep="\t", index=True)

    # Parse GMT/GMX gene sets files from a list of inputs and create a name:members dict written out as a json file
    if options.gsdb != None:
        with open(options.gsdb) as f:
            gene_sets_dbfile_list = f.read().splitlines()

    gs_data = GSEAlib.read_sets(gene_sets_dbfile_list)
    genesets = gs_data['genesets']
    genesets_descr = gs_data['descriptions']
    with open('input/raw_set_to_genes.json', 'w') as path:
        json.dump(genesets, path,  indent=2)

    # Filter gene sets to just genes in input dataset
    gs_data_subset = GSEAlib.filter_sets(genesets, input_ds.index)
    gs_data_subset_sets = gs_data_subset['genesets']
    gs_data_subset_lengths = gs_data_subset['lengths']
    passing_lengths = dict((key, value) for key, value in gs_data_subset_lengths.items(
    ) if (value >= max(options.min, 1) and value <= options.max))
    passing_sets = {key: gs_data_subset_sets[key]
                    for key in passing_lengths.keys()}
    with open('input/filtered_set_to_genes.json', 'w') as path:
        json.dump(genesets, path,  indent=2)

    # Construct GSEA Settings json file
    gsea_settings = {
        "number_of_permutations": options.nperm,
        "permutation": options.perm,
        "feature_name": "Features",
        "metric": options.rank_metric,
        "score_name": options.rank_metric,
        "algorithm": options.method,
        "exponent": options.exponent,
        "maximum_gene_set_size": options.max,
        "minimum_gene_set_size": options.min,
        "remove_gene_set_genes": True,
        "random_seed": options.seed,
        "number_of_jobs": options.cpu,
        "number_of_extreme_gene_sets_to_plot": options.nplot,
        "gene_sets_to_plot": []
    }

    with open('input/gsea_settings.json', 'w') as path:
        json.dump(gsea_settings, path,  indent=2)

    # Run GSEA
    subprocess.check_output(['gsea', 'metric-rank', 'input/gsea_settings.json', 'input/target_by_sample.tsv',
                             'input/gene_by_sample.tsv', 'input/filtered_set_to_genes.json', os.getcwd()])

    # Parse Results
    genesets_descr = pandas.DataFrame.from_dict(
        genesets_descr, orient="index", columns=["URL"])
    results = GSEAlib.result_paths(os.getcwd())
    plots = [result for result in results if "plot" in result]
    gsea_stats = pandas.read_csv(
        'set_x_statistic_x_number.tsv', sep="\t", index_col=0)
    ranked_genes = pandas.read_csv(
        'gene_x_metric_x_score.tsv', sep="\t", index_col=0)
    random_es_distribution = pandas.read_csv(
        'set_x_random_x_enrichment.tsv', sep="\t", index_col=0)

    # Add set sizes to enrichment report
    gsea_stats.insert(0, 'Size', '')
    for gs in range(len(gsea_stats)):
        gsea_stats.loc[gsea_stats.index[gs],
                       'Size'] = passing_lengths[gsea_stats.index[gs]]
    gsea_stats.to_csv(
        'set_x_statistic_x_number.tsv', sep="\t")

    # Positive Enrichment Report
    gsea_pos = gsea_stats[gsea_stats.loc[:, "Enrichment"] > 0]
    if len(gsea_pos) > 0:
        gsea_pos = genesets_descr.merge(gsea_pos, how='inner', left_index=True, right_index=True).sort_values(
            ["Enrichment"], axis=0, ascending=(False)).reset_index()
        gsea_pos.insert(1, 'Details', '')
        for gs in range(len(gsea_pos)):
            # Compute original set size, filtered set, and filtered size
            unfiltered_len = len(genesets[gsea_pos.iloc[gs]['index']])
            filtered_gs = list(set(genesets[gsea_pos.iloc[gs]['index']]) & set(
                list(ranked_genes.index.values)))
            filtered_len = len(filtered_gs)
            if "plot/" + gsea_pos.iloc[gs]['index'].lower() + ".html" in plots:
                # Get Enrichment Statistics for the set of interest
                report_set = pandas.DataFrame(gsea_pos.iloc[gs]).copy(
                    deep=True)
                report_set.rename({'index': 'Gene Set'}, axis=0, inplace=True)
                set_enrichment_score = report_set.loc['Enrichment'].values[0]
                # Only do plotting work if we need to
                heatmap_fig = GSEAlib.plot_set_heatmap(
                    input_ds, phenotypes, ranked_genes, filtered_gs, ascending=True)
                null_es_fig = GSEAlib.set_perm_indepkde_displot(
                    random_es_distribution.loc[gsea_pos.iloc[gs]['index']], set_enrichment_score)
                # Edit in the needed information to the per-set enrichment reports
                report_set.loc["Details"] = "Dataset: " + os.path.splitext(os.path.basename(options.dataset))[
                    0] + "<br>Enriched in Phenotype: \"" + str(labels[0]) + "\" of comparison " + str(labels[0]) + " vs " + str(labels[1])
                page = open(
                    "plot/" + gsea_pos.iloc[gs]['index'].lower() + ".html", 'r')
                page_str = page.read()
                leading_edge_table, leading_edge_subset = GSEAlib.get_leading_edge(
                    page_str)
                doc = dominate.document(title=gsea_pos.iloc[gs]['index'])
                doc += h3("Enrichment Details")
                doc += raw(report_set.to_html(header=False,
                                              render_links=True, escape=False, justify='left'))
                doc += raw("<br>")
                doc += h3("Enrichment Plot")
                doc += raw(page_str.replace("<!doctype html>", ""))
                doc += raw("<br>")
                doc += h3("Table: GSEA details")
                doc += raw(leading_edge_table.to_html())
                doc += p('Investigate core enrichment with ', a("MSigDB Webtools",
                                                                href='https://www.gsea-msigdb.org/gsea/msigdb/annotate.jsp?geneIdList=' + leading_edge_subset, target='_blank'), " or ", a("Query NDEx", href='https://www.ndexbio.org/iquery/?genes=' + leading_edge_subset, target='_blank'))
                doc += raw("<br>")
                doc += h3("Row Normalized Expression Heatmap for " +
                          gsea_pos.iloc[gs]['index'])  # add a title for the heatmap
                doc += raw(heatmap_fig)
                doc += raw("<br>")
                doc += h3("Random Enrichment Score Distribution for " +
                          gsea_pos.iloc[gs]['index'])  # add a title for the ES distplot
                doc += raw(null_es_fig)
                with open("plot/" + gsea_pos.iloc[gs]['index'].lower() + ".html", 'w') as f:
                    f.write(doc.render())
                # HTMLify the positive report
                gsea_pos.at[gs, "Details"] = "<a href=plot/" + \
                    gsea_pos.iloc[gs]['index'].lower(
                ) + ".html target='_blank'>Details...</a>"
        gsea_pos["index"] = gsea_pos.apply(
            lambda row: "<a href='{}' target='_blank'>{}</a>".format(row.URL, row['index']), axis=1)
        gsea_pos.drop("URL", axis=1, inplace=True)
        gsea_pos = gsea_pos.rename(
            columns={'index': 'Gene Set<br>follow link to MSigDB'})
        gsea_pos.index += 1
    gsea_pos.to_html(open('gsea_report_for_positive_enrichment.html', 'w'),
                     render_links=True, escape=False, justify='center')

    # Negative Enrichment Report
    gsea_neg = gsea_stats[gsea_stats.loc[:, "Enrichment"] < 0]
    if len(gsea_neg) > 0:
        gsea_neg = genesets_descr.merge(gsea_neg, how='inner', left_index=True, right_index=True).sort_values(
            ["Enrichment"], axis=0, ascending=(True)).reset_index()
        gsea_neg.insert(1, 'Details', '')
        for gs in range(len(gsea_neg)):
            # Compute original set size, filtered set, and filtered size
            unfiltered_len = len(genesets[gsea_neg.iloc[gs]['index']])
            filtered_gs = list(set(genesets[gsea_neg.iloc[gs]['index']]) & set(
                list(ranked_genes.index.values)))
            filtered_len = len(filtered_gs)
            if "plot/" + gsea_neg.iloc[gs]['index'].lower() + ".html" in plots:
                # Get Enrichment Statistics for the set of interest
                report_set = pandas.DataFrame(gsea_neg.iloc[gs]).copy(
                    deep=True)
                report_set.rename({'index': 'Gene Set'}, axis=0, inplace=True)
                set_enrichment_score = report_set.loc['Enrichment'].values[0]
                # Only do plotting work if we need to
                heatmap_fig = GSEAlib.plot_set_heatmap(
                    input_ds, phenotypes, ranked_genes, filtered_gs, ascending=False)
                null_es_fig = GSEAlib.set_perm_indepkde_displot(
                    random_es_distribution.loc[gsea_neg.iloc[gs]['index']], set_enrichment_score)
                # Edit in the needed information to the per-set enrichment reports
                report_set.loc["Details"] = "Dataset: " + os.path.splitext(os.path.basename(options.dataset))[
                    0] + "<br>Enriched in Phenotype: \"" + str(labels[1]) + "\" of comparison " + str(labels[0]) + " vs " + str(labels[1])
                page = open(
                    "plot/" + gsea_neg.iloc[gs]['index'].lower() + ".html", 'r')
                page_str = page.read()
                leading_edge_table, leading_edge_subset = GSEAlib.get_leading_edge(
                    page_str)
                doc = dominate.document(title=gsea_neg.iloc[gs]['index'])
                doc += h3("Enrichment Details")
                doc += raw(report_set.to_html(header=False,
                                              render_links=True, escape=False, justify='left'))
                doc += raw("<br>")
                doc += h3("Enrichment Plot")
                doc += raw(page_str.replace("<!doctype html>", ""))
                doc += raw("<br>")
                doc += h3("Table: GSEA details")
                doc += raw(leading_edge_table.to_html())
                doc += p('Investigate core enrichment with ', a("MSigDB Webtools",
                                                                href='https://www.gsea-msigdb.org/gsea/msigdb/annotate.jsp?geneIdList=' + leading_edge_subset, target='_blank'), " or ", a("Query NDEx", href='https://www.ndexbio.org/iquery/?genes=' + leading_edge_subset, target='_blank'))
                doc += raw("<br>")
                doc += h3("Row Normalized Expression Heatmap for " +
                          gsea_neg.iloc[gs]['index'])  # add a title for the heatmap
                doc += raw(heatmap_fig)
                doc += raw("<br>")
                doc += h3("Random Enrichment Score Distribution for " +
                          gsea_neg.iloc[gs]['index'])  # add a title for the ES distplot
                doc += raw(null_es_fig)
                with open("plot/" + gsea_neg.iloc[gs]['index'].lower() + ".html", 'w') as f:
                    f.write(doc.render())
                # HTMLify the negative report
                gsea_neg.at[gs, "Details"] = "<a href=plot/" + \
                    gsea_neg.iloc[gs]['index'].lower(
                ) + ".html target='_blank'>Details...</a>"
        gsea_neg["index"] = gsea_neg.apply(
            lambda row: "<a href='{}' target='_blank'>{}</a>".format(row.URL, row['index']), axis=1)
        gsea_neg.drop("URL", axis=1, inplace=True)
        gsea_neg = gsea_neg.rename(
            columns={'index': 'Gene Set<br>follow link to MSigDB'})
        gsea_neg.index += 1
    gsea_neg.to_html(open('gsea_report_for_negative_enrichment.html',
                          'w'), render_links=True, escape=False, justify='center')

    # Create Report for Dataset top markers
    dataset_markers = numpy.append(ranked_genes.sort_values(ranked_genes.columns[0], ascending=False).index.values[0:50], ranked_genes.sort_values(
        ranked_genes.columns[0], ascending=False).index.values[len(ranked_genes) - 50:len(ranked_genes)])
    heatmap_fig = GSEAlib.plot_set_heatmap(
        input_ds, phenotypes, ranked_genes, list(dataset_markers), ascending=True)
    corr_plot_fig = GSEAlib.plot_gene_rankings(ranked_genes, labels)
    global_es_distplot_fig = GSEAlib.global_es_indepkde_distplot(
        gsea_stats['Enrichment'])
    doc = dominate.document(title="Heat map and correlation plot for " +
                            os.path.splitext(os.path.basename(options.dataset))[0])
    doc += h3("Row Normalized Expression Heatmap for the top 50 features for each phenotype in " +
              os.path.splitext(os.path.basename(options.dataset))[0])  # add a title for the heatmap
    doc += raw(heatmap_fig)
    doc += raw("<br>")
    doc += h3("Ranked Gene List Correlation Profile")
    doc += raw(corr_plot_fig)
    doc += raw("<br>")
    doc += h3("Global Enrichment Score Distribution")
    doc += raw(global_es_distplot_fig)
    with open("heat_map_corr_plot.html", 'w') as f:
        f.write(doc.render())

    # Create Report Index using dominate package
    gsea_index = dominate.document(
        title="GSEA Report for Dataset " + os.path.splitext(os.path.basename(options.dataset))[0])
    gsea_index += h1("GSEA Report for Dataset " +
                     os.path.splitext(os.path.basename(options.dataset))[0])
    gsea_index += h2(str(labels[0]) + " vs. " + str(labels[1]))
    gsea_index += h3("Enrichment in phenotype: " + str(
        labels[0]) + " (" + str(sum(phenotypes['Phenotypes'] == 0)) + " samples)")
    gsea_index += ul(
        li(str(len(gsea_stats[(gsea_stats['Enrichment'] >= 0)])) + " / " + str(
            len(gsea_stats)) + " gene sets are upregulated in phenotype ",  b(str(labels[0]))),
        li(str(len(gsea_stats[(gsea_stats['Enrichment'] > 0) & (
            gsea_stats['Adjusted global p value'] < 0.25)])) + " gene sets are significant at adjusted global pValue (FDR) < 25%"),
        li(str(len(gsea_stats[(gsea_stats['Enrichment'] > 0) & (
            gsea_stats['Global p value'] < 0.01)])) + " gene sets are significantly enriched at pValue < 1%"),
        li(str(len(gsea_stats[(gsea_stats['Enrichment'] > 0) & (
            gsea_stats['Global p value'] < 0.05)])) + " gene sets are significantly enriched at pValue < 5%"),
        li(a("Detailed enrichment results in html format",
             href="gsea_report_for_positive_enrichment.html", target='_blank')),
        li(a("Guide to interpret results",
             href='http://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideFrame.html?_Interpreting_GSEA_Results', target='_blank'))
    )
    gsea_index += h3("Enrichment in phenotype: " + str(
        labels[1]) + " (" + str(sum(phenotypes['Phenotypes'] == 1)) + " samples)")
    gsea_index += ul(
        li(str(len(gsea_stats[(gsea_stats['Enrichment'] < 0)])) + " / " + str(
            len(gsea_stats)) + " gene sets are upregulated in phenotype ", b(str(labels[1]))),
        li(str(len(gsea_stats[(gsea_stats['Enrichment'] < 0) & (
            gsea_stats['Adjusted global p value'] < 0.25)])) + " gene sets are significant at adjusted global pValue (FDR) < 25%"),
        li(str(len(gsea_stats[(gsea_stats['Enrichment'] < 0) & (
            gsea_stats['Global p value'] < 0.01)])) + " gene sets are significantly enriched at pValue < 1%"),
        li(str(len(gsea_stats[(gsea_stats['Enrichment'] < 0) & (
            gsea_stats['Global p value'] < 0.05)])) + " gene sets are significantly enriched at pValue < 5%"),
        li(a("Detailed enrichment results in html format",
             href="gsea_report_for_negative_enrichment.html", target='_blank')),
        li(a("Guide to interpret results",
             href='http://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideFrame.html?_Interpreting_GSEA_Results', target='_blank'))
    )
    gsea_index += h3("Dataset details")
    if options.collapse != "none":
        gsea_index += ul(
            li("The Dataset has " + str(input_length) + " native features"),
            li("After collapsing features into gene symbols, there are: " +
               str(collapse_length) + " genes"),
            li("Collapse method: \"" + options.collapse + "\" was used to collapse features to gene symbols"))
    else:
        gsea_index += ul(
            li("The Dataset has " + str(input_length) + " features (genes)"),
            li("No probe set => gene symbol collapsing was requested, so all " +
               str(input_length) + " features were used"))
    gsea_index += h3("Gene set details")
    gsea_index += ul(
        li("Gene set size filters (min=" + str(options.min) + ", max=" + str(options.max) + ") resulted in filtering out " +
           str(len(genesets) - len(passing_sets)) + " / " + str(len(genesets)) + " gene sets"),
        li("The remaining " + str(len(passing_sets)) +
           " gene sets were used in the analysis")
    )
    gsea_index += h3("Gene markers for the " +
                     str(labels[0]) + " vs. " + str(labels[1]) + " comparison")
    gsea_index += ul(
        li("The dataset has " + str(len(ranked_genes)) + " features (genes)"),
        li("# of markers for phenotype " + str(labels[0]) + ": " + str(numpy.count_nonzero(ranked_genes.iloc[:, 0].values > 0)) + " (" + str(round(
            numpy.count_nonzero(ranked_genes.iloc[:, 0].values > 0) / len(ranked_genes) * 100, 1)) + "%) with correlation area " + str(round(simps(abs(ranked_genes.iloc[:, 0].values[ranked_genes.iloc[:, 0].values > 0])) / simps(abs(ranked_genes.iloc[:, 0].values)) * 100, 1)) + "%"),
        li("# of markers for phenotype " + str(labels[0]) + ": " + str(numpy.count_nonzero(ranked_genes.iloc[:, 0].values < 0)) + " (" + str(round(
            numpy.count_nonzero(ranked_genes.iloc[:, 0].values < 0) / len(ranked_genes) * 100, 1)) + "%) with correlation area " + str(round(simps(abs(ranked_genes.iloc[:, 0].values[ranked_genes.iloc[:, 0].values < 0])) / simps(abs(ranked_genes.iloc[:, 0].values)) * 100, 1)) + "%"),
        li(a("Detailed rank ordered gene list for all features in the dataset (.tsv file)",
             href='gene_x_metric_x_score.tsv')),
        li(a("Heat map and gene list correlation profile for all features in the dataset",
             href='heat_map_corr_plot.html'))
    )
    gsea_index += h3("Reproducibility")
    gsea_index += ul(
        li("Random seed used for permutation generation: " + str(options.seed)),
        li(a("Parameters passed to GSEA.jl (.json file)",
             href='input/gsea_settings.json'))
    )
    gsea_index += h3("Citing GSEA and MSigDB")
    gsea_index += p('To cite your use of the GSEA software please reference the following:')
    gsea_index += ul(
        li(a("Subramanian, A., Tamayo, P., et al. (2005, PNAS).",
             href='https://www.pnas.org/content/102/43/15545', target='_blank')),
        li(a("Mootha, V. K., Lindgren, C. M., et al. (2003, Nature Genetics).",
             href='http://www.nature.com/ng/journal/v34/n3/abs/ng1180.html', target='_blank'))
    )
    gsea_index += p('For use of the Molecular Signatures Database (MSigDB), to cite please reference one or more of the following',
                    br(), 'as appropriate, along with the source for the gene set as listed on the gene set page:')
    gsea_index += ul(
        li(a("Liberzon A, et al. (Bioinformatics, 2011).",
             href='https://doi.org/10.1093/bioinformatics/btr260', target='_blank')),
        li(a("Liberzon A, et al. (Cell Systems 2015).",
             href='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4707969/', target='_blank'))
    )

    with open('index.html', 'w') as f:
        f.write(gsea_index.render())

    # Zip up results
    if options.zip == True:
        gsea_files = []
        for folderName, subfolders, filenames in os.walk(os.path.relpath(os.getcwd())):
            for filename in filenames:
                # create complete filepath of file in directory
                filePath = os.path.join(folderName, filename)
                # Add file to zip
                gsea_files.append(filePath)
        with ZipFile("gsea_results.zip", "w") as gsea_zip:
            for filename in gsea_files:
                gsea_zip.write(filename)


if __name__ == '__main__':
    main()
