import os
import sys
import ast
import pandas
import numpy
import json
import math
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import plotly.express as px
# import plotly.figure_factory as ff
from scipy.stats import gaussian_kde


# Simple implementation of a GCT parser
# Accepts a GCT file and returns a Pandas Dataframe with a single index
def read_gct(gct):
    dataset = pandas.read_csv(gct, sep='\t', header=2, index_col=[
        0, 1], skip_blank_lines=True)
    dataset.index.names = ["Name", "Description"]
    dataset_descriptions = dataset.index.to_frame(index=False)
    dataset_descriptions.set_index(["Name"], inplace=True)
    dataset.index = dataset.index.droplevel(1)  # Drop gene descriptions
    return {'data': dataset, 'row_descriptions': dataset_descriptions["Description"].values, 'input_length': len(dataset.index)}


# Simple implementation of a CHIP Parser for use with ssGSEA
# Reads in a CHIP formatted file and returns a pandas dataframe containing
# the probe to gene mappings
def read_chip(chip):
    chip_df = pandas.read_csv(
        chip, sep='\t', index_col=0, skip_blank_lines=True)
    return chip_df


# Simple implementation of GSEA DEsktop's Collapse Dataset functions for use
# with ssSGEA
# Accepts an expression dataset in GCT format, a CHIP file, and a
# collapse metric and returns a pandas dataframe formatted version of the
# dataset collapsed from probe level to gene level using the specified metric.
def collapse_dataset(dataset, chip, method="sum", drop=True):
    import pandas as pd
    if isinstance(dataset, pandas.DataFrame):
        dataset = dataset
    elif isinstance(dataset, dict) == False:
        dataset = read_gct(dataset)
    if isinstance(chip, pandas.DataFrame) == False:
        chip = read_chip(chip)
    if isinstance(dataset, dict) == True:
        dataset = dataset['data']
    input_len = len(dataset.index)
    joined_df = chip.join(dataset, how='right')
    joined_df.index.name = "Name"
    joined_df.reset_index(drop=False, inplace=True)
    mappings = joined_df[["Name",
                          "Gene Symbol"]].drop_duplicates().copy().sort_values('Gene Symbol').rename(columns={'Name': 'Dataset ID(s)'})  # Save mapping details for reporting
    joined_df = joined_df.drop("Name", axis=1).dropna(subset=['Gene Symbol'])
    annotations = joined_df[["Gene Symbol",
                             "Gene Title"]].drop_duplicates().copy()  # Save gene annotations for reporting
    joined_df.drop("Gene Title", axis=1, inplace=True)
    # Do Mathematical Collapse Operations
    if method.lower() == "sum":
        collapsed_df = joined_df.groupby(["Gene Symbol"]).sum()
    if method.lower() == "mean":
        collapsed_df = joined_df.groupby(["Gene Symbol"]).mean()
    if method.lower() == "median":
        collapsed_df = joined_df.groupby(["Gene Symbol"]).median()
    if method.lower() == "max":
        collapsed_df = joined_df.groupby(["Gene Symbol"]).max()
    if method.lower() == "absmax":
        collapsed_df = joined_df.loc[joined_df.groupby(
            ['Gene Symbol']).idxmax()]
    collapsed_df.index.name = "Name"
    # Group mapping details
    mappings = pandas.DataFrame(mappings.groupby(
        'Gene Symbol', dropna=False)['Dataset ID(s)'].apply(list))
    mappings['Dataset ID(s)'] = [','.join(map(str, l))
                                 for l in mappings['Dataset ID(s)']]
    return {'data': collapsed_df, 'row_descriptions': annotations["Gene Title"].values, 'mappings': mappings, 'input_length': input_len, 'collapse_length': len(collapsed_df.index)}


# Save a GCT result to a file, ensuring the filename has the extension .gct
def write_gct(gct, file_name, check_file_extension=True):
    if check_file_extension:
        file_name = check_extension(file_name, ".gct")
    rows = str(len(gct['data']))
    columns = str(len(gct['data'].columns))
    if len(gct['row_descriptions']) != int(rows):
        sys.exit("Number of row descriptions (" +
                 len(gct['row_descriptions']) + ") not equal to number of row names (" + rows + ").")
    row_descriptions = gct['row_descriptions']
    if row_descriptions == None:
        row_descriptions = ['NA'] * int(rows)
    m = gct['data'].copy()
    m.insert(loc=0, column='Description', value=gct['row_descriptions'])
    with open(file_name, 'w') as file:
        file.write('#1.2\n' + rows + '\t' + columns + '\n')
        m.to_csv(file, sep='\t', index_label="NAME", mode='w+')
    return(file_name)


# extension e.g. '.gct'
def check_extension(file_name, extension):
    import re
    ext = re.search(extension + "$", file_name.lower())
    if ext == None:
        file_name = file_name + extension
    return file_name


# Read CLS function adapted from https://github.com/broadinstitute/gsea_python/blob/ccal-refactor/gsea/Utils.py
def read_cls(path):
    """
    Returns a Pandas Series with phenotype labels from a .cls
            path (str): Local filepath for a .cls file
            returns (pd.Series): Pandas Series with phenotype labels
    """
    lines = open(path).readlines()
    if "numeric" in lines[0]:
        labels = {0: 'Pos', 1: 'Neg'}
        phens = lines[2].strip('\n').split()
        assignments = [lines[1].strip('\n').strip("#").split()[0]
                       for i in range(len(phens))]
        return labels, pandas.concat([pandas.Series(assignments, name='Numeric'), pandas.Series(phens, name='Phenotypes')], axis=1)
    else:  # Get the categorical labels assigned to the categorical phenotypes in the CLS as 'labels'
        labels = {
            label: i for i, label in enumerate(lines[1][1:-1].split())
        }
        try:
            # get the assignments of samples to an indicator label
            assignments = lines[2][:-1].split()
            # get the unique list of indicator labels used
            unique_assignments = pandas.Series(
                assignments).drop_duplicates().tolist()
            # create an alternate list of phentype labels from the sample assignment indicator labels
            alt_labels = {
                label: i for i, label in enumerate(unique_assignments)
            }
            if labels.keys == alt_labels.keys:  # If the labels defined in the cls match the sample indicator labels exactly exactly use them to substitute the numeric phentype assignments
                phens = [labels[lab] for lab in lines[2].strip('\n').split()]
            else:  # If the sample indicator labels don't match the phenotype labels exactly, assume the indicators are ordered according to the CLS spec and substitute the numeric phennotype assigments using the order
                phens = [alt_labels[lab]
                         for lab in lines[2].strip('\n').split()]
            labels = {value: key for (key, value) in labels.items()}
            return labels, pandas.concat([pandas.Series(assignments, name='Labels'), pandas.Series(phens, name='Phenotypes')], axis=1)
        except KeyError:  # Assume phenotype row is already ints
            phens = list(map(int, lines[2].strip('\n').split()))
            labels = {value: key for (key, value) in labels.items()}
            return labels, pandas.concat([pandas.Series(phens, name='Labels'), pandas.Series(phens, name='Phenotypes')], axis=1)


# Read Maps Phenotypes to Samples, adapted from https://github.com/broadinstitute/gsea_python/blob/ccal-refactor/gsea/Utils.py
def match_phenotypes(expr, phen):
    """
    Populates the index of phen with the column names of expr
            expr (pandas.DataFrame): DataFrame from read_gct
            phen (pandas.Series): Series from read_cls
            returns (pandas.Series): phen with index set to expr columns
    """
    if len(phen) == len(expr.columns):
        phen.index = expr.columns
    else:
        #		common = set(phen['phenotypes'].index) & set(expr.columns)
        #		expr = expr[list(common)]
        #		phen = phen[list(common)]
        sys.exit(
            "The number of samples in the CLS file did not match the number of samples in the dataset.")
    return phen


# Simple GMT/GMX to Dict parser
def read_sets(gene_sets_dbfile_list):
    genesets = {}
    genesets_descr = {}
    genesets_len = {}
    for gsdb in gene_sets_dbfile_list:
        gsdb_split = gsdb.split(".")
        if gsdb_split[-1] == "gmt":
            with open(gsdb) as f:
                temp = f.read().splitlines()
                for i in range(len(temp)):
                    gs_line = temp[i].split("\t")
                    gene_set_name = gs_line[0]
                    gene_set_desc = gs_line[1]
                    gene_set_tags = gs_line[2:len(gs_line)]
                    genesets[gene_set_name] = list(set(gene_set_tags))
                    # Not used yet but should end up in reports eventually
                    genesets_descr[gene_set_name] = gene_set_desc
        elif gsdb_split[-1] == "gmx":  # is a gmx formatted file
            df_temp = pandas.read_csv(
                gsdb, sep='\t', skip_blank_lines=True).transpose().dropna(how='all')
            for i in range(len(df_temp)):
                gs_line = df_temp.iloc[i][~df_temp.iloc[i].isnull()]
                gene_set_name = gs_line.name
                gene_set_desc = gs_line[0]
                gene_set_tags = gs_line[1:len(gs_line)]
                genesets[gene_set_name] = list(set(gene_set_tags))
                # Not used yet but should end up in reports eventually
                genesets_descr[gene_set_name] = gene_set_desc
        elif gsdb_split[-1] == "json":
            with open(gsdb) as f:
                temp = json.load(f)
            gene_set_tags = {key: temp[key]['geneSymbols']
                             for key in temp.keys()}
            genesets.update(gene_set_tags)
            gene_set_desc = {key: temp[key]['msigdbURL']
                             for key in temp.keys()}
            genesets_descr.update(gene_set_desc)
        else:
            sys.exit("The Gene Set Database format was not recognised.")
    genesets_len = {key: len(value) for key, value in genesets.items()}
    return {'genesets': genesets, 'descriptions': genesets_descr, 'lengths': genesets_len}


# Restrict Gene sets to input universe
def filter_sets(genesets_dict, dataset_index):
    genesets_filtered = {}
    genesets_len_filtered = {}
    for (key, value) in genesets_dict.items():
        genesets_filtered[key] = list(
            set(genesets_dict[key]) & set(list(dataset_index.values)))
    genesets_len_filtered = {key: len(value)
                             for key, value in genesets_filtered.items()}
    return {'genesets': genesets_filtered, 'lengths': genesets_len_filtered}


# Get file paths
def result_paths(root_dir):
    file_set = set()
    for dir_, _, files in os.walk(root_dir):
        for file_name in files:
            rel_dir = os.path.relpath(dir_, root_dir)
            rel_file = os.path.join(rel_dir, file_name)
            file_set.add(rel_file)
    return list(file_set)


# Plot Heatmaps for a given set of genes
def plot_set_heatmap(input_ds, phenotypes, ranked_genes, filtered_gs, ascending):
    ranking_colorscale = go.Heatmap(z=ranked_genes, colorscale='RdBu_r', zmid=0)[
        'colorscale']  # Create a colorscale for the ranking subset
    filtered_len = len(filtered_gs)
    ranked_gs_genes = ranked_genes.loc[filtered_gs].sort_values(
        ranked_genes.columns[0], ascending=ascending)
    gs_expression = input_ds.loc[ranked_gs_genes.index].copy()
    gs_expression_norm = gs_expression.subtract(gs_expression.min(axis=1), axis=0)\
        .divide(gs_expression.max(axis=1) - gs_expression.min(axis=1), axis=0)\
        .combine_first(gs_expression)  # Row Normalize Gene Set gene expression
    # Construct plotly heatmap
    # Instantiate a plot containing slots for the main heatmap and a slot for the phenotype label bar
    layout = [[{}, {}], [{"rowspan": len(gs_expression_norm), "colspan": 1},
                         {"rowspan": len(gs_expression_norm), "colspan": 1}]]
    layout.extend([[None, None]] * (len(gs_expression_norm) - 1))
    fig = make_subplots(rows=1 + len(gs_expression_norm), cols=2, specs=layout, shared_xaxes=False,
                        row_heights=[1 / (1 + len(gs_expression_norm))] * (1 + len(gs_expression_norm)), column_widths=[len(gs_expression_norm.columns) / (1 + len(gs_expression_norm.columns)), 1 / (1 + len(gs_expression_norm.columns))], horizontal_spacing=0.01)
    # Populate the first plot slot with the phenotype label information
    # NOTE: This will cause errors if using the phenotypes['Numeric'] Structure
    fig.append_trace(go.Heatmap(z=pandas.DataFrame(phenotypes['Phenotypes']).transpose(), colorscale='spectral', showscale=False, text=pandas.DataFrame(
        phenotypes['Labels']).transpose(), x=gs_expression_norm.columns.to_list(), y=["Phenotype"], name=''), row=1, col=1)
    # Add the plot containing the normalized expression heatmap annotated with the input expression data's values
    fig.append_trace(go.Heatmap(z=gs_expression_norm, colorscale='RdBu_r', colorbar={'title': {'text': 'Row Normalized Expression', 'side': 'top'}, 'x': 1.04, 'y': .9, 'len': 200, 'lenmode': 'pixels', 'thickness': 10,
                                                                                     'orientation': 'h', 'xanchor': 'left', 'yanchor': 'bottom'}, x=gs_expression_norm.columns.to_list(), y=gs_expression_norm.index.to_list(), name="", text=gs_expression, hovertemplate="%{text}"), row=2, col=1)
    # Add a plot containing the gene rankings in the gene list
    fig.append_trace(go.Heatmap(z=ranked_gs_genes, colorscale=ranking_colorscale, colorbar={'title': {'text': ranked_gs_genes.columns.to_list()[0], 'side': 'top'}, 'x': 1.04, 'y': .9, 'len': 200, 'lenmode': 'pixels', 'thickness': 10, 'orientation': 'h', 'xanchor': 'left', 'yanchor': 'top'}, zmax=float(
        ranked_genes.max()), zmin=float(ranked_genes.min()), x=ranked_gs_genes.columns.to_list(), y=ranked_gs_genes.index.to_list(), name=""), row=2, col=2)
    # Set the plot layout parameters to fit the data dimensions
    # [ (1,1) x,y   ]  [ (1,2) x2,y2 ]
    # ⎡ (2,1) x3,y3 ⎤  ⎡ (2,2) x4,y4 ⎤
    fig = fig.update_layout(
        xaxis=dict(dtick=1, side='top', tickangle=-90, type='category', showticklabels=True, scaleanchor='x3'), yaxis=dict(),
        xaxis2=dict(), yaxis2=dict(),
        xaxis3=dict(showticklabels=False), yaxis3=dict(dtick=1, showticklabels=True),
        xaxis4=dict(dtick=1, side='top', tickangle=-90, showticklabels=True), yaxis4=dict(showticklabels=False),
        margin=dict(autoexpand=True, b=0, r=0), height=20.01 + (20 * filtered_len), width=250.01 + (22 * len(phenotypes)))
    # save the <div> into python ## Reference for output options: https://plotly.com/python-api-reference/generated/plotly.io.to_html.html
    heatmap_fig = fig.to_html(
        full_html=False, include_plotlyjs='cdn')
    # , default_height="{:.0%}".format(filtered_len / 20 if filtered_len / 20 >= 1 else 1))
    return(heatmap_fig)


# Plot Barchart of top gene ranking
def plot_gene_rankings(ranked_genes, labels):
    corplot_data = ranked_genes.sort_values(
        ranked_genes.columns[0], ascending=False)
    corplot_bar = px.bar(corplot_data, y=corplot_data.columns[0], color=corplot_data.columns[0],
                         color_continuous_scale='RdBu_r', color_continuous_midpoint=0)
    corplot_bar = corplot_bar.update_traces(marker_line_width=0)
    corplot_bar = corplot_bar.add_vline(x=numpy.where(numpy.diff(numpy.sign(corplot_data.iloc[:, 0].values)))[
        0][0] + 0.5, line_dash="dashdot", annotation_text="zero-cross", line_color="grey")
    corplot_bar = corplot_bar.add_annotation(x=0, y=0, text=str(
        labels[0]) + " (positively correlated)", font={'color': "#EF553B"}, xanchor='left', yanchor='top', showarrow=False)
    corplot_bar = corplot_bar.add_annotation(x=len(corplot_data), y=0, text=str(
        labels[1]) + " (negatively correlated)", font={'color': "#636EFA"}, xanchor='right', yanchor='bottom', showarrow=False)
    corplot_bar = corplot_bar.update_layout(width=1280,
                                            margin=dict(autoexpand=True, t=24))
    corr_plot_fig = corplot_bar.to_html(
        full_html=False, include_plotlyjs='cdn')
    return(corr_plot_fig)


# Plot Preranked Heatmap for a given set of genes
def plot_set_prerank_heatmap(input_ds, phenotypes, ranked_genes, filtered_gs, ascending):
    ranking_colorscale = go.Heatmap(z=ranked_genes, colorscale='RdBu_r', zmid=0)[
        'colorscale']  # Create a colorscale for the ranking subset
    filtered_len = len(filtered_gs)
    ranked_gs_genes = ranked_genes.loc[filtered_gs].sort_values(
        ranked_genes.columns[0], ascending=ascending)
    gs_expression = input_ds.loc[ranked_gs_genes.index].copy()
    gs_expression_norm = gs_expression.subtract(gs_expression.min(axis=1), axis=0)\
        .divide(gs_expression.max(axis=1) - gs_expression.min(axis=1), axis=0)\
        .combine_first(gs_expression)  # Row Normalize Gene Set gene expression
    # Construct plotly heatmap
    # Instantiate a plot containing slots for the main heatmap and a slot for the phenotype label bar
    layout = [[{"rowspan": len(gs_expression_norm), "colspan": 1}]]
    layout.extend([[None]] * (len(gs_expression_norm) - 1))
    fig = make_subplots(rows=len(gs_expression_norm), cols=1, specs=layout, shared_xaxes=False,
                        row_heights=[1 / len(gs_expression_norm)] * len(gs_expression_norm), column_widths=[1], horizontal_spacing=0.01)
    # Construct plotly heatmap for the ranked list
    fig.append_trace(go.Heatmap(z=ranked_gs_genes, colorscale=ranking_colorscale, colorbar={'title': {'text': ranked_gs_genes.columns.to_list()[0], 'side': 'right'}, 'x': 1.04, 'y': .9, 'len': 200, 'lenmode': 'pixels', 'thickness': 10, 'orientation': 'v', 'xanchor': 'left', 'yanchor': 'top'}, zmax=float(
        ranked_genes.max()), zmin=float(ranked_genes.min()), x=ranked_gs_genes.columns.to_list(), y=ranked_gs_genes.index.to_list(), name=""), row=1, col=1)
    fig = fig.update_layout(
        xaxis=dict(dtick=1, side='top', tickangle=-90, showticklabels=True), yaxis=dict(dtick=1, showticklabels=True),
        margin=dict(autoexpand=True, b=0, r=0), height=20.01 + (20 * filtered_len), width=250.01 + (22 * len(phenotypes)))
    # save the <div> into python ## Reference for output options: https://plotly.com/python-api-reference/generated/plotly.io.to_html.html
    heatmap_fig = fig.to_html(
        full_html=False, include_plotlyjs='cdn')
    # , default_height="{:.0%}".format(filtered_len / 20 if filtered_len / 20 >= 1 else 1))
    return(heatmap_fig)


# Plot Permutation Distplot with Indepdenent KDE
def set_perm_indepkde_displot(random_score_matrix, true_es):
    if true_es > 0:
        pos_visible = True
        neg_visible = 'legendonly'
    elif true_es < 0:
        pos_visible = 'legendonly'
        neg_visible = True
    else:
        pos_visible = True
        neg_visible = True
    pos_perm = random_score_matrix[random_score_matrix >= 0]
    neg_perm = random_score_matrix[random_score_matrix <= 0]
    xrange = numpy.arange(numpy.min(neg_perm) - 0.1,
                          numpy.max(pos_perm) + 0.1, 0.0025)
    # Compute Positive Permutation Statistics
    pos_perm_kde = gaussian_kde(pos_perm)
    pos_kde_plot = go.Scatter(x=xrange[xrange >= 0], y=pos_perm_kde.pdf(xrange)[xrange >= 0], mode='lines', line=dict(
        width=1.5, color=px.colors.qualitative.Plotly[1]), name='Perm ES Gaussian KDE (Pos)', visible=pos_visible)
    set_pos_histogram = go.Histogram(x=pos_perm, marker=dict(
        color=px.colors.qualitative.Pastel1[0]), name='Perm ES Histogram (Pos)', visible=pos_visible)
    set_pos_rug = go.Box(x=pos_perm, marker_symbol='line-ns-open',
                         marker_color=px.colors.qualitative.Plotly[1], boxpoints='all', jitter=0, name='Perm ES Rugplot (Pos)', visible=pos_visible)
    # Compute Negative Permutation Statistics
    neg_perm_kde = gaussian_kde(neg_perm)
    neg_kde_plot = go.Scatter(x=xrange[xrange <= 0], y=neg_perm_kde.pdf(xrange)[xrange <= 0], mode='lines', line=dict(
        width=1.5, color=px.colors.qualitative.Plotly[0]), name='Perm ES Gaussian KDE (Neg)', visible=neg_visible)
    set_neg_histogram = go.Histogram(x=neg_perm, marker=dict(
        color=px.colors.qualitative.Pastel1[1]), name='Perm ES Histogram (Neg)', visible=neg_visible)
    set_neg_rug = go.Box(x=neg_perm, marker_symbol='line-ns-open',
                         marker_color=px.colors.qualitative.Plotly[0], boxpoints='all', jitter=0, name='Perm ES Rugplot (Neg)', visible=neg_visible)
    # Create Plot
    set_distplot = make_subplots(rows=2, cols=1, specs=[[{"secondary_y": True}], [
                                 {}]], row_heights=[0.8, 0.2], vertical_spacing=0.075, shared_xaxes=True)
    set_distplot = set_distplot.add_trace(
        set_pos_histogram, secondary_y=False, row=1, col=1)
    set_distplot = set_distplot.add_trace(
        pos_kde_plot, secondary_y=True, row=1, col=1)
    set_distplot = set_distplot.add_trace(
        set_pos_rug, row=2, col=1)
    set_distplot = set_distplot.add_trace(
        set_neg_histogram, secondary_y=False, row=1, col=1)
    set_distplot = set_distplot.add_trace(
        neg_kde_plot, secondary_y=True, row=1, col=1)
    set_distplot = set_distplot.add_trace(
        set_neg_rug, row=2, col=1)
    set_distplot = set_distplot.add_vline(
        x=true_es, line_dash="dashdot", annotation_text="Set True ES", annotation_position="top", line_color="grey")
    set_distplot = set_distplot.update_layout(yaxis=dict(title="Number of Permuted ES"),
                                              yaxis2=dict(rangemode='tozero', title=("Permutation KDE")), bargap=0.01, xaxis=dict(title="Permutation Enrichment Scores", autorange="reversed"), margin=dict(autoexpand=True, t=24, b=0, r=0), height=800, width=1280)
    set_distplot_fig = set_distplot.to_html(
        full_html=False, include_plotlyjs='cdn')
    return(set_distplot_fig)


# Plot Permutation Distplot with Joint KDE
def set_perm_jointkde_displot(random_score_matrix, true_es):
    if true_es > 0:
        pos_visible = True
        neg_visible = 'legendonly'
    elif true_es < 0:
        pos_visible = 'legendonly'
        neg_visible = True
    else:
        pos_visible = True
        neg_visible = True
    pos_perm = random_score_matrix[random_score_matrix >= 0]
    neg_perm = random_score_matrix[random_score_matrix <= 0]
    xrange = numpy.arange(numpy.nan_to_num(numpy.min(neg_es), nan=0) - 0.1,
                          numpy.nan_to_num(numpy.max(pos_es), nan=0) + 0.1, 0.0025)
    set_kde = gaussian_kde(random_score_matrix)
    # Compute Positive Permutation Statistics
    if len(pos_es) > 0:
        pos_kde_plot = go.Scatter(x=xrange[xrange >= 0], y=set_kde.pdf(xrange)[xrange >= 0], mode='lines', line=dict(
            width=1.5, color=px.colors.qualitative.Plotly[1]), name='Perm ES Gaussian KDE (Pos)', visible=pos_visible)
        set_pos_histogram = go.Histogram(x=pos_perm, marker=dict(
            color=px.colors.qualitative.Pastel1[0]), name='Perm ES Histogram (Pos)', visible=pos_visible)
        set_pos_rug = go.Box(x=pos_perm, marker_symbol='line-ns-open',
                             marker_color=px.colors.qualitative.Plotly[1], boxpoints='all', jitter=0, name='Perm ES Rugplot (Pos)', visible=pos_visible)
    # Compute Negative Permutation Statistics
    if len(neg_es) > 0:
        neg_kde_plot = go.Scatter(x=xrange[xrange <= 0], y=set_kde.pdf(xrange)[xrange <= 0], mode='lines', line=dict(
            width=1.5, color=px.colors.qualitative.Plotly[0]), name='Perm ES Gaussian KDE (Neg)', visible=neg_visible)
        set_neg_histogram = go.Histogram(x=neg_perm, marker=dict(
            color=px.colors.qualitative.Pastel1[1]), name='Perm ES Histogram (Neg)', visible=neg_visible)
        set_neg_rug = go.Box(x=neg_perm, marker_symbol='line-ns-open',
                             marker_color=px.colors.qualitative.Plotly[0], boxpoints='all', jitter=0, name='Perm ES Rugplot (Neg)', visible=neg_visible)
    # Create Plot
    set_distplot = make_subplots(rows=2, cols=1, specs=[[{"secondary_y": True}], [
                                 {}]], row_heights=[0.8, 0.2], vertical_spacing=0.075, shared_xaxes=True)
    # Add Positive Plot Elements
    if len(pos_es) > 0:
        set_distplot = set_distplot.add_trace(
            set_pos_histogram, secondary_y=False, row=1, col=1)
        set_distplot = set_distplot.add_trace(
            pos_kde_plot, secondary_y=True, row=1, col=1)
        set_distplot = set_distplot.add_trace(
            set_pos_rug, row=2, col=1)
    # Add Negative Plot Elements
    if len(neg_es) > 0:
        set_distplot = set_distplot.add_trace(
            set_neg_histogram, secondary_y=False, row=1, col=1)
        set_distplot = set_distplot.add_trace(
            neg_kde_plot, secondary_y=True, row=1, col=1)
        set_distplot = set_distplot.add_trace(
            set_neg_rug, row=2, col=1)
    set_distplot = set_distplot.add_vline(
        x=true_es, line_dash="dashdot", annotation_text="Set True ES", annotation_position="top", line_color="grey")
    set_distplot = set_distplot.update_layout(yaxis=dict(title="Number of Permuted ES"),
                                              yaxis2=dict(rangemode='tozero', title=("Permutation KDE")), bargap=0.01, xaxis=dict(title="Permutation Enrichment Scores", autorange="reversed"), margin=dict(autoexpand=True, t=24, b=0, r=0), height=800, width=1280)
    set_distplot_fig = set_distplot.to_html(
        full_html=False, include_plotlyjs='cdn')
    return(set_distplot_fig)


# Plot Global ES Distplot with Indepdenent KDE
def global_es_indepkde_distplot(score_matrix):
    pos_es = score_matrix[score_matrix >= 0]
    neg_es = score_matrix[score_matrix <= 0]
    xrange = numpy.arange(numpy.nan_to_num(numpy.min(neg_es), nan=0) - 0.1,
                          numpy.nan_to_num(numpy.max(pos_es), nan=0) + 0.1, 0.0025)
    # Compute Positive ES Statistics
    if len(pos_es) > 0:
        pos_es_kde = gaussian_kde(pos_es)
        pos_kde_plot = go.Scatter(x=xrange[xrange >= 0], y=pos_es_kde.pdf(xrange)[xrange >= 0], mode='lines', line=dict(
            width=1.5, color=px.colors.qualitative.Plotly[1]), name='ES Gaussian KDE (Pos)')
        set_pos_histogram = go.Histogram(x=pos_es, marker=dict(
            color=px.colors.qualitative.Pastel1[0]), name='ES Histogram (Pos)')
        set_pos_rug = go.Box(x=pos_es, marker_symbol='line-ns-open',
                             marker_color=px.colors.qualitative.Plotly[1], boxpoints='all', jitter=0, name='ES Rugplot (Pos)')
    # Compute Negative ES Statistics
    if len(neg_es) > 0:
        neg_es_kde = gaussian_kde(neg_es)
        neg_kde_plot = go.Scatter(x=xrange[xrange <= 0], y=neg_es_kde.pdf(xrange)[xrange <= 0], mode='lines', line=dict(
            width=1.5, color=px.colors.qualitative.Plotly[0]), name='ES Gaussian KDE (Neg)')
        set_neg_histogram = go.Histogram(x=neg_es, marker=dict(
            color=px.colors.qualitative.Pastel1[1]), name='ES Histogram (Neg)')
        set_neg_rug = go.Box(x=neg_es, marker_symbol='line-ns-open',
                             marker_color=px.colors.qualitative.Plotly[0], boxpoints='all', jitter=0, name='ES Rugplot (Neg)')
    # Create Plot
    set_distplot = make_subplots(rows=2, cols=1, specs=[[{"secondary_y": True}], [
                                 {}]], row_heights=[0.8, 0.2], vertical_spacing=0.075, shared_xaxes=True)
    # Add Positive Plot Elements
    if len(pos_es) > 0:
        set_distplot = set_distplot.add_trace(
            set_pos_histogram, secondary_y=False, row=1, col=1)
        set_distplot = set_distplot.add_trace(
            pos_kde_plot, secondary_y=True, row=1, col=1)
        set_distplot = set_distplot.add_trace(
            set_pos_rug, row=2, col=1)
    # Add Negative Plot Elements
    if len(neg_es) > 0:
        set_distplot = set_distplot.add_trace(
            set_neg_histogram, secondary_y=False, row=1, col=1)
        set_distplot = set_distplot.add_trace(
            neg_kde_plot, secondary_y=True, row=1, col=1)
        set_distplot = set_distplot.add_trace(
            set_neg_rug, row=2, col=1)
    set_distplot = set_distplot.update_layout(yaxis=dict(title="Number of Enrichment Scores"),
                                              yaxis2=dict(rangemode='tozero', title=("ES KDE")), bargap=0.01, xaxis=dict(title="Enrichment Scores", autorange="reversed"), margin=dict(autoexpand=True, t=24, b=0, r=0), height=800, width=1280)
    set_distplot_fig = set_distplot.to_html(
        full_html=False, include_plotlyjs='cdn')
    return(set_distplot_fig)


# Plot Global ES Distplot with Joint KDE
def global_es_jointkde_distplot(score_matrix):
    pos_es = score_matrix[score_matrix >= 0]
    neg_es = score_matrix[score_matrix <= 0]
    xrange = numpy.arange(numpy.min(neg_es) - 0.1,
                          numpy.max(pos_es) + 0.1, 0.0025)
    set_kde = gaussian_kde(score_matrix)
    # Compute Positive ES Statistics
    pos_kde_plot = go.Scatter(x=xrange[xrange >= 0], y=set_kde.pdf(xrange)[xrange >= 0], mode='lines', line=dict(
        width=1.5, color=px.colors.qualitative.Plotly[1]), name='ES Gaussian KDE (Pos)')
    set_pos_histogram = go.Histogram(x=pos_es, marker=dict(
        color=px.colors.qualitative.Pastel1[0]), name='ES Histogram (Pos)')
    set_pos_rug = go.Box(x=pos_es, marker_symbol='line-ns-open',
                         marker_color=px.colors.qualitative.Plotly[1], boxpoints='all', jitter=0, name='ES Rugplot (Pos)')
    # Compute Negative ES Statistics
    neg_kde_plot = go.Scatter(x=xrange[xrange <= 0], y=set_kde.pdf(xrange)[xrange <= 0], mode='lines', line=dict(
        width=1.5, color=px.colors.qualitative.Plotly[0]), name='ES Gaussian KDE (Neg)')
    set_neg_histogram = go.Histogram(x=neg_es, marker=dict(
        color=px.colors.qualitative.Pastel1[1]), name='ES Histogram (Neg)')
    set_neg_rug = go.Box(x=neg_es, marker_symbol='line-ns-open',
                         marker_color=px.colors.qualitative.Plotly[0], boxpoints='all', jitter=0, name='ES Rugplot (Neg)')
    # Create Plot
    set_distplot = make_subplots(rows=2, cols=1, specs=[[{"secondary_y": True}], [
                                 {}]], row_heights=[0.8, 0.2], vertical_spacing=0.075, shared_xaxes=True)
    set_distplot = set_distplot.add_trace(
        set_pos_histogram, secondary_y=False, row=1, col=1)
    set_distplot = set_distplot.add_trace(
        pos_kde_plot, secondary_y=True, row=1, col=1)
    set_distplot = set_distplot.add_trace(
        set_pos_rug, row=2, col=1)
    set_distplot = set_distplot.add_trace(
        set_neg_histogram, secondary_y=False, row=1, col=1)
    set_distplot = set_distplot.add_trace(
        neg_kde_plot, secondary_y=True, row=1, col=1)
    set_distplot = set_distplot.add_trace(
        set_neg_rug, row=2, col=1)
    set_distplot = set_distplot.update_layout(yaxis=dict(title="Number of Enrichment Scores"),
                                              yaxis2=dict(rangemode='tozero', title=("ES KDE")), bargap=0.01, xaxis=dict(title="Enrichment Scores", autorange="reversed"), margin=dict(autoexpand=True, t=24, b=0, r=0), height=800, width=1280)
    set_distplot_fig = set_distplot.to_html(
        full_html=False, include_plotlyjs='cdn')
    return(set_distplot)


def get_leading_edge(page_str):
    content = ast.literal_eval(page_str.split("\n")[6].strip())
    set_members = content[0][1]['text']
    gene_list = content[0][0]['text']
    gene_list_rank = content[0][0]['x']
    gene_list_metric = content[0][0]['y']
    es_index = content[0][4]['x']
    set_es = content[0][4]['y'][0]
    gene_ranks = dict(zip(gene_list, gene_list_rank))
    gene_metrics = dict(zip(gene_list, gene_list_metric))
    set_ranks = {k: int(v) for k, v in gene_ranks.items() if k in set_members}
    set_rank_metrics = {k: v for k,
                        v in gene_metrics.items() if k in set_members}
    # Extract negative distribution
    neg_symbols = content[0][2]['text']
    neg_running_es = content[0][2]['y']
    neg_es = min(neg_running_es)
    scored_neg_genes = dict(zip(neg_symbols, neg_running_es))
    negative_set = {k: v for k, v in scored_neg_genes.items()
                    if k in set_members}
    # Extract positive distribution
    pos_symbols = content[0][3]['text']
    pos_running_es = content[0][3]['y']
    scored_pos_genes = dict(zip(pos_symbols, pos_running_es))
    # Get complete Running ES Dict
    running_es = {gene: max(
        [scored_neg_genes[gene], scored_pos_genes[gene]], key=abs) for gene in gene_list}
    set_running_es = {k: v for k, v in running_es.items() if k in set_members}
    set_le_info = pandas.DataFrame([set_ranks, set_rank_metrics, set_running_es], index=[
                                   "Rank in Gene List", "Rank Metric Score", "Running ES"]).transpose()
    set_le_info["Rank in Gene List"] = set_le_info["Rank in Gene List"].astype(
        int)
    set_le_info["Core Enrichment"] = ""
    if set_es >= 0:
        set_le_info.loc[set_le_info['Rank in Gene List'].values >=
                        es_index, 'Core Enrichment'] = 'Yes'
        set_le_info.loc[set_le_info['Rank in Gene List'].values <
                        es_index, 'Core Enrichment'] = 'No'
    else:
        set_le_info.loc[set_le_info['Rank in Gene List'].values <=
                        es_index, 'Core Enrichment'] = 'Yes'
        set_le_info.loc[set_le_info['Rank in Gene List'].values >
                        es_index, 'Core Enrichment'] = 'No'
    set_le_info = set_le_info.sort_values(
        ["Rank Metric Score"], axis=0, ascending=(False)).reset_index()
    set_le_info.rename({'index': 'Gene Symbol'}, axis=1, inplace=True)
    set_le_info.index = set_le_info.index + 1
    leading_edge = set_le_info['Gene Symbol'][set_le_info['Core Enrichment'] == 'Yes'].values
    leading_edge = ",".join(leading_edge)
    set_le_info_styled = set_le_info.style.set_table_styles(
        [{"selector": "", "props": [("border", "1px solid grey")]},
         {"selector": "tbody td", "props": [("border", "1px solid grey")]},
            {"selector": "th", "props": [("border", "1px solid grey")]}
         ]
    )
    set_le_info_styled = set_le_info_styled.apply(highlight_leading_edge, subset=[
        'Core Enrichment'], axis=1)
    return(set_le_info_styled, leading_edge)


def highlight_leading_edge(column):
    highlight = 'background-color: #e6ffe6;'
    default = ''
    return [highlight if v == 'Yes' else default for v in column]
