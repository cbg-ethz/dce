import pandas as pd

import KEGGutils as kg


def main(pathway_id, csv_fname, geneid_fname, nodes_fname, edges_fname, result_fname):
    # read DCE data and more
    df_dce = pd.read_csv(csv_fname)
    df_geneids = pd.read_csv(geneid_fname)

    # read hetnet
    df_nodes = pd.read_csv(nodes_fname, sep='\t')
    df_edges = pd.read_csv(edges_fname, sep='\t')

    # retrieve KEGG pathway genes
    pathway = kg.KEGGpathway(pathway_id=pathway_id)
    genes = {g.strip(' .')
             for e in pathway.genes.values()
             for g in e['description'].split(',')}

    # map some IDs
    ensembl2symbol = df_geneids.set_index('ENSEMBL').to_dict()['SYMBOL']
    df_dce['source'] = df_dce['source'].map(ensembl2symbol)
    df_dce['target'] = df_dce['target'].map(ensembl2symbol)

    # combine
    kegg_gene_idx = df_nodes.loc[df_nodes['name'].isin(genes), 'id'].to_list()
    df_pw = df_edges[
        (df_edges['source'].isin(kegg_gene_idx)) |
        (df_edges['target'].isin(kegg_gene_idx))
    ]

    compound_gene_interactions = ['CuG', 'CdG', 'CbG']
    df_drug = df_pw[df_pw['metaedge'].isin(compound_gene_interactions)]

    # annotate result
    df_res = df_drug.merge(
        df_nodes,
        how='left', left_on='source', right_on='id'
    ).drop(columns=['kind']).rename(columns={'name': 'compound_name'})

    df_res = df_res.merge(
        df_nodes,
        how='left', left_on='target', right_on='id'
    ).drop(columns=['kind']).rename(columns={'name': 'gene_name'})

    df_res = df_res.drop(columns=['id_x', 'id_y'])

    df_res.to_csv(result_fname, index=False)

    # further investigation
    breast_cancer_id = df_nodes[(df_nodes['kind'] == 'Disease') & (df_nodes['name'].str.contains('breast cancer'))]['id'].iloc[0]

    breast_cancer_genes = df_edges.loc[(df_edges['source'].isin([breast_cancer_id])) | (df_edges['target'].isin([breast_cancer_id])), 'target'].tolist()

    foo = df_res[~df_res['target'].isin(breast_cancer_genes)]


if __name__ == '__main__':
    main(
        snakemake.wildcards.pathway,
        snakemake.input.csv_fname,
        snakemake.input.geneid_fname,
        snakemake.input.hetnet_nodes_fname,
        snakemake.input.hetnet_edges_fname,
        snakemake.output.result_fname)
