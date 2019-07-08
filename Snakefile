configfile: 'config.yaml'
workdir: config['working_directory']


rule all:
    input:
        expand('images/pathway_comparison_{project}_{pathway}.pdf', project=config['projects'], pathway=config['pathways']),
        expand('images/group_relations_{project}_{pathway}.pdf', project=config['projects'], pathway=config['pathways'])

rule download_TCGA:
    output:
        data_dir = directory('tcga_data/{project}/GDCdata'),
        expression_file = 'tcga_data/{project}/expression_matrix.csv',
        classification_file = 'tcga_data/{project}/case_classifications.csv'
    script:
        'scripts/tcga_download.R'

rule download_KEGG:
    output:
        xml_file = 'kegg_data/{pathway}.xml',
        network_file = 'kegg_data/{pathway}.edgelist.csv',
        plot_file = 'kegg_data/{pathway}.pdf'
    script:
        'scripts/kegg_download.R'

rule compute_causal_effects:
    input:
        expression_file = 'tcga_data/{project}/expression_matrix.csv',
        classification_file = 'tcga_data/{project}/case_classifications.csv',
        network_file = 'kegg_data/{pathway}.edgelist.csv'
    output:
        normal_file = 'results/{project}/{pathway}/graph.normal.edgelist.csv',
        tumor_file = 'results/{project}/{pathway}/graph.tumor.edgelist.csv'
    script:
        'scripts/compute_causal_effects.R'

rule visualize:
    input:
        normal_file = 'results/{project}/{pathway}/graph.normal.edgelist.csv',
        tumor_file = 'results/{project}/{pathway}/graph.tumor.edgelist.csv'
    output:
        img_file = 'images/pathway_comparison_{project}_{pathway}.pdf'
    script:
        'scripts/visualization.py'

rule investigate_group_relations:
    input:
        expression_file = 'tcga_data/{project}/expression_matrix.csv',
        classification_file = 'tcga_data/{project}/case_classifications.csv',
        network_file = 'kegg_data/{pathway}.edgelist.csv'
    output:
        data_file = 'results/group_relations/{project}/{pathway}/dce_data.csv',
        plot_file = 'images/group_relations_{project}_{pathway}.pdf'
    script:
        'scripts/group_relations.R'
