from biocypher import BioCypher
from src.adapters import adapter_synthetic_pathways, adapter_synthetic_proteins

# Create an instance of BioCypher
bc = BioCypher()

# Configure the protein interaction adapter.
protein_adapter = adapter_synthetic_proteins.Adapter(
    tsv_path=adapter_synthetic_proteins.TSV_FILE_PATH_SYNTHETIC_PROTEINS,
    node_types=[adapter_synthetic_proteins.AdapterNodeType.PROTEIN],
    node_fields=[
        adapter_synthetic_proteins.AdapterProteinField.ID,
        adapter_synthetic_proteins.AdapterProteinField.PREFERRED_ID,
        adapter_synthetic_proteins.AdapterProteinField.GENE_SYMBOL,
        adapter_synthetic_proteins.AdapterProteinField.NCBI_TAX_ID,
    ],
    edge_types=[
        adapter_synthetic_proteins.AdapterEdgeType.PROTEIN_PROTEIN_INTERACTION,
        adapter_synthetic_proteins.AdapterEdgeType.BINDING,
        adapter_synthetic_proteins.AdapterEdgeType.ACTIVATION,
        adapter_synthetic_proteins.AdapterEdgeType.PHOSPHORYLATION,
        adapter_synthetic_proteins.AdapterEdgeType.UBIQUITINATION,
        adapter_synthetic_proteins.AdapterEdgeType.INHIBITION,
    ],
)

# Configure the pathway participation adapter.
pathway_adapter = adapter_synthetic_pathways.Adapter(
    tsv_path=adapter_synthetic_pathways.TSV_FILE_PATH_SYNTHETIC_PATHWAYS,
    node_types=[
        adapter_synthetic_pathways.AdapterNodeType.PROTEIN,
        adapter_synthetic_pathways.AdapterNodeType.PATHWAY,
    ],
    edge_types=[adapter_synthetic_pathways.AdapterEdgeType.PARTICIPATES_IN],
)

# Create a knowledge graph from both adapters.
for adapter in (protein_adapter, pathway_adapter):
    bc.write_nodes(adapter.get_nodes())
    bc.write_edges(adapter.get_edges())

# Generate assets for Neo4j exportation
bc.write_import_call()

# Print a summary when
bc.summary()