from biocypher import BioCypher
from src.adapters.adapter_synthetic_proteins import (
    AdapterNodeType,
    AdapterProteinField,
    AdapterEdgeType,
    Adapter,
    TSV_FILE_PATH_SYNTHETIC_PROTEINS
)

# Create an instance of BioCypher
bc = BioCypher()

# Choose the node type you want appear in the Knowledge Graph
node_types = [AdapterNodeType.PROTEIN]

# Choose protein adapter fields to include in the knowledge graph.
node_fields = [
    AdapterProteinField.ID,
    AdapterProteinField.PREFERRED_ID,
    AdapterProteinField.GENE_SYMBOL,
    AdapterProteinField.NCBI_TAX_ID,
]

# Choose the node type you want appear in the Knowledge Graph
edge_types = [
    AdapterEdgeType.PROTEIN_PROTEIN_INTERACTION,
    AdapterEdgeType.BINDING,
    AdapterEdgeType.ACTIVATION,
    AdapterEdgeType.PHOSPHORYLATION,
    AdapterEdgeType.UBIQUITINATION,
    AdapterEdgeType.INHIBITION,
]

# (there is not code here!) Choose interaction adapter fields to include in the knowledge graph.
# By default, in case of not specifying this, BioCypher will bring all the fields defined in the adapter

# Create an adapter instance
adapter = Adapter(
    tsv_path=TSV_FILE_PATH_SYNTHETIC_PROTEINS,
    node_types=node_types,
    node_fields=node_fields,
    edge_types=edge_types,
)

# Create a knowledge graph from the adapter
bc.write_nodes(adapter.get_nodes())
bc.write_edges(adapter.get_edges())

# Generate assets for Neo4j exportation
bc.write_import_call()

# Print a summary when
bc.summary()