from enum import Enum, auto
from itertools import chain
from typing import Optional, Generator
from pathlib import Path

import pandas as pd
from biocypher._logger import logger

repo_root = Path(__file__).resolve().parents[2]
TSV_FILE_PATH_SYNTHETIC_PROTEINS = repo_root / "data" / "in" / "synthetic_protein_interactions.tsv"

class AdapterNodeType(Enum):
    """
    Define types of nodes the adapter can provide.
    """

    PROTEIN = auto()

class AdapterProteinField(Enum):
    """
    Define possible fields the adapter can provide for proteins.
    """

    ID = "id"
    PREFERRED_ID = "preferred_id"
    GENE_SYMBOL = "genesymbol"
    NCBI_TAX_ID = "ncbi_tax_id"

class AdapterEdgeType(Enum):
    """
    Enum for the types of the protein adapter.
    """

    PROTEIN_PROTEIN_INTERACTION = "protein_protein_interaction"
    BINDING = "binding"
    ACTIVATION = "activation"
    PHOSPHORYLATION = "phosphorylation"
    UBIQUITINATION = "ubiquitination"
    INHIBITION = "inhibition"

class AdapterProteinProteinEdgeField(Enum):
    """
    Define possible fields the adapter can provide for protein-protein edges.
    """

    INTERACTION_TYPE = "interaction_type"
    INTERACTION_SOURCE = "interaction_source"
    IS_DIRECTED = "is_directed"
    IS_STIMULATION = "is_stimulation"
    IS_INHIBITION = "is_inhibition"
    CONSENSUS_DIRECTION = "consensus_direction"
    CONSENSUS_STIMULATION = "consensus_stimulation"
    CONSENSUS_INHIBITION = "consensus_inhibition"

class Adapter:
    def __init__(
        self,
        tsv_path: str | Path = TSV_FILE_PATH_SYNTHETIC_PROTEINS,
        node_types: Optional[list] = None,
        node_fields: Optional[list] = None,
        edge_types: Optional[list] = None,
        edge_fields: Optional[list] = None,
    ):
        self.tsv_path = tsv_path
        self._set_types_and_fields(node_types, node_fields, edge_types, edge_fields)

    def _read_tsv(self) -> pd.DataFrame:
        """
        Reads and validates the TSV file.
        Returns:
            pd.DataFrame: DataFrame containing the TSV data.
        Raises:
            FileNotFoundError: If the file does not exist.
            ValueError: If required columns are missing.
        """
        if not Path(self.tsv_path).exists():
            logger.error(f"TSV file not found: {self.tsv_path}")
            raise FileNotFoundError(f"TSV file not found: {self.tsv_path}")
        df = pd.read_table(self.tsv_path, sep="\t", header=0)
        required_columns = [
            'source', 'target', 'source_genesymbol', 'target_genesymbol',
            'ncbi_tax_id_source', 'ncbi_tax_id_target', 'type',
            'entity_type_source', 'entity_type_target',
            'is_directed', 'is_stimulation', 'is_inhibition', 'consensus_direction',
            'consensus_stimulation', 'consensus_inhibition'
        ]
        missing = [col for col in required_columns if col not in df.columns]
        if missing:
            logger.error(f"Missing columns in TSV: {missing}")
            raise ValueError(f"TSV must contain columns: {missing}")
        return df

    def get_nodes(self) -> 'Generator[tuple[str, str, dict], None, None]':
        """
        Yields node tuples for node types specified in the adapter constructor.

        Returns:
            Generator[tuple[str, str, dict], None, None]:
                Each tuple is (id, label, properties).
        """
        logger.info("Reading nodes.")
        df = self._read_tsv()
        seen_nodes = set()

        # Generator for nodes in the `source` column
        for row in df.itertuples(index=False):
            id_ = str(row.source)
            input_label = "uniprot_protein"

            if id_ in seen_nodes:
                continue

            properties = {
                'genesymbol': row.source_genesymbol,
                'ncbi_tax_id': row.ncbi_tax_id_source,
                'entity_type': row.entity_type_source,
            }

            seen_nodes.add(id_)

            yield(
                id_,
                input_label,
                properties
            )

        # Generator for nodes in the `target` column
        for row in df.itertuples(index=False):
            id_ = str(row.target)
            input_label = "uniprot_protein"

            if id_ in seen_nodes:
                continue

            properties = {
                'genesymbol': row.target_genesymbol,
                'ncbi_tax_id': row.ncbi_tax_id_target,
                'entity_type': row.entity_type_target,
            }

            seen_nodes.add(id_)

            yield(
                id_,
                input_label,
                properties
            )

    def get_edges(self) -> 'Generator[tuple[str, str, str, str, dict], None, None]':
        """
        Yields edge tuples for edge types specified in the adapter constructor.

        Returns:
            Generator[tuple[str, str, str, str, dict], None, None]:
                Each tuple is (id, source, target, type, properties).
        """
        logger.info("Generating edges.")
        df = self._read_tsv()

        for row in df.itertuples(index=False):
            source = str(row.source)
            target = str(row.target)
            type = str(row.type)

            # Include edge type in ID to avoid collisions across interaction types.
            id = f"{source}_{target}_{type}"

            properties = {
                'is_directed': row.is_directed,
                'is_stimulation': row.is_stimulation,
                'is_inhibition': row.is_inhibition,
                'consensus_direction': row.consensus_direction,
                'consensus_stimulation': row.consensus_stimulation,
                'consensus_inhibition': row.consensus_inhibition
            }

            yield (
                id,
                source,
                target,
                type,
                properties
            )

    def get_node_count(self) -> int:
        """
        Returns the number of nodes generated by the adapter.

        Returns:
            int: Number of nodes generated.
        """
        return sum(1 for _ in self.get_nodes())

    def _set_types_and_fields(self, node_types, node_fields, edge_types, edge_fields) -> None:
        """
        Sets the node and edge types and fields for the adapter.

        Args:
            node_types (Optional[list]): List of node types.
            node_fields (Optional[list]): List of node fields.
            edge_types (Optional[list]): List of edge types.
            edge_fields (Optional[list]): List of edge fields.
        """
        if node_types:
            self.node_types = node_types
        else:
            self.node_types = [type for type in AdapterNodeType]

        if node_fields:
            self.node_fields = node_fields
        else:
            self.node_fields = [
                field
                for field in chain(
                    AdapterProteinField,
                )
            ]

        if edge_types:
            self.edge_types = edge_types
        else:
            self.edge_types = [type for type in AdapterEdgeType]

        if edge_fields:
            self.edge_fields = edge_fields
        else:
            self.edge_fields = [field for field in chain()]