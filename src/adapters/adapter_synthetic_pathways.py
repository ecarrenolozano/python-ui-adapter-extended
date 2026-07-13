from enum import Enum, auto
from itertools import chain
from pathlib import Path
from typing import Generator, Optional

import pandas as pd
from biocypher._logger import logger

repo_root = Path(__file__).resolve().parents[2]
TSV_FILE_PATH_SYNTHETIC_PATHWAYS = repo_root / "data" / "in" / "synthetic_reactome_protein_pathway.tsv"


class AdapterNodeType(Enum):
    """Define types of nodes the adapter can provide."""

    PROTEIN = auto()
    PATHWAY = auto()


class AdapterProteinField(Enum):
    """Define possible fields the adapter can provide for proteins."""

    ID = "id"
    PREFERRED_ID = "preferred_id"
    GENE_SYMBOL = "genesymbol"
    NCBI_TAX_ID = "ncbi_tax_id"
    ENTITY_TYPE = "entity_type"


class AdapterPathwayField(Enum):
    """Define possible fields the adapter can provide for pathways."""

    ID = "id"
    PREFERRED_ID = "preferred_id"
    NAME = "name"
    DESCRIPTION = "description"
    SOURCE_DATABASE = "source_database"
    IS_SYNTHETIC = "is_synthetic"


class AdapterEdgeType(Enum):
    """Enum for the edge types of the pathway adapter."""

    PARTICIPATES_IN = "participates_in"


class AdapterProteinPathwayEdgeField(Enum):
    """Define possible fields the adapter can provide for protein-pathway edges."""

    EVIDENCE_CODE = "evidence_code"
    SOURCE_DATABASE = "source_database"
    IS_SYNTHETIC = "is_synthetic"


def _as_bool(value) -> bool:
    """Normalize TSV bool-like values into Python bools."""
    if isinstance(value, bool):
        return value
    if pd.isna(value):
        return False
    return str(value).strip().lower() in {"true", "1", "yes", "y", "t"}


class Adapter:
    def __init__(
        self,
        tsv_path: str | Path = TSV_FILE_PATH_SYNTHETIC_PATHWAYS,
        node_types: Optional[list] = None,
        node_fields: Optional[list] = None,
        edge_types: Optional[list] = None,
        edge_fields: Optional[list] = None,
    ):
        self.tsv_path = tsv_path
        self._set_types_and_fields(node_types, node_fields, edge_types, edge_fields)

    def _read_tsv(self) -> pd.DataFrame:
        """
        Read and validate the pathway TSV file.

        Returns:
            pd.DataFrame: DataFrame containing TSV data.

        Raises:
            FileNotFoundError: If the file does not exist.
            ValueError: If required columns are missing.
        """
        if not Path(self.tsv_path).exists():
            logger.error(f"TSV file not found: {self.tsv_path}")
            raise FileNotFoundError(f"TSV file not found: {self.tsv_path}")

        df = pd.read_table(self.tsv_path, sep="\t", header=0)
        required_columns = [
            "uniprot_id",
            "gene_symbol",
            "ncbi_tax_id",
            "entity_type_source",
            "pathway_id",
            "pathway_name",
            "pathway_description",
            "relationship_type",
            "evidence_code",
            "source_database",
            "is_synthetic",
        ]

        missing = [col for col in required_columns if col not in df.columns]
        if missing:
            logger.error(f"Missing columns in TSV: {missing}")
            raise ValueError(f"TSV must contain columns: {missing}")

        return df

    def get_nodes(self) -> Generator[tuple[str, str, dict], None, None]:
        """
        Yield node tuples for proteins and pathways.

        Returns:
            Generator[tuple[str, str, dict], None, None]:
                Each tuple is (id, label, properties).
        """
        logger.info("Reading protein and pathway nodes.")
        df = self._read_tsv()

        seen_proteins = set()
        for row in df.itertuples(index=False):
            protein_id = str(row.uniprot_id)
            if protein_id in seen_proteins:
                continue

            seen_proteins.add(protein_id)
            yield (
                protein_id,
                "uniprot_protein",
                {
                    "genesymbol": row.gene_symbol,
                    "ncbi_tax_id": row.ncbi_tax_id,
                    "entity_type": row.entity_type_source,
                },
            )

        seen_pathways = set()
        for row in df.itertuples(index=False):
            pathway_id = str(row.pathway_id)
            if pathway_id in seen_pathways:
                continue

            seen_pathways.add(pathway_id)
            yield (
                pathway_id,
                "reactome_pathway",
                {
                    "name": row.pathway_name,
                    "description": row.pathway_description,
                    "source_database": row.source_database,
                    "is_synthetic": _as_bool(row.is_synthetic),
                },
            )

    def get_edges(self) -> Generator[tuple[str, str, str, str, dict], None, None]:
        """
        Yield edge tuples for protein-pathway relationships.

        Returns:
            Generator[tuple[str, str, str, str, dict], None, None]:
                Each tuple is (id, source, target, type, properties).
        """
        logger.info("Generating protein-pathway edges.")
        df = self._read_tsv()

        for row in df.itertuples(index=False):
            source = str(row.uniprot_id)
            target = str(row.pathway_id)
            edge_type = str(row.relationship_type)
            edge_id = f"{source}_{target}_{edge_type}"

            yield (
                edge_id,
                source,
                target,
                edge_type,
                {
                    "evidence_code": row.evidence_code,
                    "source_database": row.source_database,
                    "is_synthetic": _as_bool(row.is_synthetic),
                },
            )

    def get_node_count(self) -> int:
        """Return number of nodes generated by the adapter."""
        return sum(1 for _ in self.get_nodes())

    def _set_types_and_fields(self, node_types, node_fields, edge_types, edge_fields) -> None:
        """Set node and edge types and fields for adapter initialisation."""
        if node_types:
            self.node_types = node_types
        else:
            self.node_types = [node_type for node_type in AdapterNodeType]

        if node_fields:
            self.node_fields = node_fields
        else:
            self.node_fields = [
                field
                for field in chain(
                    AdapterProteinField,
                    AdapterPathwayField,
                )
            ]

        if edge_types:
            self.edge_types = edge_types
        else:
            self.edge_types = [edge_type for edge_type in AdapterEdgeType]

        if edge_fields:
            self.edge_fields = edge_fields
        else:
            self.edge_fields = [field for field in chain(AdapterProteinPathwayEdgeField)]