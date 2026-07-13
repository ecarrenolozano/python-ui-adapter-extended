from pathlib import Path

import pytest

from src.adapters.adapter_synthetic_pathways import Adapter


def _write_tsv(path: Path, content: str) -> None:
    path.write_text(content, encoding="utf-8")


@pytest.fixture
def valid_pathway_tsv(tmp_path: Path) -> Path:
    tsv = tmp_path / "synthetic_pathways.tsv"
    _write_tsv(
        tsv,
        "\t".join(
            [
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
        )
        + "\n"
        + "\n".join(
            [
                "P1\tG1\t9606\tprotein\tR-HSA-1\tPathway A\tDesc A\tparticipates_in\tSYNTHETIC\tReactome-like dummy\tTrue",
                "P1\tG1\t9606\tprotein\tR-HSA-2\tPathway B\tDesc B\tparticipates_in\tSYNTHETIC\tReactome-like dummy\tTrue",
                "P2\tG2\t9606\tprotein\tR-HSA-1\tPathway A\tDesc A\tparticipates_in\tSYNTHETIC\tReactome-like dummy\tTrue",
            ]
        )
        + "\n",
    )
    return tsv


def test_read_tsv_raises_on_missing_required_columns(tmp_path: Path) -> None:
    tsv = tmp_path / "invalid.tsv"
    _write_tsv(tsv, "uniprot_id\tpathway_id\nP1\tR-HSA-1\n")

    adapter = Adapter(tsv_path=tsv)

    with pytest.raises(ValueError):
        next(adapter.get_nodes())


def test_get_nodes_are_deduplicated_and_labeled(valid_pathway_tsv: Path) -> None:
    adapter = Adapter(tsv_path=valid_pathway_tsv)
    nodes = list(adapter.get_nodes())

    assert len(nodes) == 4

    protein_nodes = [n for n in nodes if n[1] == "uniprot_protein"]
    pathway_nodes = [n for n in nodes if n[1] == "reactome_pathway"]

    assert len(protein_nodes) == 2
    assert len(pathway_nodes) == 2

    for node_id, label, props in protein_nodes:
        assert node_id.startswith("P")
        assert label == "uniprot_protein"
        assert {"genesymbol", "ncbi_tax_id", "entity_type"}.issubset(props.keys())

    for node_id, label, props in pathway_nodes:
        assert node_id.startswith("R-HSA-")
        assert label == "reactome_pathway"
        assert {"name", "description", "source_database", "is_synthetic"}.issubset(props.keys())


def test_get_edges_have_valid_structure(valid_pathway_tsv: Path) -> None:
    adapter = Adapter(tsv_path=valid_pathway_tsv)
    edges = list(adapter.get_edges())

    assert len(edges) == 3

    for edge_id, source, target, edge_type, props in edges:
        assert isinstance(edge_id, str)
        assert source.startswith("P")
        assert target.startswith("R-HSA-")
        assert edge_type == "participates_in"
        assert {"evidence_code", "source_database", "is_synthetic"}.issubset(props.keys())
        assert props["is_synthetic"] is True


def test_relationship_integrity_edge_endpoints_exist_in_nodes(valid_pathway_tsv: Path) -> None:
    adapter = Adapter(tsv_path=valid_pathway_tsv)
    nodes = list(adapter.get_nodes())
    edges = list(adapter.get_edges())

    node_ids = {n[0] for n in nodes}
    for _, source, target, _, _ in edges:
        assert source in node_ids
        assert target in node_ids


def test_get_node_count_matches_materialized_nodes(valid_pathway_tsv: Path) -> None:
    adapter = Adapter(tsv_path=valid_pathway_tsv)
    assert adapter.get_node_count() == len(list(adapter.get_nodes()))