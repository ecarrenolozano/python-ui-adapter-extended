from pathlib import Path

import pytest
import yaml
from biocypher import BioCypher

from src.adapters.adapter_synthetic_proteins import Adapter


def _write_tsv(path: Path, content: str) -> None:
    path.write_text(content, encoding="utf-8")


@pytest.fixture
def valid_tsv(tmp_path: Path) -> Path:
    tsv = tmp_path / "synthetic.tsv"
    _write_tsv(
        tsv,
        "\t".join(
            [
                "source",
                "target",
                "source_genesymbol",
                "target_genesymbol",
                "ncbi_tax_id_source",
                "ncbi_tax_id_target",
                "type",
                "entity_type_source",
                "entity_type_target",
                "is_directed",
                "is_stimulation",
                "is_inhibition",
                "consensus_direction",
                "consensus_stimulation",
                "consensus_inhibition",
            ]
        )
        + "\n"
        + "\n".join(
            [
                "P1\tP2\tG1\tG2\t9606\t9606\tbinding\tprotein\tprotein\tTrue\tFalse\tFalse\tTrue\tFalse\tFalse",
                "P1\tP2\tG1\tG2\t9606\t9606\tactivation\tprotein\tprotein\tTrue\tTrue\tFalse\tTrue\tTrue\tFalse",
                "P2\tP3\tG2\tG3\t9606\t9606\tphosphorylation\tprotein\tprotein\tTrue\tTrue\tFalse\tTrue\tTrue\tFalse",
            ]
        )
        + "\n",
    )
    return tsv


def test_read_tsv_raises_on_missing_required_columns(tmp_path: Path) -> None:
    tsv = tmp_path / "invalid.tsv"
    _write_tsv(tsv, "source\ttarget\nP1\tP2\n")

    adapter = Adapter(tsv_path=tsv)

    with pytest.raises(ValueError):
        next(adapter.get_nodes())


def test_get_nodes_are_deduplicated_and_schema_shaped(valid_tsv: Path) -> None:
    adapter = Adapter(tsv_path=valid_tsv)
    nodes = list(adapter.get_nodes())

    assert len(nodes) == 3

    node_ids = [n[0] for n in nodes]
    assert len(set(node_ids)) == 3

    for node_id, label, props in nodes:
        assert isinstance(node_id, str)
        assert label == "uniprot_protein"
        assert {"genesymbol", "ncbi_tax_id", "entity_type"}.issubset(props.keys())


def test_get_edges_have_unique_ids_and_valid_structure(valid_tsv: Path) -> None:
    adapter = Adapter(tsv_path=valid_tsv)
    edges = list(adapter.get_edges())

    assert len(edges) == 3

    edge_ids = [e[0] for e in edges]
    assert len(set(edge_ids)) == len(edge_ids)

    allowed_types = {"binding", "activation", "phosphorylation"}
    for edge_id, source, target, edge_type, props in edges:
        assert isinstance(edge_id, str)
        assert source.startswith("P")
        assert target.startswith("P")
        assert edge_type in allowed_types
        assert {
            "is_directed",
            "is_stimulation",
            "is_inhibition",
            "consensus_direction",
            "consensus_stimulation",
            "consensus_inhibition",
        }.issubset(props.keys())


def test_relationship_integrity_edge_endpoints_exist_in_nodes(valid_tsv: Path) -> None:
    adapter = Adapter(tsv_path=valid_tsv)
    nodes = list(adapter.get_nodes())
    edges = list(adapter.get_edges())

    node_ids = {n[0] for n in nodes}
    for _, source, target, _, _ in edges:
        assert source in node_ids
        assert target in node_ids


def test_get_node_count_matches_materialized_nodes(valid_tsv: Path) -> None:
    adapter = Adapter(tsv_path=valid_tsv)
    assert adapter.get_node_count() == len(list(adapter.get_nodes()))


def test_pipeline_output_headers_align_with_schema_config(
    valid_tsv: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    repo_root = Path(__file__).resolve().parents[1]
    schema_src = repo_root / "config" / "schema_config.yaml"
    project_dir = tmp_path / "mini_project"
    config_dir = project_dir / "config"
    config_dir.mkdir(parents=True, exist_ok=True)

    (config_dir / "schema_config.yaml").write_text(
        schema_src.read_text(encoding="utf-8"),
        encoding="utf-8",
    )

    cache_dir = (repo_root / ".cache").as_posix()
    (config_dir / "biocypher_config.yaml").write_text(
        (
            "biocypher:\n"
            "  offline: true\n"
            "  debug: false\n"
            "  schema_config_path: config/schema_config.yaml\n"
            f"  cache_directory: {cache_dir}\n\n"
            "neo4j:\n"
            "  database_name: neo4j\n"
            "  delimiter: '\\t'\n"
            "  array_delimiter: '|'\n"
            "  skip_duplicate_nodes: true\n"
            "  skip_bad_relationships: true\n"
            "  edge_labels_order: Leaves\n"
        ),
        encoding="utf-8",
    )

    monkeypatch.chdir(project_dir)
    bc = BioCypher()
    adapter = Adapter(tsv_path=valid_tsv)
    bc.write_nodes(adapter.get_nodes())
    bc.write_edges(adapter.get_edges())
    bc.write_import_call()

    out_root = project_dir / "biocypher-out"
    assert out_root.exists()
    run_dir = max(out_root.iterdir(), key=lambda p: p.stat().st_mtime)

    schema = yaml.safe_load((config_dir / "schema_config.yaml").read_text(encoding="utf-8"))
    schema_edge_labels = {
        value.get("input_label")
        for value in schema.values()
        if isinstance(value, dict)
        and value.get("represented_as") == "edge"
        and value.get("input_label")
    }

    protein_header = (run_dir / "Protein-header.csv").read_text(encoding="utf-8").splitlines()[0]
    protein_columns = set(protein_header.split("\t"))
    assert {":ID", "genesymbol", "ncbi_tax_id:long", "entity_type", "id", "preferred_id", ":LABEL"}.issubset(
        protein_columns
    )

    edge_header_files = [
        p
        for p in run_dir.glob("*-header.csv")
        if p.name != "Protein-header.csv"
    ]
    assert edge_header_files

    for edge_header_file in edge_header_files:
        edge_label = edge_header_file.stem.replace("-header", "").lower()
        assert edge_label in schema_edge_labels

        edge_columns = set(edge_header_file.read_text(encoding="utf-8").splitlines()[0].split("\t"))
        assert {
            ":START_ID",
            "id",
            "is_directed:boolean",
            "is_stimulation:boolean",
            "is_inhibition:boolean",
            "consensus_direction:boolean",
            "consensus_stimulation:boolean",
            "consensus_inhibition:boolean",
            ":END_ID",
            ":TYPE",
        }.issubset(edge_columns)
