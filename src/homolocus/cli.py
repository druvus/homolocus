"""Command-line interface: ``homolocus find`` (spec 13). A thin wrapper over the Python API.

Settings flags are generated from the configuration dataclasses, so their
types, choices and defaults live in one place. Unset flags are not passed on
and the dataclass defaults apply. Coordinates on the command line are
1-based inclusive.
"""

import argparse
import dataclasses
from dataclasses import dataclass
from typing import Any, Literal, get_args, get_origin, get_type_hints

from homolocus.config import DiamondSearch, FamilyGrouping, RegionDetection, ScoringWeights

OPTION_HELP = {
    "evalue": "pair-level E-value threshold",
    "min_identity": "pair-level minimum identity, percent",
    "sensitivity": "DIAMOND sensitivity mode: %(choices)s",
    "max_target_seqs": "DIAMOND targets per query; 0 is unlimited",
    "max_hsps": "DIAMOND HSPs per pair; 0 reports all",
    "hsp_evalue": "per-HSP E-value cutoff passed to DIAMOND",
    "threads": "DIAMOND threads; unset leaves DIAMOND's default",
    "family_min_identity": "identity linking paralogs into a family, percent",
    "family_min_coverage": "coverage both paralogs must reach, fraction",
    "min_matches": "minimum represented families per candidate (at least 2)",
    "min_reference_coverage": "minimum fraction of reference families represented",
    "min_reference_gene_coverage": "coverage at which a reference gene counts as represented",
    "min_event_coverage": "query or subject coverage of a matching event",
    "max_intervening_genes": "most genes between consecutive matching events",
    "max_inter_match_bp_gap": "largest gap between consecutive matching events, bp",
    "max_region_bp": "largest unpadded candidate span, bp",
    "max_candidate_overlap": "overlap above which a lower-ranked candidate is dropped",
    "padding_genes": "genes added on each side of a candidate",
    "padding_bp": "base pairs added on each side of a candidate",
    "edge_margin_bp": "distance from a contig end that flags a candidate as edge",
    "min_split_contribution": "share of a gene a split fragment must add",
    "max_fragment_overlap": "largest overlap between split fragments",
    "min_fusion_part_coverage": "coverage each fused part must reach",
    "content": "weight of the content component",
    "similarity": "weight of the similarity component",
    "order": "weight of the order component",
    "orientation": "weight of the orientation component",
    "spacing": "weight of the spacing component",
}
"""Help text for every configuration field, keyed by field name."""

_CONFIG_GROUPS: tuple[tuple[str, type, str], ...] = (
    ("DIAMOND search", DiamondSearch, ""),
    ("family grouping", FamilyGrouping, ""),
    ("region detection", RegionDetection, ""),
    ("scoring weights", ScoringWeights, "weight_"),
)
_DISABLE_FLAGS = {"family_min_identity": ("--no-family-grouping", "do not group paralogs")}


@dataclass(frozen=True)
class FinderSettings:
    """The four RegionFinder configurations chosen on the command line."""

    search: DiamondSearch
    family_grouping: FamilyGrouping
    detection: RegionDetection
    scoring: ScoringWeights


def build_parser() -> argparse.ArgumentParser:
    """The ``homolocus`` argument parser with its ``find`` command."""
    parser = argparse.ArgumentParser(prog="homolocus", description=__doc__.splitlines()[0])
    commands = parser.add_subparsers(dest="command", required=True)
    find = commands.add_parser("find", help="find regions homologous to a reference region")
    _add_inputs(find)
    for title, config, prefix in _CONFIG_GROUPS:
        _add_config_group(find.add_argument_group(title), config, prefix)
    return parser


def finder_settings(args: argparse.Namespace) -> FinderSettings:
    """Build the configurations from parsed flags; unset flags keep the defaults.

    Raises ValueError when a value breaks a configuration rule.
    """
    return FinderSettings(
        search=DiamondSearch(**_given(args, DiamondSearch, "")),
        family_grouping=FamilyGrouping(**_given(args, FamilyGrouping, "")),
        detection=RegionDetection(**_given(args, RegionDetection, "")),
        scoring=ScoringWeights(**_given(args, ScoringWeights, "weight_")),
    )


def _add_inputs(parser: argparse.ArgumentParser) -> None:
    inputs = parser.add_argument_group("inputs and output")
    inputs.add_argument("--reference", required=True, help="reference GenBank file")
    inputs.add_argument("--contig", help="reference record; optional for a single record")
    inputs.add_argument("--start", type=int, help="region start, 1-based inclusive")
    inputs.add_argument("--end", type=int, help="region end, 1-based inclusive")
    inputs.add_argument("--genome", required=True, help="target genome, GenBank or FASTA")
    inputs.add_argument(
        "--annotation",
        choices=("auto", "pyrodigal"),
        default="auto",
        help="target genes: supplied CDS when present, or always Pyrodigal (default: auto)",
    )
    inputs.add_argument("--output", required=True, help="output directory")


def _add_config_group(group: argparse._ArgumentGroup, config: type, prefix: str) -> None:
    hints = get_type_hints(config)
    for field in dataclasses.fields(config):
        dest = prefix + field.name
        target: argparse._ActionsContainer = group
        if dest in _DISABLE_FLAGS:
            target = group.add_mutually_exclusive_group()
            flag, text = _DISABLE_FLAGS[dest]
            target.add_argument(
                flag,
                dest=dest,
                action="store_const",
                const=None,
                default=argparse.SUPPRESS,
                help=text,
            )
        target.add_argument(
            "--" + dest.replace("_", "-"),
            dest=dest,
            default=argparse.SUPPRESS,
            help=_help(field),
            **_value_spec(hints[field.name]),
        )


def _help(field: "dataclasses.Field[Any]") -> str:
    text = OPTION_HELP[field.name]
    return text if field.default is None else f"{text} (default: {field.default})"


def _value_spec(hint: object) -> dict[str, Any]:
    """argparse ``type`` or ``choices`` for a field type; None is never typed in."""
    if get_origin(hint) is Literal:
        return {"choices": get_args(hint), "metavar": "MODE"}
    values = [arg for arg in get_args(hint) if arg is not type(None)]
    value_type = values[0] if values else hint
    return {"type": value_type, "metavar": getattr(value_type, "__name__", "VALUE").upper()}


def _given(args: argparse.Namespace, config: type, prefix: str) -> dict[str, Any]:
    names = (field.name for field in dataclasses.fields(config))
    return {name: getattr(args, prefix + name) for name in names if hasattr(args, prefix + name)}
