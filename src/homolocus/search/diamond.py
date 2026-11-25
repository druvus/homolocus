"""DIAMOND blastp: commands, tabular parsing and versions (spec 6).

This is the only module that knows about DIAMOND. It turns DIAMOND output
into Alignment models with canonical gene IDs and internal coordinates.
"""

import logging
import re
import shlex
import shutil
import subprocess
import tempfile
from collections.abc import Iterable
from pathlib import Path

from homolocus.config import DiamondSearch
from homolocus.coordinates import to_zero_based
from homolocus.models.alignment import Alignment
from homolocus.search.protein_index import ProteinIndex

OUTPUT_FIELDS = (
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "qlen",
    "slen",
)
"""Tabular columns requested from DIAMOND (spec 6.4)."""

logger = logging.getLogger(__name__)

_COLUMN = {name: position for position, name in enumerate(OUTPUT_FIELDS)}

MINIMUM_VERSION = (2, 0)

_VERSION = re.compile(r"diamond version (\d+)\.(\d+)\.(\d+)")
_SortKey = tuple[str, str, int, int]


class DiamondError(RuntimeError):
    """DIAMOND failed, or its version is unsupported."""


class DiamondNotFoundError(DiamondError):
    """The DIAMOND executable is not on PATH."""


def makedb_command(
    executable: str, fasta: Path, database: Path, config: DiamondSearch
) -> list[str]:
    """Command that builds a DIAMOND database from a protein FASTA file."""
    return [
        executable,
        "makedb",
        "--in",
        str(fasta),
        "--db",
        str(database),
        "--quiet",
        *_threads(config),
    ]


def blastp_command(
    executable: str, query: Path, database: Path, output: Path, config: DiamondSearch
) -> list[str]:
    """Command for the blastp search. Every limit is passed explicitly (spec 6.2).

    DIAMOND's per-HSP identity and coverage filters are never used.
    """
    return [
        executable,
        "blastp",
        "--query",
        str(query),
        "--db",
        str(database),
        "--out",
        str(output),
        "--outfmt",
        "6",
        *OUTPUT_FIELDS,
        "--evalue",
        str(config.hsp_evalue),
        "--max-target-seqs",
        str(config.max_target_seqs),
        "--max-hsps",
        str(config.max_hsps),
        *_sensitivity(config),
        "--quiet",
        *_threads(config),
    ]


def _sensitivity(config: DiamondSearch) -> list[str]:
    return [] if config.sensitivity == "default" else [f"--{config.sensitivity}"]


def _threads(config: DiamondSearch) -> list[str]:
    return [] if config.threads is None else ["--threads", str(config.threads)]


def parse_tabular(lines: Iterable[str], index: ProteinIndex) -> tuple[Alignment, ...]:
    """Parse DIAMOND tabular output (spec 6.4) into alignments.

    Tokens become canonical gene IDs and 1-based inclusive coordinates become
    0-based half-open. The result is sorted by query token, subject token,
    query start and subject start, so it does not depend on DIAMOND's order.
    """
    keyed: list[tuple[_SortKey, Alignment]] = []
    for number, line in enumerate(lines, start=1):
        if line.strip():
            keyed.append(_parse_line(number, line, index))
    return tuple(alignment for _, alignment in sorted(keyed, key=lambda item: item[0]))


def _parse_line(number: int, line: str, index: ProteinIndex) -> tuple[_SortKey, Alignment]:
    fields = line.rstrip("\n").split("\t")
    if len(fields) != len(OUTPUT_FIELDS):
        raise ValueError(
            f"DIAMOND output line {number}: expected {len(OUTPUT_FIELDS)} fields, "
            f"found {len(fields)}"
        )
    try:
        alignment = _alignment(fields, index)
    except KeyError as error:
        raise ValueError(f"DIAMOND output line {number}: unknown token {error}") from None
    except ValueError as error:
        raise ValueError(f"DIAMOND output line {number}: {error}") from None
    key = (fields[0], fields[1], alignment.query_start, alignment.subject_start)
    return key, alignment


def _alignment(fields: list[str], index: ProteinIndex) -> Alignment:
    def row(name: str) -> str:
        return fields[_COLUMN[name]]

    query_start, query_end = to_zero_based(int(row("qstart")), int(row("qend")))
    subject_start, subject_end = to_zero_based(int(row("sstart")), int(row("send")))
    return Alignment(
        query_id=index.gene_id_for(row("qseqid")),
        subject_id=index.gene_id_for(row("sseqid")),
        identity=float(row("pident")),
        alignment_length=int(row("length")),
        query_start=query_start,
        query_end=query_end,
        subject_start=subject_start,
        subject_end=subject_end,
        query_length=int(row("qlen")),
        subject_length=int(row("slen")),
        evalue=float(row("evalue")),
        bit_score=float(row("bitscore")),
    )


def parse_version(text: str) -> str:
    """Version string from ``diamond version`` output; 2.0 or newer is required."""
    match = _VERSION.search(text)
    if match is None:
        raise DiamondError(f"cannot read the DIAMOND version from {text.strip()!r}")
    major, minor, _ = (int(part) for part in match.groups())
    version = ".".join(match.groups())
    if (major, minor) < MINIMUM_VERSION:
        raise DiamondError(f"DIAMOND 2.0 or newer is required; found {version}")
    return version


def diamond_version(executable: str = "diamond") -> str:
    """Installed DIAMOND version. Raises if DIAMOND is missing or older than 2.0."""
    return parse_version(_run([_resolve(executable), "version"]))


def search(
    index: ProteinIndex,
    config: DiamondSearch,
    raw_output: str | Path,
    *,
    self_search: bool = False,
    executable: str = "diamond",
) -> tuple[Alignment, ...]:
    """Search reference proteins against target proteins, or against themselves.

    DIAMOND's tabular output is written unchanged to ``raw_output``; an empty
    file is written when one side has no searchable protein. In a self-search
    (``self_search=True``) hits of a protein to itself are dropped (spec 6.5).
    """
    raw = Path(raw_output)
    subjects = index.reference_fasta() if self_search else index.target_fasta()
    if not index.reference_fasta() or not subjects:
        raw.write_text("", encoding="utf-8")
        return ()
    program = _resolve(executable)
    diamond_version(program)
    with tempfile.TemporaryDirectory(prefix="homolocus-diamond-") as workdir:
        query_fasta = _write(Path(workdir) / "queries.faa", index.reference_fasta())
        subject_fasta = _write(Path(workdir) / "subjects.faa", subjects)
        database = Path(workdir) / "subjects"
        _run(makedb_command(program, subject_fasta, database, config))
        _run(blastp_command(program, query_fasta, database, raw, config))
    alignments = parse_tabular(raw.read_text(encoding="utf-8").splitlines(), index)
    if self_search:
        return tuple(a for a in alignments if a.query_id != a.subject_id)
    return alignments


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def _resolve(executable: str) -> str:
    path = shutil.which(executable)
    if path is None:
        raise DiamondNotFoundError(
            f"DIAMOND executable {executable!r} not found on PATH; install DIAMOND 2.0 or newer"
        )
    return path


def _run(command: list[str]) -> str:
    logger.debug("running: %s", shlex.join(command))
    try:
        completed = subprocess.run(command, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as error:
        raise DiamondError(
            f"DIAMOND {command[1]} failed (exit {error.returncode}): {error.stderr.strip()}"
        ) from error
    except OSError as error:
        raise DiamondError(f"could not run DIAMOND {command[1]}: {error}") from error
    return completed.stdout
