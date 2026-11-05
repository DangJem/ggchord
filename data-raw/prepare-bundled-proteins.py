"""Select complete public proteins and calculate their all-against-all BLASTp."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import tarfile

parser = argparse.ArgumentParser()
parser.add_argument("archive", type=Path)
parser.add_argument("--blastp", default="blastp")
args = parser.parse_args()
root = Path(__file__).resolve().parent.parent
out = root / "examples/protein"
out.mkdir(exist_ok=True)
ids = ["BVI_029_0020", "BVI_058_0011", "E4-10_023_0015", "ADZ16414.1"]
records = {}
members = []
with tarfile.open(args.archive) as archive:
    for name in ["emales/emales.faa", "emales/mavirus.faa"]:
        original = archive.extractfile(name).read()
        members.append({"member": name, "sha256": hashlib.sha256(original).hexdigest()})
        for record in original.decode().split(">"):
            if record and record.split()[0] in ids:
                records[record.split()[0]] = ">" + record
if set(records) != set(ids):
    raise ValueError("Complete selected protein records are missing")
fasta = out / "four-proteins.faa"
fasta.write_text("".join(records[identifier] for identifier in ids))
results = out / "all-vs-all-blastp.tsv"
arguments = ["-query", str(fasta), "-subject", str(fasta), "-outfmt", "6",
             "-evalue", "1e-5", "-max_hsps", "1", "-out", str(results)]
subprocess.run([args.blastp, *arguments], check=True)
rows = [line.split("\t") for line in results.read_text().splitlines()]
pairs = {(row[0], row[1]) for row in rows}
if pairs != {(q, s) for q in ids for s in ids}:
    raise ValueError("Selected records do not have a complete observed comparison")
provenance = {
    "url": "https://raw.githubusercontent.com/thackl/gggenomes/master/data-raw/emales.tgz",
    "archive_sha256": hashlib.sha256(args.archive.read_bytes()).hexdigest(),
    "members": members, "protein_ids": ids,
    "blastp_version": subprocess.check_output([args.blastp, "-version"], text=True).strip(),
    "blastp_arguments": ["-query", fasta.name, "-subject", fasta.name, "-outfmt", "6",
                        "-evalue", "1e-5", "-max_hsps", "1", "-out", results.name],
    "fasta_sha256": hashlib.sha256(fasta.read_bytes()).hexdigest(),
    "blastp_sha256": hashlib.sha256(results.read_bytes()).hexdigest(),
    "selection": "Three complete EMALE proteins from the cluster containing BVI_023A_0009, plus the matching complete Mavirus ADZ16414.1 protein. Bundled links retain one direction for each of six unordered pairs. All 16 computed rows are retained here."
}
(out / "sources.json").write_text(json.dumps(provenance, indent=2) + "\n")
