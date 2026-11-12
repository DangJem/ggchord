#!/usr/bin/env python3
"""Explore existing GFA1 path coordinates through optional gfapy, not a reader API.

Usage: python prototype-graph-paths.py /path/to/unchanged/LPA.gfa
Input: pinned ODGI v0.8.4 human LPA graph, commit
3f59d596cff3bbe9ddff3862e3212812254496e5, test/LPA.gfa.
Only zero-overlap, sequence-bearing paths are projected. No graph construction,
alignment, pangenome inference, all-pairs expansion or ODGI invocation occurs.
"""
import collections
import csv
import hashlib
import importlib.metadata
import json
import sys
from pathlib import Path

import gfapy

root = Path(__file__).resolve().parents[2]
folder = root / 'tools/long-term-audit'
source = Path(sys.argv[1])
expected = '9017b433f35b604bdcafd9339318f1649585bceda4ccf263f1ee1f16adf0cdf3'
assert hashlib.sha256(source.read_bytes()).hexdigest() == expected, 'Source bytes changed'
graph = gfapy.Gfa.from_file(str(source), version='gfa1')
assert all(str(link.overlap) == '0M' for link in graph.dovetails), 'Nonzero/unknown overlaps need a different projection'
nodes = {line.name: line for line in graph.segments}
assert all(line.sequence != '*' for line in nodes.values()), 'Unknown node sequence cannot be filled'
# Use the decoder's original edge orientation. Every projected transition must
# occur in the graph, either as written or in its reverse-complement traversal.
edges = set()
for link in graph.dovetails:
    fields = str(link).split('\t')
    _, left, lo, right, ro = fields[:5]
    edges.add((left, lo, right, ro))
    reverse = {'+': '-', '-': '+'}
    edges.add((right, reverse[ro], left, reverse[lo]))
selected = {'chm13__LPA__tig00000001', 'HG002__LPA__tig00000001'}
sequences = []
steps = []
all_steps = 0
repeated = {}
for path in graph.paths:
    cursor = 0
    occurrences = collections.Counter()
    previous = None
    for index, oriented in enumerate(path.segment_names, 1):
        node = oriented.line
        length = len(node.sequence)
        assert length == node.length and length > 0
        if previous is not None:
            assert (previous[0], previous[1], node.name, oriented.orient) in edges
        occurrences[node.name] += 1
        # A declared first-100-step display subset, not an altered source path.
        if path.name in selected and index <= 100:
            steps.append(dict(accver=path.name, start=cursor + 1, end=cursor + length,
                              strand=oriented.orient, node_id=node.name,
                              source_step=index, node_occurrence=occurrences[node.name]))
        cursor += length
        previous = (node.name, oriented.orient)
        all_steps += 1
    # Complete length of this supplied path, not of an inferred chromosome.
    sequences.append(dict(accver=path.name, length=cursor, sequence_unit='bp'))
    repeated[path.name] = sum(count > 1 for count in occurrences.values())
assert {row['accver'] for row in steps} == selected

def save(name, rows):
    with (folder / name).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

save('v0.27.0-graph-path-sources.csv', sequences)
save('v0.27.0-graph-path-steps.csv', steps)
# No links are emitted: a repeated node ID alone does not uniquely choose two
# occurrence anchors, or prove orthology/identity/read support.
report = dict(
    status='passed exploratory projection',
    source=dict(url='https://raw.githubusercontent.com/pangenome/odgi/3f59d596cff3bbe9ddff3862e3212812254496e5/test/LPA.gfa',
                sha256=expected, license='MIT, copyright 2019 Erik Garrison',
                description='Existing human Lipoprotein A graph; source names are opaque path IDs'),
    decoder=dict(package='gfapy', version=importlib.metadata.version('gfapy')),
    nodes=len(graph.segments), links=len(graph.dovetails), paths=len(graph.paths),
    source_steps=all_steps, selected_first_steps=len(steps),
    repeated_node_ids_by_path=repeated,
    coordinates='1-based-inclusive cumulative coordinates of original zero-overlap path steps; path strand is retained',
    scientific_pair_links_emitted=0,
    decision='Use existing sequence/feature views for an explicit path subset. No generic GFA/ODGI reader or graph Coord is justified by this prototype.',
    unsupported=['nonzero or unknown overlaps', 'missing node lengths/characters', 'GFA W/J or GFA2',
                 'native ODGI binary decoding', 'automatic occurrence pairing or orthology'],
    no_analysis_engine=True)
(folder / 'v0.27.0-graph-path-results.json').write_text(json.dumps(report, indent=2) + '\n')
print('Graph path prototype:', len(graph.paths), 'paths;', all_steps, 'source steps;', len(steps), 'selected steps; no pair links')
