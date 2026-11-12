"""Compare recorded real layouts; never substitute for a device render."""
import hashlib
import json
import math
import sys
from pathlib import Path

root = Path(__file__).resolve().parent
files = {}
comparisons = {}
endpoint_count = 0
variant = 'current' if '--current' in sys.argv else 'refined' if '--refined' in sys.argv else 'after'
for mode in (("adaptive-",) if '--adaptive-only' in sys.argv else ("", "adaptive-")):
    before_path = root / f"v0.27.0-label-size-{mode}results.json"
    after_path = root / f"v0.27.0-label-size-{mode}{variant}-results.json"
    before, after = (json.loads(path.read_text()) for path in (before_path, after_path))
    for path in (before_path, after_path):
        files[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
    for name, current in after["observations"].items():
        previous = before["observations"][name]

        def labels(case):
            fields = ("layer_id", "source_row", "actual_label", "text_x", "text_y", "anchor_x", "anchor_y", "size")
            return [tuple(row.get(field) for field in fields) for row in case["visibility"]]

        def endpoints(case):
            groups = {}
            for row in case["quality"]["leaders"]:
                groups.setdefault(row["leader_id"], []).append(row)
            return {key: ((rows[0]["x0"], rows[0]["y0"]), (rows[-1]["x1"], rows[-1]["y1"]))
                    for key, rows in groups.items()}

        assert labels(previous) == labels(current), name
        old_ends, new_ends = endpoints(previous), endpoints(current)
        assert old_ends.keys() == new_ends.keys(), name
        for key, pairs in old_ends.items():
            # JSON records mm to four decimals; compare semantic endpoint
            # preservation, without locking implementation coordinates.
            assert all(math.dist(p, q) <= .001 for p, q in zip(pairs, new_ends[key])), (name, key)
            endpoint_count += 1
        comparisons[f"{mode}{name}"] = {
            "source_labels_and_native_positions_anchors_sizes_preserved": True,
            "directed_endpoint_pairs_preserved": len(old_ends),
            "quality": current["quality"]["summary"],
        }
result = {"status": "passed recorded layout preservation comparisons",
          "directed_endpoint_pairs_preserved": endpoint_count,
          "input_json_sha256": files, "comparisons": comparisons,
          "scope": "Recorded three-source compact fixed/adaptive profiles only; not arbitrary size or ink certification"}
(root / (f"v0.27.0-label-preservation-{variant}-results.json" if variant != 'after'
         else "v0.27.0-label-preservation-results.json")).write_text(json.dumps(result, indent=2) + "\n")
print(f"Preserved labels, native positions, anchors, sizes and {endpoint_count} directed endpoint pairs")
