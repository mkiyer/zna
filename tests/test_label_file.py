"""Tests for YAML label file loading and CLI merge logic."""
import os
import tempfile
import unittest

import yaml

from zna.cli import load_label_file, merge_label_defs, build_label_defs
from zna.dtypes import LabelDef, parse_dtype


def _write_yaml(data: dict, path: str) -> None:
    with open(path, "w") as fh:
        yaml.dump(data, fh, default_flow_style=False)


class TestLoadLabelFile(unittest.TestCase):
    """Tests for load_label_file()."""

    def test_basic_load(self):
        with tempfile.NamedTemporaryFile(suffix=".yaml", mode="w", delete=False) as f:
            yaml.dump({"labels": [
                {"name": "NM", "type": "C", "description": "Edit distance"},
                {"name": "AS", "type": "i"},
            ]}, f)
            path = f.name
        try:
            defs = load_label_file(path)
            self.assertEqual(len(defs), 2)
            self.assertEqual(defs[0].name, "NM")
            self.assertEqual(defs[0].dtype.code, "C")
            self.assertEqual(defs[0].description, "Edit distance")
            self.assertEqual(defs[0].label_id, 0)
            self.assertEqual(defs[1].name, "AS")
            self.assertEqual(defs[1].dtype.code, "i")
            self.assertEqual(defs[1].description, "")
            self.assertEqual(defs[1].label_id, 1)
        finally:
            os.unlink(path)

    def test_missing_integer(self):
        with tempfile.NamedTemporaryFile(suffix=".yaml", mode="w", delete=False) as f:
            yaml.dump({"labels": [
                {"name": "NM", "type": "C", "missing": 255},
            ]}, f)
            path = f.name
        try:
            defs = load_label_file(path)
            self.assertEqual(defs[0].missing, 255)
        finally:
            os.unlink(path)

    def test_missing_float(self):
        with tempfile.NamedTemporaryFile(suffix=".yaml", mode="w", delete=False) as f:
            yaml.dump({"labels": [
                {"name": "de", "type": "f", "missing": -1.0},
            ]}, f)
            path = f.name
        try:
            defs = load_label_file(path)
            self.assertAlmostEqual(defs[0].missing, -1.0)
        finally:
            os.unlink(path)

    def test_missing_char(self):
        """Type A with single-char missing → stored as ord()."""
        with tempfile.NamedTemporaryFile(suffix=".yaml", mode="w", delete=False) as f:
            yaml.dump({"labels": [
                {"name": "tp", "type": "A", "missing": "?"},
            ]}, f)
            path = f.name
        try:
            defs = load_label_file(path)
            self.assertEqual(defs[0].missing, ord("?"))
        finally:
            os.unlink(path)

    def test_all_dtypes(self):
        """All 11 type codes should parse correctly."""
        entries = [{"name": f"T{c}", "type": c} for c in "AcCsSiIfdqQ"]
        with tempfile.NamedTemporaryFile(suffix=".yaml", mode="w", delete=False) as f:
            yaml.dump({"labels": entries}, f)
            path = f.name
        try:
            defs = load_label_file(path)
            self.assertEqual(len(defs), 11)
            for i, code in enumerate("AcCsSiIfdqQ"):
                self.assertEqual(defs[i].dtype.code, code)
        finally:
            os.unlink(path)

    def test_label_ids_sequential(self):
        entries = [{"name": f"L{i}", "type": "C"} for i in range(5)]
        with tempfile.NamedTemporaryFile(suffix=".yaml", mode="w", delete=False) as f:
            yaml.dump({"labels": entries}, f)
            path = f.name
        try:
            defs = load_label_file(path)
            for i, d in enumerate(defs):
                self.assertEqual(d.label_id, i)
        finally:
            os.unlink(path)

    # --- Error cases ---

    def test_file_not_found(self):
        with self.assertRaises(SystemExit) as ctx:
            load_label_file("/nonexistent/path.yaml")
        self.assertIn("not found", str(ctx.exception))

    def test_missing_labels_key(self):
        with tempfile.NamedTemporaryFile(suffix=".yaml", mode="w", delete=False) as f:
            yaml.dump({"definitions": []}, f)
            path = f.name
        try:
            with self.assertRaises(SystemExit) as ctx:
                load_label_file(path)
            self.assertIn("labels", str(ctx.exception))
        finally:
            os.unlink(path)

    def test_labels_not_a_list(self):
        with tempfile.NamedTemporaryFile(suffix=".yaml", mode="w", delete=False) as f:
            yaml.dump({"labels": "not a list"}, f)
            path = f.name
        try:
            with self.assertRaises(SystemExit) as ctx:
                load_label_file(path)
            self.assertIn("must be a list", str(ctx.exception))
        finally:
            os.unlink(path)

    def test_entry_not_a_mapping(self):
        with tempfile.NamedTemporaryFile(suffix=".yaml", mode="w", delete=False) as f:
            yaml.dump({"labels": ["not_a_dict"]}, f)
            path = f.name
        try:
            with self.assertRaises(SystemExit) as ctx:
                load_label_file(path)
            self.assertIn("must be a mapping", str(ctx.exception))
        finally:
            os.unlink(path)

    def test_missing_name(self):
        with tempfile.NamedTemporaryFile(suffix=".yaml", mode="w", delete=False) as f:
            yaml.dump({"labels": [{"type": "C"}]}, f)
            path = f.name
        try:
            with self.assertRaises(SystemExit) as ctx:
                load_label_file(path)
            self.assertIn("missing 'name'", str(ctx.exception))
        finally:
            os.unlink(path)

    def test_missing_type(self):
        with tempfile.NamedTemporaryFile(suffix=".yaml", mode="w", delete=False) as f:
            yaml.dump({"labels": [{"name": "NM"}]}, f)
            path = f.name
        try:
            with self.assertRaises(SystemExit) as ctx:
                load_label_file(path)
            self.assertIn("missing 'type'", str(ctx.exception))
        finally:
            os.unlink(path)

    def test_invalid_dtype(self):
        with tempfile.NamedTemporaryFile(suffix=".yaml", mode="w", delete=False) as f:
            yaml.dump({"labels": [{"name": "NM", "type": "Z"}]}, f)
            path = f.name
        try:
            with self.assertRaises(SystemExit) as ctx:
                load_label_file(path)
            self.assertIn("unknown dtype", str(ctx.exception))
        finally:
            os.unlink(path)

    def test_invalid_missing_value(self):
        with tempfile.NamedTemporaryFile(suffix=".yaml", mode="w", delete=False) as f:
            yaml.dump({"labels": [{"name": "NM", "type": "C", "missing": "not_a_num"}]}, f)
            path = f.name
        try:
            with self.assertRaises(SystemExit) as ctx:
                load_label_file(path)
            self.assertIn("invalid missing value", str(ctx.exception))
        finally:
            os.unlink(path)

    def test_tag_field(self):
        """Optional 'tag' field should be loaded into LabelDef."""
        with tempfile.NamedTemporaryFile(suffix=".yaml", mode="w", delete=False) as f:
            yaml.dump({"labels": [
                {"name": "edit_dist", "type": "C", "tag": "NM"},
                {"name": "AS", "type": "i"},
            ]}, f)
            path = f.name
        try:
            defs = load_label_file(path)
            self.assertEqual(defs[0].name, "edit_dist")
            self.assertEqual(defs[0].tag, "NM")
            self.assertEqual(defs[0].effective_tag, "NM")
            self.assertIsNone(defs[1].tag)
            self.assertEqual(defs[1].effective_tag, "AS")
        finally:
            os.unlink(path)


class TestMergeLabelDefs(unittest.TestCase):
    """Tests for merge_label_defs()."""

    def _make_yaml_defs(self) -> tuple[LabelDef, ...]:
        return (
            LabelDef(0, "NM", "Edit distance", parse_dtype("C"), missing=255),
            LabelDef(1, "AS", "Alignment score", parse_dtype("i"), missing=-1),
            LabelDef(2, "MQ", "Mapping quality", parse_dtype("C")),
        )

    def test_no_overrides(self):
        yaml_defs = self._make_yaml_defs()
        result = merge_label_defs(yaml_defs, [], [])
        self.assertEqual(len(result), 3)
        for orig, merged in zip(yaml_defs, result):
            self.assertEqual(orig.name, merged.name)
            self.assertEqual(orig.dtype.code, merged.dtype.code)
            self.assertEqual(orig.description, merged.description)
            self.assertEqual(orig.missing, merged.missing)

    def test_override_type(self):
        yaml_defs = self._make_yaml_defs()
        result = merge_label_defs(yaml_defs, ["AS:s"], [])
        # AS should now be int16 instead of int32
        as_def = [d for d in result if d.name == "AS"][0]
        self.assertEqual(as_def.dtype.code, "s")
        # Other labels unchanged
        self.assertEqual(result[0].dtype.code, "C")
        self.assertEqual(result[2].dtype.code, "C")

    def test_override_description(self):
        yaml_defs = self._make_yaml_defs()
        result = merge_label_defs(yaml_defs, [], ["NM:Number of mismatches"])
        nm_def = [d for d in result if d.name == "NM"][0]
        self.assertEqual(nm_def.description, "Number of mismatches")

    def test_append_new_label(self):
        yaml_defs = self._make_yaml_defs()
        result = merge_label_defs(yaml_defs, ["de:f"], ["de:Divergence"])
        self.assertEqual(len(result), 4)
        de_def = result[3]
        self.assertEqual(de_def.name, "de")
        self.assertEqual(de_def.dtype.code, "f")
        self.assertEqual(de_def.description, "Divergence")

    def test_override_and_append(self):
        yaml_defs = self._make_yaml_defs()
        result = merge_label_defs(yaml_defs, ["AS:s", "XS:s"], ["XS:Subopt score"])
        self.assertEqual(len(result), 4)
        # AS overridden
        as_def = [d for d in result if d.name == "AS"][0]
        self.assertEqual(as_def.dtype.code, "s")
        # XS appended
        xs_def = [d for d in result if d.name == "XS"][0]
        self.assertEqual(xs_def.dtype.code, "s")
        self.assertEqual(xs_def.description, "Subopt score")

    def test_label_ids_sequential(self):
        yaml_defs = self._make_yaml_defs()
        result = merge_label_defs(yaml_defs, ["XS:s", "de:f"], [])
        for i, d in enumerate(result):
            self.assertEqual(d.label_id, i)

    def test_missing_preserved_from_yaml(self):
        yaml_defs = self._make_yaml_defs()
        result = merge_label_defs(yaml_defs, ["NM:s"], [])
        # NM type changes but missing is preserved from YAML
        nm_def = [d for d in result if d.name == "NM"][0]
        self.assertEqual(nm_def.dtype.code, "s")
        self.assertEqual(nm_def.missing, 255)

    def test_tag_preserved_from_yaml(self):
        """Tags from YAML should be preserved through merge."""
        yaml_defs = (
            LabelDef(0, "edit_dist", "Edit distance", parse_dtype("C"), tag="NM"),
            LabelDef(1, "AS", "Alignment score", parse_dtype("i")),
        )
        result = merge_label_defs(yaml_defs, [], [])
        self.assertEqual(result[0].tag, "NM")
        self.assertEqual(result[0].effective_tag, "NM")
        self.assertIsNone(result[1].tag)

    def test_tag_override_from_cli(self):
        """CLI 3-part spec should override YAML tag."""
        yaml_defs = (
            LabelDef(0, "edit_dist", "Edit distance", parse_dtype("C"), tag="NM"),
        )
        result = merge_label_defs(yaml_defs, ["edit_dist:C:nm"], [])
        self.assertEqual(result[0].tag, "nm")

    def test_cli_tag_appended(self):
        """CLI-only labels with tags should work."""
        yaml_defs = (
            LabelDef(0, "AS", "Score", parse_dtype("i")),
        )
        result = merge_label_defs(yaml_defs, ["edit_dist:C:NM"], [])
        self.assertEqual(len(result), 2)
        self.assertEqual(result[1].name, "edit_dist")
        self.assertEqual(result[1].tag, "NM")


class TestLoadExampleFile(unittest.TestCase):
    """Test that the shipped examples/labels.yaml loads correctly."""

    def test_example_loads(self):
        example_path = os.path.join(
            os.path.dirname(__file__), "..", "examples", "labels.yaml"
        )
        if not os.path.exists(example_path):
            self.skipTest("examples/labels.yaml not found")
        defs = load_label_file(example_path)
        self.assertGreater(len(defs), 0)
        names = [d.name for d in defs]
        self.assertIn("NM", names)
        self.assertIn("AS", names)
        # IDs sequential
        for i, d in enumerate(defs):
            self.assertEqual(d.label_id, i)


if __name__ == "__main__":
    unittest.main()
