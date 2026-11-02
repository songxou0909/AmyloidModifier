import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
import zipfile

spec = importlib.util.spec_from_file_location("release_wheels", Path(__file__).parents[1] / "scripts" / "release_wheels.py")
release = importlib.util.module_from_spec(spec)
spec.loader.exec_module(release)


class WheelValidationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        (self.root / "wheels").mkdir()
        (self.root / "README.md").write_text("<!-- WHEELS:START -->\n<!-- WHEELS:END -->", encoding="utf-8")

    def wheel(self, version, source="print('test')", prefix="chimerax_amyloidmodifier"):
        path = self.root / "wheels" / f"{prefix}-{version}-py3-none-any.whl"
        with zipfile.ZipFile(path, "w") as archive:
            archive.writestr(f"chimerax_amyloidmodifier-{version}.dist-info/METADATA",
                             f"Name: ChimeraX-AmyloidModifier\nVersion: {version}\nRequires-Dist: numpy\n")
            archive.writestr("chimerax/amyloidmodifier/tool.py", source)
        return path

    def test_numeric_order_and_original_filename_case(self):
        self.wheel("1.9")
        self.wheel("1.10", prefix="ChimeraX_AmyloidModifier")
        self.assertEqual([r["version"] for r in release.catalogue(self.root)], ["1.9", "1.10"])

    def test_reject_replacing_or_removing_a_released_wheel(self):
        path = self.wheel("2.2.2")
        release.prepare(self.root)
        self.wheel("2.2.2", source="changed")
        with self.assertRaisesRegex(ValueError, "Do not replace"):
            release.catalogue(self.root)
        path.unlink()
        with self.assertRaisesRegex(ValueError, "missing"):
            release.catalogue(self.root)

    def test_reject_duplicate_or_older_version(self):
        self.wheel("2.2.2")
        release.prepare(self.root)
        duplicate = self.wheel("2.2.2.0")
        with self.assertRaises(ValueError):
            release.catalogue(self.root)
        duplicate.unlink()
        self.wheel("2.2.1")
        with self.assertRaisesRegex(ValueError, "newer version"):
            release.catalogue(self.root)

    def test_new_version_clears_old_doi_and_updates_citation(self):
        self.wheel("2.2.2")
        release.prepare(self.root)
        cff = self.root / "CITATION.cff"
        cff.write_text(cff.read_text(encoding="utf-8") + 'doi: "10.5281/zenodo.123"\n', encoding="utf-8")
        release.prepare(self.root)
        self.assertIn('doi: "10.5281/zenodo.123"', cff.read_text(encoding="utf-8"))
        self.wheel("2.3.0")
        release.prepare(self.root)
        self.assertNotIn("doi:", cff.read_text(encoding="utf-8"))
        self.assertIn('version: "2.3.0"', cff.read_text(encoding="utf-8"))
        self.assertEqual(json.loads((self.root / ".zenodo.json").read_text())["version"], "2.3.0")

    def test_reject_filename_metadata_mismatch(self):
        path = self.wheel("2.2.2")
        path.rename(self.root / "wheels" / path.name.replace("2.2.2", "2.3.0"))
        with self.assertRaisesRegex(ValueError, "disagree"):
            release.catalogue(self.root)


if __name__ == "__main__":
    unittest.main()
