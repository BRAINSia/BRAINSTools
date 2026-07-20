import re
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import cli_xml_lib as lib
import generate_tool_catalog as catalog
import reorder_cli_xml_elements as reorder


SAMPLE = """<?xml version="1.0" encoding="utf-8"?>
<executable>
  <category>Utilities.BRAINS</category>
  <title>Sample Tool</title>
  <description>Does a thing.</description>
  <acknowledgements><![CDATA[Funded by someone.]]></acknowledgements>
  <version>5.8.0</version>
  <!-- a comment that must survive -->
  <documentation-url>https://example.org</documentation-url>
  <parameters>
    <label>Input</label>
  </parameters>
</executable>
"""


class TestReorder(unittest.TestCase):
    def test_moves_acknowledgements_after_contributor_position(self):
        out = reorder.reorder_text(SAMPLE)
        order = [
            t
            for t in (
                "category",
                "title",
                "description",
                "version",
                "documentation-url",
                "acknowledgements",
            )
            if f"<{t}>" in out or f"<{t}><" in out
        ]
        positions = [out.index(f"<{t}") for t in order]
        self.assertEqual(positions, sorted(positions))
        self.assertLess(out.index("<version>"), out.index("<acknowledgements>"))

    def test_preserves_cdata_and_comments(self):
        out = reorder.reorder_text(SAMPLE)
        self.assertIn(
            "<acknowledgements><![CDATA[Funded by someone.]]></acknowledgements>", out
        )
        self.assertIn("<!-- a comment that must survive -->", out)

    def test_interleaved_comment_stays_adjacent_to_its_neighbours(self):
        # In SAMPLE, the comment sits between <version> and
        # <documentation-url>. A real reorder moves <acknowledgements> from
        # before <version> to after <documentation-url>; the comment must
        # stay pinned between <version> and <documentation-url>, not get
        # relocated to the end of the header alongside <acknowledgements>.
        out = reorder.reorder_text(SAMPLE)
        version_end = out.index("</version>") + len("</version>")
        comment_start = out.index("<!-- a comment that must survive -->")
        doc_url_start = out.index("<documentation-url>")
        ack_start = out.index("<acknowledgements>")
        self.assertLess(version_end, comment_start)
        self.assertLess(comment_start, doc_url_start)
        self.assertLess(doc_url_start, ack_start)

    def test_preserves_parameters_block_verbatim(self):
        out = reorder.reorder_text(SAMPLE)
        self.assertIn("<parameters>\n    <label>Input</label>\n  </parameters>", out)

    def test_is_idempotent(self):
        once = reorder.reorder_text(SAMPLE)
        twice = reorder.reorder_text(once)
        self.assertEqual(once, twice)

    def test_element_text_is_unchanged(self):
        before = re.findall(r"<(\w[\w-]*)>([^<]*)</\1>", SAMPLE)
        after = re.findall(r"<(\w[\w-]*)>([^<]*)</\1>", reorder.reorder_text(SAMPLE))
        self.assertEqual(sorted(before), sorted(after))

    def test_already_canonical_is_byte_identical(self):
        canonical = """<?xml version="1.0" encoding="utf-8"?>
<executable>
  <category>Utilities.BRAINS</category>
  <title>Sample Tool</title>
  <description>Does a thing.</description>
  <version>5.8.0</version>
  <!-- a comment that must survive -->
  <documentation-url>https://example.org</documentation-url>
  <acknowledgements><![CDATA[Funded by someone.]]></acknowledgements>
  <parameters>
    <label>Input</label>
  </parameters>
</executable>
"""
        self.assertEqual(reorder.reorder_text(canonical), canonical)


class TestDiscovery(unittest.TestCase):
    def test_descriptor_count(self):
        self.assertEqual(len(lib.descriptor_paths()), 69)

    def test_excludes_archive_testsuite_futurework(self):
        for p in lib.descriptor_paths():
            s = str(p)
            self.assertNotIn("ARCHIVE/", s)
            self.assertNotIn("/TestSuite/", s)
            self.assertNotIn("/future_work/", s)

    def test_published_count(self):
        published = [p for p in lib.descriptor_paths() if lib.is_published(p)]
        self.assertEqual(len(published), 67)

    def test_unbuilt_tools_are_not_published(self):
        for name in ("ComputeReflectiveCorrelationMetric", "gtractCoRegAnatomyBspline"):
            hits = [p for p in lib.descriptor_paths() if p.stem == name]
            self.assertEqual(len(hits), 1, f"{name} descriptor not found")
            self.assertFalse(lib.is_published(hits[0]), f"{name} must not be published")

    def test_built_tools_are_published(self):
        for name in (
            "BRAINSFit",
            "BRAINSStripRotation",
            "GenerateAverageLmkFile",
            "PerformMetricTest",
        ):
            hits = [p for p in lib.descriptor_paths() if p.stem == name]
            self.assertEqual(len(hits), 1, f"{name} descriptor not found")
            self.assertTrue(lib.is_published(hits[0]), f"{name} must be published")


class TestConformance(unittest.TestCase):
    def test_no_violations(self):
        violations = lib.check_all()
        self.assertEqual(violations, [], "\n".join(violations))


class TestCatalog(unittest.TestCase):
    def setUp(self):
        self.markdown = catalog.render_catalog()

    def test_includes_a_known_published_tool(self):
        self.assertIn("General Registration (BRAINS)", self.markdown)

    def test_excludes_unbuilt_tools(self):
        self.assertNotIn("Compute RC metric values", self.markdown)
        self.assertNotIn("Coregister B0 to Anatomy B-Spline", self.markdown)

    def test_groups_by_category_heading(self):
        self.assertIn("\n## Registration\n", self.markdown)

    def test_tool_count_matches_published_set(self):
        published = [p for p in lib.descriptor_paths() if lib.is_published(p)]
        self.assertEqual(self.markdown.count("\n### "), len(published))
        expected = {lib.header_text(p)["title"].strip() for p in published}
        rendered = {
            line[4:].strip()
            for line in self.markdown.splitlines()
            if line.startswith("### ")
        }
        self.assertEqual(rendered, expected)

    def test_is_deterministic(self):
        self.assertEqual(self.markdown, catalog.render_catalog())

    def test_duplicate_published_titles_do_not_crash_render(self):
        fake_paths = [lib.repo_root() / "a.xml", lib.repo_root() / "b.xml"]
        fake_headers = {
            fake_paths[0]: {
                "title": "Same Title",
                "category": "Registration",
                "description": "First.",
            },
            fake_paths[1]: {
                "title": "Same Title",
                "category": "Registration",
                "description": "Second.",
            },
        }
        real_descriptor_paths = lib.descriptor_paths
        real_is_published = lib.is_published
        real_header_text = lib.header_text
        lib.descriptor_paths = lambda: fake_paths
        lib.is_published = lambda path: True
        lib.header_text = lambda path: fake_headers[path]
        try:
            catalog.render_catalog()
        except TypeError as exc:  # pragma: no cover - regression guard
            self.fail(f"render_catalog crashed on duplicate titles: {exc}")
        finally:
            lib.descriptor_paths = real_descriptor_paths
            lib.is_published = real_is_published
            lib.header_text = real_header_text


if __name__ == "__main__":
    unittest.main()
