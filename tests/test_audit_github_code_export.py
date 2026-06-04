from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from scripts.audit_github_code_export import audit_export


class AuditGithubCodeExportTests(unittest.TestCase):
    def test_code_and_documentation_export_passes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "phase0_pipeline").mkdir()
            (root / "docs").mkdir()
            (root / "phase0_pipeline" / "module.py").write_text("print('ok')\n", encoding="utf-8")
            (root / "docs" / "README.md").write_text("# Documentation\n", encoding="utf-8")

            self.assertEqual(audit_export(root), [])

    def test_git_metadata_and_ignored_caches_are_not_publication_payload(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / ".git" / "objects").mkdir(parents=True)
            (root / "__pycache__").mkdir()
            (root / ".git" / "objects" / "object.log").write_text(
                "-----BEGIN " + "OPENSSH PRIVATE KEY-----\n", encoding="utf-8"
            )
            (root / "__pycache__" / "ignored.pyc").write_text("cache\n", encoding="utf-8")

            self.assertEqual(audit_export(root), [])

    def test_data_results_logs_and_credentials_are_blocked(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "data").mkdir()
            (root / "results").mkdir()
            (root / "data" / "gwas.tsv").write_text("SNP\n", encoding="utf-8")
            (root / "results" / "run.log").write_text("run\n", encoding="utf-8")
            (root / "key.txt").write_text(
                "-----BEGIN " + "OPENSSH PRIVATE KEY-----\n", encoding="utf-8"
            )

            failures = audit_export(root)
            self.assertTrue(any("forbidden directory" in failure for failure in failures))
            self.assertTrue(any("forbidden data/runtime suffix" in failure for failure in failures))
            self.assertTrue(any("credential material" in failure for failure in failures))


if __name__ == "__main__":
    unittest.main()
