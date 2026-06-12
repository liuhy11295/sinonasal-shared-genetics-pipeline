import subprocess
import unittest
from pathlib import Path


class RepositorySafetyTests(unittest.TestCase):
    def test_tracked_files_do_not_contain_author_private_roots(self):
        root = Path(__file__).resolve().parents[1]
        tracked = subprocess.check_output(
            ["git", "ls-files"], cwd=root, text=True
        ).splitlines()
        forbidden = (
            "/" + "home" + "/" + "lhy",
            "/" + "platform_data" + "/" + "p_user",
        )
        findings = []
        for relative in tracked:
            path = root / relative
            if not path.is_file():
                continue
            text = path.read_text(encoding="utf-8", errors="ignore")
            for marker in forbidden:
                if marker in text:
                    findings.append(f"{relative}: {marker}")
        self.assertEqual(findings, [], "\n".join(findings))

    def test_no_tracked_private_keys(self):
        root = Path(__file__).resolve().parents[1]
        tracked = subprocess.check_output(
            ["git", "ls-files"], cwd=root, text=True
        ).splitlines()
        key_names = {"id_rsa", "id_ed25519", "id_ecdsa"}
        self.assertFalse(any(Path(path).name in key_names for path in tracked))

    def test_publication_figure_scripts_use_stable_names(self):
        root = Path(__file__).resolve().parents[1]
        script_dir = root / "11_publication_outputs" / "figures" / "scripts"
        expected = {
            "build_figure4_directional_panel.R",
            "build_main_figure_components.R",
            "common.R",
            "explore_figure1_layouts.R",
            "qa_all_figures.R",
            "qa_figure1.R",
            "render_all.R",
            "render_figure1.R",
            "render_figure2.R",
            "render_figure3.R",
            "render_figure4.R",
            "render_supplementary_figures.R",
        }
        self.assertEqual(
            {path.name for path in script_dir.glob("*.R")},
            expected,
        )


if __name__ == "__main__":
    unittest.main()
