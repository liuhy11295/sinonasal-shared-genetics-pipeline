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

    def test_workflow_script_inventories_are_minimal_and_stable(self):
        root = Path(__file__).resolve().parents[1]
        expected = {
            "00_environment_setup/scripts": {
                "setup_environment.slurm",
                "install_method_tools.slurm",
                "create_input_manifest.sh",
            },
            "01_ldsc_screen/scripts": {"run_ldsc.py"},
            "02_lava_screen/scripts": {"run_lava.R", "aggregate_results.py"},
            "03_mtag_after_ldsc_lava/scripts": {
                "prepare_inputs.py",
                "run_mtag.py",
            },
            "04_screen_placo_cpassoc_genomewide/scripts": {
                "run_pair.R",
                "aggregate_results.py",
            },
            "05_coloc_abf/scripts": {"run_locus.R", "aggregate_results.py"},
            "06_coloc_susie/scripts": {"run_locus.R", "aggregate_results.py"},
            "07_magma_gene_pathway/scripts": {
                "prepare_resources.sh",
                "merge_ld_reference.sh",
                "run_pair.py",
                "aggregate_results.py",
            },
            "08_twas/scripts": {
                "select_candidate_genes.py",
                "prepare_sumstats.py",
                "download_weights.sh",
                "prepare_weights.py",
                "run_twas.sh",
                "run_twas_job.sh",
                "aggregate_results.py",
            },
            "09_lcv_mr_secondary/scripts": {
                "prepare_inputs.py",
                "run_lcv.R",
                "run_mr.R",
                "finalize_outputs.py",
            },
            "10_downstream_annotation/scripts": {
                "run_pathway_ora.R",
                "run_cell_marker_ora.py",
                "annotate_fuma_concordance.py",
                "annotate_drug_targets.py",
            },
            "11_publication_outputs/figures/scripts": {
                "common.R",
                "render_all.R",
                "render_figure1.R",
                "render_figure2.R",
                "render_figure3.R",
                "render_figure4.R",
                "render_supplementary_figures.R",
            },
            "11_publication_outputs/tables/scripts": {"build_final_tables.py"},
        }
        for relative, filenames in expected.items():
            with self.subTest(directory=relative):
                directory = root / relative
                observed = {
                    path.name
                    for path in directory.iterdir()
                    if path.is_file()
                }
                self.assertEqual(observed, filenames)


if __name__ == "__main__":
    unittest.main()
