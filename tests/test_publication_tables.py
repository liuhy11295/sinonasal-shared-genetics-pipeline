import csv
import unittest
from pathlib import Path


class PublicationTableTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source = (
            Path(__file__).resolve().parents[1]
            / "11_publication_outputs"
            / "tables"
            / "source_data"
        )

    def read(self, name):
        with (self.source / name).open(encoding="utf-8-sig", newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def test_expected_sheet_inventory_and_row_counts(self):
        expected = {
            "ST1_GWAS_and_power_QC.tsv": 61,
            "ST2_Analysis_parameters.tsv": 36,
            "ST2_Harmonization_QC.tsv": 31,
            "ST3_Final_31_pairs.tsv": 31,
            "ST4_Prioritized_genes.tsv": 123,
            "ST4_Prioritized_loci.tsv": 181,
            "ST5_Approved_drug_targets.tsv": 127,
            "ST5_Directional_evidence.tsv": 62,
            "ST5_Significant_pathways.tsv": 305,
        }
        self.assertEqual(
            {path.name for path in self.source.glob("*.tsv")},
            set(expected),
        )
        for name, count in expected.items():
            with self.subTest(name=name):
                self.assertEqual(len(self.read(name)), count)

    def test_final_pair_order_is_complete(self):
        rows = self.read("ST3_Final_31_pairs.tsv")
        self.assertEqual([row["P"] for row in rows], [f"P{i}" for i in range(1, 32)])
        self.assertEqual([int(row["pair_order"]) for row in rows], list(range(1, 32)))
        self.assertEqual(len({row["pair_id"] for row in rows}), 31)

    def test_primary_identifiers_are_populated(self):
        checks = {
            "ST1_GWAS_and_power_QC.tsv": ("trait_id",),
            "ST2_Analysis_parameters.tsv": ("section", "module_or_level"),
            "ST2_Harmonization_QC.tsv": ("trait_or_pair_id",),
            "ST3_Final_31_pairs.tsv": ("P", "pair_id", "pair_name"),
            "ST4_Prioritized_genes.tsv": ("gene_symbol",),
            "ST4_Prioritized_loci.tsv": ("P", "locus_id", "pair_id"),
            "ST5_Approved_drug_targets.tsv": ("drug_name", "target_gene_list"),
            "ST5_Directional_evidence.tsv": ("P", "pair_id", "direction"),
            "ST5_Significant_pathways.tsv": ("database", "term_id", "term_name"),
        }
        for name, fields in checks.items():
            rows = self.read(name)
            for field in fields:
                with self.subTest(name=name, field=field):
                    self.assertTrue(all(row.get(field, "").strip() for row in rows))


if __name__ == "__main__":
    unittest.main()
