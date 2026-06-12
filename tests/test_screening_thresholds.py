from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]


def load_module(path):
    spec = spec_from_file_location(path.stem, path)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class ScreeningThresholdTests(unittest.TestCase):
    def test_bh_adjust_preserves_order_and_monotonicity(self):
        module = load_module(ROOT / "01_ldsc_screen/scripts/run_ldsc_rg.py")
        adjusted = module.bh_adjust([0.01, 0.04, 0.03, 0.20])
        expected = [0.04, 0.05333333333333334, 0.05333333333333334, 0.2]
        for observed, target in zip(adjusted, expected):
            self.assertAlmostEqual(observed, target)

    def test_bh_adjust_caps_at_one(self):
        module = load_module(ROOT / "01_ldsc_screen/scripts/run_ldsc_rg.py")
        self.assertEqual(module.bh_adjust([0.8, 0.9]), [0.9, 0.9])


if __name__ == "__main__":
    unittest.main()
