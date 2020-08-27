import unittest
import sys

sys.path.append("c:\\Users\\winst\\Documents\\MEGA\\intern and work proj\\ASTAR_BII")

import nwk_to_fig


class TestNwkToFig(unittest.TestCase):
    base_path = "./test/samples/"

    def test_get_taxons_out(self):
        """
        Test that the taxons are being generated correctly from the nwk file
        """
        out = nwk_to_fig.get_taxons_out(self.base_path + "input.nwk")
        with open(self.base_path + "taxon_out.txt", "r") as f:
            expected = f.read()
        self.assertEqual(out, expected, "Taxons should be formatted in ascending order")

    def test_get_figtree(self):
        """
        Test that figtree is properly formatted
            - single apostrophe removed
            - label number changed to [&label=NUMBER]
        """
        print("test called")
        out = nwk_to_fig.get_figtree(self.base_path + "figtree_in.txt")
        with open(self.base_path + "figtree_out.txt", "r") as f:
            expected = f.read()
        with open(self.base_path + "figtree_actual_out.txt", "w") as f:
            f.write(out)
        self.assertEqual(out, expected, "figtree formatted incorrectly!")

    def test_get_setting(self):
        """
        Test that figtree settings are written correctly
        """
        out = nwk_to_fig.get_settings()
        with open(self.base_path + "settings_out.txt", "r") as f:
            expected = f.read()
        self.assertEqual(out, expected, "settings should be identical")

    def test_convert_to_nexus(self):
        """
        Runs get_taxons, get_figtree, and get_setting in order and combines the result into a final output
        """
        nwk_to_fig.convert_to_nexus(
            self.base_path + "input.nwk", self.base_path + "actual_tree"
        )
        with open(self.base_path + "actual_tree") as f:
            out = f.read()
        with open(self.base_path + "expected_tree") as f:
            expected = f.read()
        self.assertEqual(out, expected, "files should be identical")


if __name__ == "__main__":
    unittest.main()
