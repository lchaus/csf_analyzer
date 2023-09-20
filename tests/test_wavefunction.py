import unittest
from analyzer.wavefunction import Wavefunction

class TestWavefunction(unittest.TestCase):

    def setUp(self):
        # TODO 
        pass

    def tearDown(self):
        # TODO CLEAN UP
        pass

    def test_read_orbitals(self):
        # Test the read_orbitals function
        wavefunction = Wavefunction()
        wavefunction.orbitals_file = "test_orbitals.dict"  # Provide a test file
        wavefunction.read_orbitals()
        self.assertIsNotNone(wavefunction.orbitals)

    def test_extract_data(self):
        # Test the extract_data function
        wavefunction = Wavefunction()
        wavefunction.orca_output_file = "test_orca_output.txt"  # Provide a test file
        wavefunction.extract_data()
        self.assertGreater(len(wavefunction.wavefunction), 0)

    # Add more test cases for a simple reproducible example

if __name__ == '__main__':
    unittest.main()
