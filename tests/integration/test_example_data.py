"""Integration tests for pathSQE using example data from the repository"""

import pytest
import numpy as np
import os
from pathlib import Path

from pathSQE import helper


class TestExampleDataLoading:
    """Tests that verify pathSQE can load and process example data"""

    @pytest.fixture
    def si_example_dir(self):
        """Path to Si example data"""
        return Path(__file__).parent.parent.parent / "examples" / "Si_ARCS_publicData"

    @pytest.fixture
    def bz_coverage_file(self, si_example_dir):
        """Path to BZ coverage file"""
        return si_example_dir / "BZ_coverage_sorted.txt"

    def test_bz_coverage_file_exists(self, bz_coverage_file):
        """Test that example BZ coverage file exists"""
        assert bz_coverage_file.exists(), (
            f"BZ coverage file not found at {bz_coverage_file}"
        )

    def test_load_bz_coverage_data(self, bz_coverage_file):
        """Test loading BZ coverage data from example"""
        # Attempt to load with header
        try:
            data = np.genfromtxt(
                str(bz_coverage_file), skip_header=1, usecols=(0, 1), ndmin=2
            )
        except ValueError:
            # If file format is different, try alternative parsing
            with open(bz_coverage_file, "r") as f:
                lines = f.readlines()
                data_lines = [l.strip() for l in lines[1:] if l.strip()]
                assert len(data_lines) > 0, "No data found in BZ coverage file"

    def test_bz_coverage_data_quality(self, bz_coverage_file):
        """Test that BZ coverage data contains valid entries"""
        with open(bz_coverage_file, "r") as f:
            lines = f.readlines()

        # Skip header
        data_lines = lines[1:]
        assert len(data_lines) > 0, (
            "BZ coverage file should contain at least one data line"
        )

        # Check that each line has data
        for line in data_lines:
            if line.strip():  # Skip empty lines
                parts = line.split()
                assert len(parts) >= 2, f"Line should have at least 2 columns: {line}"


class TestPathSQEInputFiles:
    """Tests that verify pathSQE input files are well-formed"""

    @pytest.fixture
    def example_files(self):
        """Get all pathSQE input and define_data files recursively"""
        examples_path = Path(__file__).parent.parent.parent / "examples"
        input_files = list(examples_path.rglob("pathSQE_input.py"))
        define_files = list(examples_path.rglob("define_data.py"))
        bz_files = list(examples_path.rglob("BZ_coverage_sorted.txt"))

        return {"input": input_files, "define": define_files, "bz": bz_files}

    def test_example_input_files_exist(self, example_files):
        """Test that pathSQE_input.py files exist in examples"""
        assert len(example_files["input"]) > 0, (
            "No pathSQE_input.py files found in examples"
        )

    def test_example_define_files_exist(self, example_files):
        """Test that define_data.py files exist in examples"""
        assert len(example_files["define"]) > 0, (
            "No define_data.py files found in examples"
        )

    def test_example_bz_files_exist(self, example_files):
        """Test that BZ_coverage_sorted.txt files exist in examples"""
        assert len(example_files["bz"]) > 0, (
            "No BZ_coverage_sorted.txt files found in examples"
        )

    def test_pathsqe_input_valid_python(self, example_files):
        """Test that pathSQE_input.py files are valid Python"""
        for input_file in example_files["input"]:
            with open(input_file, "r") as f:
                code = f.read()
            try:
                compile(code, str(input_file), "exec")
            except SyntaxError as e:
                pytest.fail(f"Syntax error in {input_file}: {e}")

    def test_define_data_valid_python(self, example_files):
        """Test that define_data.py files are valid Python"""
        for data_file in example_files["define"]:
            with open(data_file, "r") as f:
                code = f.read()
            try:
                compile(code, str(data_file), "exec")
            except SyntaxError as e:
                pytest.fail(f"Syntax error in {data_file}: {e}")


class TestHelperWithExampleData:
    """Integration tests using helper functions with example data"""

    @pytest.fixture
    def si_example_dir(self):
        """Path to Si example data"""
        return Path(__file__).parent.parent.parent / "examples" / "Si_ARCS_publicData"

    def test_conv_desc_string_with_bz_coordinates(self):
        """Test conv_to_desc_string with realistic BZ coordinates"""
        # Example BZ coordinates from the code
        bz_coords = [1.0, 2.0, 3.0]
        result = helper.conv_to_desc_string(bz_coords)
        assert result == "1.0,2.0,3.0"
        assert isinstance(result, str)

    def test_conv_desc_string_with_matrix_rows(self):
        """Test conv_to_desc_string with transformation matrix rows"""
        # Example from pathSQE documentation
        # Primitive to mantid transformation matrix example
        row = np.array([-1, 1, 1])
        result = helper.conv_to_desc_string(row)
        assert result == "-1,1,1"

    def test_make_slice_desc_realistic_parameters(self):
        """Test make_slice_desc with realistic pathSQE parameters"""
        # Realistic parameters based on examples
        pathSQE_params = {
            "E bins": "0,0.5,80"  # From Si example
        }

        # Example transformation matrix from Si documentation
        prim_to_mantid = np.array([[-1, 1, 1], [1, -1, 1], [1, 1, -1]])

        q_dims_and_bins = [
            prim_to_mantid[0],
            prim_to_mantid[1],
            prim_to_mantid[2],
            [-0.5, 0.5, 0.5],
            [-0.5, 0.5, 0.5],
            [-0.5, 0.5, 0.5],
        ]

        point1 = np.array([0.0, 0.0, 0.0])  # Gamma point
        point2 = np.array([0.0, 1.0, 0.0])  # X point
        path_seg = ("Gamma", "X")  # From Si example

        result = helper.make_slice_desc(
            pathSQE_params, q_dims_and_bins, point1, point2, path_seg
        )

        # Verify expected keys are present
        assert "Name" in result
        assert "QDimension0" in result
        assert "QDimension1" in result
        assert "QDimension2" in result
        assert result["seg_start_name"] == "Gamma"
        assert result["seg_end_name"] == "X"

    def test_slice_desc_with_realistic_k_path_points(self):
        """Test make_slice_desc with multiple realistic K-path points"""
        # From Si example: [('Gamma', 'X'), ('X', 'U'), ('K', 'Gamma'), ...]
        pathSQE_params = {"E bins": "0,0.5,80"}
        prim_to_mantid = np.array([[-1, 1, 1], [1, -1, 1], [1, 1, -1]])
        q_dims_and_bins = [
            prim_to_mantid[0],
            prim_to_mantid[1],
            prim_to_mantid[2],
            [-0.5, 0.5, 0.5],
            [-0.5, 0.5, 0.5],
            [-0.5, 0.5, 0.5],
        ]

        # Examples from Si documentation
        points = {
            "Gamma": [0.0, 0.0, 0.0],
            "X": [0, 1, 0],
            "L": [0.5, 0.5, 0.5],
            "W": [0.5, 1, 0],
            "K": [0.75, 0.75, 0],
            "U": [0.25, 1, 0.25],
        }

        path_segments = [
            ("Gamma", "X"),
            ("X", "U"),
            ("K", "Gamma"),
        ]

        for start_name, end_name in path_segments:
            point1 = np.array(points[start_name])
            point2 = np.array(points[end_name])

            result = helper.make_slice_desc(
                pathSQE_params, q_dims_and_bins, point1, point2, (start_name, end_name)
            )

            assert result["seg_start_name"] == start_name
            assert result["seg_end_name"] == end_name
            assert result["good_slice"] is True
