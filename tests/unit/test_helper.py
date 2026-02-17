"""Unit tests for pathSQE.helper module"""
import pytest
import numpy as np
import os
import tempfile
from pathlib import Path

from pathSQE import helper


class TestSimpleReadPoscar:
    """Tests for simple_read_poscar function"""

    def test_read_poscar_basic(self, sample_poscar_file):
        """Test reading a basic POSCAR file"""
        cell, positions, atomic_numbers = helper.simple_read_poscar(sample_poscar_file)

        # Check cell structure
        assert len(cell) == 3
        assert len(cell[0]) == 3
        np.testing.assert_array_almost_equal(
            cell[0], [5.431, 0.0, 0.0]
        )
        np.testing.assert_array_almost_equal(
            cell[1], [0.0, 5.431, 0.0]
        )
        np.testing.assert_array_almost_equal(
            cell[2], [0.0, 0.0, 5.431]
        )

        # Check positions
        assert len(positions) == 2
        np.testing.assert_array_almost_equal(positions[0], [0.0, 0.0, 0.0])
        np.testing.assert_array_almost_equal(positions[1], [0.25, 0.25, 0.25])

        # Check atomic numbers (Si has atomic number 14)
        assert atomic_numbers == [14, 14]

    def test_read_poscar_multiple_species(self, temp_dir):
        """Test reading a POSCAR with multiple species"""
        poscar_content = """Mixed Si and Ge
1.0
5.431 0.0 0.0
0.0 5.431 0.0
0.0 0.0 5.431
Si Ge
1 1
Direct
0.0 0.0 0.0
0.25 0.25 0.25
"""
        poscar_path = os.path.join(temp_dir, "POSCAR_mixed")
        with open(poscar_path, "w") as f:
            f.write(poscar_content)

        cell, positions, atomic_numbers = helper.simple_read_poscar(poscar_path)

        assert len(positions) == 2
        assert atomic_numbers[0] == 14  # Si
        assert atomic_numbers[1] == 32  # Ge

    def test_read_poscar_scaling(self, temp_dir):
        """Test that lattice parameter scaling is correctly applied"""
        poscar_content = """Scaled lattice
2.0
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
Si
1
Direct
0.5 0.5 0.5
"""
        poscar_path = os.path.join(temp_dir, "POSCAR_scaled")
        with open(poscar_path, "w") as f:
            f.write(poscar_content)

        cell, positions, atomic_numbers = helper.simple_read_poscar(poscar_path)

        # Lattice vectors should be scaled by alat=2.0
        np.testing.assert_array_almost_equal(cell[0], [2.0, 0.0, 0.0])
        np.testing.assert_array_almost_equal(cell[1], [0.0, 2.0, 0.0])
        np.testing.assert_array_almost_equal(cell[2], [0.0, 0.0, 2.0])

    def test_read_poscar_with_comments(self, temp_dir):
        """Test that POSCAR files with inline comments are handled correctly"""
        poscar_content = """Test with comments
1.0
5.431 0.0 0.0  ! x lattice vector
0.0 5.431 0.0  ! y lattice vector
0.0 0.0 5.431  ! z lattice vector
Si  ! Element
2   ! Number of atoms
Direct
0.0 0.0 0.0    ! Position 1
0.25 0.25 0.25 ! Position 2
"""
        poscar_path = os.path.join(temp_dir, "POSCAR_comments")
        with open(poscar_path, "w") as f:
            f.write(poscar_content)

        cell, positions, atomic_numbers = helper.simple_read_poscar(poscar_path)

        assert len(positions) == 2
        assert len(cell) == 3
        assert atomic_numbers == [14, 14]

    def test_read_poscar_nonexistent_file(self):
        """Test that reading a nonexistent file raises an error"""
        with pytest.raises(FileNotFoundError):
            helper.simple_read_poscar("/nonexistent/path/POSCAR")

    def test_read_poscar_invalid_coordinates(self, temp_dir):
        """Test that non-Direct coordinate format raises an error"""
        poscar_content = """Invalid format
1.0
5.431 0.0 0.0
0.0 5.431 0.0
0.0 0.0 5.431
Si
2
Cartesian
0.0 0.0 0.0
1.357 1.357 1.357
"""
        poscar_path = os.path.join(temp_dir, "POSCAR_cartesian")
        with open(poscar_path, "w") as f:
            f.write(poscar_content)

        with pytest.raises(ValueError, match="can only deal with 'direct' POSCARs"):
            helper.simple_read_poscar(poscar_path)


class TestConvToDescString:
    """Tests for conv_to_desc_string function"""

    def test_numpy_array(self, sample_numpy_array):
        """Test conversion of numpy array"""
        result = helper.conv_to_desc_string(sample_numpy_array)
        assert result == "1.0,2.0,3.0,4.5"

    def test_list(self, sample_list):
        """Test conversion of list"""
        result = helper.conv_to_desc_string(sample_list)
        assert result == "1.0,2.0,3.0,4.5"

    def test_tuple(self, sample_tuple):
        """Test conversion of tuple"""
        result = helper.conv_to_desc_string(sample_tuple)
        assert result == "1.0,2.0,3.0,4.5"

    def test_integers(self):
        """Test conversion with integers"""
        result = helper.conv_to_desc_string([1, 2, 3])
        assert result == "1,2,3"

    def test_mixed_types(self):
        """Test conversion with mixed int and float"""
        result = helper.conv_to_desc_string([1, 2.5, 3])
        assert result == "1,2.5,3"

    def test_single_element(self):
        """Test conversion with single element"""
        result = helper.conv_to_desc_string([5.0])
        assert result == "5.0"

    def test_negative_values(self):
        """Test conversion with negative values"""
        result = helper.conv_to_desc_string([-1.0, -2.5, 3.0])
        assert result == "-1.0,-2.5,3.0"

    def test_large_array(self):
        """Test conversion with large array"""
        large_list = list(range(100))
        result = helper.conv_to_desc_string(large_list)
        expected = ",".join(str(i) for i in range(100))
        assert result == expected

    def test_zero_values(self):
        """Test conversion with zeros"""
        result = helper.conv_to_desc_string([0, 0.0, 0])
        assert result == "0,0.0,0"


class TestMakeSliceDesc:
    """Tests for make_slice_desc function"""

    def test_make_slice_desc_basic(self):
        """Test basic slice description creation"""
        pathSQE_params = {
            "E bins": "-10,0.1,10"
        }
        q_dims_and_bins = [
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1],
            [-0.5, 0.5, 0.5],
            [-0.5, 0.5, 0.5],
            [-0.5, 0.5, 0.5],
        ]
        point1 = np.array([0.0, 0.0, 0.0])
        point2 = np.array([1.0, 0.0, 0.0])
        path_seg = ("Gamma", "X")

        result = helper.make_slice_desc(
            pathSQE_params, q_dims_and_bins, point1, point2, path_seg
        )

        assert "Name" in result
        assert "QDimension0" in result
        assert "QDimension1" in result
        assert "QDimension2" in result
        assert result["seg_start_name"] == "Gamma"
        assert result["seg_end_name"] == "X"
        assert result["good_slice"] is True

    def test_make_slice_desc_inv_angstrom_ratio(self):
        """Test that inv_angstrom_ratio is computed correctly"""
        pathSQE_params = {"E bins": "-10,0.1,10"}
        q_dims_and_bins = [
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1],
            [-0.5, 0.5, 0.5],
            [-0.5, 0.5, 0.5],
            [-0.5, 0.5, 0.5],
        ]
        point1 = np.array([0.0, 0.0, 0.0])
        point2 = np.array([0.5, 0.0, 0.0])
        path_seg = ("Gamma", "X")

        result = helper.make_slice_desc(
            pathSQE_params, q_dims_and_bins, point1, point2, path_seg
        )

        # Distance between points is 0.5, so ratio should be 0.5/0.5 = 1.0
        assert result["inv_angstrom_ratio"] == 1.0

    def test_make_slice_desc_different_segments(self):
        """Test slice description with different path segments"""
        pathSQE_params = {"E bins": "-10,0.1,10"}
        q_dims_and_bins = [
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1],
            [-0.5, 0.5, 0.5],
            [-0.5, 0.5, 0.5],
            [-0.5, 0.5, 0.5],
        ]
        point1 = np.array([0.5, 0.5, 0.0])
        point2 = np.array([0.0, 0.0, 0.0])
        path_seg = ("M", "Gamma")

        result = helper.make_slice_desc(
            pathSQE_params, q_dims_and_bins, point1, point2, path_seg
        )

        assert result["seg_start_name"] == "M"
        assert result["seg_end_name"] == "Gamma"
        assert result["good_slice"] is True
