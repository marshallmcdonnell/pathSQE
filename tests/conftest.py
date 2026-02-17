"""Pytest configuration and fixtures for pathSQE tests"""
import pytest
import numpy as np
import tempfile
import os
from pathlib import Path


@pytest.fixture
def temp_dir():
    """Fixture providing a temporary directory for test files"""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield tmpdir


@pytest.fixture
def sample_poscar_content():
    """Fixture providing a sample POSCAR file content"""
    return """Cubic Si
1.0
5.431 0.0 0.0
0.0 5.431 0.0
0.0 0.0 5.431
Si
2
Direct
0.0 0.0 0.0
0.25 0.25 0.25
"""


@pytest.fixture
def sample_poscar_file(temp_dir, sample_poscar_content):
    """Fixture providing a temporary POSCAR file"""
    poscar_path = os.path.join(temp_dir, "POSCAR")
    with open(poscar_path, "w") as f:
        f.write(sample_poscar_content)
    return poscar_path


@pytest.fixture
def sample_numpy_array():
    """Fixture providing a sample numpy array for testing"""
    return np.array([1.0, 2.0, 3.0, 4.5])


@pytest.fixture
def sample_list():
    """Fixture providing a sample list for testing"""
    return [1.0, 2.0, 3.0, 4.5]


@pytest.fixture
def sample_tuple():
    """Fixture providing a sample tuple for testing"""
    return (1.0, 2.0, 3.0, 4.5)
