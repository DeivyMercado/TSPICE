import pytest
import tspice

def test_import():
    """Test that the package can be imported."""
    assert tspice is not None

def test_version():
    """Test that the package has a version."""
    assert hasattr(tspice, '__version__')
    assert isinstance(tspice.__version__, str)

def test_classes_exist():
    """Test that main classes exist."""
    assert hasattr(tspice, 'Body')
    assert hasattr(tspice, 'BodyResponse')

def test_functions_exist():
    """Test that main functions exist."""
    assert hasattr(tspice, 'initialize')
