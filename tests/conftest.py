"""
Pytest configuration and shared fixtures for nf-scrnaseq tests
"""

import pytest
from pathlib import Path


@pytest.fixture(scope="session")
def project_root():
    """Get project root directory"""
    return Path(__file__).parent.parent


@pytest.fixture(scope="session")
def modules_dir(project_root):
    """Get modules directory"""
    return project_root / 'modules' / 'local'


@pytest.fixture(scope="session")
def test_data_dir(project_root):
    """Get test data directory"""
    return project_root / 'test_data' / '10x_sample'


def pytest_configure(config):
    """Configure pytest with custom markers"""
    config.addinivalue_line(
        "markers", "slow: marks tests as slow (deselect with '-m \"not slow\"')"
    )
    config.addinivalue_line(
        "markers", "integration: marks tests as integration tests"
    )
    config.addinivalue_line(
        "markers", "unit: marks tests as unit tests"
    )
