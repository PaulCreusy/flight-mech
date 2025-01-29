import pytest
import subprocess
from pathlib import Path
from unittest.mock import patch
import matplotlib.pyplot as plt
import os

def pytest_collect_file(parent, file_path: Path):
    """Make pytest collect all Python files in examples/ as tests."""
    if file_path.suffix == ".py" and "conftest.py" not in file_path.name:
        return ExampleFile.from_parent(parent, path=file_path)

class ExampleFile(pytest.File):
    def collect(self):
        yield ExampleItem.from_parent(self, name=self.path.name)

class ExampleItem(pytest.Item):
    def runtest(self):
        """Run the script as a standalone subprocess."""
        result = subprocess.run(
            ["python", str(self.path)],
            capture_output=True,
            text=True,
            # Add the env variable to trigger patch in subprocess
            env={**os.environ, "PYTEST_RUNNING": "1"}
        )
        if result.returncode != 0:
            raise RuntimeError(
                f"Example {self.path} failed:\n{result.stdout}\n{result.stderr}")

    def reportinfo(self):
        return self.path, 0, f"running example: {self.path}"

@pytest.fixture(autouse=True)
def mock_show():
    """Automatically mock plt.show() in all tests to prevent blocking execution."""
    with patch.object(plt, "show"):
        yield
