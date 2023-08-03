import pytest

import yaml
from pathlib import Path
from typing import List

from troute.config import Config

TEST_DIR = Path(__file__).parent
ROOT_TEST_DIR = TEST_DIR / "../../../test"


def config_files() -> List[Path]:
    files = list(ROOT_TEST_DIR.glob("*/*.yaml"))
    return files


@pytest.mark.parametrize("file", config_files())
def test_naive_deserialization(file: Path):
    data = yaml.load(file.read_text(), Loader=yaml.Loader)
    Config(**data)
