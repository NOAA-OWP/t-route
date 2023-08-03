from pathlib import Path
from pydantic import (
    FilePath as PydanticFilePath,
    DirectoryPath as PydanticDirectoryPath,
)

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pydantic.typing import CallableGenerator

from ._utils import strict_set


class FilePath(PydanticFilePath):
    """
    Coerce into a `pathlib.Path` type. If strict mode is enabled (see _utils.use_strict), will raise
    if file does not exist.
    """

    @classmethod
    def __get_validators__(cls) -> "CallableGenerator":
        yield cls.validate

    @classmethod
    def validate(cls, value: Path) -> Path:
        if strict_set():
            for validator in PydanticFilePath.__get_validators__():
                value = validator(value)
            return value
        else:
            return Path(value)


class DirectoryPath(PydanticDirectoryPath):
    """
    Coerce into a `pathlib.Path` type. If strict mode is enabled (see _utils.use_strict), will raise
    if directory does not exist.
    """

    @classmethod
    def __get_validators__(cls) -> "CallableGenerator":
        yield cls.validate

    @classmethod
    def validate(cls, value: Path) -> Path:
        if strict_set():
            for validator in PydanticDirectoryPath.__get_validators__():
                value = validator(value)
            return value
        else:
            return Path(value)
