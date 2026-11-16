"""Records shared by documentation parsing and repository checks."""

from dataclasses import dataclass


@dataclass(frozen=True, order=True)
class Issue:
    path: str
    line: int
    message: str

    def __str__(self) -> str:
        return f"{self.path}:{self.line}: {self.message}"


@dataclass(frozen=True)
class VisibleLine:
    number: int
    raw: str
    masked: str
    spans: tuple[str, ...]


@dataclass(frozen=True)
class Fence:
    tag: str
    start: int
    body: tuple[tuple[int, str], ...]
