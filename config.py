from __future__ import annotations

import json
from pathlib import Path
from typing import Any


def load_config(path: str | Path) -> dict[str, Any]:
    """Load a JSON configuration file."""
    with Path(path).open("r", encoding="utf-8") as fin:
        return json.load(fin)


def read_range(spec: dict[str, float | int] | list[float]) -> list[float]:
    """Read either a literal list or a {min,max,num} scan spec."""
    if isinstance(spec, list):
        if not spec:
            raise ValueError("range list must not be empty")
        return [float(x) for x in spec]

    xmin = float(spec["min"])
    xmax = float(spec["max"])
    num = int(spec["num"])
    if num < 1:
        raise ValueError("range num must be >= 1")
    if num == 1:
        return [xmin]
    return [xmin + (xmax - xmin) * i / (num - 1) for i in range(num)]


def section(config: dict[str, Any], name: str) -> dict[str, Any]:
    if name not in config or not isinstance(config[name], dict):
        raise KeyError(f'config missing object "{name}"')
    return config[name]
