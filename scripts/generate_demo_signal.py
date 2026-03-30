\
from __future__ import annotations

from pathlib import Path

import numpy as np


class DemoSignalGenerator:
    def __init__(self, samples: int, dt_seconds: float, random_seed: int) -> None:
        self._samples = samples
        self._dt_seconds = dt_seconds
        self._random_seed = random_seed

    def generate(self) -> np.ndarray:
        rng = np.random.default_rng(self._random_seed)
        time_axis = np.arange(self._samples, dtype=np.float64) * self._dt_seconds

        useful_signal = (
            3.5 * np.sin(2.0 * np.pi * 50.0 * time_axis)
            + 2.2 * np.sin(2.0 * np.pi * 120.0 * time_axis + 0.4)
            + 1.8 * np.sin(2.0 * np.pi * 400.0 * time_axis + 1.1)
        )
        noise = 0.8 * rng.normal(size=self._samples)
        return (useful_signal + noise).astype(np.float32)

    def save(self, output_path: Path) -> None:
        signal = self.generate()
        output_path.parent.mkdir(parents=True, exist_ok=True)
        signal.tofile(output_path)


def main() -> None:
    generator = DemoSignalGenerator(
        samples=300_000,
        dt_seconds=0.0001,
        random_seed=42,
    )
    generator.save(Path("examples/demo_signal_float32.bin"))
    print("Demo signal generated: examples/demo_signal_float32.bin")


if __name__ == "__main__":
    main()
