\
from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, List

import matplotlib.pyplot as plt
import numpy as np


@dataclass(frozen=True)
class SignalConfig:
    samples: int
    dt_seconds: float
    dtype: str
    endianness: str
    min_peak_ratio: float
    neighbor_bins: int
    top_k: int
    use_hann_window: bool


@dataclass(frozen=True)
class PeakInfo:
    bin_index: int
    frequency_hz: float
    amplitude: float


@dataclass(frozen=True)
class AnalysisSummary:
    detected_frequencies_hz: List[float]
    noise_floor: float
    max_amplitude: float
    kept_bins: List[int]


class BinarySignalLoader:
    _DTYPE_SIZES = {
        "float32": 4,
        "float64": 8,
        "int16": 2,
        "int32": 4,
        "int64": 8,
    }

    _NUMPY_DTYPES = {
        "float32": np.float32,
        "float64": np.float64,
        "int16": np.int16,
        "int32": np.int32,
        "int64": np.int64,
    }

    def __init__(self, expected_samples: int) -> None:
        self._expected_samples = expected_samples

    def infer_dtype(self, file_path: Path) -> str:
        if file_path.suffix.lower() == ".npy":
            array = np.load(file_path, allow_pickle=False)
            return str(array.dtype)

        file_size = file_path.stat().st_size
        matches = [
            dtype_name
            for dtype_name, item_size in self._DTYPE_SIZES.items()
            if file_size == self._expected_samples * item_size
        ]

        if len(matches) == 1:
            return matches[0]

        supported = ", ".join(sorted(self._DTYPE_SIZES))
        raise ValueError(
            "Unable to infer dtype from file size. "
            f"File size: {file_size} bytes, expected samples: {self._expected_samples}, "
            f"matched dtypes: {matches or 'none'}. "
            f"Use --dtype explicitly. Supported dtypes: {supported}."
        )

    def load(self, file_path: Path, dtype_name: str, endianness: str) -> np.ndarray:
        if file_path.suffix.lower() == ".npy":
            signal = np.load(file_path, allow_pickle=False)

            signal = np.asarray(signal)

            if signal.ndim != 1:
                signal = signal.reshape(-1)

            if signal.size != self._expected_samples:
                raise ValueError(
                    f"Expected {self._expected_samples} samples, but got {signal.size} "
                    f"from NPY file {file_path.name}."
                )

            return signal.astype(np.float64, copy=False)

        resolved_dtype_name = self.infer_dtype(file_path) if dtype_name == "auto" else dtype_name
        if resolved_dtype_name not in self._NUMPY_DTYPES:
            supported = ", ".join(sorted(self._NUMPY_DTYPES))
            raise ValueError(f"Unsupported dtype: {resolved_dtype_name}. Supported: {supported}.")

        numpy_dtype = np.dtype(self._NUMPY_DTYPES[resolved_dtype_name])
        if endianness == "little":
            numpy_dtype = numpy_dtype.newbyteorder("<")
        elif endianness == "big":
            numpy_dtype = numpy_dtype.newbyteorder(">")
        else:
            raise ValueError("Endianness must be either 'little' or 'big'.")

        signal = np.fromfile(file_path, dtype=numpy_dtype)
        if signal.size != self._expected_samples:
            raise ValueError(
                f"Expected {self._expected_samples} samples, but got {signal.size}. "
                "Check the dtype or the file integrity."
            )

        return signal.astype(np.float64, copy=False)
class FourierSignalAnalyzer:
    def __init__(self, config: SignalConfig) -> None:
        self._config = config

    def analyze(self, signal: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[PeakInfo], AnalysisSummary]:
        centered_signal = signal - np.mean(signal)
        window = self._build_window(centered_signal.size)
        windowed_signal = centered_signal * window

        spectrum = np.fft.rfft(windowed_signal)
        frequencies = np.fft.rfftfreq(centered_signal.size, d=self._config.dt_seconds)

        amplitude = self._compute_amplitude(windowed_signal.size, spectrum, window)
        noise_floor = self._estimate_noise_floor(amplitude)
        peaks = self._detect_peaks(frequencies, amplitude, noise_floor)
        filtered_spectrum = self._filter_spectrum(spectrum, peaks)
        filtered_signal = np.fft.irfft(filtered_spectrum, n=centered_signal.size) + np.mean(signal)

        summary = AnalysisSummary(
            detected_frequencies_hz=[peak.frequency_hz for peak in peaks],
            noise_floor=float(noise_floor),
            max_amplitude=float(np.max(amplitude)),
            kept_bins=sorted(self._build_kept_bin_set(peaks)),
        )
        return frequencies, amplitude, filtered_signal, peaks, summary

    def _build_window(self, size: int) -> np.ndarray:
        if self._config.use_hann_window:
            return np.hanning(size)
        return np.ones(size, dtype=np.float64)

    def _compute_amplitude(
        self,
        sample_count: int,
        spectrum: np.ndarray,
        window: np.ndarray,
    ) -> np.ndarray:
        coherent_gain = float(np.mean(window))
        amplitude = 2.0 * np.abs(spectrum) / (sample_count * coherent_gain)
        amplitude[0] = np.abs(spectrum[0]) / (sample_count * coherent_gain)
        if sample_count % 2 == 0:
            amplitude[-1] = np.abs(spectrum[-1]) / (sample_count * coherent_gain)
        return amplitude

    def _estimate_noise_floor(self, amplitude: np.ndarray) -> float:
        usable = amplitude[1:]
        if usable.size == 0:
            return 0.0
        return float(np.median(usable))

    def _detect_peaks(
        self,
        frequencies: np.ndarray,
        amplitude: np.ndarray,
        noise_floor: float,
    ) -> list[PeakInfo]:
        if amplitude.size < 3:
            return []

        minimum_threshold = max(noise_floor * self._config.min_peak_ratio, 1e-12)
        local_maxima: list[PeakInfo] = []

        for index in range(1, amplitude.size - 1):
            is_local_maximum = amplitude[index] >= amplitude[index - 1] and amplitude[index] >= amplitude[index + 1]
            is_positive_frequency = frequencies[index] > 0.0
            if not (is_local_maximum and is_positive_frequency):
                continue

            refined_frequency = self._quadratic_peak_frequency(frequencies, amplitude, index)
            local_maxima.append(
                PeakInfo(
                    bin_index=index,
                    frequency_hz=float(refined_frequency),
                    amplitude=float(amplitude[index]),
                )
            )

        if not local_maxima:
            return []

        local_maxima.sort(key=lambda item: item.amplitude, reverse=True)

        split_index = self._find_signal_noise_split(local_maxima)
        selected = local_maxima[:split_index]

        if not selected:
            selected = [peak for peak in local_maxima if peak.amplitude >= minimum_threshold]

        selected = [peak for peak in selected if peak.amplitude >= minimum_threshold]
        return selected[: self._config.top_k]

    def _find_signal_noise_split(self, local_maxima: list[PeakInfo]) -> int:
        if len(local_maxima) == 1:
            return 1

        for index in range(len(local_maxima) - 1):
            current_amplitude = local_maxima[index].amplitude
            next_amplitude = local_maxima[index + 1].amplitude
            if next_amplitude <= 0.0:
                return index + 1
            if current_amplitude >= self._config.min_peak_ratio * next_amplitude:
                return index + 1

        return 0

    def _quadratic_peak_frequency(
        self,
        frequencies: np.ndarray,
        amplitude: np.ndarray,
        index: int,
    ) -> float:
        left = amplitude[index - 1]
        center = amplitude[index]
        right = amplitude[index + 1]
        denominator = left - 2.0 * center + right

        if abs(denominator) < 1e-18:
            return float(frequencies[index])

        delta = 0.5 * (left - right) / denominator
        bin_width = frequencies[1] - frequencies[0]
        return float(frequencies[index] + delta * bin_width)

    def _build_kept_bin_set(self, peaks: Iterable[PeakInfo]) -> set[int]:
        kept_bins: set[int] = {0}
        for peak in peaks:
            start = max(0, peak.bin_index - self._config.neighbor_bins)
            finish = peak.bin_index + self._config.neighbor_bins + 1
            kept_bins.update(range(start, finish))
        return kept_bins

    def _filter_spectrum(self, spectrum: np.ndarray, peaks: list[PeakInfo]) -> np.ndarray:
        kept_bins = self._build_kept_bin_set(peaks)
        filtered = np.zeros_like(spectrum)
        for bin_index in kept_bins:
            if 0 <= bin_index < spectrum.size:
                filtered[bin_index] = spectrum[bin_index]
        return filtered


class SignalPlotter:
    def __init__(self, output_dir: Path) -> None:
        self._output_dir = output_dir
        self._output_dir.mkdir(parents=True, exist_ok=True)

    def save_signal_plot(self, time_ms: np.ndarray, raw_signal: np.ndarray, filtered_signal: np.ndarray) -> Path:
        figure = plt.figure(figsize=(14, 6))
        axes = figure.add_subplot(111)
        axes.plot(time_ms, raw_signal, linewidth=0.8, alpha=0.7, label="Raw signal")
        axes.plot(time_ms, filtered_signal, linewidth=1.2, label="Filtered signal")
        axes.set_title("Signal in time domain")
        axes.set_xlabel("Time, ms")
        axes.set_ylabel("Amplitude")
        axes.grid(True)
        axes.legend()

        output_path = self._output_dir / "signal_time_domain.png"
        figure.tight_layout()
        figure.savefig(output_path, dpi=160)
        plt.close(figure)
        return output_path

    def save_spectrum_plot(self, frequencies: np.ndarray, amplitude: np.ndarray, peaks: list[PeakInfo]) -> Path:
        figure = plt.figure(figsize=(14, 6))
        axes = figure.add_subplot(111)
        axes.plot(frequencies, amplitude, linewidth=0.8)
        for peak in peaks:
            axes.axvline(peak.frequency_hz, linestyle="--", linewidth=1.0)
        axes.set_title("Amplitude spectrum")
        axes.set_xlabel("Frequency, Hz")
        axes.set_ylabel("Amplitude")
        axes.grid(True)

        output_path = self._output_dir / "spectrum.png"
        figure.tight_layout()
        figure.savefig(output_path, dpi=160)
        plt.close(figure)
        return output_path

    def save_filtered_only_plot(self, time_ms: np.ndarray, filtered_signal: np.ndarray) -> Path:
        figure = plt.figure(figsize=(14, 6))
        axes = figure.add_subplot(111)
        axes.plot(time_ms, filtered_signal, linewidth=1.0)
        axes.set_title("Filtered signal only")
        axes.set_xlabel("Time, ms")
        axes.set_ylabel("Amplitude")
        axes.grid(True)

        output_path = self._output_dir / "filtered_signal.png"
        figure.tight_layout()
        figure.savefig(output_path, dpi=160)
        plt.close(figure)
        return output_path


class FourierFilteringApp:
    def __init__(self, config: SignalConfig) -> None:
        self._config = config
        self._loader = BinarySignalLoader(expected_samples=config.samples)
        self._analyzer = FourierSignalAnalyzer(config=config)

    def run(self, input_path: Path, output_dir: Path) -> AnalysisSummary:
        signal = self._loader.load(
            file_path=input_path,
            dtype_name=self._config.dtype,
            endianness=self._config.endianness,
        )

        frequencies, amplitude, filtered_signal, peaks, summary = self._analyzer.analyze(signal=signal)

        time_ms = np.arange(signal.size, dtype=np.float64) * self._config.dt_seconds * 1000.0
        plotter = SignalPlotter(output_dir=output_dir)
        plotter.save_signal_plot(time_ms=time_ms, raw_signal=signal, filtered_signal=filtered_signal)
        plotter.save_spectrum_plot(frequencies=frequencies, amplitude=amplitude, peaks=peaks)
        plotter.save_filtered_only_plot(time_ms=time_ms, filtered_signal=filtered_signal)

        self._save_numeric_outputs(
            output_dir=output_dir,
            time_ms=time_ms,
            filtered_signal=filtered_signal,
            frequencies=frequencies,
            amplitude=amplitude,
            peaks=peaks,
            summary=summary,
        )

        return summary

    def _save_numeric_outputs(
        self,
        output_dir: Path,
        time_ms: np.ndarray,
        filtered_signal: np.ndarray,
        frequencies: np.ndarray,
        amplitude: np.ndarray,
        peaks: list[PeakInfo],
        summary: AnalysisSummary,
    ) -> None:
        output_dir.mkdir(parents=True, exist_ok=True)

        filtered_data = np.column_stack((time_ms, filtered_signal))
        np.savetxt(
            output_dir / "filtered_signal.csv",
            filtered_data,
            delimiter=",",
            header="time_ms,filtered_signal",
            comments="",
        )

        spectrum_data = np.column_stack((frequencies, amplitude))
        np.savetxt(
            output_dir / "spectrum.csv",
            spectrum_data,
            delimiter=",",
            header="frequency_hz,amplitude",
            comments="",
        )

        with (output_dir / "detected_frequencies.txt").open("w", encoding="utf-8") as stream:
            if not peaks:
                stream.write("No peaks were detected.\n")
            else:
                stream.write("Detected useful frequencies (Hz):\n")
                for peak in peaks:
                    stream.write(
                        f"- {peak.frequency_hz:.6f} Hz (bin={peak.bin_index}, amplitude={peak.amplitude:.6f})\n"
                    )

        with (output_dir / "summary.json").open("w", encoding="utf-8") as stream:
            json.dump(asdict(summary), stream, indent=2)


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Solve the Fourier filtering task for a binary signal file.")
    parser.add_argument("--input", required=True, help="Path to the binary file with the noisy signal.")
    parser.add_argument("--output", default="output", help="Directory where all generated files will be stored.")
    parser.add_argument("--samples", type=int, default=100_000, help="Expected number of samples.")
    parser.add_argument(
        "--dt-ms",
        type=float,
        default=0.1,
        help="Sampling interval in milliseconds.",
    )
    parser.add_argument(
        "--dtype",
        default="auto",
        choices=["auto", "float32", "float64", "int16", "int32", "int64"],
        help="Input binary data type. Use 'auto' to infer from file size.",
    )
    parser.add_argument(
        "--endianness",
        default="little",
        choices=["little", "big"],
        help="Byte order of the binary file.",
    )
    parser.add_argument(
        "--min-peak-ratio",
        type=float,
        default=2.0,
        help="Minimum ratio between a peak amplitude and the estimated noise floor.",
    )
    parser.add_argument(
        "--neighbor-bins",
        type=int,
        default=2,
        help="How many bins around each detected peak should be preserved during inverse FFT.",
    )
    parser.add_argument(
        "--top-k",
        type=int,
        default=20,
        help="Maximum number of strongest frequency peaks to keep in the report.",
    )
    parser.add_argument(
        "--no-hann-window",
        action="store_true",
        help="Disable Hann window before FFT.",
    )
    return parser


def main() -> None:
    parser = build_argument_parser()
    args = parser.parse_args()

    config = SignalConfig(
        samples=args.samples,
        dt_seconds=args.dt_ms / 1000.0,
        dtype=args.dtype,
        endianness=args.endianness,
        min_peak_ratio=args.min_peak_ratio,
        neighbor_bins=args.neighbor_bins,
        top_k=args.top_k,
        use_hann_window=not args.no_hann_window,
    )
    app = FourierFilteringApp(config=config)
    summary = app.run(input_path=Path(args.input), output_dir=Path(args.output))

    print("Detected useful frequencies (Hz):")
    if not summary.detected_frequencies_hz:
        print("  No peaks were detected.")
    else:
        for frequency in summary.detected_frequencies_hz:
            print(f"  {frequency:.6f}")
    print(f"Noise floor: {summary.noise_floor:.6f}")
    print(f"Max amplitude: {summary.max_amplitude:.6f}")


if __name__ == "__main__":
    main()
