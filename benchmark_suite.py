#!/usr/bin/env python3
"""
Automated benchmarking and regression detection for HOTPANTS.
Tracks performance across multiple workloads and detects regressions.
"""
import json
import subprocess
import tempfile
import time
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path
from typing import Optional

import numpy as np
from astropy.io import fits
from scipy.ndimage import gaussian_filter


@dataclass
class BenchmarkWorkload:
    """Configuration for a single benchmark workload."""
    name: str
    nx: int
    ny: int
    n_stars: int
    nrx: int
    nry: int
    use_fft: bool = True
    description: str = ""


@dataclass
class BenchmarkResult:
    """Result of a single benchmark run."""
    workload_name: str
    timestamp: str
    elapsed_seconds: float
    image_size: int  # pixels
    num_regions: int
    use_fft: bool
    git_sha: Optional[str] = None
    git_branch: Optional[str] = None

    def to_dict(self):
        return asdict(self)


class HotpantsBenchmark:
    """Orchestrates HOTPANTS benchmarking and tracks results."""

    def __init__(self, hotpants_bin: Path, results_dir: Path = None):
        self.hotpants_bin = hotpants_bin
        self.results_dir = results_dir or Path("/tmp/hotpants_benchmarks")
        self.results_dir.mkdir(parents=True, exist_ok=True)
        self.results_file = self.results_dir / "benchmark_results.jsonl"

        # Standard workloads for regression testing
        self.standard_workloads = [
            BenchmarkWorkload(
                name="small_single",
                nx=256, ny=256, n_stars=50,
                nrx=1, nry=1,
                description="Small image, single region"
            ),
            BenchmarkWorkload(
                name="medium_single",
                nx=512, ny=512, n_stars=120,
                nrx=1, nry=1,
                description="Medium image, single region (baseline)"
            ),
            BenchmarkWorkload(
                name="large_single",
                nx=1024, ny=1024, n_stars=250,
                nrx=1, nry=1,
                description="Large image, single region"
            ),
            BenchmarkWorkload(
                name="medium_tiled",
                nx=512, ny=512, n_stars=120,
                nrx=2, nry=2,
                description="Medium image, 4 regions"
            ),
            BenchmarkWorkload(
                name="large_tiled",
                nx=1024, ny=1024, n_stars=250,
                nrx=4, nry=4,
                description="Large image, 16 regions"
            ),
        ]

    def make_synthetic_images(self, tmpdir: Path, nx: int, ny: int, n_stars: int):
        """Generate synthetic template and science images."""
        rng = np.random.default_rng(0xC0FFEE)
        margin = 40
        x_stars = rng.uniform(margin, nx - margin, n_stars)
        y_stars = rng.uniform(margin, ny - margin, n_stars)
        fluxes = rng.uniform(8_000, 45_000, n_stars)

        noiseless = np.zeros((ny, nx), dtype=np.float64)
        for x, y, f in zip(x_stars, y_stars, fluxes):
            xi, yi = int(round(x)), int(round(y))
            noiseless[yi, xi] += f

        bg = 1000.0
        gain = 1.0
        rdnoise = 5.0
        psf_template = 1.5
        psf_science = np.sqrt(psf_template**2 + 1.5**2)

        # Template
        template_conv = gaussian_filter(noiseless, sigma=psf_template)
        template_conv += bg
        template_noise = np.sqrt(np.maximum(template_conv, 0) / gain + rdnoise**2) * rng.standard_normal((ny, nx))
        template = (template_conv + template_noise).astype(np.float32)

        # Science
        science_conv = gaussian_filter(noiseless * 1.05, sigma=psf_science)
        science_conv += bg * 0.98
        science_noise = np.sqrt(np.maximum(science_conv, 0) / gain + rdnoise**2) * rng.standard_normal((ny, nx))
        science = (science_conv + science_noise).astype(np.float32)

        template_path = tmpdir / "template.fits"
        science_path = tmpdir / "science.fits"
        diff_path = tmpdir / "diff.fits"

        fits.PrimaryHDU(template).writeto(str(template_path), overwrite=True)
        fits.PrimaryHDU(science).writeto(str(science_path), overwrite=True)

        return template_path, science_path, diff_path

    def run_workload(self, workload: BenchmarkWorkload, num_trials: int = 1) -> Optional[BenchmarkResult]:
        """Run a benchmark workload and return result."""
        print(f"\n{'='*70}")
        print(f"Benchmark: {workload.name}")
        print(f"  Image size: {workload.nx}×{workload.ny} ({workload.nx*workload.ny:,} pixels)")
        print(f"  Stars: {workload.n_stars}")
        print(f"  Regions: {workload.nrx}×{workload.nry}")
        print(f"  FFT: {'enabled' if workload.use_fft else 'disabled'}")
        print(f"{'='*70}")

        with tempfile.TemporaryDirectory(prefix="hotpants_bench_") as tmpdir:
            tmpdir = Path(tmpdir)
            template, science, diff = self.make_synthetic_images(tmpdir, workload.nx, workload.ny, workload.n_stars)

            times = []
            for trial in range(num_trials):
                cmd = [
                    str(self.hotpants_bin),
                    "-inim", str(science),
                    "-tmplim", str(template),
                    "-outim", str(diff),
                    "-tu", "60000", "-tl", "-200",
                    "-iu", "60000", "-il", "-200",
                    "-tg", "1.0", "-ig", "1.0",
                    "-tr", "5.0", "-ir", "5.0",
                    "-r", "8",
                    "-nsx", "5", "-nsy", "5",
                    "-nss", "3",
                    "-ft", "8",
                    "-nrx", str(workload.nrx),
                    "-nry", str(workload.nry),
                    "-v", "0",
                ]

                if workload.use_fft:
                    cmd.append("-fft")

                start = time.perf_counter()
                result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
                elapsed = time.perf_counter() - start

                if result.returncode != 0:
                    print(f"  ERROR: {result.stderr}")
                    return None

                times.append(elapsed)
                print(f"  Trial {trial + 1}/{num_trials}: {elapsed:.2f}s")

            avg_time = np.mean(times)
            std_time = np.std(times) if len(times) > 1 else 0.0

            print(f"  Average: {avg_time:.2f}s (±{std_time:.2f}s)")

            result = BenchmarkResult(
                workload_name=workload.name,
                timestamp=datetime.now().isoformat(),
                elapsed_seconds=avg_time,
                image_size=workload.nx * workload.ny,
                num_regions=workload.nrx * workload.nry,
                use_fft=workload.use_fft,
                git_sha=self._get_git_sha(),
                git_branch=self._get_git_branch(),
            )

            return result

    def _get_git_sha(self) -> Optional[str]:
        """Get current git commit SHA."""
        try:
            result = subprocess.run(
                ["git", "rev-parse", "HEAD"],
                capture_output=True,
                text=True,
                cwd=self.hotpants_bin.parent.parent,
                timeout=5
            )
            return result.stdout.strip() if result.returncode == 0 else None
        except Exception:
            return None

    def _get_git_branch(self) -> Optional[str]:
        """Get current git branch."""
        try:
            result = subprocess.run(
                ["git", "rev-parse", "--abbrev-ref", "HEAD"],
                capture_output=True,
                text=True,
                cwd=self.hotpants_bin.parent.parent,
                timeout=5
            )
            return result.stdout.strip() if result.returncode == 0 else None
        except Exception:
            return None

    def save_result(self, result: BenchmarkResult):
        """Append result to benchmark log."""
        with open(self.results_file, "a") as f:
            f.write(json.dumps(result.to_dict()) + "\n")

    def load_results(self) -> list[BenchmarkResult]:
        """Load all recorded benchmark results."""
        results = []
        if not self.results_file.exists():
            return results

        with open(self.results_file) as f:
            for line in f:
                try:
                    data = json.loads(line)
                    results.append(BenchmarkResult(**data))
                except json.JSONDecodeError:
                    pass
        return results

    def detect_regressions(self, threshold_pct: float = 10.0) -> dict:
        """
        Compare current results against historical average.
        Returns dict of workload_name -> regression info.
        """
        results = self.load_results()
        if not results:
            return {}

        # Group by workload and branch
        by_workload = {}
        for result in results:
            key = (result.workload_name, result.git_branch or "unknown")
            if key not in by_workload:
                by_workload[key] = []
            by_workload[key].append(result)

        regressions = {}
        for (workload_name, branch), workload_results in by_workload.items():
            if len(workload_results) < 2:
                continue

            # Exclude latest run when computing baseline
            baseline = workload_results[:-1]
            latest = workload_results[-1]

            baseline_avg = np.mean([r.elapsed_seconds for r in baseline])
            regression_pct = ((latest.elapsed_seconds - baseline_avg) / baseline_avg) * 100

            if abs(regression_pct) > threshold_pct:
                regressions[workload_name] = {
                    "branch": branch,
                    "baseline_avg": baseline_avg,
                    "latest": latest.elapsed_seconds,
                    "regression_pct": regression_pct,
                    "num_baseline_runs": len(baseline),
                }

        return regressions

    def generate_report(self, output_file: Optional[Path] = None) -> str:
        """Generate a comprehensive performance report."""
        results = self.load_results()
        if not results:
            return "No benchmark results found."

        # Group by workload
        by_workload = {}
        for result in results:
            if result.workload_name not in by_workload:
                by_workload[result.workload_name] = []
            by_workload[result.workload_name].append(result)

        report = []
        report.append("# HOTPANTS Performance Benchmark Report\n")
        report.append(f"Generated: {datetime.now().isoformat()}\n\n")

        # Regression detection
        regressions = self.detect_regressions()
        if regressions:
            report.append("## ⚠️ Performance Regressions Detected\n")
            for workload_name, info in regressions.items():
                pct = info["regression_pct"]
                trend = "**SLOWER**" if pct > 0 else "**FASTER**"
                report.append(f"\n### {workload_name} [{trend}]\n")
                report.append(f"- Baseline avg: {info['baseline_avg']:.2f}s ({info['num_baseline_runs']} runs)\n")
                report.append(f"- Latest: {info['latest']:.2f}s\n")
                report.append(f"- Change: {pct:+.1f}%\n")
        else:
            report.append("## ✓ No regressions detected\n\n")

        # Per-workload summary
        report.append("## Workload Performance Summary\n\n")
        report.append("| Workload | Latest (s) | Avg (s) | Std | Runs | Status |\n")
        report.append("|----------|------------|--------|-----|------|--------|\n")

        for workload_name in sorted(by_workload.keys()):
            results_for_workload = by_workload[workload_name]
            times = [r.elapsed_seconds for r in results_for_workload]
            latest = results_for_workload[-1].elapsed_seconds
            avg = np.mean(times)
            std = np.std(times)
            n = len(times)

            status = "🔴 REGRESSION" if workload_name in regressions else "🟢 OK"
            report.append(f"| {workload_name} | {latest:.2f} | {avg:.2f} | {std:.2f} | {n} | {status} |\n")

        # Detailed history
        report.append("\n## Detailed History\n\n")
        for workload_name in sorted(by_workload.keys()):
            report.append(f"### {workload_name}\n\n")
            report.append("| Timestamp | Time (s) | Branch | SHA |\n")
            report.append("|-----------|----------|--------|-----|\n")
            for result in sorted(by_workload[workload_name], key=lambda r: r.timestamp):
                timestamp = result.timestamp[:19]  # YYYY-MM-DDTHH:MM:SS
                branch = result.git_branch or "unknown"
                sha = (result.git_sha or "unknown")[:8]
                report.append(f"| {timestamp} | {result.elapsed_seconds:.2f} | {branch} | {sha} |\n")
            report.append("\n")

        report_text = "".join(report)

        if output_file:
            output_file.write_text(report_text)

        return report_text

    def run_suite(self, num_trials: int = 1, workloads: Optional[list[BenchmarkWorkload]] = None):
        """Run full benchmark suite."""
        if workloads is None:
            workloads = self.standard_workloads

        print(f"\n{'='*70}")
        print("HOTPANTS BENCHMARK SUITE")
        print(f"{'='*70}\n")

        results_collected = []
        for workload in workloads:
            result = self.run_workload(workload, num_trials=num_trials)
            if result:
                self.save_result(result)
                results_collected.append(result)

        # Generate and display report
        report = self.generate_report()
        print("\n" + "="*70)
        print("REPORT")
        print("="*70)
        print(report)

        # Save report to file
        report_file = self.results_dir / "latest_report.md"
        self.generate_report(report_file)
        print(f"\nReport saved to: {report_file}")

        # Check for regressions and exit with appropriate code
        regressions = self.detect_regressions()
        if regressions:
            print(f"\n⚠️ {len(regressions)} regression(s) detected!")
            return 1
        else:
            print("\n✓ All benchmarks passed!")
            return 0


def main():
    import argparse

    parser = argparse.ArgumentParser(description="HOTPANTS benchmark suite")
    parser.add_argument(
        "--hotpants",
        type=Path,
        default=Path("/home/user/hotpants/build/hotpants"),
        help="Path to hotpants binary"
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        help="Directory for benchmark results"
    )
    parser.add_argument(
        "--trials",
        type=int,
        default=1,
        help="Number of trials per workload"
    )
    parser.add_argument(
        "--regression-threshold",
        type=float,
        default=10.0,
        help="Regression detection threshold (percent)"
    )
    parser.add_argument(
        "--report-only",
        action="store_true",
        help="Only generate report from existing results"
    )

    args = parser.parse_args()

    if not args.hotpants.exists():
        print(f"Error: {args.hotpants} not found")
        return 1

    bench = HotpantsBenchmark(args.hotpants, args.results_dir)

    if args.report_only:
        report = bench.generate_report()
        print(report)
        return 0
    else:
        return bench.run_suite(num_trials=args.trials)


if __name__ == "__main__":
    exit(main())
