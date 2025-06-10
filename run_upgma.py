#!/usr/bin/env python3
"""
run_upgma.py

Cross-platform wrapper to build and run the UPGMA container (Docker or Podman) in one of two modes:

  1) Sequences mode (default): uses the built-in DNA sequences hardcoded in upgma.py
     → `python run_upgma.py --mode sequences`

  2) Matrix mode: supply your own distance-matrix file
     → `python run_upgma.py --mode matrix --file distances.csv`

This script will:

  1. Detect whether 'docker' or 'podman' is installed.
  2. Build (or rebuild) the image named `upgma:latest` from the current folder.
  3. Run the container, mounting ./output and (if requested) your distances file.
  4. In “matrix” mode, it passes `python upgma.py --distances /app/distances.csv`.
     In “sequences” mode, it simply relies on the Dockerfile's default CMD.

Usage Examples:
  # 1) Sequences mode (no external file):
  python run_upgma.py --mode sequences

  # 2) Matrix mode (distance-matrix CSV/TSV):
  python run_upgma.py --mode matrix --file ./distances.csv
"""

import argparse
import os
import shutil
import subprocess
import sys

# ───────────── CLI PARSER ─────────────
def parse_args():
    p = argparse.ArgumentParser(
        prog="run_upgma.py",
        description="Build-and-run UPGMA container (Docker or Podman)."
    )
    p.add_argument(
        "--mode",
        choices=["sequences", "matrix"],
        default="sequences",
        help=(
            "'sequences' to use hard-coded DNA in the container (default), "
            "or 'matrix' to supply your own distance matrix."
        )
    )
    p.add_argument(
        "--file", "-f",
        metavar="PATH",
        default="distances.csv",
        help=(
            "Path to your distance-matrix file (CSV). "
            "Required if --mode matrix. (default: distances.csv)"
        )
    )
    p.add_argument(
        "--image", "-i",
        default="upgma:latest",
        help="Container image name to build/run (default: upgma:latest)."
    )
    return p.parse_args()

# ───────────── HELPER: FIND CONTAINER ENGINE ─────────────
def find_container_engine():
    """
    Detect and return either 'docker' or 'podman' if:
      - the CLI tool exists
      - the service/daemon is actually responding

    Returns:
        str: either 'docker' or 'podman'

    Raises:
        SystemExit: if neither engine is usable
    """
    def is_usable(cmd):
        try:
            subprocess.run([cmd, "info"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
            return True
        except (subprocess.CalledProcessError, FileNotFoundError):
            return False

    if shutil.which("docker") and is_usable("docker"):
        return "docker"
    if shutil.which("podman") and is_usable("podman"):
        return "podman"

    print("Error: Neither Docker nor Podman is installed and running.", file=sys.stderr)
    sys.exit(1)

# ───────────── MAIN ─────────────
if __name__ == "__main__":
    args = parse_args()

    engine = find_container_engine()
    image_name = args.image

    # 1) Build (or rebuild) the image from the current directory
    print(f"▶ Building image '{image_name}' using {engine}...")
    build_cmd = [engine, "build", "-t", image_name, "."]
    try:
        subprocess.run(build_cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error: '{engine} build' returned nonzero exit code {e.returncode}.", file=sys.stderr)
        sys.exit(e.returncode)

    # 2) Ensure output folder exists (host side)
    host_output = os.path.abspath("output")
    os.makedirs(host_output, exist_ok=True)

    # 3) Build volume bindings
    #
    #    - Always bind‐mount host ./output → container /app/output
    #    - If in matrix mode, also bind‐mount host distances file → /app/distances.csv
    binds = []
    host_output_posix = host_output.replace("\\", "/")
    binds.append(f"{host_output_posix}:/app/output")

    matrix_arg = []
    if args.mode == "matrix":
        if not args.file:
            print("Error: --mode matrix requires --file <distance_matrix.csv>", file=sys.stderr)
            sys.exit(1)
        if not os.path.isfile(args.file):
            print(f"Error: distance matrix file '{args.file}' does not exist.", file=sys.stderr)
            sys.exit(1)

        host_matrix = os.path.abspath(args.file)
        host_matrix_posix = host_matrix.replace("\\", "/")
        binds.append(f"{host_matrix_posix}:/app/distances.csv")
        matrix_arg = ["python", "upgma.py", "--distances", "/app/distances.csv"]

    # 4) Construct the container run command:
    #
    #    For sequences mode:
    #      <engine> run --rm -v host_output:/app/output image_name
    #
    #    For matrix mode:
    #      <engine> run --rm -v host_output:/app/output -v host_matrix:/app/distances.csv image_name \
    #                   python upgma.py --distances /app/distances.csv
    #
    cmd = [engine, "run", "--rm"]
    for b in binds:
        cmd += ["-v", b]
    cmd.append(image_name)
    cmd += matrix_arg

    # 5) Print & execute
    print("\n▶ Running command:")
    print(" ".join(cmd))
    print("\n==========================\n")
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error: container run returned nonzero exit code {e.returncode}", file=sys.stderr)
        sys.exit(e.returncode)

    print("\n▷ Finished! Check 'output/dendrogram.png' (and console logs) for results.")
