"""Shared loading helpers for the figure scripts."""

import glob
import os

import numpy as np

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
VNA_DIR = os.path.join(REPO, "data", "vna")
RAMAN_DIR = os.path.join(REPO, "data", "raman")
OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")

AIR_FILE = "ses06062026_bat0001_s11s21_empty.csv"


def load_export(path):
    """Return a dict with frequency (Hz), |S21| (dB), S21 phase (deg), Re and Im."""
    with open(path, "r") as fh:
        header = fh.readline().rstrip("\n").rstrip(",").split(",")
    col = {name.strip(): i for i, name in enumerate(header)}
    data = np.genfromtxt(path, delimiter=",", skip_header=1)
    mag = data[:, col["S21_Magnitude (linear)"]]
    with np.errstate(divide="ignore"):
        s21_db = 20.0 * np.log10(mag)
    return dict(
        f=data[:, col["Frequency"]],
        db=s21_db,
        phase=data[:, col["S21_Phase"]],
        re=data[:, col["S21_Real"]],
        im=data[:, col["S21_Imaginary"]],
    )


def load_all(data_dir=VNA_DIR):
    """Air baseline plus the five PE and five sand loadings (|S21| in dB)."""
    air = load_export(os.path.join(data_dir, AIR_FILE))
    pe = np.array([load_export(p)["db"]
                   for p in sorted(glob.glob(os.path.join(data_dir, "*_PE_*.csv")))])
    sand = np.array([load_export(p)["db"]
                     for p in sorted(glob.glob(os.path.join(data_dir, "*_sand_*.csv")))])
    if pe.shape[0] != 5 or sand.shape[0] != 5:
        raise SystemExit(f"Expected 5 PE and 5 sand files in {data_dir}")
    return air, pe, sand


def out_path(name):
    os.makedirs(OUT_DIR, exist_ok=True)
    return os.path.join(OUT_DIR, name)
