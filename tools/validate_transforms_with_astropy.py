#!/usr/bin/env python3
"""
Cross-check astroforces frame transforms against Astropy.

This script validates two paths:
1) "gmst_approx" ECI<->ECEF GMST-spin path
2) GCRF<->ITRF path using CIP/EOP inputs

Requires:
  pip install astropy pyerfa numpy
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

try:
    from astropy import coordinates as coord
    from astropy import units as u
    from astropy.time import Time
    from astropy.utils import iers
    import erfa
except ImportError as exc:
    print(
        "error: missing dependency for Astropy validation. Install with: pip install astropy pyerfa numpy",
        file=sys.stderr,
    )
    raise SystemExit(3) from exc


@dataclass
class Vec3:
    x: float
    y: float
    z: float

    def as_np(self) -> np.ndarray:
        return np.array([self.x, self.y, self.z], dtype=float)


def parse_vec(line: str) -> Vec3:
    m = re.search(r"=([+-]?\d.*)$", line.strip())
    if not m:
        raise RuntimeError(f"cannot parse vector line: {line}")
    parts = [float(x) for x in m.group(1).split(",")]
    if len(parts) != 3:
        raise RuntimeError(f"expected 3 components, got {len(parts)} in line: {line}")
    return Vec3(parts[0], parts[1], parts[2])


def run_cli(cli: Path, args: list[str]) -> list[str]:
    out = subprocess.check_output([str(cli), *args], text=True)
    # Keep only key=value data records; ignore logger output lines.
    return [ln for ln in out.splitlines() if ln.strip() and "=" in ln]


def norm(v: np.ndarray) -> float:
    return float(np.linalg.norm(v))


def simple_case(cli: Path, t_utc: Time, r_gcrf_m: Vec3, v_gcrf_mps: Vec3) -> tuple[float, float]:
    lines = run_cli(
        cli,
        [
            "gmst_approx",
            f"{t_utc.unix:.6f}",
            f"{r_gcrf_m.x:.16e}",
            f"{r_gcrf_m.y:.16e}",
            f"{r_gcrf_m.z:.16e}",
            f"{v_gcrf_mps.x:.16e}",
            f"{v_gcrf_mps.y:.16e}",
            f"{v_gcrf_mps.z:.16e}",
        ],
    )
    r_ecef_cpp = parse_vec(lines[0]).as_np()
    v_ecef_cpp = parse_vec(lines[1]).as_np()

    gcrs = coord.GCRS(
        x=r_gcrf_m.x * u.m,
        y=r_gcrf_m.y * u.m,
        z=r_gcrf_m.z * u.m,
        v_x=v_gcrf_mps.x * u.m / u.s,
        v_y=v_gcrf_mps.y * u.m / u.s,
        v_z=v_gcrf_mps.z * u.m / u.s,
        representation_type="cartesian",
        differential_type="cartesian",
        obstime=t_utc,
    )
    itrs = gcrs.transform_to(coord.ITRS(obstime=t_utc))
    r_ecef_ast = itrs.cartesian.xyz.to_value(u.m)
    v_ecef_ast = itrs.cartesian.differentials["s"].d_xyz.to_value(u.m / u.s)
    return norm(r_ecef_cpp - r_ecef_ast), norm(v_ecef_cpp - v_ecef_ast)


def rigorous_case(cli: Path, t_utc: Time, r_gcrf_m: Vec3, v_gcrf_mps: Vec3) -> tuple[float, float]:
    t_tt = t_utc.tt
    jd_tt = float(t_tt.jd)
    jd_utc = float(t_utc.jd)
    x, y, s = erfa.xys00a(t_tt.jd1, t_tt.jd2)
    tab = iers.IERS_Auto.open()
    xp, yp = tab.pm_xy(t_utc)
    dut1 = float(t_utc.delta_ut1_utc)

    lines = run_cli(
        cli,
        [
            "gcrf_to_itrf",
            f"{jd_utc:.16f}",
            f"{jd_tt:.16f}",
            f"{x:.16e}",
            f"{y:.16e}",
            f"{s:.16e}",
            f"{xp.to_value(u.rad):.16e}",
            f"{yp.to_value(u.rad):.16e}",
            f"{dut1:.16e}",
            "0.0",
            "0.0",
            "0.0",
            f"{r_gcrf_m.x:.16e}",
            f"{r_gcrf_m.y:.16e}",
            f"{r_gcrf_m.z:.16e}",
            f"{v_gcrf_mps.x:.16e}",
            f"{v_gcrf_mps.y:.16e}",
            f"{v_gcrf_mps.z:.16e}",
        ],
    )
    r_itrf_cpp = parse_vec(lines[0]).as_np()
    v_itrf_cpp = parse_vec(lines[1]).as_np()

    gcrs = coord.GCRS(
        x=r_gcrf_m.x * u.m,
        y=r_gcrf_m.y * u.m,
        z=r_gcrf_m.z * u.m,
        v_x=v_gcrf_mps.x * u.m / u.s,
        v_y=v_gcrf_mps.y * u.m / u.s,
        v_z=v_gcrf_mps.z * u.m / u.s,
        representation_type="cartesian",
        differential_type="cartesian",
        obstime=t_utc,
    )
    itrs = gcrs.transform_to(coord.ITRS(obstime=t_utc))
    r_itrf_ast = itrs.cartesian.xyz.to_value(u.m)
    v_itrf_ast = itrs.cartesian.differentials["s"].d_xyz.to_value(u.m / u.s)
    return norm(r_itrf_cpp - r_itrf_ast), norm(v_itrf_cpp - v_itrf_ast)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cli", default="build/macos-debug/frame_transform_cli")
    parser.add_argument("--unix-utc", type=float, default=1.0e9)
    args = parser.parse_args()

    cli = Path(args.cli)
    if not cli.exists():
        print(f"error: CLI not found: {cli}", file=sys.stderr)
        return 2

    iers.conf.auto_max_age = None

    t_utc = Time(args.unix_utc, format="unix", scale="utc")
    r = Vec3(7000e3, -1200e3, 2500e3)
    v = Vec3(1200.0, 6800.0, -900.0)

    simple_pos_err_m, simple_vel_err_mps = simple_case(cli, t_utc, r, v)
    rig_pos_err_m, rig_vel_err_mps = rigorous_case(cli, t_utc, r, v)

    print("=== Astropy Transform Cross-Check ===")
    print(f"epoch_utc_unix_s={args.unix_utc:.6f}")
    print(f"simple_pos_err_m={simple_pos_err_m:.6e}")
    print(f"simple_vel_err_mps={simple_vel_err_mps:.6e}")
    print(f"rigorous_pos_err_m={rig_pos_err_m:.6e}")
    print(f"rigorous_vel_err_mps={rig_vel_err_mps:.6e}")
    print("note=machine precision vs astropy requires matching model assumptions and EOP/CIP inputs")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
