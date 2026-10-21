"""
hematite_coords.py
==================
Companion Python module for the Hematite (alpha-Fe2O3) Crystal & Spin-Wave (CPW) /
SAW (IDT) Coordinate Transformation Calculator.

Provides the exact same 3x3 SO(3) rotation matrices, Miller-Bravais conversions,
5-frame transformations (Lab, Mount, Sub, Wave, Cryst), easy-plane magnetic
equilibrium (m_hat, l_hat, B_perp_c), magnet rotation sweeps, and inverse target
angle solvers as the interactive web application (`js/math_core.js`).

Coordinate Systems (5 linked orthonormal right-handed frames):
  1. Lab   (F_Lab)   : (X_L, Y_L, Z_L) - Horizontal magnet plane (X_L, Y_L), vertical Z_L.
                       B_Lab(phi_M) = B0 * (cos(phi_M), sin(phi_M), 0)^T.
  2. Mount (F_Mount) : Sample holder frame at psi_S = 0:
                       - 'IP'  (horizontal): z_m = Z_L, x_m = X_L, y_m = Y_L
                       - 'OOP' (vertical)  : z_m = X_L, x_m = Y_L, y_m = Z_L
  3. Sub   (F_Sub)   : (x_s, y_s, z_s) - Substrate edges (x_s, y_s) and surface normal z_s = n,
                       rotated by psi_S around z_s in the holder, and linked to F_Cryst via
                       cut angles (theta_cut, phi_cut) and optional wafer-flat offset alpha_flat.
  4. Wave  (F_Wave)  : (k_hat, t_hat, n_hat) - Microstructure (IDT / CPW) frame rotated by
                       theta_dev around n_hat from x_s:
                       k_hat = propagation direction, t_hat = n_hat x k_hat (transverse / fingers).
  5. Cryst (F_Cryst) : (e_a, e_m, e_c) - Orthonormal right-handed crystal basis of alpha-Fe2O3:
                       e_a || [11-20] (C2 axis),
                       e_m || [1-100] (perp to C2 in basal plane),
                       e_c || [0001]  (C3 trigonal axis, hard magnetic axis / DMI vector D).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, List, Literal, Optional, Sequence, Tuple, Union

import numpy as np

# Room-temperature lattice constants of hematite (alpha-Fe2O3, space group R-3c, No. 167)
A_LATTICE_DEFAULT: float = 5.0356  # Angstrom
C_LATTICE_DEFAULT: float = 13.7489  # Angstrom
THETA_R_CUT_DEG: float = math.degrees(
    math.atan(C_LATTICE_DEFAULT / (math.sqrt(3.0) * A_LATTICE_DEFAULT))
)  # ~57.60955 deg

FrameName = Literal["Lab", "Mount", "Sub", "Wave", "Cryst"]
MountMode = Literal["OOP", "IP"]
CutPreset = Literal["c-cut", "a-cut", "m-cut", "r-cut", "custom"]

CUT_PRESETS: Dict[str, Dict[str, float]] = {
    "c-cut": {"theta_cut_deg": 0.0, "phi_cut_deg": 0.0},
    "a-cut": {"theta_cut_deg": 90.0, "phi_cut_deg": 0.0},
    "m-cut": {"theta_cut_deg": 90.0, "phi_cut_deg": 90.0},
    "r-cut": {"theta_cut_deg": THETA_R_CUT_DEG, "phi_cut_deg": 90.0},
}


def wrap_deg_360(deg: float) -> float:
    val = float(deg) % 360.0
    if val < 0.0:
        val += 360.0
    return 0.0 if abs(val) < 1e-12 or abs(val - 360.0) < 1e-12 else val


def wrap_deg_180(deg: float) -> float:
    w = wrap_deg_360(deg + 180.0) - 180.0
    return 180.0 if abs(w + 180.0) < 1e-12 else w


def normalize3(v: Sequence[float], fallback: Tuple[float, float, float] = (1.0, 0.0, 0.0)) -> np.ndarray:
    arr = np.asarray(v, dtype=float).reshape(3)
    n = float(np.linalg.norm(arr))
    if n < 1e-12:
        return np.asarray(fallback, dtype=float)
    return arr / n


def angle_between_deg(a: Sequence[float], b: Sequence[float]) -> float:
    va = np.asarray(a, dtype=float).reshape(3)
    vb = np.asarray(b, dtype=float).reshape(3)
    na = float(np.linalg.norm(va))
    nb = float(np.linalg.norm(vb))
    if na < 1e-12 or nb < 1e-12:
        return 0.0
    c = float(np.clip(np.dot(va, vb) / (na * nb), -1.0, 1.0))
    return math.degrees(math.acos(c))


def rot_x(angle_deg: float) -> np.ndarray:
    r = math.radians(angle_deg)
    c, s = math.cos(r), math.sin(r)
    return np.array([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]], dtype=float)


def rot_y(angle_deg: float) -> np.ndarray:
    r = math.radians(angle_deg)
    c, s = math.cos(r), math.sin(r)
    return np.array([[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]], dtype=float)


def rot_z(angle_deg: float) -> np.ndarray:
    r = math.radians(angle_deg)
    c, s = math.cos(r), math.sin(r)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]], dtype=float)


# ============================================================================
# Miller-Bravais Crystallography in F_Cryst = (e_a, e_m, e_c)
# ============================================================================

def miller_bravais_direction_to_crystal(
    uvtw: Sequence[float],
    a: float = A_LATTICE_DEFAULT,
    c: float = C_LATTICE_DEFAULT,
    normalize_vec: bool = True,
) -> np.ndarray:
    """
    Converts 4-index Miller-Bravais direction [u, v, t, w] (where t = -(u+v))
    to a 3D vector in the orthonormal crystal basis F_Cryst = (e_a, e_m, e_c).
    """
    u, v, _, w = uvtw
    vec = np.array(
        [
            1.5 * a * (u + v),
            0.5 * math.sqrt(3.0) * a * (u - v),
            c * w,
        ],
        dtype=float,
    )
    return normalize3(vec) if normalize_vec else vec


def crystal_to_miller_bravais_direction(
    vec_cryst: Sequence[float],
    a: float = A_LATTICE_DEFAULT,
    c: float = C_LATTICE_DEFAULT,
) -> np.ndarray:
    """
    Converts a 3D vector in F_Cryst = (e_a, e_m, e_c) into Miller-Bravais direction
    indices [u, v, t, w].
    """
    xa, xm, xc = np.asarray(vec_cryst, dtype=float).reshape(3)
    u_plus_v = (2.0 * xa) / (3.0 * a)
    u_minus_v = (2.0 * xm) / (math.sqrt(3.0) * a)
    u = 0.5 * (u_plus_v + u_minus_v)
    v = 0.5 * (u_plus_v - u_minus_v)
    t = -(u + v)
    w = xc / c
    return np.array([u, v, t, w], dtype=float)


def miller_bravais_plane_to_crystal_normal(
    hkil: Sequence[float],
    a: float = A_LATTICE_DEFAULT,
    c: float = C_LATTICE_DEFAULT,
    normalize_vec: bool = True,
) -> np.ndarray:
    """
    Converts 4-index Miller-Bravais plane (h, k, i, l) (where i = -(h+k))
    to the reciprocal lattice vector / surface normal in F_Cryst = (e_a, e_m, e_c).
    """
    h, k, _, l = hkil
    g = np.array(
        [
            (h + k) / a,
            (h - k) / (math.sqrt(3.0) * a),
            l / c,
        ],
        dtype=float,
    )
    return normalize3(g) if normalize_vec else g


def crystal_normal_to_miller_bravais_plane(
    vec_cryst: Sequence[float],
    a: float = A_LATTICE_DEFAULT,
    c: float = C_LATTICE_DEFAULT,
) -> np.ndarray:
    """
    Converts a surface normal vector in F_Cryst = (e_a, e_m, e_c) into Miller-Bravais
    plane indices (h, k, i, l).
    """
    ga, gm, gc = np.asarray(vec_cryst, dtype=float).reshape(3)
    h_plus_k = ga * a
    h_minus_k = gm * math.sqrt(3.0) * a
    h = 0.5 * (h_plus_k + h_minus_k)
    k = 0.5 * (h_plus_k - h_minus_k)
    i = -(h + k)
    l = gc * c
    return np.array([h, k, i, l], dtype=float)


# ============================================================================
# Core 3x3 Rotation Matrices Linking the 5 Frames
# ============================================================================

def r_cryst_from_sub_base(theta_cut_deg: float, phi_cut_deg: float) -> np.ndarray:
    """
    Base rotation matrix R_{Cryst <- Sub}^{(0)}(theta_cut_deg, phi_cut_deg)
    whose columns are the unrotated substrate axes (x_s^(0), y_s^(0), z_s^(0))
    expressed in F_Cryst = (e_a, e_m, e_c).
    """
    th = math.radians(theta_cut_deg)
    ph = math.radians(phi_cut_deg)
    st, ct = math.sin(th), math.cos(th)
    sp, cp = math.sin(ph), math.cos(ph)

    if abs(st) < 1e-10:
        if ct >= 0.0:
            return np.column_stack(([cp, sp, 0.0], [-sp, cp, 0.0], [0.0, 0.0, 1.0]))
        return np.column_stack(([cp, sp, 0.0], [sp, -cp, 0.0], [0.0, 0.0, -1.0]))

    xs0 = np.array([-sp, cp, 0.0], dtype=float)
    ys0 = np.array([-ct * cp, -ct * sp, st], dtype=float)
    zs0 = np.array([st * cp, st * sp, ct], dtype=float)
    return np.column_stack((xs0, ys0, zs0))


def r_cryst_from_sub(
    theta_cut_deg: float, phi_cut_deg: float, alpha_flat_deg: float = 0.0
) -> np.ndarray:
    """
    Full rotation matrix R_{Cryst <- Sub}(theta_cut_deg, phi_cut_deg, alpha_flat_deg).
    """
    return r_cryst_from_sub_base(theta_cut_deg, phi_cut_deg) @ rot_z(alpha_flat_deg)


def r_lab_from_mount(mount_mode: MountMode = "OOP") -> np.ndarray:
    """
    Mount-to-Lab rotation matrix R_{Lab <- Mount}(mount_mode) at psi_S = 0.
    """
    if mount_mode == "IP":
        return np.eye(3, dtype=float)
    # OOP configuration: z_s = X_L, x_s(0) = Y_L, y_s(0) = Z_L
    return np.array(
        [
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ],
        dtype=float,
    )


def r_lab_from_sub(mount_mode: MountMode = "OOP", psi_s_deg: float = 0.0) -> np.ndarray:
    """
    Substrate-to-Lab rotation matrix R_{Lab <- Sub}(mount_mode, psi_s_deg).
    """
    return r_lab_from_mount(mount_mode) @ rot_z(psi_s_deg)


def r_sub_from_wave(theta_dev_deg: float = 0.0) -> np.ndarray:
    """
    Wave-to-Substrate rotation matrix R_{Sub <- Wave}(theta_dev_deg).
    """
    return rot_z(theta_dev_deg)


# ============================================================================
# High-Level Setup Class for Evaluation Scripts
# ============================================================================

@dataclass
class HematiteSetup:
    name: str = "Setup A"
    mount_mode: MountMode = "OOP"
    cut_preset: CutPreset = "r-cut"
    theta_cut_deg: float = THETA_R_CUT_DEG
    phi_cut_deg: float = 90.0
    alpha_flat_deg: float = 0.0
    device_type: str = "CPW"
    theta_dev_deg: float = 0.0
    psi_s_deg: float = 0.0
    phi_m_deg: float = 90.0
    b0_tesla: float = 0.5
    a_lattice: float = A_LATTICE_DEFAULT
    c_lattice: float = C_LATTICE_DEFAULT

    def resolve_cut_angles(self) -> Tuple[float, float]:
        if self.cut_preset in CUT_PRESETS and self.cut_preset != "custom":
            if self.cut_preset == "r-cut":
                th = math.degrees(
                    math.atan(self.c_lattice / (math.sqrt(3.0) * self.a_lattice))
                )
                return th, CUT_PRESETS["r-cut"]["phi_cut_deg"]
            preset = CUT_PRESETS[self.cut_preset]
            return preset["theta_cut_deg"], preset["phi_cut_deg"]
        return float(self.theta_cut_deg), float(self.phi_cut_deg)

    def frame_matrices(self) -> Dict[str, np.ndarray]:
        th_cut, ph_cut = self.resolve_cut_angles()
        R_Lab_Mount = r_lab_from_mount(self.mount_mode)
        R_Mount_Sub = rot_z(self.psi_s_deg)
        R_Lab_Sub = R_Lab_Mount @ R_Mount_Sub
        R_Sub_Wave = r_sub_from_wave(self.theta_dev_deg)
        R_Cryst_Sub = r_cryst_from_sub(th_cut, ph_cut, self.alpha_flat_deg)

        R_Sub_from = {
            "Sub": np.eye(3, dtype=float),
            "Mount": R_Mount_Sub.T,
            "Lab": R_Lab_Sub.T,
            "Wave": R_Sub_Wave,
            "Cryst": R_Cryst_Sub.T,
        }
        R_to_from_Sub = {
            "Sub": np.eye(3, dtype=float),
            "Mount": R_Mount_Sub,
            "Lab": R_Lab_Sub,
            "Wave": R_Sub_Wave.T,
            "Cryst": R_Cryst_Sub,
        }
        return {
            "R_Lab_Mount": R_Lab_Mount,
            "R_Mount_Sub": R_Mount_Sub,
            "R_Lab_Sub": R_Lab_Sub,
            "R_Sub_Wave": R_Sub_Wave,
            "R_Cryst_Sub": R_Cryst_Sub,
            "R_Sub_from": R_Sub_from,
            "R_to_from_Sub": R_to_from_Sub,
        }

    def get_rotation_matrix(self, to_frame: FrameName, from_frame: FrameName) -> np.ndarray:
        """
        Returns the 3x3 rotation matrix R_{to_frame <- from_frame}.
        """
        if to_frame == from_frame:
            return np.eye(3, dtype=float)
        mats = self.frame_matrices()
        return mats["R_to_from_Sub"][to_frame] @ mats["R_Sub_from"][from_frame]

    def transform_vector(
        self, vec: Sequence[float], from_frame: FrameName, to_frame: FrameName
    ) -> np.ndarray:
        R = self.get_rotation_matrix(to_frame, from_frame)
        return R @ np.asarray(vec, dtype=float).reshape(3)

    def vector_in_all_frames(
        self, vec: Sequence[float], from_frame: FrameName
    ) -> Dict[FrameName, np.ndarray]:
        mats = self.frame_matrices()
        v_sub = mats["R_Sub_from"][from_frame] @ np.asarray(vec, dtype=float).reshape(3)
        return {
            frame: mats["R_to_from_Sub"][frame] @ v_sub
            for frame in ("Lab", "Mount", "Sub", "Wave", "Cryst")
        }

    def evaluate(self, phi_m_deg: Optional[float] = None) -> Dict[str, object]:
        """
        Evaluates all vectors, projections, and easy-plane equilibrium state.
        """
        phi_m = self.phi_m_deg if phi_m_deg is None else float(phi_m_deg)
        phi_rad = math.radians(phi_m)
        b0 = float(self.b0_tesla)

        b_hat_lab = np.array([math.cos(phi_rad), math.sin(phi_rad), 0.0], dtype=float)
        b_vec_lab = b0 * b_hat_lab

        b_hat = self.vector_in_all_frames(b_hat_lab, "Lab")
        B = self.vector_in_all_frames(b_vec_lab, "Lab")

        k_hat = self.vector_in_all_frames([1.0, 0.0, 0.0], "Wave")
        t_hat = self.vector_in_all_frames([0.0, 1.0, 0.0], "Wave")
        n_hat = self.vector_in_all_frames([0.0, 0.0, 1.0], "Sub")
        xs_hat = self.vector_in_all_frames([1.0, 0.0, 0.0], "Sub")
        ys_hat = self.vector_in_all_frames([0.0, 1.0, 0.0], "Sub")
        ea_hat = self.vector_in_all_frames([1.0, 0.0, 0.0], "Cryst")
        em_hat = self.vector_in_all_frames([0.0, 1.0, 0.0], "Cryst")
        ec_hat = self.vector_in_all_frames([0.0, 0.0, 1.0], "Cryst")

        # Sample-plane projection:
        b_in_sample_ratio = float(math.hypot(b_hat["Sub"][0], b_hat["Sub"][1]))
        b_perp_sample_ratio = float(b_hat["Sub"][2])
        phi_b_sub_deg = (
            wrap_deg_360(math.degrees(math.atan2(b_hat["Sub"][1], b_hat["Sub"][0])))
            if b_in_sample_ratio > 1e-9
            else 0.0
        )
        phi_b_wave_deg = (
            wrap_deg_180(math.degrees(math.atan2(b_hat["Wave"][1], b_hat["Wave"][0])))
            if b_in_sample_ratio > 1e-9
            else 0.0
        )

        # Easy-plane (0001) projection & equilibrium (m_hat || B_easy, l_hat = m_hat x e_c):
        b_cryst = b_hat["Cryst"]
        b_easy_ratio = float(math.hypot(b_cryst[0], b_cryst[1]))
        b_parallel_c_ratio = float(b_cryst[2])
        is_singular = b_easy_ratio < 1e-7

        if not is_singular:
            m_hat_cryst = np.array([b_cryst[0] / b_easy_ratio, b_cryst[1] / b_easy_ratio, 0.0])
        else:
            m_hat_cryst = np.array([1.0, 0.0, 0.0])
        l_hat_cryst = np.array([m_hat_cryst[1], -m_hat_cryst[0], 0.0])

        m_hat = self.vector_in_all_frames(m_hat_cryst, "Cryst")
        l_hat = self.vector_in_all_frames(l_hat_cryst, "Cryst")
        b_easy_proj = self.vector_in_all_frames(
            [b_cryst[0] * b0, b_cryst[1] * b0, 0.0], "Cryst"
        )

        angle_b_k_deg = angle_between_deg(b_hat["Sub"], k_hat["Sub"])
        angle_b_n_deg = angle_between_deg(b_hat["Sub"], n_hat["Sub"])
        angle_b_c_deg = angle_between_deg(b_hat["Cryst"], [0.0, 0.0, 1.0])
        angle_m_k_deg = angle_between_deg(m_hat["Sub"], k_hat["Sub"])
        angle_l_k_deg = angle_between_deg(l_hat["Sub"], k_hat["Sub"])

        return {
            "vectors": {
                "B": B,
                "b_hat": b_hat,
                "b_easy_proj": b_easy_proj,
                "m_hat": m_hat,
                "l_hat": l_hat,
                "k_hat": k_hat,
                "t_hat": t_hat,
                "n_hat": n_hat,
                "xs_hat": xs_hat,
                "ys_hat": ys_hat,
                "ea_hat": ea_hat,
                "em_hat": em_hat,
                "ec_hat": ec_hat,
            },
            "metrics": {
                "b0_tesla": b0,
                "b_in_sample_ratio": b_in_sample_ratio,
                "b_perp_sample_ratio": b_perp_sample_ratio,
                "b_easy_ratio": b_easy_ratio,
                "b_parallel_c_ratio": b_parallel_c_ratio,
                "b_easy_tesla": b_easy_ratio * b0,
                "b_parallel_c_tesla": b_parallel_c_ratio * b0,
                "is_easy_plane_singular": is_singular,
                "phi_b_sub_deg": phi_b_sub_deg,
                "phi_b_wave_deg": phi_b_wave_deg,
                "angle_b_k_deg": angle_b_k_deg,
                "angle_b_n_deg": angle_b_n_deg,
                "angle_b_c_deg": angle_b_c_deg,
                "angle_m_k_deg": angle_m_k_deg,
                "angle_l_k_deg": angle_l_k_deg,
                "b_dot_k": float(np.dot(b_hat["Sub"], k_hat["Sub"])),
                "m_dot_k": float(np.dot(m_hat["Sub"], k_hat["Sub"])),
                "l_dot_k": float(np.dot(l_hat["Sub"], k_hat["Sub"])),
                "cpw_torque_m": float(np.linalg.norm(np.cross(k_hat["Sub"], m_hat["Sub"]))),
                "cpw_torque_l": float(np.linalg.norm(np.cross(k_hat["Sub"], l_hat["Sub"]))),
            },
        }

    def sweep_magnet(self, phi_m_array_deg: Sequence[float]) -> Dict[str, np.ndarray]:
        """
        Evaluates a magnet rotation sweep over an array of magnet angles phi_M (in degrees).
        """
        phi_arr = np.asarray(phi_m_array_deg, dtype=float)
        n = len(phi_arr)
        out = {
            "phi_m_deg": phi_arr,
            "b_easy_ratio": np.zeros(n),
            "b_in_sample_ratio": np.zeros(n),
            "angle_b_k_deg": np.zeros(n),
            "angle_m_k_deg": np.zeros(n),
            "angle_l_k_deg": np.zeros(n),
            "b_dot_k": np.zeros(n),
            "m_dot_k": np.zeros(n),
            "l_dot_k": np.zeros(n),
            "cpw_torque_m": np.zeros(n),
            "cpw_torque_l": np.zeros(n),
        }
        for idx, phi in enumerate(phi_arr):
            m = self.evaluate(phi_m_deg=float(phi))["metrics"]
            for k in out:
                if k != "phi_m_deg":
                    out[k][idx] = m[k]
        return out

    def solve_target_angles(
        self, target_vec_sub: Sequence[float]
    ) -> Dict[str, List[Dict[str, Union[str, float]]]]:
        """
        Inverse solver returning the (psi_S, phi_M) angles required to align B
        with a desired unit vector `target_vec_sub` expressed in F_Sub.
        """
        u = normalize3(target_vec_sub)
        ux, uy, uz = float(u[0]), float(u[1]), float(u[2])
        xy_norm = math.hypot(ux, uy)

        phi1 = math.degrees(math.acos(max(-1.0, min(1.0, uz))))
        psi1 = (
            wrap_deg_360(math.degrees(math.atan2(-uy, ux)))
            if xy_norm > 1e-9
            else wrap_deg_360(self.psi_s_deg)
        )
        phi2 = wrap_deg_360(360.0 - phi1)
        psi2 = wrap_deg_360(psi1 + 180.0)

        ip_err_deg = math.degrees(math.acos(max(-1.0, min(1.0, xy_norm))))
        az_sub_deg = (
            wrap_deg_360(math.degrees(math.atan2(uy, ux))) if xy_norm > 1e-9 else 0.0
        )
        cur_psi = wrap_deg_360(self.psi_s_deg)

        return {
            "oop_solutions": [
                {
                    "mount_mode": "OOP",
                    "psi_s_deg": psi1,
                    "phi_m_deg": phi1,
                    "error_deg": 0.0,
                },
                {
                    "mount_mode": "OOP",
                    "psi_s_deg": psi2,
                    "phi_m_deg": phi2,
                    "error_deg": 0.0,
                },
            ],
            "ip_solutions": [
                {
                    "mount_mode": "IP",
                    "psi_s_deg": cur_psi,
                    "phi_m_deg": wrap_deg_360(cur_psi + az_sub_deg),
                    "error_deg": ip_err_deg,
                },
                {
                    "mount_mode": "IP",
                    "psi_s_deg": wrap_deg_360(-az_sub_deg),
                    "phi_m_deg": 0.0,
                    "error_deg": ip_err_deg,
                },
            ],
        }

    def solve_sweep_landmarks_for_fixed_psi(self) -> Dict[str, float]:
        """
        For the currently fixed manual sample mounting angle psi_S, computes the
        exact magnet angles phi_M during a live magnet sweep where:
          - B is closest to +k_hat (`phi_m_b_parallel_k`)
          - B is closest to +t_hat (`phi_m_b_perp_k`)
          - B lies 100% inside the (0001) Easy-Plane (`phi_m_in_easy_plane`)
          - B lies 100% inside the Sample Plane (`phi_m_in_sample_plane`)
          - Weak FM moment m_hat aligns with +k_easy (`phi_m_m_align_k`)
          - Neel vector l_hat aligns with +k_easy (`phi_m_l_align_k`)
        """
        st = self.evaluate()
        mats = self.frame_matrices()
        R_Lab_Sub = mats["R_Lab_Sub"]
        R_Sub_Cryst = mats["R_Sub_from"]["Cryst"]

        def phi_max_dot(v_lab: np.ndarray) -> float:
            h = math.hypot(float(v_lab[0]), float(v_lab[1]))
            return wrap_deg_360(math.degrees(math.atan2(float(v_lab[1]), float(v_lab[0])))) if h > 1e-9 else 0.0

        def phi_align_easy(u_easy_cryst: Sequence[float]) -> float:
            u_norm = normalize3([u_easy_cryst[0], u_easy_cryst[1], 0.0])
            v_perp_cryst = np.array([-u_norm[1], u_norm[0], 0.0])
            v_perp_lab = R_Lab_Sub @ (R_Sub_Cryst @ v_perp_cryst)
            u_easy_lab = R_Lab_Sub @ (R_Sub_Cryst @ u_norm)
            if math.hypot(float(v_perp_lab[0]), float(v_perp_lab[1])) < 1e-9:
                return phi_max_dot(u_easy_lab)
            c1 = wrap_deg_360(math.degrees(math.atan2(-float(v_perp_lab[0]), float(v_perp_lab[1]))))
            c2 = wrap_deg_360(c1 + 180.0)
            r1 = math.radians(c1)
            dot1 = math.cos(r1) * float(u_easy_lab[0]) + math.sin(r1) * float(u_easy_lab[1])
            return c1 if dot1 >= 0.0 else c2

        k_lab = st["vectors"]["k_hat"]["Lab"]
        t_lab = st["vectors"]["t_hat"]["Lab"]
        ec_lab = st["vectors"]["ec_hat"]["Lab"]
        n_lab = st["vectors"]["n_hat"]["Lab"]
        k_cryst = st["vectors"]["k_hat"]["Cryst"]
        k_easy = normalize3([k_cryst[0], k_cryst[1], 0.0])

        phi_easy = (
            wrap_deg_360(math.degrees(math.atan2(-float(ec_lab[0]), float(ec_lab[1]))))
            if math.hypot(float(ec_lab[0]), float(ec_lab[1])) > 1e-9
            else 0.0
        )
        if phi_easy >= 180.0:
            phi_easy -= 180.0
        phi_samp = (
            wrap_deg_360(math.degrees(math.atan2(-float(n_lab[0]), float(n_lab[1]))))
            if math.hypot(float(n_lab[0]), float(n_lab[1])) > 1e-9
            else 90.0
        )
        if phi_samp >= 180.0:
            phi_samp -= 180.0

        return {
            "phi_m_b_parallel_k": phi_max_dot(k_lab),
            "phi_m_b_perp_k": phi_max_dot(t_lab),
            "phi_m_in_easy_plane": phi_easy,
            "phi_m_in_sample_plane": phi_samp,
            "phi_m_m_align_k": phi_align_easy(k_easy),
            "phi_m_l_align_k": phi_align_easy([-k_easy[1], k_easy[0], 0.0]),
        }


if __name__ == "__main__":
    print(f"Hematite r-cut (1-102) angle theta_r = {THETA_R_CUT_DEG:.5f} deg")
    demo = HematiteSetup(
        name="r-cut OOP Demo",
        mount_mode="OOP",
        cut_preset="r-cut",
        theta_dev_deg=35.0,
        psi_s_deg=25.0,
        phi_m_deg=70.0,
        b0_tesla=0.5,
    )
    res = demo.evaluate()
    print("B in F_Sub  :", np.round(res["vectors"]["B"]["Sub"], 5))
    print("B in F_Wave :", np.round(res["vectors"]["B"]["Wave"], 5))
    print("B in F_Cryst:", np.round(res["vectors"]["B"]["Cryst"], 5))
    print("m_hat in F_Wave:", np.round(res["vectors"]["m_hat"]["Wave"], 5))
    print("l_hat in F_Wave:", np.round(res["vectors"]["l_hat"]["Wave"], 5))
    print("Metrics:", {k: round(v, 4) if isinstance(v, float) else v for k, v in res["metrics"].items()})
