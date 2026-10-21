/**
 * math_core.js
 * Core 3D linear algebra, Miller-Bravais crystallography, 5-frame coordinate
 * transformations, easy-plane equilibrium physics, and inverse angle solver
 * for spin-wave (CPW) and SAW (IDT) experiments on hematite (alpha-Fe2O3).
 *
 * Coordinate Frames:
 *   1. Lab   (F_Lab)   : (X_L, Y_L, Z_L) - Horizontal magnet plane (X_L, Y_L), vertical Z_L.
 *                        B_Lab(phi_M) = B0 * (cos(phi_M), sin(phi_M), 0)^T.
 *   2. Mount (F_Mount) : Sample holder frame at psi_S = 0.
 *                        - IP  (horizontal): z_m = Z_L, x_m = X_L, y_m = Y_L (I_3)
 *                        - OOP (vertical)  : z_m = X_L, x_m = Y_L, y_m = Z_L
 *   3. Sub   (F_Sub)   : (x_s, y_s, z_s) - Substrate edges (x_s, y_s) and surface normal z_s = n,
 *                        rotated by psi_S around z_s in the holder, and linked to F_Cryst via
 *                        cut angles (theta_cut, phi_cut) and optional wafer-flat offset alpha_flat.
 *   4. Wave  (F_Wave)  : (k_hat, t_hat, n_hat) - Microstructure (IDT / CPW) frame rotated by
 *                        theta_dev around n_hat from x_s:
 *                        k_hat = propagation direction, t_hat = n_hat x k_hat (transverse / fingers).
 *   5. Cryst (F_Cryst) : (e_a, e_m, e_c) - Orthonormal right-handed crystal basis of alpha-Fe2O3:
 *                        e_a || [11-20] (C2 axis), e_m || [1-100] (perp to C2 in basal plane),
 *                        e_c || [0001]  (C3 axis, hard magnetic axis / DMI vector D).
 */

export const HEMATITE_CONSTANTS = {
    a: 5.0356, // Angstrom (alpha-Fe2O3 room-temperature lattice constant a)
    c: 13.7489, // Angstrom (alpha-Fe2O3 room-temperature lattice constant c)
    spaceGroup: 'R-3c (No. 167)',
    get thetaRCutDeg() {
        return (Math.atan(this.c / (Math.sqrt(3) * this.a)) * 180) / Math.PI;
    }
};

export const DEG2RAD = Math.PI / 180;
export const RAD2DEG = 180 / Math.PI;

export function degToRad(deg) {
    return deg * DEG2RAD;
}

export function radToDeg(rad) {
    return rad * RAD2DEG;
}

export function wrapDeg360(deg) {
    const m = deg % 360;
    return m < 0 ? m + 360 : (Object.is(m, -0) ? 0 : m);
}

export function wrapDeg180(deg) {
    const w = wrapDeg360(deg + 180) - 180;
    return w === -180 ? 180 : w;
}

// ============================================================================
// 1. 3D Vector & 3x3 Matrix Algebra (Exact SO(3) Operations)
// ============================================================================

export function vec3(x = 0, y = 0, z = 0) {
    return [x, y, z];
}

export function add3(a, b) {
    return [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
}

export function sub3(a, b) {
    return [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
}

export function scale3(v, s) {
    return [v[0] * s, v[1] * s, v[2] * s];
}

export function dot3(a, b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

export function cross3(a, b) {
    return [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0]
    ];
}

export function norm3(v) {
    return Math.hypot(v[0], v[1], v[2]);
}

export function normalize3(v, fallback = [1, 0, 0]) {
    const n = norm3(v);
    if (n < 1e-12) return [...fallback];
    return [v[0] / n, v[1] / n, v[2] / n];
}

export function angleBetweenDeg(a, b) {
    const na = norm3(a);
    const nb = norm3(b);
    if (na < 1e-12 || nb < 1e-12) return 0;
    const c = Math.max(-1, Math.min(1, dot3(a, b) / (na * nb)));
    return Math.acos(c) * RAD2DEG;
}

export function eye3() {
    return [
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1]
    ];
}

export function matFromCols(c0, c1, c2) {
    return [
        [c0[0], c1[0], c2[0]],
        [c0[1], c1[1], c2[1]],
        [c0[2], c1[2], c2[2]]
    ];
}

export function matCols(M) {
    return [
        [M[0][0], M[1][0], M[2][0]],
        [M[0][1], M[1][1], M[2][1]],
        [M[0][2], M[1][2], M[2][2]]
    ];
}

export function transpose3(M) {
    return [
        [M[0][0], M[1][0], M[2][0]],
        [M[0][1], M[1][1], M[2][1]],
        [M[0][2], M[1][2], M[2][2]]
    ];
}

export function matVec3(M, v) {
    return [
        M[0][0] * v[0] + M[0][1] * v[1] + M[0][2] * v[2],
        M[1][0] * v[0] + M[1][1] * v[1] + M[1][2] * v[2],
        M[2][0] * v[0] + M[2][1] * v[1] + M[2][2] * v[2]
    ];
}

export function matMul3(A, B) {
    const out = [
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0]
    ];
    for (let i = 0; i < 3; i++) {
        for (let j = 0; j < 3; j++) {
            out[i][j] = A[i][0] * B[0][j] + A[i][1] * B[1][j] + A[i][2] * B[2][j];
        }
    }
    return out;
}

export function rotX(angleDeg) {
    const r = degToRad(angleDeg);
    const c = Math.cos(r);
    const s = Math.sin(r);
    return [
        [1, 0, 0],
        [0, c, -s],
        [0, s, c]
    ];
}

export function rotY(angleDeg) {
    const r = degToRad(angleDeg);
    const c = Math.cos(r);
    const s = Math.sin(r);
    return [
        [c, 0, s],
        [0, 1, 0],
        [-s, 0, c]
    ];
}

export function rotZ(angleDeg) {
    const r = degToRad(angleDeg);
    const c = Math.cos(r);
    const s = Math.sin(r);
    return [
        [c, -s, 0],
        [s, c, 0],
        [0, 0, 1]
    ];
}

// ============================================================================
// 2. Miller-Bravais Crystallography in F_Cryst = (e_a, e_m, e_c)
// ============================================================================
// Basis definition:
//   e_a = [1, 0, 0]^T || [11-20]  (a1 + a2 direction, C2 axis)
//   e_m = [0, 1, 0]^T || [1-100]  (a1 - a2 direction, mirror-plane basal axis)
//   e_c = [0, 0, 1]^T || [0001]   (c-axis, C3 trigonal axis)
//
// Direct lattice vectors:
//   a1 = a * ( 1/2 e_a + sqrt(3)/2 e_m )
//   a2 = a * ( 1/2 e_a - sqrt(3)/2 e_m )
//   a3 = -(a1 + a2) = -a * e_a
//   c  = c * e_c
//
// For 4-index Miller-Bravais direction [u, v, t, w] with t = -(u + v):
//   R_[uvtw] = u*a1 + v*a2 + t*a3 + w*c
//            = (3/2 * a * (u + v)) e_a + (sqrt(3)/2 * a * (u - v)) e_m + (c * w) e_c
//
// For 4-index Miller-Bravais plane (h, k, i, l) with i = -(h + k):
//   Reciprocal vector G_(hkil) in (e_a, e_m, e_c):
//   G_(hkil) = ((h + k)/a) e_a + ((h - k)/(sqrt(3)*a)) e_m + (l/c) e_c
//   Satisfies G_(hkil) . R_[uvtw] = h*u + k*v + i*t + l*w and |G|^2 = 1/d_{hkil}^2.

export function millerBravaisDirectionToCrystal(
    uvtw,
    a = HEMATITE_CONSTANTS.a,
    c = HEMATITE_CONSTANTS.c,
    normalizeVec = true
) {
    const [u, v, , w] = uvtw;
    const vec = [
        1.5 * a * (u + v),
        0.5 * Math.sqrt(3) * a * (u - v),
        c * w
    ];
    return normalizeVec ? normalize3(vec) : vec;
}

export function crystalToMillerBravaisDirection(
    vecCryst,
    a = HEMATITE_CONSTANTS.a,
    c = HEMATITE_CONSTANTS.c
) {
    const [xa, xm, xc] = vecCryst;
    const uPlusV = (2 * xa) / (3 * a);
    const uMinusV = (2 * xm) / (Math.sqrt(3) * a);
    const u = 0.5 * (uPlusV + uMinusV);
    const v = 0.5 * (uPlusV - uMinusV);
    const t = -(u + v);
    const w = xc / c;
    return [u, v, t, w];
}

export function millerBravaisPlaneToCrystalNormal(
    hkil,
    a = HEMATITE_CONSTANTS.a,
    c = HEMATITE_CONSTANTS.c,
    normalizeVec = true
) {
    const [h, k, , l] = hkil;
    const g = [
        (h + k) / a,
        (h - k) / (Math.sqrt(3) * a),
        l / c
    ];
    return normalizeVec ? normalize3(g) : g;
}

export function crystalNormalToMillerBravaisPlane(
    vecCryst,
    a = HEMATITE_CONSTANTS.a,
    c = HEMATITE_CONSTANTS.c
) {
    const [ga, gm, gc] = vecCryst;
    const hPlusK = ga * a;
    const hMinusK = gm * Math.sqrt(3) * a;
    const h = 0.5 * (hPlusK + hMinusK);
    const k = 0.5 * (hPlusK - hMinusK);
    const i = -(h + k);
    const l = gc * c;
    return [h, k, i, l];
}

/**
 * Formats a 4-index Miller-Bravais array [u, v, t, w] into nice readable string
 * by scaling to small integers if close, or 2-decimal numbers.
 */
export function formatMillerBravais(indices, isPlane = false) {
    const { finalVals, useInts } = _scaleMillerBravais(indices);
    if (!finalVals) return isPlane ? '(0 0 0 0)' : '[0 0 0 0]';

    const parts = finalVals.map((val) => {
        if (useInts) {
            const r = Math.round(val);
            return r < 0 ? `${Math.abs(r)}\u0305` : `${r}`;
        }
        const rounded = Math.abs(val) < 0.005 ? 0 : val;
        return rounded.toFixed(2);
    });

    if (useInts && finalVals.every((v) => Math.abs(Math.round(v)) < 10)) {
        const compact = finalVals
            .map((val) => {
                const r = Math.round(val);
                return r < 0 ? `${Math.abs(r)}\u0305` : `${r}`;
            })
            .join('');
        return isPlane ? `(${compact})` : `[${compact}]`;
    }
    return isPlane ? `(${parts.join(' ')})` : `[${parts.join(' ')}]`;
}

/**
 * Formats a 4-index Miller-Bravais array [u, v, t, w] into a KaTeX LaTeX string
 * with properly centered \bar{n} overbars over negative indices.
 */
export function formatMillerBravaisLatex(indices, isPlane = false) {
    const { finalVals, useInts } = _scaleMillerBravais(indices);
    if (!finalVals) return isPlane ? '(0000)' : '[0000]';

    if (useInts && finalVals.every((v) => Math.abs(Math.round(v)) < 10)) {
        const compact = finalVals
            .map((val) => {
                const r = Math.round(val);
                return r < 0 ? `\\bar{${Math.abs(r)}}` : `${r}`;
            })
            .join('');
        return isPlane ? `(${compact})` : `[${compact}]`;
    }

    const parts = finalVals.map((val) => {
        if (useInts) {
            const r = Math.round(val);
            return r < 0 ? `\\overline{${Math.abs(r)}}` : `${r}`;
        }
        const rounded = Math.abs(val) < 0.005 ? 0 : val;
        return rounded < 0 ? `\\overline{${Math.abs(rounded).toFixed(1)}}` : rounded.toFixed(1);
    });
    return isPlane ? `(${parts.join('\\,')})` : `[${parts.join('\\,')}]`;
}

function _gcdInt(a, b) {
    let x = Math.abs(Math.round(a));
    let y = Math.abs(Math.round(b));
    while (y !== 0) {
        const t = y;
        y = x % y;
        x = t;
    }
    return x;
}

function _scaleMillerBravais(indices) {
    const maxAbs = Math.max(...indices.map(Math.abs));
    if (maxAbs < 1e-9) return { finalVals: null, useInts: true };
    const scaled = indices.map((x) => x / maxAbs);

    let bestInts = [0, 0, 0, 0];
    let bestScore = Infinity;

    for (let m = 1; m <= 6; m++) {
        let u = Math.round(scaled[0] * m);
        let v = Math.round(scaled[1] * m);
        let t = -(u + v);
        let w = Math.round(scaled[3] * m);

        if (u === 0 && v === 0 && w === 0) continue;

        const g = _gcdInt(_gcdInt(u, v), _gcdInt(t, w)) || 1;
        u = Math.round(u / g);
        v = Math.round(v / g);
        t = -(u + v);
        w = Math.round(w / g);

        const candLen = Math.hypot(u, v, t, w);
        const origLen = Math.hypot(scaled[0], scaled[1], scaled[2], scaled[3]);
        const dirErr =
            1 -
            Math.abs(
                (u * scaled[0] + v * scaled[1] + t * scaled[2] + w * scaled[3]) /
                    (candLen * origLen || 1)
            );
        // Prefer exact small integers; slightly penalize larger multipliers
        const score = dirErr + m * 0.002;
        if (score < bestScore) {
            bestScore = score;
            bestInts = [u, v, t, w];
        }
        if (dirErr < 0.0015) break;
    }

    return { finalVals: bestInts, useInts: true };
}

// ============================================================================
// 3. Cut Presets & Substrate-to-Crystal Orientation (F_Sub <-> F_Cryst)
// ============================================================================

export const CUT_PRESETS = {
    'c-cut': {
        id: 'c-cut',
        label: 'c-cut (0001)',
        shortLabel: 'c-cut',
        planeHKIL: [0, 0, 0, 1],
        thetaCutDeg: 0,
        phiCutDeg: 0,
        xsLabel: '[112\u03050] (e_a)',
        ysLabel: '[11\u030500] (e_m)',
        zsLabel: '(0001) || [0001] (e_c)',
        description: 'Basalebene (0001) liegt exakt in der Probenoberfläche (Easy-Plane = Probenebene, n̂ = ê_c).'
    },
    'a-cut': {
        id: 'a-cut',
        label: 'a-cut (112\u03050)',
        shortLabel: 'a-cut',
        planeHKIL: [1, 1, -2, 0],
        thetaCutDeg: 90,
        phiCutDeg: 0,
        xsLabel: '[11\u030500] (ê_m, in Easy-Plane)',
        ysLabel: '[0001] (ê_c, c-Achse)',
        zsLabel: '(112\u03050) || [112\u03050] (ê_a)',
        description: 'Probenoberfläche steht senkrecht zur Easy-Plane; c-Achse [0001] (ŷ_s) und [11\u030500] (x̂_s) liegen in der Probenebene.'
    },
    'm-cut': {
        id: 'm-cut',
        label: 'm-cut (101\u03050) / (11\u030500)',
        shortLabel: 'm-cut',
        planeHKIL: [1, -1, 0, 0],
        thetaCutDeg: 90,
        phiCutDeg: 90,
        xsLabel: '[1\u03051\u030520] (-ê_a, in Easy-Plane)',
        ysLabel: '[0001] (ê_c, c-Achse)',
        zsLabel: '(11\u030500) || [11\u030500] (ê_m)',
        description: 'Prismenebene senkrecht zur Easy-Plane; c-Achse [0001] (ŷ_s) und C2-Achse [1\u03051\u030520] (x̂_s) liegen in der Probenebene.'
    },
    'r-cut': {
        id: 'r-cut',
        label: 'r-cut (11\u030502)',
        shortLabel: 'r-cut',
        planeHKIL: [1, -1, 0, 2],
        get thetaCutDeg() {
            return HEMATITE_CONSTANTS.thetaRCutDeg;
        },
        phiCutDeg: 90,
        xsLabel: '[1\u03051\u030520] (-ê_a, in Easy-Plane)',
        ysLabel: '[1\u0305101] (proj. +c-Achse)',
        zsLabel: '(11\u030502) Normale (57.61° zu c)',
        description: 'Natürliche Spalt-/Wachstumsebene (11\u030502): Probennormale um θ_r = arctan(c/(√3 a)) ≈ 57.61° gegen die c-Achse verkippt.'
    },
    'custom': {
        id: 'custom',
        label: 'Benutzerdefiniert (Custom Cut)',
        shortLabel: 'Custom',
        planeHKIL: null,
        thetaCutDeg: 45,
        phiCutDeg: 0,
        xsLabel: 'x̂_s^(0) (Easy-Plane Schnittlinie)',
        ysLabel: 'ŷ_s^(0) (proj. +c-Richtung)',
        zsLabel: 'n̂(θ_cut, φ_cut)',
        description: 'Frei wählbarer Schnittwinkel (θ_cut, φ_cut) oder Miller-Bravais-Ebene (hkil).'
    }
};

/**
 * Computes the base rotation matrix R_{Cryst <- Sub}^{(0)}(thetaCutDeg, phiCutDeg)
 * whose columns are the unrotated substrate axes (x_s^(0), y_s^(0), z_s^(0))
 * expressed in F_Cryst = (e_a, e_m, e_c).
 *
 * Convention:
 *   - Surface normal in F_Cryst:
 *       z_s^(0) = n = (sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta))^T
 *   - For theta_cut == 0 (c-cut (0001)):
 *       x_s^(0) = (cos(phi), sin(phi), 0)^T  (at phi=0: e_a || [11-20])
 *       y_s^(0) = (-sin(phi), cos(phi), 0)^T (at phi=0: e_m || [1-100])
 *       z_s^(0) = (0, 0, 1)^T                (e_c || [0001])
 *   - For theta_cut == 180 ((000-1)):
 *       x_s^(0) = (cos(phi), sin(phi), 0)^T
 *       y_s^(0) = (sin(phi), -cos(phi), 0)^T
 *       z_s^(0) = (0, 0, -1)^T
 *   - For 0 < theta_cut < 180 (inclined or vertical cuts: a-cut, m-cut, r-cut, custom):
 *       x_s^(0) = (-sin(phi), cos(phi), 0)^T  (line of nodes: intersection of substrate & (0001) easy plane)
 *       y_s^(0) = (-cos(theta)*cos(phi), -cos(theta)*sin(phi), sin(theta))^T (in-plane projection of +e_c)
 *       z_s^(0) = (sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta))^T   (outward surface normal n)
 *
 * Notice that (x_s^(0), y_s^(0), z_s^(0)) is strictly orthonormal and right-handed (det = +1) everywhere!
 */
export function rCrystFromSubBase(thetaCutDeg, phiCutDeg) {
    const th = degToRad(thetaCutDeg);
    const ph = degToRad(phiCutDeg);
    const st = Math.sin(th);
    const ct = Math.cos(th);
    const sp = Math.sin(ph);
    const cp = Math.cos(ph);

    if (Math.abs(st) < 1e-10) {
        if (ct >= 0) {
            // theta_cut = 0 deg (c-cut)
            return matFromCols([cp, sp, 0], [-sp, cp, 0], [0, 0, 1]);
        } else {
            // theta_cut = 180 deg
            return matFromCols([cp, sp, 0], [sp, -cp, 0], [0, 0, -1]);
        }
    }

    const xs0 = [-sp, cp, 0];
    const ys0 = [-ct * cp, -ct * sp, st];
    const zs0 = [st * cp, st * sp, ct];
    return matFromCols(xs0, ys0, zs0);
}

/**
 * Full rotation matrix R_{Cryst <- Sub}(thetaCutDeg, phiCutDeg, alphaFlatDeg)
 * including the in-plane wafer-flat / edge reference offset alpha_flat.
 */
export function rCrystFromSub(thetaCutDeg, phiCutDeg, alphaFlatDeg = 0) {
    const R0 = rCrystFromSubBase(thetaCutDeg, phiCutDeg);
    const RFlat = rotZ(alphaFlatDeg);
    return matMul3(R0, RFlat);
}

// ============================================================================
// 4. Sample Holder & Wave Matrices (F_Lab <-> F_Mount <-> F_Sub <-> F_Wave)
// ============================================================================

/**
 * Mount-to-Lab rotation matrix R_{Lab <- Mount}(mountMode) at psi_S = 0.
 *   - 'IP' (horizontal in magnet plane):
 *       z_s = Z_L = (0,0,1)^T, x_s(0) = X_L = (1,0,0)^T, y_s(0) = Y_L = (0,1,0)^T
 *   - 'OOP' (vertical holder, normal along magnet zero axis X_L):
 *       z_s = X_L = (1,0,0)^T, x_s(0) = Y_L = (0,1,0)^T, y_s(0) = Z_L = (0,0,1)^T
 */
export function rLabFromMount(mountMode = 'OOP') {
    if (mountMode === 'IP') {
        return eye3();
    }
    // OOP configuration: columns in F_Lab are (Y_L, Z_L, X_L) -> det = +1
    return [
        [0, 0, 1],
        [1, 0, 0],
        [0, 1, 0]
    ];
}

/**
 * Substrate-to-Lab rotation matrix R_{Lab <- Sub}(mountMode, psiSDeg)
 * including sample rotation psi_S around its surface normal z_s.
 */
export function rLabFromSub(mountMode = 'OOP', psiSDeg = 0) {
    const RMount = rLabFromMount(mountMode);
    const RPsi = rotZ(psiSDeg);
    return matMul3(RMount, RPsi);
}

/**
 * Wave-to-Substrate rotation matrix R_{Sub <- Wave}(thetaDevDeg)
 * where k_hat = cos(theta_dev) x_s + sin(theta_dev) y_s,
 *       t_hat = -sin(theta_dev) x_s + cos(theta_dev) y_s,
 *       n_hat = z_s.
 */
export function rSubFromWave(thetaDevDeg = 0) {
    return rotZ(thetaDevDeg);
}

// ============================================================================
// 5. Complete 5-Frame Transformation Engine & Easy-Plane Equilibrium
// ============================================================================

export const FRAME_LIST = ['Lab', 'Mount', 'Sub', 'Wave', 'Cryst'];

export const FRAME_METADATA = {
    Lab: {
        id: 'Lab',
        name: 'Labor-System (ℱ_Lab)',
        shortName: 'ℱ_Lab',
        axesLabels: ['X̂_L (φ_M=0°)', 'Ŷ_L (φ_M=90°)', 'Ẑ_L (Vertikal)'],
        axesShort: ['X_L', 'Y_L', 'Z_L']
    },
    Mount: {
        id: 'Mount',
        name: 'Halterungs-System (ℱ_Mount)',
        shortName: 'ℱ_Mount',
        axesLabels: ['x̂_m (ψ_S=0°)', 'ŷ_m', 'ẑ_m = n̂'],
        axesShort: ['x_m', 'y_m', 'z_m']
    },
    Sub: {
        id: 'Sub',
        name: 'Substrat-System (ℱ_Sub)',
        shortName: 'ℱ_Sub',
        axesLabels: ['x̂_s (Probenkante)', 'ŷ_s (Probenkante)', 'ẑ_s = n̂ (Normale)'],
        axesShort: ['x_s', 'y_s', 'z_s']
    },
    Wave: {
        id: 'Wave',
        name: 'Wellen-System (ℱ_Wave)',
        shortName: 'ℱ_Wave',
        axesLabels: ['k̂ (Ausbreitung)', 't̂ = n̂×k̂ (Finger/Leiter)', 'n̂ (Normale)'],
        axesShort: ['k', 't', 'n']
    },
    Cryst: {
        id: 'Cryst',
        name: 'Kristall-System (ℱ_Cryst)',
        shortName: 'ℱ_Cryst',
        axesLabels: ['ê_a ∥ [112\u03050]', 'ê_m ∥ [11\u030500]', 'ê_c ∥ [0001]'],
        axesShort: ['e_a', 'e_m', 'e_c']
    }
};

export function createDefaultSetup(overrides = {}) {
    return {
        name: 'c-cut (0001)',
        mountMode: 'IP', // 'OOP' | 'IP'
        cutPreset: 'c-cut', // 'c-cut' | 'a-cut' | 'm-cut' | 'r-cut' | 'custom'
        thetaCutDeg: 0,
        phiCutDeg: 0,
        alphaFlatDeg: 0,
        deviceType: 'CPW', // 'CPW' (Spin Waves) | 'IDT' (SAW)
        thetaDevDeg: 0,
        psiSDeg: 0,
        phiMDeg: 0,
        b0Tesla: 0.5,
        aLattice: HEMATITE_CONSTANTS.a,
        cLattice: HEMATITE_CONSTANTS.c,
        ...overrides
    };
}

/**
 * Resolves effective (thetaCutDeg, phiCutDeg) from a setup object.
 */
export function resolveCutAngles(setup) {
    const preset = CUT_PRESETS[setup.cutPreset];
    if (preset && setup.cutPreset !== 'custom') {
        if (setup.cutPreset === 'r-cut') {
            const a = setup.aLattice || HEMATITE_CONSTANTS.a;
            const c = setup.cLattice || HEMATITE_CONSTANTS.c;
            return {
                thetaCutDeg: Math.atan(c / (Math.sqrt(3) * a)) * RAD2DEG,
                phiCutDeg: preset.phiCutDeg
            };
        }
        return {
            thetaCutDeg: preset.thetaCutDeg,
            phiCutDeg: preset.phiCutDeg
        };
    }
    return {
        thetaCutDeg: Number(setup.thetaCutDeg ?? 0),
        phiCutDeg: Number(setup.phiCutDeg ?? 0)
    };
}

/**
 * Builds all 3x3 rotation matrices R_{To <- From} between all 5 frames.
 */
export function computeFrameMatrices(setup) {
    const { thetaCutDeg, phiCutDeg } = resolveCutAngles(setup);
    const alphaFlatDeg = Number(setup.alphaFlatDeg ?? 0);
    const psiSDeg = Number(setup.psiSDeg ?? 0);
    const thetaDevDeg = Number(setup.thetaDevDeg ?? 0);
    const mountMode = setup.mountMode || 'OOP';

    // Core links relative to F_Sub:
    const R_Lab_Mount = rLabFromMount(mountMode);
    const R_Mount_Sub = rotZ(psiSDeg);
    const R_Lab_Sub = matMul3(R_Lab_Mount, R_Mount_Sub);
    const R_Sub_Wave = rSubFromWave(thetaDevDeg);
    const R_Cryst_Sub = rCrystFromSub(thetaCutDeg, phiCutDeg, alphaFlatDeg);

    // Express every frame's basis in F_Sub (R_{Sub <- Frame}):
    const R_Sub_from = {
        Sub: eye3(),
        Mount: transpose3(R_Mount_Sub),
        Lab: transpose3(R_Lab_Sub),
        Wave: R_Sub_Wave,
        Cryst: transpose3(R_Cryst_Sub)
    };

    // Express F_Sub in every frame (R_{Frame <- Sub}):
    const R_to_from_Sub = {
        Sub: eye3(),
        Mount: R_Mount_Sub,
        Lab: R_Lab_Sub,
        Wave: transpose3(R_Sub_Wave),
        Cryst: R_Cryst_Sub
    };

    /**
     * Returns 3x3 rotation matrix R_{toFrame <- fromFrame}
     */
    function getRotationMatrix(toFrame, fromFrame) {
        if (toFrame === fromFrame) return eye3();
        return matMul3(R_to_from_Sub[toFrame], R_Sub_from[fromFrame]);
    }

    /**
     * Transforms a 3D vector from `fromFrame` into all 5 frames.
     */
    function vectorInAllFrames(vec, fromFrame) {
        const vSub = matVec3(R_Sub_from[fromFrame], vec);
        return {
            Lab: matVec3(R_to_from_Sub.Lab, vSub),
            Mount: matVec3(R_to_from_Sub.Mount, vSub),
            Sub: vSub,
            Wave: matVec3(R_to_from_Sub.Wave, vSub),
            Cryst: matVec3(R_to_from_Sub.Cryst, vSub)
        };
    }

    return {
        thetaCutDeg,
        phiCutDeg,
        R_Lab_Mount,
        R_Mount_Sub,
        R_Lab_Sub,
        R_Sub_Wave,
        R_Cryst_Sub,
        R_Sub_from,
        R_to_from_Sub,
        getRotationMatrix,
        vectorInAllFrames
    };
}

/**
 * Evaluates all geometric and magnetic quantities for a given setup state.
 */
export function evaluateSetup(setup) {
    const matrices = computeFrameMatrices(setup);
    const { thetaCutDeg, phiCutDeg, vectorInAllFrames } = matrices;

    const phiMDeg = Number(setup.phiMDeg ?? 0);
    const psiSDeg = Number(setup.psiSDeg ?? 0);
    const thetaDevDeg = Number(setup.thetaDevDeg ?? 0);
    const alphaFlatDeg = Number(setup.alphaFlatDeg ?? 0);
    const b0 = Number(setup.b0Tesla ?? 1.0);
    const a = Number(setup.aLattice ?? HEMATITE_CONSTANTS.a);
    const c = Number(setup.cLattice ?? HEMATITE_CONSTANTS.c);

    // 1. Unit vectors of key physical directions in all frames:
    const phiMRad = degToRad(phiMDeg);
    const bHatLab = [Math.cos(phiMRad), Math.sin(phiMRad), 0];
    const bVecLab = scale3(bHatLab, b0);

    const bHat = vectorInAllFrames(bHatLab, 'Lab');
    const B = vectorInAllFrames(bVecLab, 'Lab');

    const kHat = vectorInAllFrames([1, 0, 0], 'Wave');
    const tHat = vectorInAllFrames([0, 1, 0], 'Wave');
    const nHat = vectorInAllFrames([0, 0, 1], 'Sub');

    const xsHat = vectorInAllFrames([1, 0, 0], 'Sub');
    const ysHat = vectorInAllFrames([0, 1, 0], 'Sub');

    const eaHat = vectorInAllFrames([1, 0, 0], 'Cryst');
    const emHat = vectorInAllFrames([0, 1, 0], 'Cryst');
    const ecHat = vectorInAllFrames([0, 0, 1], 'Cryst');

    // Line of nodes (intersection of sample plane and (0001) easy plane):
    const nodesCross = cross3(ecHat.Sub, nHat.Sub);
    const nodesNorm = norm3(nodesCross);
    const nodesSub = nodesNorm > 1e-9 ? scale3(nodesCross, 1 / nodesNorm) : [1, 0, 0];
    const nodesHat = vectorInAllFrames(nodesSub, 'Sub');

    // 2. Sample-plane projection of B (perpendicular to nHat):
    const bParallelSubVec = [bHat.Sub[0], bHat.Sub[1], 0];
    const bPerpSampleRatio = bHat.Sub[2]; // B . nHat / B0
    const bInSampleRatio = Math.hypot(bHat.Sub[0], bHat.Sub[1]); // |B_||s| / B0
    const bSampleHatSub = bInSampleRatio > 1e-9
        ? [bHat.Sub[0] / bInSampleRatio, bHat.Sub[1] / bInSampleRatio, 0]
        : [1, 0, 0];
    const bSampleProjVec = vectorInAllFrames(scale3(bParallelSubVec, b0), 'Sub');
    const bSampleHat = vectorInAllFrames(bSampleHatSub, 'Sub');

    const phiBSubDeg = bInSampleRatio > 1e-9
        ? wrapDeg360(Math.atan2(bHat.Sub[1], bHat.Sub[0]) * RAD2DEG)
        : 0;
    const phiBWaveDeg = bInSampleRatio > 1e-9
        ? wrapDeg180(Math.atan2(bHat.Wave[1], bHat.Wave[0]) * RAD2DEG)
        : 0;
    const elevationSampleDeg = Math.asin(Math.max(-1, Math.min(1, bPerpSampleRatio))) * RAD2DEG;

    // 3. Basal Easy-Plane (0001) projection of B (perpendicular to e_c):
    const bCryst = bHat.Cryst;
    const bParallelCRatio = bCryst[2]; // B . e_c / B0
    const bEasyRatio = Math.hypot(bCryst[0], bCryst[1]); // B_perp_c / B0
    const bEasyVecCryst = [bCryst[0] * b0, bCryst[1] * b0, 0];
    const bEasyProjVec = vectorInAllFrames(bEasyVecCryst, 'Cryst');

    // 4. Easy-Plane Equilibrium: Weak FM Moment m_hat || B_easy and Neel vector l_hat = m_hat x e_c
    const isEasyPlaneSingular = bEasyRatio < 1e-7;
    const mHatCryst = !isEasyPlaneSingular
        ? [bCryst[0] / bEasyRatio, bCryst[1] / bEasyRatio, 0]
        : [1, 0, 0]; // fallback along e_a when B || c
    // l_hat = m_hat x e_c in F_Cryst: (m_a, m_m, 0) x (0, 0, 1) = (m_m, -m_a, 0)
    const lHatCryst = [mHatCryst[1], -mHatCryst[0], 0];

    const mHat = vectorInAllFrames(mHatCryst, 'Cryst');
    const lHat = vectorInAllFrames(lHatCryst, 'Cryst');

    // Azimuthal angles in the basal (0001) easy plane (measured from e_a || [11-20] toward e_m || [1-100]):
    const phiEasyMDeg = wrapDeg360(Math.atan2(mHatCryst[1], mHatCryst[0]) * RAD2DEG);
    const phiEasyLDeg = wrapDeg360(Math.atan2(lHatCryst[1], lHatCryst[0]) * RAD2DEG);

    // 5. Key Physical Angles & Coupling Projections:
    const angleBKDeg = angleBetweenDeg(bHat.Sub, kHat.Sub);
    const angleBNDeg = angleBetweenDeg(bHat.Sub, nHat.Sub);
    const angleBCDeg = angleBetweenDeg(bHat.Cryst, [0, 0, 1]);
    const angleMKDeg = angleBetweenDeg(mHat.Sub, kHat.Sub);
    const angleLKDeg = angleBetweenDeg(lHat.Sub, kHat.Sub);
    const angleKCDeg = angleBetweenDeg(kHat.Cryst, [0, 0, 1]);
    const angleNCDeg = angleBetweenDeg(nHat.Cryst, [0, 0, 1]);

    // Signed in-plane angles of m_hat and l_hat projections on the sample relative to k_hat:
    const mInSampleNorm = Math.hypot(mHat.Wave[0], mHat.Wave[1]);
    const lInSampleNorm = Math.hypot(lHat.Wave[0], lHat.Wave[1]);
    const phiMWaveInPlaneDeg = mInSampleNorm > 1e-9
        ? wrapDeg180(Math.atan2(mHat.Wave[1], mHat.Wave[0]) * RAD2DEG)
        : 0;
    const phiLWaveInPlaneDeg = lInSampleNorm > 1e-9
        ? wrapDeg180(Math.atan2(lHat.Wave[1], lHat.Wave[0]) * RAD2DEG)
        : 0;

    // CPW / SAW Coupling Geometry Indicators:
    // For a CPW with signal line along t_hat, the in-plane RF Oersted field is along k_hat (h_rf || k_hat).
    // The RF torque on m_hat from h_rf || k_hat has magnitude |k_hat x m_hat| = sin(angle(m, k)),
    // and on l_hat has |k_hat x l_hat| = sin(angle(l, k)).
    const cpwTorqueM = norm3(cross3(kHat.Sub, mHat.Sub));
    const cpwTorqueL = norm3(cross3(kHat.Sub, lHat.Sub));

    // Miller-Bravais indices of key sample/wave directions:
    const mbXs = crystalToMillerBravaisDirection(xsHat.Cryst, a, c);
    const mbYs = crystalToMillerBravaisDirection(ysHat.Cryst, a, c);
    const mbZsDir = crystalToMillerBravaisDirection(nHat.Cryst, a, c);
    const mbZsPlane = crystalNormalToMillerBravaisPlane(nHat.Cryst, a, c);
    const mbK = crystalToMillerBravaisDirection(kHat.Cryst, a, c);
    const mbT = crystalToMillerBravaisDirection(tHat.Cryst, a, c);
    const mbB = crystalToMillerBravaisDirection(bHat.Cryst, a, c);
    const mbM = crystalToMillerBravaisDirection(mHat.Cryst, a, c);
    const mbL = crystalToMillerBravaisDirection(lHat.Cryst, a, c);

    return {
        setup: {
            ...setup,
            thetaCutDeg,
            phiCutDeg,
            alphaFlatDeg,
            psiSDeg,
            phiMDeg,
            thetaDevDeg,
            b0Tesla: b0,
            aLattice: a,
            cLattice: c
        },
        matrices,
        vectors: {
            B,
            bHat,
            bSampleProjVec,
            bSampleHat,
            bEasyProjVec,
            mHat,
            lHat,
            kHat,
            tHat,
            nHat,
            xsHat,
            ysHat,
            eaHat,
            emHat,
            ecHat,
            nodesHat
        },
        metrics: {
            b0Tesla: b0,
            bInSampleTesla: bInSampleRatio * b0,
            bPerpSampleTesla: bPerpSampleRatio * b0,
            bInSampleRatio,
            bPerpSampleRatio,
            bEasyTesla: bEasyRatio * b0,
            bParallelCTesla: bParallelCRatio * b0,
            bEasyRatio,
            bParallelCRatio,
            isEasyPlaneSingular,
            phiBSubDeg,
            phiBWaveDeg,
            elevationSampleDeg,
            phiEasyMDeg,
            phiEasyLDeg,
            angleBKDeg,
            angleBNDeg,
            angleBCDeg,
            angleMKDeg,
            angleLKDeg,
            angleKCDeg,
            angleNCDeg,
            phiMWaveInPlaneDeg,
            phiLWaveInPlaneDeg,
            bDotK: dot3(bHat.Sub, kHat.Sub),
            bDotT: dot3(bHat.Sub, tHat.Sub),
            bDotN: dot3(bHat.Sub, nHat.Sub),
            mDotK: dot3(mHat.Sub, kHat.Sub),
            mDotT: dot3(mHat.Sub, tHat.Sub),
            mDotN: dot3(mHat.Sub, nHat.Sub),
            lDotK: dot3(lHat.Sub, kHat.Sub),
            lDotT: dot3(lHat.Sub, tHat.Sub),
            lDotN: dot3(lHat.Sub, nHat.Sub),
            cpwTorqueM,
            cpwTorqueL
        },
        millerBravais: {
            xs: { raw: mbXs, label: formatMillerBravais(mbXs, false), latex: formatMillerBravaisLatex(mbXs, false) },
            ys: { raw: mbYs, label: formatMillerBravais(mbYs, false), latex: formatMillerBravaisLatex(mbYs, false) },
            zsDir: { raw: mbZsDir, label: formatMillerBravais(mbZsDir, false), latex: formatMillerBravaisLatex(mbZsDir, false) },
            zsPlane: { raw: mbZsPlane, label: formatMillerBravais(mbZsPlane, true), latex: formatMillerBravaisLatex(mbZsPlane, true) },
            k: { raw: mbK, label: formatMillerBravais(mbK, false), latex: formatMillerBravaisLatex(mbK, false) },
            t: { raw: mbT, label: formatMillerBravais(mbT, false), latex: formatMillerBravaisLatex(mbT, false) },
            B: { raw: mbB, label: formatMillerBravais(mbB, false), latex: formatMillerBravaisLatex(mbB, false) },
            m: { raw: mbM, label: formatMillerBravais(mbM, false), latex: formatMillerBravaisLatex(mbM, false) },
            l: { raw: mbL, label: formatMillerBravais(mbL, false), latex: formatMillerBravaisLatex(mbL, false) }
        }
    };
}

/**
 * Computes the Miller-Bravais direction [uvtw] (both plain label and KaTeX latex)
 * for any in-plane device angle thetaDevDeg on the sample surface.
 */
export function computeSampleAngleMillerBravais(setup, angleDeg) {
    const { thetaCutDeg, phiCutDeg } = resolveCutAngles(setup);
    const alphaFlatDeg = Number(setup.alphaFlatDeg ?? 0);
    const a = Number(setup.aLattice || HEMATITE_CONSTANTS.a);
    const c = Number(setup.cLattice || HEMATITE_CONSTANTS.c);
    const R_Cryst_Sub = rCrystFromSub(thetaCutDeg, phiCutDeg, alphaFlatDeg);
    const rad = degToRad(angleDeg);
    const vSub = [Math.cos(rad), Math.sin(rad), 0];
    const vCryst = matVec3(R_Cryst_Sub, vSub);
    const raw = crystalToMillerBravaisDirection(vCryst, a, c);
    return {
        raw,
        label: formatMillerBravais(raw, false),
        latex: formatMillerBravaisLatex(raw, false)
    };
}

// ============================================================================
// 6. Magnet Rotation Sweep Engine (0° ... 360° phi_M)
// ============================================================================

export function computeMagnetSweep(setup, stepDeg = 2) {
    const phiMArray = [];
    const bEasyRatio = [];
    const bInSampleRatio = [];
    const angleBKDeg = [];
    const angleMKDeg = [];
    const angleLKDeg = [];
    const bDotK = [];
    const mDotK = [];
    const lDotK = [];
    const cpwTorqueM = [];
    const cpwTorqueL = [];

    for (let phi = 0; phi <= 360 + 1e-9; phi += stepDeg) {
        const st = evaluateSetup({ ...setup, phiMDeg: phi });
        phiMArray.push(phi);
        bEasyRatio.push(st.metrics.bEasyRatio);
        bInSampleRatio.push(st.metrics.bInSampleRatio);
        angleBKDeg.push(st.metrics.angleBKDeg);
        angleMKDeg.push(st.metrics.angleMKDeg);
        angleLKDeg.push(st.metrics.angleLKDeg);
        bDotK.push(st.metrics.bDotK);
        mDotK.push(st.metrics.mDotK);
        lDotK.push(st.metrics.lDotK);
        cpwTorqueM.push(st.metrics.cpwTorqueM);
        cpwTorqueL.push(st.metrics.cpwTorqueL);
    }

    return {
        phiMArray,
        bEasyRatio,
        bInSampleRatio,
        angleBKDeg,
        angleMKDeg,
        angleLKDeg,
        bDotK,
        mDotK,
        lDotK,
        cpwTorqueM,
        cpwTorqueL
    };
}

// ============================================================================
// 7. Inverse Target Angle Finder (Exact Analytical 2-DOF & 1-DOF Solvers)
// ============================================================================

export const TARGET_PRESETS = [
    {
        id: 'B_parallel_k',
        label: 'B ∥ +k̂ (Parallel zur Wellenausbreitung, BVW / longitudinal)',
        shortLabel: 'B ∥ +k̂',
        category: 'Wave'
    },
    {
        id: 'B_antiparallel_k',
        label: 'B ∥ −k̂ (Antiparallel zur Wellenausbreitung)',
        shortLabel: 'B ∥ −k̂',
        category: 'Wave'
    },
    {
        id: 'B_parallel_t',
        label: 'B ∥ +t̂ (B ⊥ k̂ in Probenebene, DE / transversal parallel zu Fingern/CPW)',
        shortLabel: 'B ⊥ k̂ (in-plane)',
        category: 'Wave'
    },
    {
        id: 'B_parallel_n',
        label: 'B ∥ +n̂ (Perfekt Out-of-Plane entlang Probennormale)',
        shortLabel: 'B ∥ +n̂ (OOP)',
        category: 'Sample'
    },
    {
        id: 'B_easy_and_sample',
        label: 'B in Easy-Plane (0001) UND in Probenebene (Schnittlinie ê_c × n̂)',
        shortLabel: 'B ∈ Easy ∩ Sample',
        category: 'EasyPlane'
    },
    {
        id: 'm_align_k',
        label: 'B in Easy-Plane (0001) mit m̂ maximal parallel zu k̂ (l̂ ⊥ k̂_proj)',
        shortLabel: 'm̂ ∥ k̂_easy (B ∈ Easy)',
        category: 'EasyPlane'
    },
    {
        id: 'l_align_k',
        label: 'B in Easy-Plane (0001) mit Néel-Vektor l̂ maximal parallel zu k̂',
        shortLabel: 'l̂ ∥ k̂_easy (B ∈ Easy)',
        category: 'EasyPlane'
    },
    {
        id: 'B_parallel_ea',
        label: 'B ∥ +ê_a ([112\u03050] C2-Kristallachse in Easy-Plane)',
        shortLabel: 'B ∥ [112\u03050]',
        category: 'Crystal'
    },
    {
        id: 'B_parallel_em',
        label: 'B ∥ +ê_m ([11\u030500] Kristallachse in Easy-Plane)',
        shortLabel: 'B ∥ [11\u030500]',
        category: 'Crystal'
    },
    {
        id: 'B_parallel_ec',
        label: 'B ∥ +ê_c ([0001] harte c-Achse, B_⊥c = 0)',
        shortLabel: 'B ∥ [0001] (c)',
        category: 'Crystal'
    },
    {
        id: 'custom_mb',
        label: 'Benutzerdefinierte Kristallrichtung [u v t w]',
        shortLabel: 'Custom [uvtw]',
        category: 'Custom'
    }
];

/**
 * Resolves the target unit vector u_Sub in F_Sub for a given target preset or custom [uvtw].
 */
export function resolveTargetVectorSub(setup, targetId, customUVTW = [1, 1, -2, 0]) {
    const state = evaluateSetup(setup);
    const { kHat, tHat, nHat, eaHat, emHat, ecHat, nodesHat } = state.vectors;
    const { R_Sub_from } = state.matrices;

    switch (targetId) {
        case 'B_parallel_k':
            return { uSub: kHat.Sub, description: 'B ∥ +k̂ (entlang Wellenvektor)' };
        case 'B_antiparallel_k':
            return { uSub: scale3(kHat.Sub, -1), description: 'B ∥ −k̂ (entgegengesetzt zu k̂)' };
        case 'B_parallel_t':
            return { uSub: tHat.Sub, description: 'B ∥ +t̂ (transversal in Probenebene, ⊥ k̂)' };
        case 'B_parallel_n':
            return { uSub: nHat.Sub, description: 'B ∥ +n̂ (entlang Probennormale)' };
        case 'B_easy_and_sample':
            return {
                uSub: nodesHat.Sub,
                description: 'B entlang Schnittlinie aus Easy-Plane (0001) und Probenoberfläche'
            };
        case 'm_align_k': {
            // Project k_hat onto the basal (0001) easy plane:
            const kCryst = kHat.Cryst;
            const kEasyCryst = normalize3([kCryst[0], kCryst[1], 0], [1, 0, 0]);
            const uSub = matVec3(R_Sub_from.Cryst, kEasyCryst);
            return {
                uSub,
                description: 'B in Easy-Plane (0001), sodass m̂ entlang der Easy-Plane-Projektion von k̂ zeigt'
            };
        }
        case 'l_align_k': {
            // We want l_hat = m_hat x e_c || k_easy => m_hat = e_c x k_easy = (-k_m, k_a, 0) in F_Cryst:
            const kCryst = kHat.Cryst;
            const kEasyCryst = normalize3([kCryst[0], kCryst[1], 0], [1, 0, 0]);
            const mTargetCryst = [-kEasyCryst[1], kEasyCryst[0], 0];
            const uSub = matVec3(R_Sub_from.Cryst, mTargetCryst);
            return {
                uSub,
                description: 'B in Easy-Plane (0001), sodass der Néel-Vektor l̂ entlang der Projektion von k̂ zeigt'
            };
        }
        case 'B_parallel_ea':
            return { uSub: eaHat.Sub, description: 'B ∥ [112\u03050] (ê_a)' };
        case 'B_parallel_em':
            return { uSub: emHat.Sub, description: 'B ∥ [11\u030500] (ê_m)' };
        case 'B_parallel_ec':
            return { uSub: ecHat.Sub, description: 'B ∥ [0001] (c-Achse ê_c)' };
        case 'custom_mb': {
            const vCryst = millerBravaisDirectionToCrystal(
                customUVTW,
                setup.aLattice || HEMATITE_CONSTANTS.a,
                setup.cLattice || HEMATITE_CONSTANTS.c,
                true
            );
            const uSub = matVec3(R_Sub_from.Cryst, vCryst);
            return {
                uSub,
                description: `B ∥ ${formatMillerBravais(customUVTW, false)}`
            };
        }
        default:
            return { uSub: kHat.Sub, description: 'B ∥ +k̂' };
    }
}

/**
 * Computes the key magnet angles phi_M during a live magnet sweep at FIXED sample mounting psi_S.
 * Returns exact phi_M values for:
 *   - B closest to +k_hat (longitudinal / BVW)
 *   - B closest to +t_hat (transverse / DE, B perp k_hat)
 *   - B 100% in the (0001) Easy-Plane (B . e_c = 0)
 *   - B 100% in the Sample Plane (B . n_hat = 0)
 *   - Weak FM moment m_hat aligned with +k_easy in the (0001) plane
 *   - Neel vector l_hat aligned with +k_easy in the (0001) plane
 */
export function solveSweepLandmarksForFixedPsi(setup) {
    const st = evaluateSetup(setup);
    const { R_Lab_Sub, R_Sub_from } = st.matrices;
    const toLab = (vSub) => matVec3(R_Lab_Sub, vSub);

    // Helper: maximize b_Lab(phi_M) . vLab => phi_M = atan2(vY, vX)
    const phiMaxDot = (vLab) => {
        const h = Math.hypot(vLab[0], vLab[1]);
        const phi = h > 1e-9 ? wrapDeg360(Math.atan2(vLab[1], vLab[0]) * RAD2DEG) : 0;
        const errDeg = Math.acos(Math.max(-1, Math.min(1, h))) * RAD2DEG;
        return { phiMDeg: phi, oppositePhiMDeg: wrapDeg360(phi + 180), errorDeg: errDeg };
    };

    // Helper: find phi_M where b_Easy(phi_M) points along target direction uEasyCryst = [ua, um, 0] in F_Cryst
    const phiAlignEasyDir = (uEasyCryst) => {
        const uNorm = normalize3([uEasyCryst[0], uEasyCryst[1], 0], [1, 0, 0]);
        // Orthogonal vector in (0001) basal plane:
        const vPerpCryst = [-uNorm[1], uNorm[0], 0];
        const vPerpLab = toLab(matVec3(R_Sub_from.Cryst, vPerpCryst));
        const uEasyLab = toLab(matVec3(R_Sub_from.Cryst, uNorm));

        const hPerp = Math.hypot(vPerpLab[0], vPerpLab[1]);
        if (hPerp < 1e-9) {
            return phiMaxDot(uEasyLab);
        }
        // b_Lab(phi) . vPerpLab = 0 => cos(phi)*vX + sin(phi)*vY = 0
        const cand1 = wrapDeg360(Math.atan2(-vPerpLab[0], vPerpLab[1]) * RAD2DEG);
        const cand2 = wrapDeg360(cand1 + 180);
        // Pick the candidate that has positive dot product with uEasyLab:
        const r1 = degToRad(cand1);
        const dot1 = Math.cos(r1) * uEasyLab[0] + Math.sin(r1) * uEasyLab[1];
        const bestPhi = dot1 >= 0 ? cand1 : cand2;
        const oppPhi = dot1 >= 0 ? cand2 : cand1;
        const stBest = evaluateSetup({ ...setup, phiMDeg: bestPhi });
        const mBest = stBest.vectors.mHat.Cryst;
        const alignErr = angleBetweenDeg(mBest, uNorm);
        return {
            phiMDeg: bestPhi,
            oppositePhiMDeg: oppPhi,
            errorDeg: alignErr,
            bEasyRatioAtPhi: stBest.metrics.bEasyRatio
        };
    };

    const kLab = st.vectors.kHat.Lab;
    const tLab = st.vectors.tHat.Lab;
    const ecLab = st.vectors.ecHat.Lab;
    const nLab = st.vectors.nHat.Lab;

    // 1. B || +k_hat and B || +t_hat
    const bParallelK = phiMaxDot(kLab);
    const bPerpK = phiMaxDot(tLab);

    // 2. B in (0001) Easy-Plane (b_Lab . ecLab = 0) -> choose primary branch in [0°, 180°)
    const ecHoriz = Math.hypot(ecLab[0], ecLab[1]);
    let phiEasy1 = ecHoriz > 1e-9 ? wrapDeg360(Math.atan2(-ecLab[0], ecLab[1]) * RAD2DEG) : 0;
    if (phiEasy1 >= 180) phiEasy1 -= 180;
    const phiEasy2 = wrapDeg360(phiEasy1 + 180);

    // 3. B in Sample Plane (b_Lab . nLab = 0) -> choose primary branch in [0°, 180°)
    const nHoriz = Math.hypot(nLab[0], nLab[1]);
    let phiSample1 = nHoriz > 1e-9 ? wrapDeg360(Math.atan2(-nLab[0], nLab[1]) * RAD2DEG) : 90;
    if (phiSample1 >= 180) phiSample1 -= 180;
    const phiSample2 = wrapDeg360(phiSample1 + 180);

    // 4. m_hat || k_easy and l_hat || k_easy in (0001) Easy-Plane
    const kCryst = st.vectors.kHat.Cryst;
    const kEasyCryst = normalize3([kCryst[0], kCryst[1], 0], [1, 0, 0]);
    const mAlignK = phiAlignEasyDir(kEasyCryst);
    // For l_hat = m_hat x e_c || k_easy, we need m_hat || e_c x k_easy = [-k_m, k_a, 0]
    const lAlignK = phiAlignEasyDir([-kEasyCryst[1], kEasyCryst[0], 0]);

    // What is on the sample at phi_M = 0° and phi_M = 90°?
    const stPhi0 = evaluateSetup({ ...setup, phiMDeg: 0 });
    const stPhi90 = evaluateSetup({ ...setup, phiMDeg: 90 });

    return {
        bParallelK,
        bPerpK,
        bInEasyPlane: { phiM1: phiEasy1, phiM2: phiEasy2 },
        bInSamplePlane: { phiM1: phiSample1, phiM2: phiSample2 },
        mAlignK,
        lAlignK,
        atPhi0: {
            phiBSubDeg: stPhi0.metrics.phiBSubDeg,
            phiBWaveDeg: stPhi0.metrics.phiBWaveDeg,
            bInSampleRatio: stPhi0.metrics.bInSampleRatio,
            bPerpSampleRatio: stPhi0.metrics.bPerpSampleRatio
        },
        atPhi90: {
            phiBSubDeg: stPhi90.metrics.phiBSubDeg,
            phiBWaveDeg: stPhi90.metrics.phiBWaveDeg,
            bInSampleRatio: stPhi90.metrics.bInSampleRatio,
            bPerpSampleRatio: stPhi90.metrics.bPerpSampleRatio
        }
    };
}

/**
 * Solves for (psi_S, phi_M) to align B with target unit vector uSub in F_Sub.
 * Returns:
 *   - oopSolutions: 2 exact (psi_S, phi_M) pairs for OOP mount (0.00° error for ANY 3D vector!)
 *   - ipSolutions:  exact or closest-in-plane (psi_S, phi_M) pairs for IP mount
 *   - fixedPsiSweep: for the current setup (mountMode & current psi_S), finds the optimal phi_M
 *                    that maximizes B . u_target AND the exact phi_M where B crosses the (0001) easy plane.
 */
export function solveTargetAngles(setup, targetId = 'B_parallel_k', customUVTW = [1, 1, -2, 0]) {
    const { uSub, description } = resolveTargetVectorSub(setup, targetId, customUVTW);
    const u = normalize3(uSub);
    const [ux, uy, uz] = u;
    const xyNorm = Math.hypot(ux, uy);

    // 1. Exact OOP Mount Solutions:
    //    b_Sub^OOP(phi_M, psi_S) = (sin(phi_M)*cos(psi_S), -sin(phi_M)*sin(psi_S), cos(phi_M))^T
    const phi1 = Math.acos(Math.max(-1, Math.min(1, uz))) * RAD2DEG;
    const psi1 = xyNorm > 1e-9 ? wrapDeg360(Math.atan2(-uy, ux) * RAD2DEG) : wrapDeg360(setup.psiSDeg || 0);
    const phi2 = wrapDeg360(360 - phi1);
    const psi2 = wrapDeg360(psi1 + 180);

    const oopSolutions = [
        {
            branchLabel: 'OOP Hauptlösung (sin φ_M ≥ 0)',
            mountMode: 'OOP',
            psiSDeg: psi1,
            phiMDeg: phi1,
            errorDeg: 0
        },
        {
            branchLabel: 'OOP Komplementärlösung (+180° Symmetrie)',
            mountMode: 'OOP',
            psiSDeg: psi2,
            phiMDeg: phi2,
            errorDeg: 0
        }
    ];

    // 2. IP Mount Solutions:
    //    b_Sub^IP(phi_M, psi_S) = (cos(phi_M - psi_S), sin(phi_M - psi_S), 0)^T
    //    Closest in-plane unit vector is (ux, uy, 0) / xyNorm, residual out-of-plane error = asin(|uz|).
    const ipErrorDeg = Math.acos(Math.max(-1, Math.min(1, xyNorm))) * RAD2DEG;
    const targetAzimuthSubDeg = xyNorm > 1e-9 ? wrapDeg360(Math.atan2(uy, ux) * RAD2DEG) : 0;
    const currentPsi = wrapDeg360(setup.psiSDeg || 0);

    const ipSolutions = [
        {
            branchLabel: `IP Lösung (bei aktuellem ψ_S = ${currentPsi.toFixed(1)}°)`,
            mountMode: 'IP',
            psiSDeg: currentPsi,
            phiMDeg: wrapDeg360(currentPsi + targetAzimuthSubDeg),
            errorDeg: ipErrorDeg
        },
        {
            branchLabel: 'IP Lösung (Magnet bei φ_M = 0°)',
            mountMode: 'IP',
            psiSDeg: wrapDeg360(-targetAzimuthSubDeg),
            phiMDeg: 0,
            errorDeg: ipErrorDeg
        }
    ];

    // 3. Single-Axis Magnet Sweep (Fixed current mountMode & current psi_S):
    //    Transform uSub and e_c into F_Lab for the current mountMode and psi_S:
    const stCurrent = evaluateSetup(setup);
    const uLab = matVec3(stCurrent.matrices.R_Lab_Sub, u);
    const uLabHoriz = Math.hypot(uLab[0], uLab[1]);
    const bestPhiMForCurrentPsi = uLabHoriz > 1e-9
        ? wrapDeg360(Math.atan2(uLab[1], uLab[0]) * RAD2DEG)
        : 0;
    const fixedPsiErrorDeg = Math.acos(Math.max(-1, Math.min(1, uLabHoriz))) * RAD2DEG;

    // Also compute the exact phi_M where B crosses the (0001) easy plane (B . e_c = 0) for current psi_S:
    const ecLab = stCurrent.vectors.ecHat.Lab;
    const ecHoriz = Math.hypot(ecLab[0], ecLab[1]);
    const easyCrossingPhiM1 = ecHoriz > 1e-9
        ? wrapDeg360(Math.atan2(-ecLab[0], ecLab[1]) * RAD2DEG)
        : 0;
    const easyCrossingPhiM2 = wrapDeg360(easyCrossingPhiM1 + 180);

    return {
        targetId,
        description,
        uSub: u,
        oopSolutions,
        ipSolutions,
        fixedPsiSweep: {
            mountMode: setup.mountMode,
            psiSDeg: currentPsi,
            bestPhiMDeg: bestPhiMForCurrentPsi,
            oppositePhiMDeg: wrapDeg360(bestPhiMForCurrentPsi + 180),
            errorDeg: fixedPsiErrorDeg,
            easyCrossingPhiM1,
            easyCrossingPhiM2
        }
    };
}

// ============================================================================
// 7b. Bidirectional Angle Perspective Translators (Wave / Sub / Cryst -> Lab)
// ============================================================================

/**
 * Computes the manual sample mounting angle psi_S for intuitive experimental
 * mounting descriptions in both IP (horizontal) and OOP (vertical) holders.
 */
export function computeMountingPresetPsiS(setup, presetId, customPhiMRefDeg = 0) {
    const st = evaluateSetup(setup);
    const thetaDev = Number(setup.thetaDevDeg || 0);
    const nodesSub = st.vectors.nodesHat.Sub;
    const phiNodesSubDeg = Math.atan2(nodesSub[1], nodesSub[0]) * RAD2DEG;

    switch (presetId) {
        case 'k_ref':
            // IP: at phi_M = 0°, B || +k_hat; OOP: at phi_M = 90°, B || +k_hat
            return wrapDeg360(-thetaDev);
        case 't_ref':
            // IP: at phi_M = 0°, B || +t_hat (perp k_hat); OOP: at phi_M = 90°, B || +t_hat
            return wrapDeg360(-(thetaDev + 90));
        case 'xs_ref':
            // +x_s at phi_M = 0° (IP) or phi_M = 90° (OOP)
            return 0;
        case 'ys_ref':
            // +y_s (proj. c-axis) at phi_M = 0° (IP) or phi_M = 90° (OOP)
            return 270;
        case 'nodes_ref':
            // Easy-Plane intersection line at phi_M = 0° (IP) or phi_M = 90° (OOP)
            return wrapDeg360(-phiNodesSubDeg);
        case 'calib_ip_k':
            // Calibration in IP mount: user specifies at which magnet angle phi_M^ref B is parallel to +k_hat
            return wrapDeg360(Number(customPhiMRefDeg || 0) - thetaDev);
        default:
            return wrapDeg360(setup.psiSDeg || 0);
    }
}

/**
 * Given any desired 3D unit vector uSub in F_Sub, returns the hardware angles
 * (psiSDeg, phiMDeg, mountMode) that realize this field direction.
 * If setup.mountMode === 'IP' and uSub has a non-zero z-component (OOP),
 * `autoSwitchOOP = true` switches to 'OOP' so that every 3D angle can be reached seamlessly.
 */
export function hardwareAnglesFromSubVector(setup, uSub, autoSwitchOOP = true) {
    const u = normalize3(uSub);
    const [ux, uy, uz] = u;
    const xyNorm = Math.hypot(ux, uy);

    const needOOP = Math.abs(uz) > 1e-4;
    const targetMount = needOOP && autoSwitchOOP ? 'OOP' : setup.mountMode;

    if (targetMount === 'OOP') {
        const phiM = Math.acos(Math.max(-1, Math.min(1, uz))) * RAD2DEG;
        const psiS = xyNorm > 1e-9
            ? wrapDeg360(Math.atan2(-uy, ux) * RAD2DEG)
            : wrapDeg360(setup.psiSDeg || 0);
        return { mountMode: 'OOP', psiSDeg: psiS, phiMDeg: phiM };
    } else {
        // IP mode: B_Sub = (cos(phi_M - psi_S), sin(phi_M - psi_S), 0)^T
        const azDeg = xyNorm > 1e-9 ? wrapDeg360(Math.atan2(uy, ux) * RAD2DEG) : 0;
        const psiS = wrapDeg360(setup.psiSDeg || 0);
        const phiM = wrapDeg360(psiS + azDeg);
        return { mountMode: 'IP', psiSDeg: psiS, phiMDeg: phiM };
    }
}

/**
 * Perspective B (Wave Frame F_Wave):
 * Sets B from in-plane wave azimuth phiBWaveDeg (angle to k_hat in sample plane)
 * and out-of-plane elevation elevationSampleDeg (angle from sample surface toward n_hat).
 */
export function hardwareAnglesFromWavePerspective(
    setup,
    phiBWaveDeg,
    elevationSampleDeg = 0,
    autoSwitchOOP = true
) {
    const elRad = degToRad(Math.max(-90, Math.min(90, elevationSampleDeg)));
    const azRad = degToRad(phiBWaveDeg);
    const cel = Math.cos(elRad);
    const sel = Math.sin(elRad);
    const uWave = [cel * Math.cos(azRad), cel * Math.sin(azRad), sel];
    const mats = computeFrameMatrices(setup);
    const uSub = matVec3(mats.R_Sub_from.Wave, uWave);
    return hardwareAnglesFromSubVector(setup, uSub, autoSwitchOOP);
}

/**
 * Perspective C (Substrate Frame F_Sub):
 * Sets B from in-plane substrate azimuth phiBSubDeg (angle to x_s)
 * and polar angle angleBNDeg = <(B, n_hat) = 90° - elevationSampleDeg.
 */
export function hardwareAnglesFromSubPerspective(
    setup,
    phiBSubDeg,
    angleBNDeg = 90,
    autoSwitchOOP = true
) {
    const thRad = degToRad(Math.max(0, Math.min(180, angleBNDeg)));
    const azRad = degToRad(phiBSubDeg);
    const st = Math.sin(thRad);
    const ct = Math.cos(thRad);
    const uSub = [st * Math.cos(azRad), st * Math.sin(azRad), ct];
    return hardwareAnglesFromSubVector(setup, uSub, autoSwitchOOP);
}

/**
 * Perspective D (Crystal / Easy-Plane Frame F_Cryst):
 * Sets B from basal (0001) easy-plane azimuth phiEasyMDeg (angle of m_hat / B_perp_c
 * relative to e_a || [11-20]) and polar angle angleBCDeg = <(B, e_c) (90° = B purely in Easy-Plane).
 */
export function hardwareAnglesFromCrystalPerspective(
    setup,
    phiEasyMDeg,
    angleBCDeg = 90,
    autoSwitchOOP = true
) {
    const thRad = degToRad(Math.max(0, Math.min(180, angleBCDeg)));
    const azRad = degToRad(phiEasyMDeg);
    const st = Math.sin(thRad);
    const ct = Math.cos(thRad);
    const uCryst = [st * Math.cos(azRad), st * Math.sin(azRad), ct];
    const mats = computeFrameMatrices(setup);
    const uSub = matVec3(mats.R_Sub_from.Cryst, uCryst);
    return hardwareAnglesFromSubVector(setup, uSub, autoSwitchOOP);
}

// ============================================================================
// 8. Python Snippet Generator for Evaluation Scripts
// ============================================================================

export function generatePythonSnippet(setupA, setupB = null, compareActive = false) {
    const { thetaCutDeg: thA, phiCutDeg: phA } = resolveCutAngles(setupA);
    const lines = [
        `# Automatisch generiertes Setup aus dem Hematite Crystal & Spin-Wave/SAW Calculator`,
        `import numpy as np`,
        `from hematite_coords import HematiteSetup`,
        ``,
        `setup_a = HematiteSetup(`,
        `    name=${JSON.stringify(setupA.name || 'Setup A')},`,
        `    mount_mode=${JSON.stringify(setupA.mountMode)},  # 'OOP' (vertikal) oder 'IP' (horizontal)`,
        `    cut_preset=${JSON.stringify(setupA.cutPreset)},`,
        `    theta_cut_deg=${Number(thA).toFixed(4)},`,
        `    phi_cut_deg=${Number(phA).toFixed(4)},`,
        `    alpha_flat_deg=${Number(setupA.alphaFlatDeg || 0).toFixed(4)},`,
        `    device_type=${JSON.stringify(setupA.deviceType || 'CPW')},`,
        `    theta_dev_deg=${Number(setupA.thetaDevDeg || 0).toFixed(4)},`,
        `    psi_s_deg=${Number(setupA.psiSDeg || 0).toFixed(4)},`,
        `    phi_m_deg=${Number(setupA.phiMDeg || 0).toFixed(4)},`,
        `    b0_tesla=${Number(setupA.b0Tesla || 0.5).toFixed(4)},`,
        `)`,
        ``,
        `state_a = setup_a.evaluate()`,
        `print("B in F_Sub (T):", np.round(state_a["vectors"]["B"]["Sub"], 4))`,
        `print("B in F_Wave (k, t, n) (T):", np.round(state_a["vectors"]["B"]["Wave"], 4))`,
        `print("B in F_Cryst (e_a, e_m, e_c) (T):", np.round(state_a["vectors"]["B"]["Cryst"], 4))`,
        `print(f"B_perp_c / B0 = {state_a['metrics']['b_easy_ratio']:.4f}, <(m, k) = {state_a['metrics']['angle_m_k_deg']:.2f} deg, <(l, k) = {state_a['metrics']['angle_l_k_deg']:.2f} deg")`,
        ``,
        `# 0..360 deg Magnet-Sweep fuer Auswerteskripte:`,
        `phi_m_deg = np.linspace(0, 360, 361)`,
        `sweep_a = setup_a.sweep_magnet(phi_m_deg)`
    ];

    if (compareActive && setupB) {
        const { thetaCutDeg: thB, phiCutDeg: phB } = resolveCutAngles(setupB);
        lines.push(
            ``,
            `setup_b = HematiteSetup(`,
            `    name=${JSON.stringify(setupB.name || 'Setup B')},`,
            `    mount_mode=${JSON.stringify(setupB.mountMode)},`,
            `    cut_preset=${JSON.stringify(setupB.cutPreset)},`,
            `    theta_cut_deg=${Number(thB).toFixed(4)},`,
            `    phi_cut_deg=${Number(phB).toFixed(4)},`,
            `    alpha_flat_deg=${Number(setupB.alphaFlatDeg || 0).toFixed(4)},`,
            `    device_type=${JSON.stringify(setupB.deviceType || 'IDT')},`,
            `    theta_dev_deg=${Number(setupB.thetaDevDeg || 0).toFixed(4)},`,
            `    psi_s_deg=${Number(setupB.psiSDeg || 0).toFixed(4)},`,
            `    phi_m_deg=${Number(setupB.phiMDeg || 0).toFixed(4)},`,
            `    b0_tesla=${Number(setupB.b0Tesla || 0.5).toFixed(4)},`,
            `)`,
            `state_b = setup_b.evaluate()`,
            `sweep_b = setup_b.sweep_magnet(phi_m_deg)`
        );
    }

    return lines.join('\n');
}
