// --- LATTICE CONSTANTS ---
// These are needed for the transformations.
const A_LATTICE = 1.0;
const C_LATTICE = 1.365 * 2; // c/a ratio for Hematite is ~2.73

/**
 * Converts Miller-Bravais direction indices [h, k, i, l] to a cartesian vector [x, y, z]
 * in the crystal's cartesian system.
 * This assumes a standard orientation where a1 is along x and c is along z.
 * @param {number[]} mb - Miller-Bravais direction indices [h, k, i, l].
 * @returns {number[]} Cartesian vector [x, y, z].
 */
export function millerBravaisToCartesian(mb) {
    const [h, k, _, l] = mb;
    // Use basis a1 = (1,0), a2 = (-1/2, sqrt(3)/2) which gives the common
    // conversion used in many crystallography references.
    // x = (3/2) * a * h
    // y = (sqrt(3)/2) * a * (h + 2*k)
    const x = A_LATTICE * (3 / 2) * h;
    const y = A_LATTICE * (Math.sqrt(3) / 2) * (h + 2 * k);
    const z = C_LATTICE * l;
    return [x, y, z];
}

/**
 * Converts a cartesian vector [x, y, z] in the crystal system to Miller-Bravais direction indices [h, k, i, l].
 * @param {number[]} cart - Cartesian vector [x, y, z].
 * @returns {number[]} Miller-Bravais direction indices [h, k, i, l].
 */
export function cartesianToMillerBravais(cart) {
    const [x, y, z] = cart;
    // inverse of the above conversion
    const h = (2 * x) / (3 * A_LATTICE);
    const k = ((2 * y) / (A_LATTICE * Math.sqrt(3)) - h) / 2;
    const l = z / C_LATTICE;
    const i = -(h + k);
    return [h, k, i, l];
}

/**
 * Converts Miller-Bravais plane indices (h, k, i, l) to a Miller-Bravais direction [u, v, t, w]
 * that is normal to the plane.
 * @param {number[]} plane - Miller-Bravais plane indices (h, k, i, l).
 * @returns {number[]} Miller-Bravais direction indices [u, v, t, w].
 */
export function millerBravaisPlaneToDirection(plane) {
    const [h, k, i, l] = plane;
    const c_a_ratio = C_LATTICE / A_LATTICE;
    // Formula for hexagonal systems to convert a plane (h,k,i,l) to a direction [u,v,t,w]
    // w = l * (3/2) * (a/c)^2 = l * (3/2) * (1 / (c/a)^2)
    const w = l * (3/2) * (1 / (c_a_ratio * c_a_ratio));
    return [h, k, i, w];
}

/**
 * Transforms a cartesian vector from the lab system to the crystal system.
 * @param {number[]} labVector - The vector [x, y, z] in the lab system.
 * @param {number[][]} transformationMatrix - The 3x3 matrix to transform from crystal to lab.
 * @returns {number[]} The cartesian vector [x, y, z] in the crystal system.
 */
export function labToCrystal(labVector, transformationMatrix) {
    // A general transformation matrix is not necessarily orthonormal.
    // Therefore, its inverse is not its transpose. A full inversion is required.
    const M_inv = invertMatrix(transformationMatrix);
    if (!M_inv) {
        console.error("Matrix is singular and cannot be inverted.");
        return [0, 0, 0];
    }
    return multiplyMatrixVector(M_inv, labVector);
}

/**
 * Transforms a cartesian vector from the crystal system to the lab system.
 * @param {number[]} crystalVector - The vector [x, y, z] in the crystal system.
 * @param {number[][]} transformationMatrix - The 3x3 matrix to transform from crystal to lab.
 * @returns {number[]} The cartesian vector [x, y, z] in the lab system.
 */
export function crystalToLab(crystalVector, transformationMatrix) {
    return multiplyMatrixVector(transformationMatrix, crystalVector);
}


// --- VECTOR AND MATRIX HELPERS ---

/**
 * Normalizes a 3D vector.
 * @param {number[]} v - The vector [x, y, z].
 * @returns {number[]} The normalized vector.
 */
export function normalize(v) {
    const len = Math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    if (len === 0) return [0, 0, 0];
    return [v[0] / len, v[1] / len, v[2] / len];
}

/**
 * Calculates the cross product of two 3D vectors.
 * @param {number[]} a - The first vector [x, y, z].
 * @param {number[]} b - The second vector [x, y, z].
 * @returns {number[]} The resulting vector.
 */
export function crossProduct(a, b) {
    return [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0]
    ];
}

/**
 * Multiplies a 3x3 matrix with a 3x1 vector.
 * @param {number[][]} m - The 3x3 matrix.
 * @param {number[]} v - The 3x1 vector.
 * @returns {number[]} The resulting 3x1 vector.
 */
function multiplyMatrixVector(m, v) {
    return [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2]
    ];
}

/**
 * Transposes a 3x3 matrix.
 * @param {number[][]} m - The 3x3 matrix.
 * @returns {number[][]} The transposed matrix.
 */
export function transposeMatrix(m) {
    return [
        [m[0][0], m[1][0], m[2][0]],
        [m[0][1], m[1][1], m[2][1]],
        [m[0][2], m[1][2], m[2][2]]
    ];
}

/**
 * Inverts a 3x3 matrix.
 * @param {number[][]} m - The 3x3 matrix.
 * @returns {number[][] | null} The inverted matrix, or null if the matrix is not invertible.
 */
function invertMatrix(m) {
    const m00 = m[0][0], m01 = m[0][1], m02 = m[0][2];
    const m10 = m[1][0], m11 = m[1][1], m12 = m[1][2];
    const m20 = m[2][0], m21 = m[2][1], m22 = m[2][2];

    const det = m00 * (m11 * m22 - m21 * m12) -
        m01 * (m10 * m22 - m12 * m20) +
        m02 * (m10 * m21 - m11 * m20);

    if (det === 0) {
        return null; // Matrix is not invertible
    }

    const invDet = 1.0 / det;

    return [
        [(m11 * m22 - m21 * m12) * invDet, (m02 * m21 - m01 * m22) * invDet, (m01 * m12 - m02 * m11) * invDet],
        [(m12 * m20 - m10 * m22) * invDet, (m00 * m22 - m02 * m20) * invDet, (m02 * m10 - m00 * m12) * invDet],
        [(m10 * m21 - m20 * m11) * invDet, (m20 * m01 - m00 * m21) * invDet, (m00 * m11 - m10 * m01) * invDet]
    ];
}