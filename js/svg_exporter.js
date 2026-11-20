/**
 * Figure export: one uniform schematic of the sample and the magnetic field.
 *
 * The sample always sits in the same pose and at the same scale (fixed camera in the sample frame
 * F_Sub), so figures made from different setups line up in a thesis. Drawn: slab, emitter (and
 * receiver) shaped as in the 3D view and turned by theta_dev, the wave vector between them, the
 * field arrow with its in-plane and out-of-plane angle, and a small triad with the Miller-Bravais
 * directions of the sample edges and the surface normal. Optionally a sweep: the field turning in
 * the magnet's rotation plane (path of the tip, faint steps), or B0 running over a range.
 * The figure is cropped tightly around what is drawn.
 *
 * Geometry lives in sample units (top surface at z = 0); the SVG viewBox is in pt, so font sizes
 * and line widths in the style are print sizes when the figure is placed at its natural size.
 */

const DEG = Math.PI / 180;
const PT_PER_MM = 72 / 25.4;

// Camera in F_Sub. Oblique: tilted about the front edge only, so that edge lies level, and the depth
// sheared to the right (cabinet projection) so the slab reads as 3D. Corner: low view onto the
// front-left corner; top: looking down the normal.
const VIEWS = {
    oblique: { azDeg: -90, elDeg: 34, shear: 0.32 },
    corner: { azDeg: -66, elDeg: 28 },
    top: { azDeg: -90, elDeg: 90 }
};

const HX = 1.5; // half edges of the sample: long along the 0° edge x_s
const HY = 1;
const FIELD_LEN = 0.86;
const DEV_DIST = 1; // emitter / receiver centres at -/+ DEV_DIST along k
const DEV_SCALE = 1; // microstructures as in the 3D view
const TRIAD_SCALE = 0.75; // direction triad, labels included, a bit smaller than the rest
const MIN_HALF_GAP = 0.36; // free half-distance between the antennas, room for k
const K_GAP = 0.06; // k starts this far from the emitter
const MOMENT_LEN = 0.8;
const R_DEV = 0.36; // arc radii: antenna angle, in-plane field angle, elevation
const R_PHI = 0.48;
const R_ELEV = 0.6;
const TINY_DEG = 0.05;
const LIGHT_DIR = [(-0.55) / Math.hypot(-0.55, -0.45, 0.7), (-0.45) / Math.hypot(-0.55, -0.45, 0.7), 0.7 / Math.hypot(-0.55, -0.45, 0.7)];

// ---------------------------------------------------------------------------
// Styles: everything that should look the same in every figure of a thesis
// ---------------------------------------------------------------------------

const SYMBOLS = { field: 'B', wave: 'k', moment: 'm', phi: 'φ', theta: 'θ', dev: 'α', sweep: '' };
const DEFAULT_EFFECTS = { extrude: false, angleFill: false, waveLine: false };

export const BUILTIN_FIGURE_STYLES = {
    Standard: {
        widthMm: 65,
        fontFamily: "'Helvetica Neue', Helvetica, Arial, sans-serif",
        fontSizePt: 9,
        vectorStyle: 'bold',
        sampleThickness: 0.1,
        arrowHeadPt: 6,
        colors: {
            sample: '#f8fafc',
            antenna: '#e6a825',
            field: '#2563a8',
            wave: '#c43728',
            angle: '#1e7a46',
            axes: '#1e293b',
            moment: '#7c3aed',
            text: '#111827',
            background: 'none'
        },
        lines: { outline: 0.72, field: 1.5, wave: 1.2, arc: 0.8, guide: 0.55, axes: 0.8 },
        effects: { ...DEFAULT_EFFECTS },
        symbols: { ...SYMBOLS }
    },
    'Paper (Serif)': {
        widthMm: 65,
        fontFamily: "Georgia, 'Times New Roman', Times, serif",
        fontSizePt: 9.5,
        vectorStyle: 'bold',
        sampleThickness: 0.1,
        arrowHeadPt: 6.2,
        colors: {
            sample: '#ffffff',
            antenna: '#e2a422',
            field: '#2563a8',
            wave: '#be3126',
            angle: '#1e7a46',
            axes: '#111111',
            moment: '#7c3aed',
            text: '#111111',
            background: 'none'
        },
        lines: { outline: 0.8, field: 1.5, wave: 1.2, arc: 0.85, guide: 0.55, axes: 0.85 },
        effects: { ...DEFAULT_EFFECTS },
        symbols: { ...SYMBOLS }
    },
    'Farbig (Flat)': {
        widthMm: 65,
        fontFamily: "'Helvetica Neue', Helvetica, Arial, sans-serif",
        fontSizePt: 9,
        vectorStyle: 'bold',
        sampleThickness: 0.11,
        arrowHeadPt: 6,
        colors: {
            sample: '#c6e0e7',
            antenna: '#eab12a',
            field: '#19325c',
            wave: '#d13b27',
            angle: '#334155',
            axes: '#334155',
            moment: '#7c3aed',
            text: '#0f172a',
            background: 'none'
        },
        lines: { outline: 0.55, field: 1.5, wave: 1.15, arc: 0.75, guide: 0.55, axes: 0.75 },
        effects: { extrude: true, angleFill: true, waveLine: false },
        symbols: { ...SYMBOLS }
    },
    'Schwarz-Weiß': {
        widthMm: 65,
        fontFamily: "Georgia, 'Times New Roman', Times, serif",
        fontSizePt: 9.5,
        vectorStyle: 'bold',
        sampleThickness: 0.1,
        arrowHeadPt: 6,
        colors: {
            sample: '#ffffff',
            antenna: '#9ca3af',
            field: '#111111',
            wave: '#111111',
            angle: '#374151',
            axes: '#111111',
            moment: '#52525b',
            text: '#000000',
            background: 'none'
        },
        lines: { outline: 0.75, field: 1.45, wave: 1.05, arc: 0.75, guide: 0.5, axes: 0.75 },
        effects: { ...DEFAULT_EFFECTS },
        symbols: { ...SYMBOLS }
    },
    Folie: {
        widthMm: 90,
        fontFamily: "'Helvetica Neue', Helvetica, Arial, sans-serif",
        fontSizePt: 14,
        vectorStyle: 'bold',
        sampleThickness: 0.11,
        arrowHeadPt: 10,
        colors: {
            sample: '#e2eef3',
            antenna: '#eca924',
            field: '#1d4ed8',
            wave: '#dc2626',
            angle: '#15803d',
            axes: '#1e293b',
            moment: '#7c3aed',
            text: '#0f172a',
            background: '#ffffff'
        },
        lines: { outline: 0.9, field: 2.5, wave: 2.0, arc: 1.2, guide: 0.9, axes: 1.2 },
        effects: { extrude: false, angleFill: true, waveLine: false },
        symbols: { ...SYMBOLS }
    }
};

export const DEFAULT_FIGURE_OPTIONS = {
    view: 'oblique', // 'oblique' (front edge level, depth sheared) | 'corner' | 'top'
    angleRef: 'wave', // in-plane field angle measured from k ('wave') or from the 0° edge ('edge')
    angleLabels: 'value', // 'value' | 'symbol' | 'both'
    showAntenna: true,
    devices: 'pair', // 'pair' | 'single' | 'IDT-CPW' | 'CPW-IDT' (emitter first)
    showWave: true,
    showAxes: true,
    showCAxis: true,
    showMoment: false,
    latex: false,
    sweep: 'none', // 'none' | 'angle' (magnet turns the field) | 'field' (B0 over a range)
    sweepFrom: '', // phi_M in ° or B0 in mT; empty: the current value (angle) or 0 (field)
    sweepTo: '', // B0 sweep end in mT; empty: the current B0
    sweepSpan: '', // angle sweep: turn in ° from the start, negative = other sense; empty: none yet
    sweepStep: '0', // faint field arrows every … ° along the turn; 0: none
    panel: '' // sub-figure letter in the top left corner, e.g. '(a)'; empty: none
};

/** Fills gaps of a (possibly older or partial) saved style from the standard style. */
export function completeFigureStyle(style) {
    const base = BUILTIN_FIGURE_STYLES.Standard;
    const s = style || {};
    return {
        ...base,
        ...s,
        colors: { ...base.colors, ...(s.colors || {}) },
        lines: { ...base.lines, ...(s.lines || {}) },
        effects: { ...base.effects, ...(s.effects || {}) },
        symbols: { ...base.symbols, ...(s.symbols || {}) }
    };
}

// ---------------------------------------------------------------------------
// Antennas: the same shapes as in the 3D view, in the wave frame [along k, along t]
// ---------------------------------------------------------------------------

const DEVICE_SETS = {
    pair: (own) => [own, own],
    single: (own) => [own],
    'IDT-CPW': () => ['IDT', 'CPW'],
    'CPW-IDT': () => ['CPW', 'IDT']
};

/** Interdigital transducer (two combs) or shorted coplanar waveguide; fingers / lines along t. */
function microstructure(type) {
    if (type === 'IDT') {
        const w = 0.042;
        const x0 = -5.5 * w;
        const yOuter = 0.28;
        const yBus = yOuter - 1.35 * w;
        const yTipTop = -yBus + w;
        const yTipBot = yBus - w;
        const comb = (y0, yB, yTip, fingers) => [
            [x0, y0],
            [-x0, y0],
            [-x0, yB],
            ...fingers.flatMap((i) => [[x0 + (i + 1) * w, yB], [x0 + (i + 1) * w, yTip], [x0 + i * w, yTip], [x0 + i * w, yB]]),
            [x0, yB]
        ].filter((p, i, all) => i === 0 || p[0] !== all[i - 1][0] || p[1] !== all[i - 1][1]);
        return [comb(yOuter, yBus, yTipTop, [8, 4, 0]), comb(-yOuter, -yBus, yTipBot, [10, 6, 2])];
    }
    const w = 0.056;
    const yBot = -0.35;
    const yTop = 0.35;
    const yIn = yTop - w;
    return [[
        [-2.5 * w, yTop], [2.5 * w, yTop], [2.5 * w, yBot], [1.5 * w, yBot], [1.5 * w, yIn], [0.5 * w, yIn],
        [0.5 * w, yBot], [-0.5 * w, yBot], [-0.5 * w, yIn], [-1.5 * w, yIn], [-1.5 * w, yBot], [-2.5 * w, yBot]
    ]];
}
const halfLength = (type) => (type === 'IDT' ? 5.5 * 0.042 : 2.5 * 0.056);

// ---------------------------------------------------------------------------
// Small vector + colour helpers
// ---------------------------------------------------------------------------

const add = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const mul = (a, s) => [a[0] * s, a[1] * s, a[2] * s];
const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const len = (a) => Math.hypot(a[0], a[1], a[2]);
const unit = (a) => mul(a, 1 / (len(a) || 1));
const clamp = (v, lo, hi) => Math.max(lo, Math.min(hi, v));
const angleDeg = (a, b) => Math.acos(clamp(dot(unit(a), unit(b)), -1, 1)) / DEG;

function shade(hex, f) {
    const m = /^#?([0-9a-f]{6})$/i.exec(hex || '');
    if (!m) return hex;
    const n = parseInt(m[1], 16);
    const ch = [(n >> 16) & 255, (n >> 8) & 255, n & 255].map((c) => Math.round(clamp(f < 0 ? c * (1 + f) : c + (255 - c) * f, 0, 255)));
    return '#' + ch.map((c) => c.toString(16).padStart(2, '0')).join('');
}

function hexLum(hex) {
    const m = /^#?([0-9a-f]{6})$/i.exec(hex || '');
    if (!m) return 0.5;
    const n = parseInt(m[1], 16);
    return (0.299 * ((n >> 16) & 255) + 0.587 * ((n >> 8) & 255) + 0.114 * (n & 255)) / 255;
}

const esc = (s) => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;');
const n2 = (v) => Number(v.toFixed(2));
const fmtAngle = (deg) => `${Number(deg.toFixed(1))}°`;
const OVERLINE = '̅';

// Advance widths in em (Arial / Times metrics), good enough to place labels and draw index bars
const EM_SANS = { digit: 0.556, other: 0.56, '[': 0.278, ']': 0.278, '(': 0.333, ')': 0.333, '°': 0.4, '.': 0.278, ',': 0.278, ' ': 0.278, '=': 0.584, '-': 0.333, '−': 0.584, Δ: 0.68 };
const EM_SERIF = { digit: 0.5, other: 0.5, '[': 0.333, ']': 0.333, '(': 0.333, ')': 0.333, '°': 0.4, '.': 0.25, ',': 0.25, ' ': 0.25, '=': 0.564, '-': 0.333, '−': 0.564, Δ: 0.64 };
const isSerif = (family) => /serif|times|georgia|palatino|garamond|cmu|computer modern|latin modern/i.test(family) && !/sans/i.test(family);

// ---------------------------------------------------------------------------
// Labels: "φ_B", "θ_{dev}" -> italic symbol + subscript; optional LaTeX output
// ---------------------------------------------------------------------------

const GREEK_TEX = {
    α: '\\alpha', β: '\\beta', γ: '\\gamma', δ: '\\delta', ε: '\\varepsilon', ζ: '\\zeta', η: '\\eta', θ: '\\theta',
    ϑ: '\\vartheta', κ: '\\kappa', λ: '\\lambda', μ: '\\mu', ν: '\\nu', ξ: '\\xi', π: '\\pi', ρ: '\\rho', σ: '\\sigma',
    τ: '\\tau', φ: '\\varphi', ϕ: '\\phi', χ: '\\chi', ψ: '\\psi', ω: '\\omega', Δ: '\\Delta', Θ: '\\Theta', Φ: '\\Phi', Ψ: '\\Psi', Ω: '\\Omega'
};

function parseSymbol(str) {
    const m = /^([^_]*)(?:_\{?([^}]*)\}?)?$/.exec(String(str));
    return { main: m ? m[1] : String(str), sub: m && m[2] ? m[2] : '' };
}

const texGreek = (s) => [...s].map((c) => (GREEK_TEX[c] ? GREEK_TEX[c] + ' ' : c)).join('').trim();

function symbolTex(str, vectorStyle) {
    const { main, sub: s } = parseSymbol(str);
    let body = texGreek(main);
    if (vectorStyle === 'bold') body = `\\mathbf{${body}}`;
    else if (vectorStyle === 'hat') body = `\\hat{${body}}`;
    else if (vectorStyle === 'arrow') body = `\\vec{${body}}`;
    if (s) body += s.length > 1 ? `_{\\mathrm{${s}}}` : `_{${texGreek(s)}}`;
    return body;
}

const VECTOR_MARK = { hat: '̂', arrow: '⃗' };

/** SVG tspans for a symbol; bold vectors render upright (like LaTeX \mathbf), Δ/∇ stay upright. */
function symbolSpans(str, fs, vectorStyle) {
    const { main, sub: s } = parseSymbol(str);
    const m = /^([Δ∇]+)(.*)$/.exec(main);
    const prefix = m ? m[1] : '';
    const varName = m ? m[2] : main;
    const bold = vectorStyle === 'bold';
    const varAttr = bold ? ' font-style="normal" font-weight="bold"' : ' font-style="italic"';
    let out = '';
    if (prefix) out += `<tspan font-style="normal">${esc(prefix)}</tspan>`;
    if (varName || !prefix) {
        out += `<tspan${varAttr}>${esc(varName + (VECTOR_MARK[vectorStyle] || ''))}</tspan>`;
    }
    if (s) {
        const italic = s.length === 1 ? ' font-style="italic"' : '';
        out += `<tspan dy="${n2(fs * 0.28)}" font-size="${n2(fs * 0.72)}"${italic}>${esc(s)}</tspan><tspan dy="${n2(-fs * 0.28)}">​</tspan>`;
    }
    return out;
}

// ---------------------------------------------------------------------------
// Projection + 2D geometry
// ---------------------------------------------------------------------------

function makeCamera(view) {
    const { azDeg, elDeg, shear = 0 } = VIEWS[view] || VIEWS.oblique;
    const az = azDeg * DEG;
    const el = elDeg * DEG;
    const right = [-Math.sin(az), Math.cos(az), 0];
    const up = [-Math.sin(el) * Math.cos(az), -Math.sin(el) * Math.sin(az), Math.cos(el)];
    // shear: depth y moves right on screen; toViewer is the matching (oblique) projection direction
    const v0 = [Math.cos(el) * Math.cos(az), Math.cos(el) * Math.sin(az), Math.sin(el)];
    const a = (-shear * v0[1]) / (1 + shear * right[1]);
    const vv = [v0[0] + a * right[0], v0[1] + a * right[1], v0[2]];
    const vl = Math.hypot(...vv);
    const toViewer = vv.map((c) => c / vl);
    const proj = (p) => [dot(p, right) + shear * p[1], dot(p, up)];
    return { right, up, toViewer, proj };
}

/** Convex hull of 2D points (monotone chain). */
function convexHull(points) {
    const p = [...points].sort((a, b) => a[0] - b[0] || a[1] - b[1]);
    const cross = (o, a, b) => (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0]);
    const half = (list) => {
        const h = [];
        for (const q of list) {
            while (h.length >= 2 && cross(h[h.length - 2], h[h.length - 1], q) <= 0) h.pop();
            h.push(q);
        }
        h.pop();
        return h;
    };
    return [...half(p), ...half([...p].reverse())];
}

/** Separating-axis test for two convex polygons. */
function convexOverlap(a, b) {
    for (const poly of [a, b]) {
        for (let i = 0; i < poly.length; i++) {
            const p = poly[i];
            const q = poly[(i + 1) % poly.length];
            const ax = [q[1] - p[1], p[0] - q[0]];
            const span = (pts) => pts.map((r) => r[0] * ax[0] + r[1] * ax[1]);
            const sa = span(a);
            const sb = span(b);
            if (Math.max(...sa) < Math.min(...sb) || Math.max(...sb) < Math.min(...sa)) return false;
        }
    }
    return true;
}

/** Area of a convex polygon inside the box {l, r, t, b} (Sutherland-Hodgman clip). */
function areaInBox(poly, box) {
    const edges = [
        (p) => p[0] - box.l,
        (p) => box.r - p[0],
        (p) => p[1] - box.t,
        (p) => box.b - p[1],
    ];
    let pts = poly;
    for (const inside of edges) {
        const out = [];
        pts.forEach((p, i) => {
            const q = pts[(i + 1) % pts.length];
            const dp = inside(p);
            const dq = inside(q);
            if (dp >= 0) out.push(p);
            if (dp >= 0 !== dq >= 0) {
                const t = dp / (dp - dq);
                out.push([p[0] + (q[0] - p[0]) * t, p[1] + (q[1] - p[1]) * t]);
            }
        });
        pts = out;
        if (!pts.length) return 0;
    }
    let a = 0;
    pts.forEach((p, i) => {
        const q = pts[(i + 1) % pts.length];
        a += p[0] * q[1] - q[0] * p[1];
    });
    return Math.abs(a) / 2;
}

/** Does segment a-b touch the box {l, r, t, b}? (Liang-Barsky clip) */
function segHitsBox(a, b, box) {
    let t0 = 0;
    let t1 = 1;
    const dx = b[0] - a[0];
    const dy = b[1] - a[1];
    const p = [-dx, dx, -dy, dy];
    const q = [a[0] - box.l, box.r - a[0], a[1] - box.t, box.b - a[1]];
    for (let i = 0; i < 4; i++) {
        if (p[i] === 0) {
            if (q[i] < 0) return false;
        } else {
            const t = q[i] / p[i];
            if (p[i] < 0) {
                if (t > t1) return false;
                if (t > t0) t0 = t;
            } else {
                if (t < t0) return false;
                if (t < t1) t1 = t;
            }
        }
    }
    return true;
}

/** Do segments a-b and c-d cross (inner points only, touching ends do not count)? */
function segsCross(a, b, c, d) {
    const cr = (p, q, r) => (q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0]);
    return cr(c, d, a) * cr(c, d, b) < -1e-9 && cr(a, b, c) * cr(a, b, d) < -1e-9;
}

/** Point of the box {l, r, t, b} closest to p */
const nearestInBox = (p, b) => [Math.max(b.l, Math.min(b.r, p[0])), Math.max(b.t, Math.min(b.b, p[1]))];
const distToBox = (p, b) => {
    const q = nearestInBox(p, b);
    return Math.hypot(q[0] - p[0], q[1] - p[1]);
};

const boxPoly = (b) => [[b.l, b.t], [b.r, b.t], [b.r, b.b], [b.l, b.b]];
const grow = (b, m) => ({ l: b.l - m, r: b.r + m, t: b.t - m, b: b.b + m });
const unionBox = (a, b) => ({ l: Math.min(a.l, b.l), r: Math.max(a.r, b.r), t: Math.min(a.t, b.t), b: Math.max(a.b, b.b) });
const boxesTouch = (a, b) => a.l < b.r && b.l < a.r && a.t < b.b && b.t < a.b;
function boundsOf(pts) {
    const xs = pts.map((p) => p[0]);
    const ys = pts.map((p) => p[1]);
    return { l: Math.min(...xs), r: Math.max(...xs), t: Math.min(...ys), b: Math.max(...ys) };
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

/**
 * @param state   result of evaluateSetup()
 * @param options see DEFAULT_FIGURE_OPTIONS
 * @param styleIn see BUILTIN_FIGURE_STYLES (partial styles are completed)
 * @returns SVG document string
 */
const num = (x, fallback) => (x === '' || x == null || !Number.isFinite(Number(x)) ? fallback : Number(x));

function resolveSweepSymbol(st, isOopSweep) {
    const custom = String(st.symbols?.sweep || '').trim();
    if (custom && custom !== 'Δφ' && custom !== 'Δθ') return custom;
    const base = isOopSweep ? st.symbols?.theta || 'θ' : st.symbols?.phi || 'φ';
    return `Δ${base}`;
}

// Sweep and drawn field, the same for the figure and its caption: the magnet turns the field in the
// lab xy plane, and the drawn field is where the turn starts
function fieldSetup(state, o) {
    let bAt = null;
    let swFrom = 0;
    let swSpan = 0;
    let isOopSweep = state.setup?.mountMode === 'OOP';
    if (o.sweep === 'angle' && state.matrices) {
        const e1 = state.matrices.vectorInAllFrames([1, 0, 0], 'Lab').Sub;
        const e2 = state.matrices.vectorInAllFrames([0, 1, 0], 'Lab').Sub;
        bAt = (phi) => unit(add(mul(e1, Math.cos(phi * DEG)), mul(e2, Math.sin(phi * DEG))));
        const nRot = cross(e1, e2);
        isOopSweep = Math.abs(nRot[2]) < 0.5;
        swFrom = num(o.sweepFrom, Number(state.setup.phiMDeg) || 0);
        swSpan = Math.max(-360, Math.min(360, num(o.sweepSpan, 0)));
    }
    const sweepAngle = Math.abs(swSpan) > TINY_DEG;
    const b0mT = (Number(state.setup.b0Tesla) || 0) * 1000;
    const fieldRange = o.sweep === 'field' ? [num(o.sweepFrom, 0), num(o.sweepTo, b0mT)] : null;
    const bHat = sweepAngle ? bAt(swFrom) : unit(state.vectors.bHat.Sub);
    return { bAt, swFrom, swSpan, sweepAngle, isOopSweep, b0mT, fieldRange, bHat };
}

export function generateFigureSVG(state, options = {}, styleIn = null) {
    const o = { ...DEFAULT_FIGURE_OPTIONS, ...options };
    const st = completeFigureStyle(styleIn);
    const C = st.colors;
    const L = st.lines;
    const E = st.effects;
    const fs = Number(st.fontSizePt) || 9;
    const latex = !!o.latex;
    const cam = makeCamera(o.view);
    const isTop = o.view === 'top';
    const thick = isTop ? 0 : Math.max(0, Number(st.sampleThickness) || 0);
    const EM = isSerif(st.fontFamily) ? EM_SERIF : EM_SANS;

    const v = state.vectors;
    const kHat = v.kHat.Sub;
    const tHat = v.tHat.Sub;

    // ---- sweep: the magnet turns the field in the lab xy plane; the drawn field is where it starts
    const { bAt, swFrom, swSpan, sweepAngle, isOopSweep, fieldRange, bHat } = fieldSetup(state, o);
    const sweepSym = resolveSweepSymbol(st, isOopSweep);
    const xEdge = [1, 0, 0];
    const nrm = [0, 0, 1];
    const origin = [0, 0, 0];

    // ---- scale: the drawn sample is widthMm wide in every figure; the crop follows the content
    const corners = (z) => [[-HX, -HY, z], [HX, -HY, z], [HX, HY, z], [-HX, HY, z]];
    const slab = [...corners(0), ...corners(-thick)];
    const slabX = slab.map((p) => cam.proj(p)[0]);
    const scale = ((Number(st.widthMm) || 55) * PT_PER_MM) / (Math.max(...slabX) - Math.min(...slabX));
    const toSvg = (p) => {
        const [x, y] = cam.proj(p);
        return [x * scale, -y * scale];
    };
    const pt = (p) => toSvg(p).map(n2).join(',');
    const pathOf = (pts, close = false) => 'M' + pts.map(pt).join('L') + (close ? 'Z' : '');

    const layers = { sample: [], antenna: [], guides: [], arcs: [], sweep: [], arrows: [], labels: [], axes: [] };

    // Everything drawn leaves its screen points (for the crop) and, where a label would look bad
    // on top of it, segments or areas with a penalty weight
    const ink = [];
    const lines = [];
    const areas = [];
    // Segments carry the id of their drawn object, so a label crossing an arc or an arrowhead
    // pays the weight once, not once per segment
    let lineId = 0;
    const addLine = (pts, w, closed = false) => {
        ink.push(...pts);
        const id = lineId++;
        for (let i = 0; i + 1 < pts.length; i++) lines.push([pts[i], pts[i + 1], w, id]);
        if (closed) lines.push([pts[pts.length - 1], pts[0], w, id]);
        return id;
    };
    // Drawn arcs (screen points): an angle value must sit closer to its own arc than to any other
    const marks = [];
    const thin = new Set(); // lines a value may be seen across (the sweep rim), though not sit on
    const W_LINE = 12;
    const W_GUIDE = 8;
    const W_EDGE = 3;
    const W_ANTENNA = 10;

    // ---- text metrics (estimates from Arial / Times widths)
    const charW = (c, size = fs) => size * (/\d/.test(c) ? EM.digit : EM[c] ?? EM.other);
    const textWidth = (plain) => [...String(plain)].filter((c) => c !== OVERLINE).reduce((s, c) => s + charW(c), 0);
    const symbolWidth = (str) => {
        const { main, sub: s } = parseSymbol(str);
        return textWidth(main) * 1.08 + [...s].reduce((w, c) => w + charW(c, fs * 0.72), 0);
    };

    // Miller-Bravais label "[1̅1̅20]": one centred chunk per character, bars drawn as short lines
    // (text-decoration overlines sit too high and too heavy)
    const millerLabel = (label, latexLabel) => {
        const width = textWidth(label);
        if (latex) return { content: esc(`$${latexLabel}$`), width };
        if (/\s/.test(label)) return { content: esc(label), width };
        const chars = [];
        [...label].forEach((c) => {
            if (c === OVERLINE && chars.length) chars[chars.length - 1].bar = true;
            else chars.push({ c });
        });
        const content = (left, y) => {
            let x = left;
            let spans = '';
            let bars = '';
            chars.forEach(({ c, bar }) => {
                const w = charW(c);
                const cx = x + w / 2;
                spans += `<tspan x="${n2(cx)}">${esc(c)}</tspan>`;
                if (bar) bars += `M${n2(cx - fs * 0.18)},${n2(y - fs * 0.79)}H${n2(cx + fs * 0.18)}`;
                x += w;
            });
            const text = `<text y="${n2(y)}" text-anchor="middle" fill="${C.text}" font-size="${fs}">${spans}</text>`;
            return bars ? `<g>${text}<path d="${bars}" stroke="${C.text}" stroke-width="${n2(fs * 0.06)}" stroke-linecap="butt" fill="none"/></g>` : text;
        };
        return { content, width };
    };

    // Label boxes: just beyond a screen point, pushed along a screen direction (unit, y down),
    // or centred on a point
    const labelBox = (c, w) => {
        if (c.center) {
            const y = c.sy + fs * 0.28;
            return { l: c.sx - w / 2, r: c.sx + w / 2, t: y - fs * 0.76, b: y + fs * 0.2, x: c.sx, y, anchor: 'middle' };
        }
        const x = c.sx + c.dx * c.gap;
        const y = c.sy + c.dy * c.gap + fs * 0.35 + c.dy * fs * 0.35;
        const anchor = c.dx > 0.35 ? 'start' : c.dx < -0.35 ? 'end' : 'middle';
        const l = anchor === 'start' ? x : anchor === 'end' ? x - w : x - w / 2;
        return { l, r: l + w, t: y - fs * 0.76, b: y + fs * 0.2, x, y, anchor };
    };
    // Spot for a w-wide label whose box clears the screen point (ex, ey) by gap in direction (ux, uy):
    // the label sits on the extension of that direction, at the same distance whatever its width
    const beyond = (ex, ey, ux, uy, w, gap, pen = 0) => {
        const h = fs * 0.96;
        const t = Math.min(Math.abs(ux) > 1e-6 ? w / 2 / Math.abs(ux) : Infinity, Math.abs(uy) > 1e-6 ? h / 2 / Math.abs(uy) : Infinity);
        return { sx: ex + ux * (t + gap), sy: ey + uy * (t + gap), center: true, pen };
    };
    const turn = ([ux, uy], deg) => [ux * Math.cos(deg * DEG) - uy * Math.sin(deg * DEG), ux * Math.sin(deg * DEG) + uy * Math.cos(deg * DEG)];

    // Labels are collected and placed at the end, most important first. Each takes the cheapest of
    // its candidate spots: no other label, no line through it, little antenna or sample edge under
    // it, close to what it labels, and not far outside the drawing (that would make it bigger)
    const labelReqs = [];
    // An angle value also names its arc (own: line id and screen points, anchor: middle of the arc):
    // the way from the arc to the value crosses no arrow, and no other arc is closer to it
    const wantLabel = (prio, cands, content, color, w, arcOf = null) => labelReqs.push({ prio, cands, content, color, w, arcOf });
    const arcCost = (q, b) => {
        const { id, pts, anchor } = q.arcOf;
        let cost = 0;
        const sight = [anchor, nearestInBox(anchor, b)];
        const hit = new Map();
        for (const [a, z, w, lid] of lines) {
            if (lid === id || hit.has(lid) || w < W_GUIDE || !segsCross(sight[0], sight[1], a, z)) continue;
            hit.set(lid, w >= W_LINE && !thin.has(lid) ? 15 : 4);
        }
        hit.forEach((c) => (cost += c));
        const near = (ps) => Math.min(...ps.map((p) => distToBox(p, b)));
        const own = near(pts);
        const other = Math.min(Infinity, ...marks.filter((mk) => mk.id !== id && (mk.angle || !q.arcOf.angle)).map((mk) => near(mk.pts)));
        if (other < own + fs * 0.3) cost += 3 + Math.max(0, own - other) / fs;
        return cost;
    };
    const placeLabels = (core, taken) => {
        const placed = taken ? [taken] : [];
        const m = fs * 0.15;
        // the spot at an arrow tip is kept for that arrow's name
        const tips = labelReqs.filter((q) => q.cands[0]?.tip).map((q) => ({ q, box: grow(labelBox(q.cands[0], q.w), m) }));
        const names = []; // placed arrow names: a value right beside one reads as part of it
        labelReqs
            .sort((a, b) => a.prio - b.prio)
            .forEach((q) => {
                let best = null;
                let bestCost = Infinity;
                q.cands.forEach((c, i) => {
                    const b = labelBox(c, q.w);
                    const e = grow(b, m);
                    let cost = i * 0.02 + (c.pen || 0);
                    if (placed.some((p) => boxesTouch(e, p))) cost += 100;
                    tips.forEach((t) => t.q !== q && boxesTouch(e, t.box) && (cost += c.leader ? 20 : 6)); // a value set off belongs nowhere near a tip
                    const hit = new Map();
                    // a value on its own arc hides the arc it names
                    for (const [a, z, w, id] of lines) if (!hit.has(id) && segHitsBox(a, z, e)) hit.set(id, id === q.arcOf?.id ? W_LINE : w);
                    hit.forEach((w) => (cost += w));
                    // a label grazing an antenna is fine, one covering it is not
                    const boxArea = (e.r - e.l) * (e.b - e.t);
                    for (const [hull, w] of areas) {
                        const a = areaInBox(hull, e);
                        if (a > 0) cost += w * (0.2 + (3 * a) / boxArea);
                    }
                    const over = Math.max(0, core.l - b.l) + Math.max(0, b.r - core.r) + Math.max(0, core.t - b.t) + Math.max(0, b.b - core.b);
                    cost += (1.5 * over) / fs;
                    if (q.arcOf) cost += arcCost(q, b) + names.filter((nb) => boxesTouch(b, grow(nb, fs * 0.6))).length * 4;
                    if (cost < bestCost) {
                        bestCost = cost;
                        best = b;
                        best.leader = c.leader ? c.from || q.arcOf.anchor : null;
                    }
                });
                placed.push(best);
                if (q.cands[0]?.tip) names.push(best);
                // a value set off from its arc gets a thin line to it
                if (best.leader) {
                    const [ax, ay] = best.leader;
                    const [bx, by] = nearestInBox(best.leader, grow(best, fs * 0.12));
                    const d = Math.hypot(bx - ax, by - ay) || 1;
                    const s = [ax + ((bx - ax) / d) * fs * 0.15, ay + ((by - ay) / d) * fs * 0.15];
                    addLine([s, [bx, by]], W_GUIDE);
                    layers.arcs.push(`<path d="M${n2(s[0])},${n2(s[1])}L${n2(bx)},${n2(by)}" stroke="${q.color === C.text ? C.angle : q.color}" stroke-width="${n2(L.guide * 0.7)}" fill="none"/>`);
                }
                layers.labels.push(
                    typeof q.content === 'function'
                        ? q.content(best.l, best.y)
                        : `<text x="${n2(best.x)}" y="${n2(best.y)}" text-anchor="${best.anchor}" fill="${q.color}" font-size="${fs}">${q.content}</text>`
                );
            });
        return placed;
    };

    const screenDir = (from, to) => {
        const a = toSvg(from);
        const b = toSvg(to);
        const d = Math.hypot(b[0] - a[0], b[1] - a[1]) || 1;
        return [(b[0] - a[0]) / d, (b[1] - a[1]) / d];
    };

    // Slim arrowhead with a swept-back base; the shaft ends inside it
    const headPath = (tip, ux, uy, hl, color) => {
        const hw = hl * 0.29;
        const bx = tip[0] - ux * hl;
        const by = tip[1] - uy * hl;
        const nx = tip[0] - ux * hl * 0.78;
        const ny = tip[1] - uy * hl * 0.78;
        return `<path d="M${n2(tip[0])},${n2(tip[1])}L${n2(bx - uy * hw)},${n2(by + ux * hw)}L${n2(nx)},${n2(ny)}L${n2(bx + uy * hw)},${n2(by - ux * hw)}Z" fill="${color}"/>`;
    };

    // Arrow in 3D, head drawn in screen space so it keeps its size when foreshortened
    const arrow = (from, to, color, width, { dashed = false, opacity = 1, w = W_LINE } = {}) => {
        const a = toSvg(from);
        const b = toSvg(to);
        const dx = b[0] - a[0];
        const dy = b[1] - a[1];
        const d = Math.hypot(dx, dy);
        if (d < 0.5) return '';
        const ux = dx / d;
        const uy = dy / d;
        const hl = Math.min(st.arrowHeadPt * (0.7 + width / 3.5), d * 0.6);
        const hw = hl * 0.29;
        addLine([a, b], w);
        addLine([b, [b[0] - ux * hl - uy * hw, b[1] - uy * hl + ux * hw], [b[0] - ux * hl + uy * hw, b[1] - uy * hl - ux * hw]], w, true);
        const dash = dashed ? ` stroke-dasharray="${n2(width * 1.6)},${n2(width * 1.4)}" stroke-linecap="butt"` : '';
        const svg =
            `<path d="M${n2(a[0])},${n2(a[1])}L${n2(b[0] - ux * hl * 0.7)},${n2(b[1] - uy * hl * 0.7)}" stroke="${color}" stroke-width="${n2(width)}"${dash} fill="none"/>` +
            headPath(b, ux, uy, hl, color);
        return opacity < 1 ? `<g opacity="${opacity}">${svg}</g>` : svg;
    };

    const guideLine = (from, to, color = C.angle) => {
        addLine([toSvg(from), toSvg(to)], W_GUIDE);
        return `<path d="${pathOf([from, to])}" stroke="${color}" stroke-width="${n2(L.guide)}" stroke-dasharray="${n2(L.guide * 1.5)},${n2(L.guide * 2.6)}" fill="none"/>`;
    };

    // Arc from direction a to direction b (3D unit vectors) around the point c
    const arc = (a, b, r, c = origin) => {
        const ang = angleDeg(a, b);
        const ua = unit(a);
        let perp = sub(unit(b), mul(ua, dot(ua, unit(b))));
        if (len(perp) < 1e-9) perp = unit([-ua[1], ua[0], 0]);
        perp = unit(perp);
        const dirAt = (f) => add(mul(ua, Math.cos(ang * f * DEG)), mul(perp, Math.sin(ang * f * DEG)));
        const steps = Math.max(8, Math.ceil(ang / 3));
        const pts = Array.from({ length: steps + 1 }, (_, i) => add(c, mul(dirAt(i / steps), r)));
        return { ang, pts, mid: dirAt(0.5), dirAt, r, c };
    };

    const angleText = (symbol, deg, mode = o.angleLabels) => {
        if (latex) {
            const s = symbolTex(symbol, 'none');
            const val = `${Number(deg.toFixed(1))}^\\circ`;
            return esc(mode === 'symbol' ? `$${s}$` : mode === 'both' ? `$${s} = ${val}$` : `$${val}$`);
        }
        if (mode === 'symbol') return symbolSpans(symbol, fs, 'none');
        if (mode === 'both') return `${symbolSpans(symbol, fs, 'none')} = ${fmtAngle(deg)}`;
        return fmtAngle(deg);
    };
    const angleWidth = (symbol, deg, mode = o.angleLabels) =>
        mode === 'symbol' ? symbolWidth(symbol) : mode === 'both' ? symbolWidth(symbol) + textWidth(` = ${fmtAngle(deg)}`) : textWidth(fmtAngle(deg));
    const vectorText = (symbol) => (latex ? esc(`$${symbolTex(symbol, st.vectorStyle)}$`) : symbolSpans(symbol, fs, st.vectorStyle));

    // Arc with its value right at it: just outside the middle of the arc, a bit further along or
    // out if that spot is taken, or inside the wedge. Never far away.
    let arcCount = 0;
    const drawArc = (a, b, r, symbol, c = origin, color = C.angle) => {
        const A = arc(a, b, r, c);
        if (A.ang < TINY_DEG) return;
        const scr = A.pts.map(toSvg);
        const id = addLine(scr, W_GUIDE);
        marks.push({ id, pts: scr, angle: true });
        const isBelow = A.pts.some((p) => p[2] < -1e-6);
        if (E.angleFill) {
            layers.arcs.push(`<path d="${pathOf([c, ...A.pts], true)}" fill="${color}" fill-opacity="${isBelow ? 0.04 : 0.09}" stroke="none"/>`);
        }
        const below = isBelow ? ` stroke-dasharray="${n2(L.arc * 2.5)},${n2(L.arc * 2)}"` : '';
        layers.arcs.push(`<path d="${pathOf(A.pts)}" stroke="${color}" stroke-width="${n2(L.arc)}"${below} fill="none"/>`);
        const w = angleWidth(symbol, A.ang);
        const cands = [];
        for (const [gap, pen] of [[0.2, 0], [0.55, 0.8], [1, 2.5]]) {
            for (const f of [0.5, 0.35, 0.65, 0.2, 0.8]) {
                const d = A.dirAt(f);
                const [sx, sy] = toSvg(add(c, mul(d, r)));
                const [ux, uy] = screenDir(c, add(c, d));
                cands.push({ ...beyond(sx, sy, ux, uy, w, gap * fs, pen + Math.abs(f - 0.5) * 1.5), from: [sx, sy], leader: gap >= 1 });
            }
        }
        const [ix, iy] = toSvg(add(c, mul(A.mid, r * 0.55)));
        cands.push({ sx: ix, sy: iy, center: true, pen: 0.6 });
        // flat wedges leave no room inside: then beside an end of the arc, outside the wedge
        const [mx, my] = toSvg(add(c, mul(A.mid, r)));
        const awayRays = [0, 1].map((f) => {
            const [sx, sy] = toSvg(add(c, mul(A.dirAt(f), r)));
            const [rx, ry] = screenDir(c, add(c, A.dirAt(f)));
            const s = (mx - sx) * -ry + (my - sy) * rx > 0 ? -1 : 1; // normal of the ray, away from the arc
            return { sx, sy, nx: -ry * s, ny: rx * s, s };
        });
        for (const { sx, sy, nx, ny, s } of awayRays) {
            for (const [g, pen] of [[0.25, 1.2], [0.7, 1.65]]) {
                cands.push(
                    { ...beyond(sx, sy, nx, ny, w, fs * g, pen), from: [sx, sy] },
                    { ...beyond(sx, sy, ...turn([nx, ny], s * -35), w, fs * g, pen + 0.2), from: [sx, sy] }
                );
            }
        }
        // crowded: further out, with a line to the arc (radially outward or away from wedge ends)
        for (const [gap, pen] of [[1.6, 2.2], [2.5, 3]]) {
            for (const f of [0.5, 0.3, 0.7]) {
                const d = A.dirAt(f);
                const [sx, sy] = toSvg(add(c, mul(d, r)));
                const [ux, uy] = screenDir(c, add(c, d));
                cands.push({ ...beyond(sx, sy, ux, uy, w, gap * fs, pen + Math.abs(f - 0.5) * 1.5), from: [sx, sy], leader: true });
            }
            for (const { sx, sy, nx, ny, s } of awayRays) {
                cands.push(
                    { ...beyond(sx, sy, nx, ny, w, gap * fs, pen + 0.4), from: [sx, sy], leader: true },
                    { ...beyond(sx, sy, ...turn([nx, ny], s * -30), w, gap * fs, pen + 0.6), from: [sx, sy], leader: true }
                );
            }
        }
        wantLabel(1 + arcCount++ * 0.1, cands, angleText(symbol, A.ang), color, w, { id, pts: scr, anchor: [mx, my], angle: true });
    };

    // Label at the free end of an arrow, or beside its shaft
    const arrowLabelCands = (from, to, w) => {
        const a = toSvg(from);
        const b = toSvg(to);
        const u = screenDir(from, to);
        const cands = [[0, 0], [45, 0.5], [-45, 0.5], [90, 1], [-90, 1]].map(([deg, pen]) => {
            const [ux, uy] = turn(u, deg);
            return { ...beyond(b[0], b[1], ux, uy, w, fs * 0.25, pen), tip: deg === 0 };
        });
        for (const f of [0.8, 0.6, 0.3]) {
            const sx = a[0] + (b[0] - a[0]) * f;
            const sy = a[1] + (b[1] - a[1]) * f;
            // long labels beside a slanted shaft need more room to clear it
            for (const [g, pen] of w > fs * 4 ? [[0.3, 0.8], [0.9, 1.6], [1.6, 2.6]] : [[0.3, 0.8]]) for (const s of [1, -1]) cands.push(beyond(sx, sy, -u[1] * s, u[0] * s, w, fs * g, pen + (0.8 - f)));
        }
        return cands;
    };

    // ---- sample slab: crisp pure-flat planes with light-aware ink outlines
    const top = corners(0);
    const lightSample = hexLum(C.sample) > 0.88;
    const sampleStroke = lightSample ? C.axes : shade(C.sample, -0.48);
    if (thick > 0) {
        const bot = corners(-thick);
        const sides = [
            { n: [0, -1, 0], q: [top[0], top[1], bot[1], bot[0]] },
            { n: [1, 0, 0], q: [top[1], top[2], bot[2], bot[1]] },
            { n: [0, 1, 0], q: [top[2], top[3], bot[3], bot[2]] },
            { n: [-1, 0, 0], q: [top[3], top[0], bot[0], bot[3]] }
        ];
        sides
            .filter((f) => dot(f.n, cam.toViewer) > 1e-6)
            .forEach((f) => {
                addLine(f.q.map(toSvg), W_EDGE, true);
                const lit = dot(f.n, LIGHT_DIR) > 0;
                const fill = shade(C.sample, lightSample ? (lit ? -0.07 : -0.15) : lit ? -0.16 : -0.28);
                layers.sample.push(`<path d="${pathOf(f.q, true)}" fill="${fill}" stroke="${sampleStroke}" stroke-width="${n2(L.outline)}" stroke-linejoin="round"/>`);
            });
    }
    addLine(top.map(toSvg), W_EDGE, true);
    layers.sample.push(`<path d="${pathOf(top, true)}" fill="${C.sample}" stroke="${sampleStroke}" stroke-width="${n2(L.outline)}" stroke-linejoin="round"/>`);

    // ---- emitter at -k and (optionally) receiver at +k, shaped as in the 3D view. A turned pair
    // that would not fit moves closer together (keeping room for k) and shrinks only as needed.
    const wave = (kk, tt) => add(mul(kHat, kk), mul(tHat, tt));
    const ownType = state.setup.deviceType === 'IDT' ? 'IDT' : 'CPW';
    const devTypes = o.showAntenna ? (DEVICE_SETS[o.devices] || DEVICE_SETS.pair)(ownType) : [];
    const shapes = devTypes.map(microstructure);
    const sideOf = (i) => (i === 0 ? -1 : 1);
    // receiver mirrored to the emitter, as on our samples: along k, and an antenna (symmetric along k)
    // also across k, so its shorting bar (+t in the shape) faces the front edge, the emitter's the back
    const back = tHat[1] < 0 ? -1 : 1;
    const tSide = (i) => back * (i === 1 && shapes.length > 1 && devTypes[1] === 'CPW' ? -1 : 1);
    const devPoint = (i, dist, s) => ([pk, pt]) => wave(sideOf(i) * (dist - pk * s), tSide(i) * pt * s);
    const overflow = (dist, s) =>
        Math.max(0, ...shapes.flatMap((polys, i) => polys.flat().map((q) => {
            const p = devPoint(i, dist, s)(q);
            return Math.max(Math.abs(p[0]) / (HX - 0.08), Math.abs(p[1]) / (HY - 0.08));
        })));
    const hkMax = Math.max(0, ...devTypes.map(halfLength));
    let devDist = DEV_DIST;
    let devScale = DEV_SCALE;
    if (shapes.length) {
        let best = null;
        for (let dist = DEV_DIST; dist > MIN_HALF_GAP; dist -= 0.02) {
            let s = Math.min(DEV_SCALE, (dist - MIN_HALF_GAP) / hkMax);
            while (s > 0.3 && overflow(dist, s) > 1) s *= 0.98;
            const score = s / DEV_SCALE + (0.5 * dist) / DEV_DIST; // full size first, then far apart
            if (!best || score > best.score + 1e-6) best = { dist, s, score };
        }
        devDist = best.dist;
        devScale = best.s;
    }
    const hDev = E.extrude && !isTop && thick > 0 ? Math.min(thick * 0.2, 0.018) : 0;
    const devPolys = [];
    const devTopPolys = [];
    shapes.forEach((polys, i) => {
        const pts3 = polys.map((poly) => poly.map(devPoint(i, devDist, devScale)));
        const lifted = hDev > 0 ? pts3.map((poly) => poly.map((p) => [p[0], p[1], hDev])) : pts3;
        const screen = lifted.flat().map(toSvg);
        ink.push(...screen);
        areas.push([convexHull(screen), W_ANTENNA]);
        devPolys.push(...pts3);
        devTopPolys.push(...lifted);
    });
    if (devPolys.length) {
        if (hDev > 0) {
            const walls = [];
            devPolys.forEach((poly) => {
                let area2 = 0;
                for (let j = 0; j < poly.length; j++) {
                    const a = poly[j];
                    const b = poly[(j + 1) % poly.length];
                    area2 += a[0] * b[1] - b[0] * a[1];
                }
                const ccw = area2 >= 0 ? 1 : -1;
                for (let j = 0; j < poly.length; j++) {
                    const p = poly[j];
                    const q = poly[(j + 1) % poly.length];
                    const dx = q[0] - p[0];
                    const dy = q[1] - p[1];
                    const edgeLen = Math.hypot(dx, dy);
                    if (edgeLen < 1e-9) continue;
                    const n = [ccw * (dy / edgeLen), -ccw * (dx / edgeLen), 0];
                    if (dot(n, cam.toViewer) <= 1e-6) continue;
                    const cen = [(p[0] + q[0]) / 2, (p[1] + q[1]) / 2, hDev / 2];
                    const quad = [[p[0], p[1], hDev], [q[0], q[1], hDev], [q[0], q[1], 0], [p[0], p[1], 0]];
                    const fill = shade(C.antenna, dot(n, LIGHT_DIR) > 0 ? -0.18 : -0.34);
                    walls.push({ depth: dot(cen, cam.toViewer), quad, fill });
                }
            });
            walls.sort((a, b) => a.depth - b.depth);
            const wallSw = n2(Math.max(0.15, L.outline * 0.35));
            walls.forEach((w) => {
                layers.antenna.push(`<path d="${pathOf(w.quad, true)}" fill="${w.fill}" stroke="${w.fill}" stroke-width="${wallSw}" stroke-linejoin="round"/>`);
            });
        }
        const antStroke = shade(C.antenna, -0.42);
        const antSw = n2(Math.max(0.35, L.outline * 0.65));
        layers.antenna.push(
            `<path d="${devTopPolys.map((q) => pathOf(q, true)).join('')}" fill="${C.antenna}" stroke="${antStroke}" stroke-width="${antSw}" stroke-linejoin="round"/>`
        );
    }

    // k runs from the emitter to the receiver (or as far without one)
    const kLim = Math.min(HX / Math.max(Math.abs(kHat[0]), 1e-9), HY / Math.max(Math.abs(kHat[1]), 1e-9)) - 0.2;
    const kFrom = shapes.length ? -devDist + halfLength(devTypes[0]) * devScale + K_GAP : -Math.min(0.8, kLim);
    const kTo = shapes.length > 1 ? devDist - halfLength(devTypes[1]) * devScale - K_GAP : shapes.length ? devDist : Math.min(0.8, kLim);

    // Optional clean flat vector sine wave on the surface along k
    if (E.waveLine && o.showWave) {
        const wk0 = shapes.length ? -devDist + halfLength(devTypes[0]) * devScale + 0.03 : kFrom;
        const wk1 = shapes.length > 1 ? devDist - halfLength(devTypes[1]) * devScale - 0.03 : kTo;
        if (wk1 - wk0 > 0.25) {
            const amp = Math.min(0.11 * devScale, 0.11);
            const pts = Array.from({ length: 73 }, (_, i) => {
                const u = i / 72;
                return wave(wk0 + (wk1 - wk0) * u, amp * Math.sin(Math.PI * u) * Math.sin(2 * Math.PI * 4.5 * u));
            });
            layers.sample.push(`<path d="${pathOf(pts)}" stroke="${C.wave}" stroke-opacity="0.42" stroke-width="${n2(L.arc)}" fill="none"/>`);
        }
    }
    const P0 = mul(kHat, (kFrom + kTo) / 2);
    const rs = Math.max(0.55, Math.min(1, (kTo - kFrom) / 2 / 0.75)); // in-plane arcs shrink with a short k
    const rDev = R_DEV * rs;
    // A turning field is as long as the antennas allow where it runs near the surface; its turn is an
    // arc on the shaft, outside the angle arcs
    let fl = FIELD_LEN;
    if (sweepAngle) {
        let maxK = 0;
        for (let i = 0; i <= 72; i++) {
            const d = bAt(swFrom + (swSpan * i) / 72);
            if (Math.abs(d[2]) < 0.25) maxK = Math.max(maxK, Math.abs(dot(d, kHat)));
        }
        if (maxK > 1e-6) fl = Math.max(0.55, Math.min(FIELD_LEN, ((kTo - kFrom) / 2 - 0.04) / maxK));
    }
    // a full turn rings the angle arcs: wider, so their values fit inside it
    const rSweep = (Math.abs(swSpan) >= 360 - 1e-6 ? 0.85 : 0.7) * fl;
    const rPhi = Math.min(R_PHI * rs, (sweepAngle ? 0.4 : 0.6) * fl);
    const rElev = Math.min(R_ELEV, (sweepAngle ? 0.52 : 0.75) * fl);

    // ---- angles: antenna turn (0° edge -> k), field in plane, field out of plane
    const devDeg = angleDeg(xEdge, kHat);
    const showDev = devDeg > TINY_DEG && (shapes.length > 0 || o.showWave);
    const bn = dot(bHat, nrm);
    const bPlane = [bHat[0], bHat[1], 0];
    const pDir = len(bPlane) > 1e-6 ? unit(bPlane) : null;
    const elevDeg = Math.asin(Math.min(1, Math.abs(bn))) / DEG;
    const showElev = elevDeg > TINY_DEG;
    const refDir = o.angleRef === 'edge' ? xEdge : kHat;
    const rayLen = (r) => r + 0.12; // short reference rays: just past their arc

    // The field starts on the k arrow, halfway between the antennas; a field (nearly) along k
    // starts beside it instead, so the two arrows do not cover each other
    // Field arrow stays above the surface: pointing away from it, or (B·n < 0) coming down onto it.
    // A turning field always starts at its foot (its tip runs round); below the surface it is dashed.
    const down = !sweepAngle && bn < -1e-6;
    const alongK = o.showWave && pDir && (isTop || elevDeg < 35) && Math.abs(dot(pDir, kHat)) > Math.cos(12 * DEG);
    const foot = alongK
        ? add(P0, add(mul(tHat, (tHat[1] >= 0 ? 1 : -1) * 0.3), mul(pDir, sweepAngle ? 0 : (-fl * Math.cos(elevDeg * DEG) * (down ? -1 : 1)) / 2)))
        : P0;
    const fTail = down ? add(foot, mul(bHat, -fl)) : foot;
    const fHead = down ? foot : add(foot, mul(bHat, fl));
    const fOuter = down ? fTail : fHead; // the end away from the surface
    const outerDir = down ? mul(bHat, -1) : bHat;
    const shadowDir = pDir && down ? mul(pDir, -1) : pDir;

    // Reference ray of the field angles; along k the k arrow itself is the ray
    if (pDir ? angleDeg(refDir, pDir) > TINY_DEG : !isTop) {
        const kShaftFree = o.showWave && !alongK && kTo - (kFrom + kTo) / 2 > rPhi + 0.05;
        const r = pDir ? rPhi : rElev;
        if (o.angleRef === 'edge') layers.guides.push(guideLine(foot, add(foot, mul(xEdge, rayLen(r)))));
        else if (!kShaftFree) layers.guides.push(guideLine(foot, add(foot, mul(kHat, rayLen(r)))));
    }

    // Antenna turn at the tail of k, between the k arrow and the 0° edge direction. Beyond 90° that
    // wedge would lie on the emitter: then it goes to the field, on the side away from it.
    if (showDev) {
        if (devDeg <= 90 + 1e-6) {
            const T = mul(kHat, kFrom);
            layers.guides.push(guideLine(T, add(T, mul(xEdge, rayLen(rDev)))));
            if (!o.showWave) layers.guides.push(guideLine(T, add(T, mul(kHat, rayLen(rDev)))));
            drawArc(xEdge, kHat, rDev, st.symbols.dev, T);
        } else {
            const busy = pDir && !alongK ? [pDir, shadowDir, arc(refDir, pDir, 1).mid] : [];
            const devMid = arc(xEdge, kHat, 1).mid;
            const clearance = (m) => Math.min(180, ...busy.map((b) => angleDeg(m, b)));
            const flip = clearance(mul(devMid, -1)) > clearance(devMid);
            const edgeRay = flip ? mul(xEdge, -1) : xEdge;
            layers.guides.push(guideLine(P0, add(P0, mul(edgeRay, rayLen(rDev)))));
            if (flip || !o.showWave) layers.guides.push(guideLine(P0, add(P0, mul(kHat, flip ? -rayLen(rDev) : rayLen(rDev)))));
            drawArc(edgeRay, flip ? mul(kHat, -1) : kHat, rDev, st.symbols.dev, P0);
        }
    }

    if (pDir) {
        // in-plane angle to k (or to the 0° edge); a tilted field gets its surface projection and vertical drop-line
        if (showElev && !isTop) {
            const projPt = add(foot, mul(shadowDir, fl * Math.cos(elevDeg * DEG)));
            const [ppx, ppy] = toSvg(projPt);
            layers.guides.push(guideLine(foot, projPt, C.angle));
            layers.guides.push(guideLine(fOuter, projPt, C.angle));
            layers.guides.push(`<circle cx="${n2(ppx)}" cy="${n2(ppy)}" r="${n2(Math.max(1.0, L.guide * 1.7))}" fill="${C.angle}"/>`);
            if (down) layers.guides.push(guideLine(foot, add(foot, mul(pDir, rayLen(rPhi))))); // ray for the in-plane angle
            drawArc(shadowDir, outerDir, rElev, st.symbols.theta, foot);
        }
        drawArc(refDir, pDir, rPhi, st.symbols.phi, foot);
    } else if (!isTop) {
        drawArc(refDir, outerDir, rElev, st.symbols.theta, foot); // field along the normal
    }

    // ---- arrows: k first (it is long), the field and m on top
    if (o.showWave) {
        const a = mul(kHat, kFrom);
        const b = mul(kHat, kTo);
        layers.arrows.push(arrow(a, b, C.wave, L.wave));
        wantLabel(5, arrowLabelCands(a, b, symbolWidth(st.symbols.wave)), vectorText(st.symbols.wave), C.wave, symbolWidth(st.symbols.wave));
    }
    // ---- angle sweep: a lightly filled wedge from the drawn field to where it ends, its rim an
    // arc across the shaft with the head at the end; faint field arrows at the end and (optionally)
    // on the way. Below the surface it is fainter and dashed. Seen edge-on (top view, field turning
    // out of the surface) the wedge would be a line: then the turn goes into the field label.
    let sweepNote = '';
    let sweepNoteW = 0;
    if (sweepAngle) {
        const sgn = Math.sign(swSpan);
        const span = Math.abs(swSpan);
        const full = span >= 360 - 1e-6;
        const step = Math.abs(num(o.sweepStep, 0));
        const ax = cam.proj(add(foot, bAt(swFrom)));
        const ay = cam.proj(add(foot, bAt(swFrom + 90)));
        const o2 = cam.proj(foot);
        const flat = Math.abs((ax[0] - o2[0]) * (ay[1] - o2[1]) - (ax[1] - o2[1]) * (ay[0] - o2[0])) < 0.12;
        if (flat) {
            const mode = o.angleLabels === 'value' ? 'both' : o.angleLabels;
            sweepNote = `, ${angleText(sweepSym, span, mode)}`;
            sweepNoteW = textWidth(', ') + angleWidth(sweepSym, span, mode);
        } else {
            const ghosts = [];
            if (step >= 1) for (let a = step; a < span - 1e-6; a += step) ghosts.push(a);
            if (!full) ghosts.push(span);
            ghosts.forEach((a) => {
                const d = bAt(swFrom + sgn * a);
                layers.sweep.push(arrow(foot, add(foot, mul(d, fl)), C.field, L.field * 0.75, { dashed: d[2] < -1e-6, opacity: 0.35, w: 4 }));
            });
            const n = Math.max(12, Math.ceil(span / 3));
            const at = (deg) => add(foot, mul(bAt(swFrom + sgn * deg), rSweep));
            const below = (p) => p[2] < -1e-6;
            // runs of points above / below the surface, interpolated cleanly across z = 0
            const runsOf = (path) => {
                const runs = [];
                let curBelow = below(path[0]);
                let run = [path[0]];
                for (let i = 1; i < path.length; i++) {
                    const b = below(path[i]);
                    if (b !== curBelow) {
                        const p = path[i - 1];
                        const q = path[i];
                        const dz = q[2] - p[2];
                        const t = Math.abs(dz) > 1e-9 ? clamp(-p[2] / dz, 0, 1) : 0.5;
                        const crossPt = add(p, mul(sub(q, p), t));
                        run.push(crossPt);
                        runs.push({ pts: run, below: curBelow });
                        run = [crossPt, q];
                        curBelow = b;
                    } else {
                        run.push(path[i]);
                    }
                }
                runs.push({ pts: run, below: curBelow });
                return runs;
            };
            const wedge = Array.from({ length: n + 1 }, (_, i) => at((span * i) / n));
            const [cx, cy] = toSvg(foot);
            runsOf(wedge).forEach((r) =>
                layers.sweep.push(`<path d="${pathOf([foot, ...r.pts], true)}" fill="${C.field}" fill-opacity="${r.below ? 0.035 : 0.08}" stroke="none"/>`)
            );
            // the rim starts on the shaft; a full turn stops short of it, so its head reads
            const end = full ? span - 16 : span;
            const rim = Array.from({ length: n + 1 }, (_, i) => at((end * i) / n));
            const scr = rim.map(toSvg);
            const id = addLine(scr, W_LINE); // text on the rim would hide where the turn runs
            marks.push({ id, pts: scr });
            thin.add(id);
            const sw = L.arc * 1.05;
            const hlRim = st.arrowHeadPt * 0.72;
            const trimEnd2D = (pts, trimDist) => {
                if (trimDist <= 0 || pts.length < 2) return pts;
                const out = pts.map((p) => [...p]);
                let rem = trimDist;
                while (out.length >= 2 && rem > 1e-6) {
                    const a = out[out.length - 2];
                    const b = out[out.length - 1];
                    const d = Math.hypot(b[0] - a[0], b[1] - a[1]);
                    if (d > rem + 1e-6) {
                        const f = (d - rem) / d;
                        out[out.length - 1] = [a[0] + (b[0] - a[0]) * f, a[1] + (b[1] - a[1]) * f];
                        rem = 0;
                    } else {
                        rem -= d;
                        out.pop();
                    }
                }
                return out;
            };
            const rimRuns = runsOf(rim);
            const tip = scr[scr.length - 1];
            const prev = scr[Math.max(0, scr.length - 3)];
            const dd = Math.hypot(tip[0] - prev[0], tip[1] - prev[1]) || 1;
            const ux = (tip[0] - prev[0]) / dd;
            const uy = (tip[1] - prev[1]) / dd;
            const tipTrim = full ? 0 : Math.min(0.6, dd * 0.35);
            const adjTip = [tip[0] - ux * tipTrim, tip[1] - uy * tipTrim];
            rimRuns.forEach((r, idx) => {
                let sPts = r.pts.map(toSvg);
                if (idx === rimRuns.length - 1) sPts = trimEnd2D(sPts, tipTrim + hlRim * 0.65);
                if (sPts.length >= 2) {
                    layers.sweep.push(
                        `<path d="M${sPts.map((p) => p.map(n2).join(',')).join('L')}" stroke="${C.field}" stroke-width="${n2(sw)}"${r.below ? ` stroke-dasharray="${n2(sw * 2.5)},${n2(sw * 2)}"` : ''} fill="none"/>`
                    );
                }
            });
            layers.sweep.push(headPath(adjTip, ux, uy, hlRim, C.field));
            // the turn just outside the middle of the rim, or further out with a line to it
            const w = angleWidth(sweepSym, span);
            const cands = [];
            for (const f of [0.5, 0.35, 0.65, 0.2, 0.8]) {
                const p = scr[Math.round(f * n)];
                const d = Math.hypot(p[0] - cx, p[1] - cy) || 1;
                for (const [g, pen, leader] of [[0.25, 0, false], [0.7, 0.8, false], [1.8, 2, true], [2.8, 3, true]]) {
                    cands.push({ ...beyond(p[0], p[1], (p[0] - cx) / d, (p[1] - cy) / d, w, g * fs, pen + Math.abs(f - 0.5) * 1.5), from: p, leader });
                }
                // beside the rim: a wedge seen nearly edge-on points its outside at the field tip
                const q = scr[Math.min(n, Math.round(f * n) + 1)];
                const r = scr[Math.max(0, Math.round(f * n) - 1)];
                const t = Math.hypot(q[0] - r[0], q[1] - r[1]) || 1;
                for (const sd of [1, -1]) {
                    cands.push({ ...beyond(p[0], p[1], (-(q[1] - r[1]) / t) * sd, ((q[0] - r[0]) / t) * sd, w, 0.3 * fs, 1 + Math.abs(f - 0.5) * 1.5), from: p });
                }
            }
            wantLabel(3, cands, angleText(sweepSym, span), C.field, w, { id, pts: scr, anchor: scr[Math.round(n / 2)] });
        }
    }

    // field strength sweep: the range in the label; through zero also the reversed field, faint
    let rangeNote = '';
    let rangeNoteW = 0;
    if (fieldRange) {
        const [lo, hi] = fieldRange;
        const fmt = (x) => `${x < 0 ? '−' : ''}${Number(Math.abs(x).toFixed(1))}`;
        const plain = ` = ${fmt(lo)}…${fmt(hi)} mT`;
        rangeNote = latex ? '' : esc(plain);
        rangeNoteW = textWidth(plain);
        if (Math.min(lo, hi) < 0 && Math.max(lo, hi) > 0 && (pDir || !isTop)) {
            const r = mul(bHat, -1);
            const rDown = !sweepAngle && dot(r, nrm) < -1e-6;
            layers.sweep.push(
                rDown
                    ? arrow(add(foot, mul(r, -fl)), foot, C.field, L.field * 0.75, { opacity: 0.35, w: 4 })
                    : arrow(foot, add(foot, mul(r, fl)), C.field, L.field * 0.75, { opacity: 0.35, w: 4 })
            );
        }
    }

    const fieldLabel = latex && fieldRange
        ? esc(`$${symbolTex(st.symbols.field, st.vectorStyle)} = ${fieldRange.map((x) => Number(x.toFixed(1))).join('\\ldots ')}\\,\\mathrm{mT}$`)
        : vectorText(st.symbols.field) + rangeNote;
    const fieldLabelW = symbolWidth(st.symbols.field) + rangeNoteW;
    if (isTop && !pDir) {
        // Top view cannot show the tilt: the field along the normal becomes ⊙ / ⊗
        const [cx, cy] = toSvg(foot);
        const r = fs * 0.62;
        const q = r * 0.6;
        const mark = down
            ? `<path d="M${n2(cx - q)},${n2(cy - q)}L${n2(cx + q)},${n2(cy + q)}M${n2(cx - q)},${n2(cy + q)}L${n2(cx + q)},${n2(cy - q)}" stroke="${C.field}" stroke-width="${n2(L.field * 0.8)}"/>`
            : `<circle cx="${n2(cx)}" cy="${n2(cy)}" r="${n2(L.field * 1.3)}" fill="${C.field}"/>`;
        layers.arrows.push(`<circle cx="${n2(cx)}" cy="${n2(cy)}" r="${n2(r)}" fill="${C.sample}" stroke="${C.field}" stroke-width="${n2(L.field * 0.8)}"/>`, mark);
        addLine(Array.from({ length: 13 }, (_, i) => [cx + r * Math.cos(i * 30 * DEG), cy + r * Math.sin(i * 30 * DEG)]), W_LINE);
        const w = fieldLabelW + sweepNoteW;
        const cands = [[0.71, -0.71], [-0.71, -0.71], [0.71, 0.71], [-0.71, 0.71], [1, 0], [-1, 0]].map(([ux, uy], i) => beyond(cx, cy, ux, uy, w, r + fs * 0.15, i * 0.3));
        wantLabel(4, cands, fieldLabel + sweepNote, C.field, w);
    } else {
        layers.arrows.push(arrow(fTail, fHead, C.field, L.field, { dashed: !down && bn < -1e-6 }));
        if (!down) {
            const [fx, fy] = toSvg(foot);
            layers.arrows.push(`<circle cx="${n2(fx)}" cy="${n2(fy)}" r="${n2(Math.max(1.15, L.field * 0.9))}" fill="${C.field}" stroke="${C.sample}" stroke-width="0.6"/>`);
        }
        // top view: the tilt goes into the label
        const tiltMode = o.angleLabels === 'value' ? 'both' : o.angleLabels;
        const tilt = isTop && showElev ? `, ${angleText(st.symbols.theta, elevDeg, tiltMode)}` : '';
        const w = fieldLabelW + (tilt ? textWidth(', ') + angleWidth(st.symbols.theta, elevDeg, tiltMode) : 0) + sweepNoteW;
        wantLabel(0.5, arrowLabelCands(foot, fOuter, w), fieldLabel + tilt + sweepNote, C.field, w); // the name first: values keep off it
    }
    // m belongs to the current field, not to a turning one
    if (o.showMoment && !sweepAngle && !state.metrics.isEasyPlaneSingular) {
        const mHat = unit(v.mHat.Sub);
        const mDown = dot(mHat, nrm) < -1e-6;
        const mOuter = add(foot, mul(mHat, mDown ? -MOMENT_LEN : MOMENT_LEN));
        layers.arrows.push(mDown ? arrow(mOuter, foot, C.moment, L.wave) : arrow(foot, mOuter, C.moment, L.wave));
        const w = symbolWidth(st.symbols.moment);
        wantLabel(6, arrowLabelCands(foot, mOuter, w), vectorText(st.symbols.moment), C.moment, w);
    }

    // ---- crystal directions of the sample edges and the normal, in a free corner
    let core = boundsOf(ink);
    let triadBox = null; // the other labels keep off it
    if (o.showAxes) {
        const mb = state.millerBravais;
        const axLen = fs * 3;
        const gap = fs * 0.25;
        const ring = fs * 0.36;
        // normal as a direction where [hkil] ∥ (hkil) (c, a, m); otherwise (r-cut) the exact plane
        const norm4 = (p) => p.map((x) => x / Math.max(...p.map(Math.abs)));
        const sameIdx = (p, q) => [1, -1].some((s) => norm4(p).every((x, i) => Math.abs(x - s * norm4(q)[i]) < 1e-3));
        const nLabel = sameIdx(mb.zsDir.raw, mb.zsPlane.raw) ? mb.zsDir : mb.zsPlane;
        const axes = [
            { d: nrm, ...millerLabel(nLabel.label, nLabel.latex) },
            { d: xEdge, ...millerLabel(mb.xs.label, mb.xs.latex) },
            { d: [0, 1, 0], ...millerLabel(mb.ys.label, mb.ys.latex) }
        ];
        const cHat = unit(v.ecHat.Sub);
        if (o.showCAxis && axes.every((a) => Math.abs(dot(a.d, cHat)) < 0.999)) {
            const up = dot(cHat, nrm) >= 0;
            axes.push({ d: up ? cHat : mul(cHat, -1), ...millerLabel(up ? '[0001]' : `[0001${OVERLINE}]`, up ? '[0001]' : '[000\\bar{1}]'), dashed: true });
        }
        const o0 = cam.proj(origin);
        axes.forEach((a) => {
            const p = cam.proj(a.d);
            const dx = p[0] - o0[0];
            const dy = -(p[1] - o0[1]);
            const d = Math.hypot(dx, dy);
            a.toward = d < 0.15; // axis towards the viewer (top view): dot in a ring
            a.u = d ? [dx / d, dy / d] : [0, 0];
            a.ex = a.u[0] * axLen * Math.min(1, d);
            a.ey = a.u[1] * axLen * Math.min(1, d);
            // label on the extension of its axis; the ring's label down left
            a.cands = (ox, oy) =>
                a.toward
                    ? [[-0.71, 0.71], [0, 1], [-1, 0], [0.71, 0.71]].map(([ux, uy], i) => beyond(ox, oy, ux, uy, a.width, ring + fs * 0.3, i))
                    : [[0, 0], [30, 1], [-30, 1], [60, 2], [-60, 2], [90, 3], [-90, 3], [120, 4], [-120, 4]].map(([deg, pen]) => {
                          const [ux, uy] = turn(a.u, deg);
                          return beyond(ox + a.ex, oy + a.ey, ux, uy, a.width, gap, pen);
                      });
        });

        // Label spots chosen together (axes close in direction would push each other aside one by
        // one): no two labels touching, none on an axis, each as close to its own axis end as possible
        const axSegs = axes.filter((a) => !a.toward).map((a) => [[0, 0], [a.ex, a.ey]]);
        const opts = axes.map((a) => a.cands(0, 0).map((cand) => ({ cand, box: labelBox(cand, a.width) })));
        let pick = opts.map(() => 0);
        let pickCost = Infinity;
        const choose = (i, sel, cost) => {
            if (cost >= pickCost) return;
            if (i === axes.length) {
                pickCost = cost;
                pick = [...sel];
                return;
            }
            opts[i].forEach(({ cand, box }, j) => {
                const e = grow(box, fs * 0.15);
                let add = cand.pen;
                sel.forEach((k, m) => boxesTouch(e, opts[m][k].box) && (add += 100));
                axSegs.forEach(([p, q]) => segHitsBox(p, q, e) && (add += 12));
                sel.push(j);
                choose(i + 1, sel, cost + add);
                sel.pop();
            });
        };
        choose(0, [], 0);

        // Footprint relative to the triad origin: axes plus the chosen labels
        let box = { l: -ring, r: ring, t: -ring, b: ring };
        axes.forEach((a, i) => {
            box = unionBox(box, { l: a.ex, r: a.ex, t: a.ey, b: a.ey });
            box = unionBox(box, opts[i][pick[i]].box);
        });

        // The triad is drawn smaller than the rest (k), so it fits into a free corner of the drawing
        const k = latex ? 1 : TRIAD_SCALE; // LaTeX sets its own text size, the labels would outgrow the box
        box = { l: box.l * k, r: box.r * k, t: box.t * k, b: box.b * k };

        // Cheapest free spot: not on the sample or any arrow, smallest overall figure, then bottom left
        const hull = convexHull(slab.map(toSvg));
        const step = fs * 0.4;
        let best = null;
        for (let i = -40; i <= 40; i++) {
            for (let j = -40; j <= 40; j++) {
                const ox = core.l - box.l + i * step;
                const oy = core.b - box.b - j * step;
                const r = { l: ox + box.l, r: ox + box.r, t: oy + box.t, b: oy + box.b };
                const e = grow(r, fs * 0.35);
                if (convexOverlap(boxPoly(e), hull)) continue;
                if (lines.some(([a, z]) => segHitsBox(a, z, e))) continue;
                const u = unionBox(core, r);
                const cost = (u.r - u.l) * (u.b - u.t) + (r.l - core.l + core.b - r.b) * fs * 0.05;
                if (!best || cost < best.cost) best = { cost, ox, oy };
            }
        }
        const { ox, oy } = best || { ox: core.l - box.l, oy: core.b - box.b + fs };

        // Drawn around (0, 0) at full size and scaled as one group; lines keep their stroke width
        const g = [];
        const sw = L.axes / k;
        const at = ([x, y]) => [ox + x * k, oy + y * k];
        axes.forEach((a, i) => {
            if (a.toward) {
                g.push(`<circle cx="0" cy="0" r="${n2(ring)}" fill="none" stroke="${C.axes}" stroke-width="${n2(sw)}"/>`);
                g.push(`<circle cx="0" cy="0" r="${n2(sw * 1.2)}" fill="${C.axes}"/>`);
            } else {
                const [ux, uy] = a.u;
                const hl = fs * 0.52;
                addLine([at([0, 0]), at([a.ex, a.ey])], W_LINE);
                g.push(
                    `<path d="M0,0L${n2(a.ex - ux * hl * 0.7)},${n2(a.ey - uy * hl * 0.7)}" stroke="${C.axes}" stroke-width="${n2(sw)}" fill="none"${a.dashed ? ` stroke-dasharray="${n2(sw * 3)},${n2(sw * 2)}"` : ''}/>`,
                    headPath([a.ex, a.ey], ux, uy, hl, C.axes)
                );
            }
            const b = opts[i][pick[i]].box;
            g.push(typeof a.content === 'function' ? a.content(b.l, b.y) : `<text x="${n2(b.x)}" y="${n2(b.y)}" text-anchor="${b.anchor}" fill="${C.text}" font-size="${fs}">${a.content}</text>`);
        });
        if (axes.some((a) => !a.toward)) g.push(`<circle cx="0" cy="0" r="${n2(sw * 1.05)}" fill="${C.axes}"/>`);
        layers.axes.push(`<g transform="translate(${n2(ox)},${n2(oy)}) scale(${n2(k)})">${g.join('')}</g>`);
        triadBox = { l: ox + box.l, r: ox + box.r, t: oy + box.t, b: oy + box.b };
        ink.push([triadBox.l, triadBox.t], [triadBox.r, triadBox.b]);
        core = boundsOf(ink);
    }

    const placed = placeLabels(core, triadBox);

    // ---- tight crop around everything drawn; the panel letter in its top left corner
    const frame = grow(placed.reduce(unionBox, core), fs * 0.3);
    const panel = String(o.panel || '').trim();
    if (panel) {
        const w = textWidth(panel) * 1.06;
        const at = (t) => ({ l: frame.l, r: frame.l + w, t, b: t + fs * 0.96 });
        let b = at(frame.t);
        const taken = (box) => ink.some(([x, y]) => x > box.l && x < box.r && y > box.t && y < box.b) || placed.some((p) => boxesTouch(box, p));
        if (taken(grow(b, fs * 0.3))) {
            frame.t -= fs * 1.2;
            b = at(frame.t);
        }
        layers.labels.push(`<text x="${n2(b.l)}" y="${n2(b.t + fs * 0.76)}" fill="${C.text}" font-size="${fs}" font-weight="bold">${esc(panel)}</text>`);
    }
    const W = n2(frame.r - frame.l);
    const H = n2(frame.b - frame.t);
    const X = n2(frame.l);
    const Y = n2(frame.t);

    const s = state.setup;
    const sweepDesc = sweepAngle
        ? ` · Sweep ${sweepSym} (φ_M=${fmtAngle(swFrom)}…${fmtAngle(swFrom + swSpan)})`
        : fieldRange
          ? ` · Sweep B_0=${fieldRange.join('…')} mT`
          : '';
    const desc = `${s.name || s.cutPreset} · ${s.deviceType} θ_dev=${fmtAngle(s.thetaDevDeg)} · Einbau ${s.mountMode} ψ_S=${fmtAngle(s.psiSDeg)} · φ_M=${fmtAngle(sweepAngle ? swFrom : s.phiMDeg)}${sweepDesc}`;
    const layer = (id, label, items) =>
        items.length ? `  <g id="${id}" inkscape:groupmode="layer" inkscape:label="${label}">\n    ${items.join('\n    ')}\n  </g>\n` : '';
    const bg = C.background && C.background !== 'none' ? `  <rect id="background" x="${X}" y="${Y}" width="${W}" height="${H}" fill="${C.background}"/>\n` : '';

    return (
        `<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n` +
        `<svg xmlns="http://www.w3.org/2000/svg" xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape" ` +
        `width="${n2(W / PT_PER_MM)}mm" height="${n2(H / PT_PER_MM)}mm" viewBox="${X} ${Y} ${W} ${H}" ` +
        `font-family="${esc(st.fontFamily)}" stroke-linecap="round">\n` +
        `  <title>Probe und Magnetfeld</title>\n  <desc>${esc(desc)}</desc>\n` +
        bg +
        layer('sample', 'Probe', layers.sample) +
        layer('antenna', 'Antenne', layers.antenna) +
        layer('sweep', 'Sweep', layers.sweep) +
        layer('guides', 'Hilfslinien', layers.guides) +
        layer('arcs', 'Winkel', layers.arcs) +
        layer('arrows', 'Pfeile', layers.arrows) +
        layer('axes', 'Kristallachsen', layers.axes) +
        layer('labels', 'Beschriftung', layers.labels) +
        `</svg>\n`
    );
}

/**
 * Caption for the figure (LaTeX), with the angles as drawn: what is on it, the wave and field
 * directions, the sweep. The figure options and style decide what is mentioned and the symbols.
 */
export function figureCaption(state, options = {}, styleIn = null) {
    const o = { ...DEFAULT_FIGURE_OPTIONS, ...options };
    const st = completeFigureStyle(styleIn);
    const sym = (k) => symbolTex(st.symbols[k], 'none');
    const vec = (k) => symbolTex(st.symbols[k], st.vectorStyle);
    const deg = (x) => `${Number(x.toFixed(1))}^\\circ`;
    const mT = (x) => `${Number(x.toFixed(1))}\\,\\mathrm{mT}`;
    const s = state.setup;
    const kHat = state.vectors.kHat.Sub;
    const { swSpan, sweepAngle, isOopSweep, b0mT, fieldRange, bHat } = fieldSetup(state, o);
    const sweepTex = symbolTex(resolveSweepSymbol(st, isOopSweep), 'none');
    const xEdge = [1, 0, 0];
    const edge = `$${state.millerBravais.xs.latex}$`;
    const parts = [];

    const cut = s.cutPreset === 'custom' ? 'eigener Schnitt' : s.cutPreset;
    const name = { IDT: 'IDT', CPW: 'Antenne' };
    const own = s.deviceType === 'IDT' ? 'IDT' : 'CPW';
    const devs = o.showAntenna ? (DEVICE_SETS[o.devices] || DEVICE_SETS.pair)(own) : [];
    const two = { IDT: 'zwei IDTs', CPW: 'zwei Antennen' };
    const devText =
        devs.length > 1
            ? devs[0] === devs[1]
                ? ` mit ${two[devs[0]]} als Sender und Empfänger`
                : ` mit Sender (${name[devs[0]]}) und Empfänger (${name[devs[1]]})`
            : devs.length
              ? ` mit Sender (${name[devs[0]]})`
              : '';
    parts.push(`Probe (${cut}, Einbau ${s.mountMode})${devText}.`);

    const devDeg = angleDeg(xEdge, kHat);
    if (o.showWave) {
        const from = devs.length > 1 ? 'läuft vom Sender zum Empfänger' : devs.length ? 'läuft vom Sender weg' : 'läuft';
        const dir = devDeg > TINY_DEG ? `, um $${sym('dev')} = ${deg(devDeg)}$ gegen die 0°-Kante ${edge} gedreht` : ` entlang der 0°-Kante ${edge}`;
        parts.push(`Die Welle $${vec('wave')}$ ${from}${dir}.`);
    }

    const strength = fieldRange ? `$B_0 = ${fieldRange.map((x) => Number(x.toFixed(1))).join('\\ldots ')}\\,\\mathrm{mT}$` : `$B_0 = ${mT(b0mT)}$`;
    const bn = bHat[2];
    const bPlane = [bHat[0], bHat[1], 0];
    const elev = Math.asin(Math.min(1, Math.abs(bn))) / DEG;
    let field;
    if (len(bPlane) < 1e-6) field = `steht senkrecht ${bn < 0 ? 'in die Probe hinein' : 'auf der Probe'}`;
    else {
        const ref = o.angleRef === 'edge' ? xEdge : kHat;
        const phi = angleDeg(ref, unit(bPlane));
        const refText = o.angleRef === 'edge' ? 'zur 0°-Kante' : `zu $${vec('wave')}$`;
        field = phi > TINY_DEG ? `liegt in der Ebene unter $${sym('phi')} = ${deg(phi)}$ ${refText}` : `zeigt in der Ebene entlang ${o.angleRef === 'edge' ? 'der 0°-Kante' : `$${vec('wave')}$`}`;
        if (elev > TINY_DEG) field += ` und ist um $${sym('theta')} = ${deg(elev)}$ ${bn < 0 ? 'unter die' : 'aus der'} Probenebene gekippt`;
    }
    parts.push(`Das Feld $${vec('field')}$ (${strength}${sweepAngle ? ', zu Beginn des Sweeps' : ''}) ${field}.`);

    if (sweepAngle) {
        const span = Math.abs(swSpan);
        parts.push(
            span >= 360 - 1e-6
                ? `Der Bogen zeigt eine volle Drehung des Feldes in der Drehebene des Magneten, die Pfeilspitze den Drehsinn.`
                : `Der Bogen zeigt, wie sich das Feld in der Drehebene des Magneten um $${sweepTex} = ${deg(span)}$ dreht; die Pfeilspitze zeigt den Drehsinn.`
        );
    } else if (fieldRange && Math.min(...fieldRange) < 0 && Math.max(...fieldRange) > 0) {
        parts.push('Der blasse Pfeil zeigt die umgekehrte Feldrichtung.');
    }
    if (o.showMoment && !sweepAngle && !state.metrics.isEasyPlaneSingular) parts.push(`$${vec('moment')}$ zeigt die Gleichgewichtslage.`);
    if (o.showAxes) parts.push('Daneben die Kristallrichtungen der Probenkanten und der Normalen (Miller-Bravais).');
    return parts.join(' ');
}
