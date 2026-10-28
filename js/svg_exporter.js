/**
 * Figure export: one uniform schematic of the sample and the magnetic field.
 *
 * The sample always sits in the same pose and at the same scale (fixed camera in the sample frame
 * F_Sub), so figures made from different setups line up in a thesis. Drawn: slab, the two antennas
 * turned by theta_dev, the wave vector, the field arrow with its in-plane and out-of-plane angle,
 * and a small triad with the Miller-Bravais indices of the sample edges and the surface normal.
 * The figure is cropped tightly around what is drawn.
 *
 * Geometry lives in sample units (top surface at z = 0); the SVG viewBox is in pt, so font sizes
 * and line widths in the style are print sizes when the figure is placed at its natural size.
 */

const DEG = Math.PI / 180;
const PT_PER_MM = 72 / 25.4;

// Camera in F_Sub. Oblique: low view onto the front-left corner; top: looking down the normal.
const VIEWS = {
    oblique: { azDeg: -66, elDeg: 28 },
    top: { azDeg: -90, elDeg: 90 }
};

const HX = 1.5; // half edges of the sample: long along the 0° edge x_s
const HY = 0.8;
const FIELD_LEN = 0.86;
const ANT_IN = 0.92; // antennas sit at |k| = ANT_IN … ANT_OUT, |t| <= ANT_T (shrunk to fit the sample)
const ANT_OUT = 1.28;
const ANT_T = 0.4;
const MOMENT_LEN = 0.8;
const R_DEV = 0.3; // arc radii: antenna angle, in-plane field angle, elevation
const R_PHI = 0.44;
const R_ELEV = 0.64;
const TINY_DEG = 0.05;

// ---------------------------------------------------------------------------
// Styles: everything that should look the same in every figure of a thesis
// ---------------------------------------------------------------------------

const SYMBOLS = { field: 'B', wave: 'k', moment: 'm', phi: 'φ', theta: 'θ', dev: 'α' };

export const BUILTIN_FIGURE_STYLES = {
    Standard: {
        widthMm: 65,
        fontFamily: 'Arial, Helvetica, sans-serif',
        fontSizePt: 9,
        vectorStyle: 'bold',
        sampleThickness: 0.08,
        arrowHeadPt: 6,
        colors: {
            sample: '#6cb8b4',
            antenna: '#d9a82e',
            field: '#111111',
            wave: '#c2410c',
            angle: '#1f2937',
            axes: '#374151',
            moment: '#7c3aed',
            text: '#111111',
            background: 'none'
        },
        lines: { outline: 0.5, field: 1.4, wave: 1.1, arc: 0.7, guide: 0.6, axes: 0.8 },
        symbols: { ...SYMBOLS }
    },
    'Schwarz-Weiß': {
        widthMm: 65,
        fontFamily: 'Times New Roman, Times, serif',
        fontSizePt: 10,
        vectorStyle: 'bold',
        sampleThickness: 0.08,
        arrowHeadPt: 6,
        colors: {
            sample: '#d4d4d4',
            antenna: '#7a7a7a',
            field: '#000000',
            wave: '#000000',
            angle: '#000000',
            axes: '#000000',
            moment: '#555555',
            text: '#000000',
            background: 'none'
        },
        lines: { outline: 0.5, field: 1.3, wave: 0.9, arc: 0.6, guide: 0.5, axes: 0.7 },
        symbols: { ...SYMBOLS }
    },
    Folie: {
        widthMm: 90,
        fontFamily: 'Arial, Helvetica, sans-serif',
        fontSizePt: 14,
        vectorStyle: 'bold',
        sampleThickness: 0.1,
        arrowHeadPt: 10,
        colors: {
            sample: '#4fa3a5',
            antenna: '#e0a526',
            field: '#1d3a8a',
            wave: '#d9480f',
            angle: '#111827',
            axes: '#374151',
            moment: '#7c3aed',
            text: '#111827',
            background: '#ffffff'
        },
        lines: { outline: 0.8, field: 2.6, wave: 2, arc: 1.2, guide: 1, axes: 1.3 },
        symbols: { ...SYMBOLS }
    }
};

export const DEFAULT_FIGURE_OPTIONS = {
    view: 'oblique', // 'oblique' | 'top'
    angleRef: 'wave', // in-plane field angle measured from k ('wave') or from the 0° edge ('edge')
    angleLabels: 'value', // 'value' | 'symbol' | 'both'
    showAntenna: true,
    showWave: true,
    showAxes: true,
    showCAxis: true,
    showMoment: false,
    latex: false
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
        symbols: { ...base.symbols, ...(s.symbols || {}) }
    };
}

// ---------------------------------------------------------------------------
// Small vector + colour helpers
// ---------------------------------------------------------------------------

const add = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const mul = (a, s) => [a[0] * s, a[1] * s, a[2] * s];
const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const len = (a) => Math.hypot(a[0], a[1], a[2]);
const unit = (a) => mul(a, 1 / (len(a) || 1));
const angleDeg = (a, b) => Math.acos(Math.max(-1, Math.min(1, dot(unit(a), unit(b))))) / DEG;

function shade(hex, f) {
    const m = /^#?([0-9a-f]{6})$/i.exec(hex || '');
    if (!m) return hex;
    const n = parseInt(m[1], 16);
    const ch = [(n >> 16) & 255, (n >> 8) & 255, n & 255].map((c) => Math.round(f < 0 ? c * (1 + f) : c + (255 - c) * f));
    return '#' + ch.map((c) => c.toString(16).padStart(2, '0')).join('');
}

const esc = (s) => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;');
const n2 = (v) => Number(v.toFixed(2));
const fmtAngle = (deg) => `${Number(deg.toFixed(1))}°`;
const OVERLINE = '̅';

// Advance widths in em (Arial / Times metrics), good enough to place labels and draw index bars
const EM_SANS = { digit: 0.556, other: 0.56, '[': 0.278, ']': 0.278, '(': 0.333, ')': 0.333, '°': 0.4, '.': 0.278, ',': 0.278, ' ': 0.278, '=': 0.584, '-': 0.333, '−': 0.584 };
const EM_SERIF = { digit: 0.5, other: 0.5, '[': 0.333, ']': 0.333, '(': 0.333, ')': 0.333, '°': 0.4, '.': 0.25, ',': 0.25, ' ': 0.25, '=': 0.564, '-': 0.333, '−': 0.564 };
const isSerif = (family) => /serif|times|georgia|palatino|garamond|cmu|computer modern|latin modern/i.test(family) && !/sans/i.test(family);

// ---------------------------------------------------------------------------
// Labels: "φ_B", "θ_{dev}" -> italic symbol + subscript; optional LaTeX output
// ---------------------------------------------------------------------------

const GREEK_TEX = {
    α: '\\alpha', β: '\\beta', γ: '\\gamma', δ: '\\delta', ε: '\\varepsilon', ζ: '\\zeta', η: '\\eta', θ: '\\theta',
    ϑ: '\\vartheta', κ: '\\kappa', λ: '\\lambda', μ: '\\mu', ν: '\\nu', ξ: '\\xi', π: '\\pi', ρ: '\\rho', σ: '\\sigma',
    τ: '\\tau', φ: '\\varphi', ϕ: '\\phi', χ: '\\chi', ψ: '\\psi', ω: '\\omega', Θ: '\\Theta', Φ: '\\Phi', Ψ: '\\Psi', Ω: '\\Omega'
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

/** SVG tspans for a symbol; vectorStyle only for vectors (B, k, m), plain italic otherwise. */
function symbolSpans(str, fs, vectorStyle) {
    const { main, sub: s } = parseSymbol(str);
    const weight = vectorStyle === 'bold' ? ' font-weight="bold"' : '';
    let out = `<tspan font-style="italic"${weight}>${esc(main + (VECTOR_MARK[vectorStyle] || ''))}</tspan>`;
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
    const { azDeg, elDeg } = VIEWS[view] || VIEWS.oblique;
    const az = azDeg * DEG;
    const el = elDeg * DEG;
    const right = [-Math.sin(az), Math.cos(az), 0];
    const up = [-Math.sin(el) * Math.cos(az), -Math.sin(el) * Math.sin(az), Math.cos(el)];
    const toViewer = [Math.cos(el) * Math.cos(az), Math.cos(el) * Math.sin(az), Math.sin(el)];
    return { right, up, toViewer, proj: (p) => [dot(p, right), dot(p, up)] };
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
export function generateFigureSVG(state, options = {}, styleIn = null) {
    const o = { ...DEFAULT_FIGURE_OPTIONS, ...options };
    const st = completeFigureStyle(styleIn);
    const C = st.colors;
    const L = st.lines;
    const fs = Number(st.fontSizePt) || 9;
    const latex = !!o.latex;
    const cam = makeCamera(o.view);
    const isTop = o.view === 'top';
    const thick = isTop ? 0 : Math.max(0, Number(st.sampleThickness) || 0);
    const EM = isSerif(st.fontFamily) ? EM_SERIF : EM_SANS;

    const v = state.vectors;
    const kHat = v.kHat.Sub;
    const tHat = v.tHat.Sub;
    const bHat = unit(v.bHat.Sub);
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

    const layers = { sample: [], antenna: [], guides: [], arcs: [], arrows: [], labels: [], axes: [] };

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
    };
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
                if (bar) bars += `M${n2(cx - fs * 0.2)},${n2(y - fs * 0.84)}H${n2(cx + fs * 0.2)}`;
                x += w;
            });
            const text = `<text y="${n2(y)}" text-anchor="middle" fill="${C.text}" font-size="${fs}">${spans}</text>`;
            return bars ? `<g>${text}<path d="${bars}" stroke="${C.text}" stroke-width="${n2(fs * 0.065)}" stroke-linecap="butt" fill="none"/></g>` : text;
        };
        return { content, width };
    };

    // Label box just beyond a screen point, pushed away along a screen direction (unit, y down)
    const labelBox = (c, w) => {
        const x = c.sx + c.dx * c.gap;
        const y = c.sy + c.dy * c.gap + fs * 0.35 + c.dy * fs * 0.35;
        const anchor = c.dx > 0.35 ? 'start' : c.dx < -0.35 ? 'end' : 'middle';
        const l = anchor === 'start' ? x : anchor === 'end' ? x - w : x - w / 2;
        return { l, r: l + w, t: y - fs * 0.76, b: y + fs * 0.2, x, y, anchor };
    };

    // Labels are collected and placed at the end, most important first. Each takes the cheapest of
    // its candidate spots: no other label, no line through it, little antenna or sample edge under
    // it, and not far outside the drawing (that would make the figure bigger)
    const labelReqs = [];
    const wantLabel = (prio, cands, content, color, w) => labelReqs.push({ prio, cands, content, color, w });
    const around = (sx, sy, dx, dy, gap, turns) =>
        turns.map((t) => {
            const c = Math.cos(t * DEG);
            const s = Math.sin(t * DEG);
            return { sx, sy, dx: dx * c - dy * s, dy: dx * s + dy * c, gap };
        });
    const placeLabels = (core) => {
        const placed = [];
        const m = fs * 0.15;
        labelReqs
            .sort((a, b) => a.prio - b.prio)
            .forEach((q) => {
                let best = null;
                let bestCost = Infinity;
                q.cands.forEach((c, i) => {
                    const b = labelBox(c, q.w);
                    const e = grow(b, m);
                    let cost = i * 0.02;
                    if (placed.some((p) => boxesTouch(e, p))) cost += 100;
                    const hit = new Map();
                    for (const [a, z, w, id] of lines) if (!hit.has(id) && segHitsBox(a, z, e)) hit.set(id, w);
                    hit.forEach((w) => (cost += w));
                    // a label grazing an antenna is fine, one covering it is not
                    const boxArea = (e.r - e.l) * (e.b - e.t);
                    for (const [hull, w] of areas) {
                        const a = areaInBox(hull, e);
                        if (a > 0) cost += w * (0.2 + (3 * a) / boxArea);
                    }
                    const over = Math.max(0, core.l - b.l) + Math.max(0, b.r - core.r) + Math.max(0, core.t - b.t) + Math.max(0, b.b - core.b);
                    cost += (1.5 * over) / fs;
                    if (cost < bestCost) {
                        bestCost = cost;
                        best = b;
                    }
                });
                placed.push(best);
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

    // Arrow in 3D, head drawn in screen space so it keeps its size when foreshortened
    const arrow = (from, to, color, width) => {
        const a = toSvg(from);
        const b = toSvg(to);
        const dx = b[0] - a[0];
        const dy = b[1] - a[1];
        const d = Math.hypot(dx, dy);
        if (d < 0.5) return '';
        const ux = dx / d;
        const uy = dy / d;
        const hl = Math.min(st.arrowHeadPt * (0.6 + width / 3.5), d * 0.6);
        const hw = hl * 0.42;
        const bx = b[0] - ux * hl * 0.82;
        const by = b[1] - uy * hl * 0.82;
        const h1 = [b[0] - ux * hl - uy * hw, b[1] - uy * hl + ux * hw];
        const h2 = [b[0] - ux * hl + uy * hw, b[1] - uy * hl - ux * hw];
        addLine([a, b], W_LINE);
        addLine([b, h1, h2], W_LINE, true);
        return (
            `<path d="M${n2(a[0])},${n2(a[1])}L${n2(bx)},${n2(by)}" stroke="${color}" stroke-width="${n2(width)}" fill="none"/>` +
            `<path d="M${n2(b[0])},${n2(b[1])}L${n2(h1[0])},${n2(h1[1])}L${n2(h2[0])},${n2(h2[1])}Z" fill="${color}"/>`
        );
    };

    const guideLine = (from, to, color = C.angle) => {
        addLine([toSvg(from), toSvg(to)], W_GUIDE);
        return `<path d="${pathOf([from, to])}" stroke="${color}" stroke-width="${n2(L.guide)}" stroke-dasharray="${n2(L.guide * 1.5)},${n2(L.guide * 2.6)}" fill="none"/>`;
    };

    // Arc from direction a to direction b (3D unit vectors) around the origin
    const arc = (a, b, r) => {
        const ang = angleDeg(a, b);
        const ua = unit(a);
        let perp = sub(unit(b), mul(ua, dot(ua, unit(b))));
        if (len(perp) < 1e-9) perp = unit([-ua[1], ua[0], 0]);
        perp = unit(perp);
        const dirAt = (f) => add(mul(ua, Math.cos(ang * f * DEG)), mul(perp, Math.sin(ang * f * DEG)));
        const steps = Math.max(8, Math.ceil(ang / 3));
        const pts = Array.from({ length: steps + 1 }, (_, i) => mul(dirAt(i / steps), r));
        return { ang, pts, mid: dirAt(0.5), dirAt, r };
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

    // Arc with its value outside (or inside) the arc, wherever it is free
    let arcCount = 0;
    const drawArc = (a, b, r, symbol, color = C.angle) => {
        const A = arc(a, b, r);
        if (A.ang < TINY_DEG) return;
        addLine(A.pts.map(toSvg), W_GUIDE);
        layers.arcs.push(`<path d="${pathOf(A.pts)}" stroke="${color}" stroke-width="${n2(L.arc)}" fill="none"/>`);
        const cands = [];
        for (const gap of [fs * 0.3, fs * 0.9, fs * 1.6, fs * 2.4, fs * 3.2, -fs * 0.3]) {
            for (const f of [0.5, 0.35, 0.65, 0.2, 0.8]) {
                const d = A.dirAt(f);
                const [sx, sy] = toSvg(mul(d, A.r));
                const [dx, dy] = screenDir(origin, d);
                if (gap < 0) cands.push({ sx, sy, dx: -dx, dy: -dy, gap: -gap });
                else cands.push(...around(sx, sy, dx, dy, gap, [0, 35, -35]));
            }
        }
        wantLabel(10 + arcCount++, cands, angleText(symbol, A.ang), C.text, angleWidth(symbol, A.ang));
    };

    // ---- sample slab: visible side faces (auto-shaded), then the top face
    const top = corners(0);
    const sampleStroke = shade(C.sample, -0.45);
    if (thick > 0) {
        const bot = corners(-thick);
        const sides = [
            { n: [0, -1, 0], q: [top[0], top[1], bot[1], bot[0]] },
            { n: [1, 0, 0], q: [top[1], top[2], bot[2], bot[1]] },
            { n: [0, 1, 0], q: [top[2], top[3], bot[3], bot[2]] },
            { n: [-1, 0, 0], q: [top[3], top[0], bot[0], bot[3]] }
        ];
        sides.forEach((f, i) => {
            if (dot(f.n, cam.toViewer) <= 1e-6) return;
            addLine(f.q.map(toSvg), W_EDGE, true);
            layers.sample.push(`<path d="${pathOf(f.q, true)}" fill="${shade(C.sample, i % 2 ? -0.28 : -0.16)}" stroke="${sampleStroke}" stroke-width="${n2(L.outline)}" stroke-linejoin="round"/>`);
        });
    }
    addLine(top.map(toSvg), W_EDGE, true);
    layers.sample.push(`<path d="${pathOf(top, true)}" fill="${C.sample}" stroke="${sampleStroke}" stroke-width="${n2(L.outline)}" stroke-linejoin="round"/>`);

    // ---- antennas in the wave frame (k, t): emitter and receiver either side of the centre,
    // shrunk as a whole when the turned pair would not fit on the sample
    const wave = (kk, tt) => add(mul(kHat, kk), mul(tHat, tt));
    const fit = Math.min(
        1,
        ...[ANT_T, -ANT_T].map((tt) => {
            const p = wave(ANT_OUT, tt);
            return Math.min((HX - 0.1) / Math.max(Math.abs(p[0]), 1e-9), (HY - 0.1) / Math.max(Math.abs(p[1]), 1e-9));
        })
    );
    const antIn = ANT_IN * fit;
    const antOut = ANT_OUT * fit;
    const antT = ANT_T * fit;
    const waveLen = Math.min(0.68, antIn * 0.85);
    if (o.showAntenna) {
        const polys = [];
        const isIdt = state.setup.deviceType === 'IDT';
        const toWave = (list) => list.map(([kk, tt]) => wave(kk, tt));
        for (const side of [-1, 1]) {
            const k0 = Math.min(antIn * side, antOut * side);
            const k1 = Math.max(antIn * side, antOut * side);
            const sidePolys = [];
            if (isIdt) {
                // two interleaved combs: bus bar along k, fingers along t
                const bar = 0.05 * fit;
                const gapT = 0.06 * fit;
                const fingers = 6;
                const pitch = (k1 - k0) / fingers;
                const fw = pitch * 0.5;
                for (const sgn of [1, -1]) {
                    const T = (t) => sgn * t;
                    const tip = -(antT - bar - gapT);
                    const pts = [[k0, T(antT)], [k1, T(antT)], [k1, T(antT - bar)]];
                    for (let i = fingers - 1; i >= 0; i--) {
                        if ((i % 2 === 0) !== (sgn > 0)) continue;
                        const kc = k0 + pitch * (i + 0.5);
                        pts.push([kc + fw / 2, T(antT - bar)], [kc + fw / 2, T(tip)], [kc - fw / 2, T(tip)], [kc - fw / 2, T(antT - bar)]);
                    }
                    pts.push([k0, T(antT - bar)]);
                    sidePolys.push(toWave(pts));
                }
            } else {
                // coplanar waveguide: ground - signal - ground along t, widening into probe pads at +t
                const kc = (k0 + k1) / 2;
                const s = 0.022 * fit;
                const g = 0.024 * fit;
                const G = 0.095 * fit;
                const ps = 0.038 * fit;
                const pg = 0.03 * fit;
                const tPad = antT - 0.12 * fit;
                const tTaper = tPad - 0.06 * fit;
                sidePolys.push(toWave([[kc - s, -antT], [kc + s, -antT], [kc + s, tTaper], [kc + ps, tPad], [kc + ps, antT], [kc - ps, antT], [kc - ps, tPad], [kc - s, tTaper]]));
                for (const sg of [1, -1]) {
                    const inL = s + g;
                    const inP = ps + pg;
                    const out = s + g + G;
                    sidePolys.push(toWave([[kc + sg * inL, -antT], [kc + sg * out, -antT], [kc + sg * out, antT], [kc + sg * inP, antT], [kc + sg * inP, tPad], [kc + sg * inL, tTaper]]));
                }
            }
            const screen = sidePolys.flat().map(toSvg);
            ink.push(...screen);
            areas.push([convexHull(screen), W_ANTENNA]);
            polys.push(...sidePolys);
        }
        const stroke = shade(C.antenna, -0.4);
        layers.antenna.push(
            `<path d="${polys.map((q) => pathOf(q, true)).join('')}" fill="${C.antenna}" stroke="${stroke}" stroke-width="${n2(L.outline * 0.6)}" stroke-linejoin="round"/>`
        );
    }

    // ---- angles: antenna turn (0° edge -> k), field in plane, field out of plane
    const devDeg = angleDeg(xEdge, kHat);
    const showDev = devDeg > TINY_DEG && (o.showAntenna || o.showWave);
    const bn = dot(bHat, nrm);
    const bPlane = [bHat[0], bHat[1], 0];
    const pDir = len(bPlane) > 1e-6 ? unit(bPlane) : null;
    const elevDeg = Math.asin(Math.min(1, Math.abs(bn))) / DEG;
    const showElev = elevDeg > TINY_DEG;
    const refDir = o.angleRef === 'edge' ? xEdge : kHat;
    const rayLen = (r) => r + 0.14; // short reference rays: just past their arc

    // Field arrow stays above the surface: pointing away from it, or (B·n < 0) coming down onto it
    const down = bn < -1e-6;
    const fTail = down ? mul(bHat, -FIELD_LEN) : origin;
    const fHead = down ? origin : mul(bHat, FIELD_LEN);
    const fOuter = down ? fTail : fHead; // the end away from the sample centre
    const outerDir = down ? mul(bHat, -1) : bHat;
    const shadowDir = pDir && down ? mul(pDir, -1) : pDir;

    if (pDir && o.angleRef === 'edge') layers.guides.push(guideLine(origin, mul(xEdge, rayLen(R_PHI))));
    if (pDir && o.angleRef === 'wave' && !o.showWave) layers.guides.push(guideLine(origin, mul(kHat, rayLen(R_PHI))));

    // Antenna angle edge -> k, or the equal vertical angle (-edge -> -k) when the field crowds that side
    if (showDev) {
        const busy = pDir ? [pDir, shadowDir, arc(refDir, pDir, 1).mid] : [];
        const devMid = arc(xEdge, kHat, 1).mid;
        const clearance = (m) => Math.min(180, ...busy.map((b) => angleDeg(m, b)));
        const flip = clearance(devMid) < 40 && clearance(mul(devMid, -1)) > clearance(devMid);
        const edgeRay = flip ? mul(xEdge, -1) : xEdge;
        if (!(pDir && o.angleRef === 'edge' && !flip)) layers.guides.push(guideLine(origin, mul(edgeRay, rayLen(R_DEV))));
        if (flip || !o.showWave) layers.guides.push(guideLine(origin, mul(kHat, flip ? -rayLen(R_DEV) : rayLen(R_DEV))));
        drawArc(edgeRay, flip ? mul(kHat, -1) : kHat, R_DEV, st.symbols.dev);
    }

    if (pDir) {
        // in-plane angle to k (or to the 0° edge); a tilted field gets its shadow on the surface
        if (showElev && !isTop) {
            const foot = mul(shadowDir, FIELD_LEN * Math.cos(elevDeg * DEG));
            layers.guides.push(guideLine(origin, foot));
            if (down) layers.guides.push(guideLine(origin, mul(pDir, rayLen(R_PHI)))); // ray for the in-plane angle
            drawArc(shadowDir, outerDir, R_ELEV, st.symbols.theta);
        }
        drawArc(refDir, pDir, R_PHI, st.symbols.phi);
    } else if (!isTop) {
        drawArc(refDir, outerDir, R_ELEV, st.symbols.theta); // field along the normal
    }

    // ---- arrows: field first, k and m on top (k lies on B when the field runs along it)
    const fieldLabel = vectorText(st.symbols.field);
    if (isTop && !pDir) {
        // Top view cannot show the tilt: the field along the normal becomes ⊙ / ⊗
        const [cx, cy] = toSvg(origin);
        const r = fs * 0.75;
        const q = r * 0.6;
        const mark = down
            ? `<path d="M${n2(cx - q)},${n2(cy - q)}L${n2(cx + q)},${n2(cy + q)}M${n2(cx - q)},${n2(cy + q)}L${n2(cx + q)},${n2(cy - q)}" stroke="${C.field}" stroke-width="${n2(L.field * 0.8)}"/>`
            : `<circle cx="${n2(cx)}" cy="${n2(cy)}" r="${n2(L.field * 1.3)}" fill="${C.field}"/>`;
        layers.arrows.push(`<circle cx="${n2(cx)}" cy="${n2(cy)}" r="${n2(r)}" fill="none" stroke="${C.field}" stroke-width="${n2(L.field * 0.8)}"/>`, mark);
        addLine(Array.from({ length: 13 }, (_, i) => [cx + r * Math.cos((i * 30) * DEG), cy + r * Math.sin((i * 30) * DEG)]), W_LINE);
        wantLabel(1, around(cx + r * 0.7, cy - r * 0.7, 0.7, -0.7, fs * 0.3, [0, 90, -90, 180]), fieldLabel, C.field, symbolWidth(st.symbols.field));
    } else {
        layers.arrows.push(arrow(fTail, fHead, C.field, L.field));
        const p = toSvg(fOuter);
        const [dx, dy] = screenDir(origin, outerDir);
        // top view: the tilt goes into the label
        const tiltMode = o.angleLabels === 'value' ? 'both' : o.angleLabels;
        const tilt = isTop && showElev ? `, ${angleText(st.symbols.theta, elevDeg, tiltMode)}` : '';
        const w = symbolWidth(st.symbols.field) + (tilt ? textWidth(', ') + angleWidth(st.symbols.theta, elevDeg, tiltMode) : 0);
        // beyond the outer end, or beside the shaft when that end is crowded
        const o0 = toSvg(origin);
        const beside = [];
        for (const f of [0.75, 0.5, 0.3, 0.1]) {
            const sx = o0[0] + (p[0] - o0[0]) * f;
            const sy = o0[1] + (p[1] - o0[1]) * f;
            beside.push({ sx, sy, dx: -dy, dy: dx, gap: fs * 0.55 }, { sx, sy, dx: dy, dy: -dx, gap: fs * 0.55 });
        }
        wantLabel(1, [...around(p[0], p[1], dx, dy, fs * 0.5, [0, 40, -40, 80, -80, 130, -130]), ...beside], fieldLabel + tilt, C.field, w);
    }
    if (o.showMoment && !state.metrics.isEasyPlaneSingular) {
        const mHat = unit(v.mHat.Sub);
        const mDown = dot(mHat, nrm) < -1e-6;
        const mFrom = mDown ? mul(mHat, -MOMENT_LEN) : origin;
        const mTo = mDown ? origin : mul(mHat, MOMENT_LEN);
        layers.arrows.push(arrow(mFrom, mTo, C.moment, L.wave));
        const [dx, dy] = screenDir(origin, mDown ? mul(mHat, -1) : mHat);
        const p = toSvg(mDown ? mFrom : mTo);
        wantLabel(3, around(p[0], p[1], dx, dy, fs * 0.5, [0, 40, -40, 80, -80, 130, -130]), vectorText(st.symbols.moment), C.moment, symbolWidth(st.symbols.moment));
    }
    if (o.showWave) {
        layers.arrows.push(arrow(origin, mul(kHat, waveLen), C.wave, L.wave));
        const tip = toSvg(mul(kHat, waveLen));
        const [dx, dy] = screenDir(origin, kHat);
        const cands = [];
        for (const gap of [fs * 0.45, fs * 1.1]) cands.push(...around(tip[0], tip[1], dx, dy, gap, [0, 50, -50, 90, -90, 135, -135]));
        wantLabel(2, cands, vectorText(st.symbols.wave), C.wave, symbolWidth(st.symbols.wave));
    }

    // ---- triad: sample edges and normal with Miller-Bravais indices, in the free corner bottom left
    let core = boundsOf(ink);
    if (o.showAxes) {
        const mb = state.millerBravais;
        const axLen = fs * 2.4;
        const o0 = cam.proj(origin);
        const screenOf = (d) => {
            const p = cam.proj(d);
            const dx = p[0] - o0[0];
            const dy = -(p[1] - o0[1]);
            return { ux: dx / (Math.hypot(dx, dy) || 1), uy: dy / (Math.hypot(dx, dy) || 1), d: Math.hypot(dx, dy) };
        };
        const axes = [
            { d: nrm, ...millerLabel(mb.zsPlane.label, mb.zsPlane.latex) },
            { d: xEdge, ...millerLabel(mb.xs.label, mb.xs.latex) },
            { d: [0, 1, 0], ...millerLabel(mb.ys.label, mb.ys.latex) }
        ];
        const cHat = unit(v.ecHat.Sub);
        if (o.showCAxis && axes.every((a) => Math.abs(dot(a.d, cHat)) < 0.999)) {
            const up = dot(cHat, nrm) >= 0;
            axes.push({ d: up ? cHat : mul(cHat, -1), ...millerLabel(up ? '[0001]' : `[0001${OVERLINE}]`, up ? '[0001]' : '[000\\bar{1}]'), dashed: true });
        }
        axes.forEach((a) => {
            const s = screenOf(a.d);
            a.toward = s.d < 0.15; // axis towards the viewer (top view): dot in a circle
            a.ux = s.ux;
            a.uy = s.uy;
            a.ex = s.ux * axLen * Math.min(1, s.d);
            a.ey = s.uy * axLen * Math.min(1, s.d);
            a.cand0 = a.toward ? { sx: -fs * 0.2, sy: fs * 0.2, dx: -0.75, dy: 0.66, gap: fs * 0.45 } : { sx: a.ex, sy: a.ey, dx: s.ux, dy: s.uy, gap: fs * 0.35 };
        });

        // Footprint relative to the triad origin: axes plus each label at its first spot
        let box = { l: -fs * 0.5, r: fs * 0.5, t: -fs * 0.5, b: fs * 0.5 };
        axes.forEach((a) => {
            box = unionBox(box, { l: a.ex, r: a.ex, t: a.ey, b: a.ey });
            box = unionBox(box, labelBox(a.cand0, a.width));
        });

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

        axes.forEach((a) => {
            const turns = [0, 35, -35, 70, -70, 110, -110];
            if (a.toward) {
                layers.axes.push(`<circle cx="${n2(ox)}" cy="${n2(oy)}" r="${n2(fs * 0.38)}" fill="none" stroke="${C.axes}" stroke-width="${n2(L.axes)}"/>`);
                layers.axes.push(`<circle cx="${n2(ox)}" cy="${n2(oy)}" r="${n2(L.axes * 1.1)}" fill="${C.axes}"/>`);
                ink.push([ox - fs * 0.4, oy - fs * 0.4], [ox + fs * 0.4, oy + fs * 0.4]);
                wantLabel(0, around(ox + a.cand0.sx, oy + a.cand0.sy, a.cand0.dx, a.cand0.dy, a.cand0.gap, [0, 40, -40]), a.content, C.text, a.width);
                return;
            }
            const ex = ox + a.ex;
            const ey = oy + a.ey;
            const hl = fs * 0.55;
            const { ux, uy } = a;
            addLine([[ox, oy], [ex, ey]], W_LINE);
            layers.axes.push(
                `<path d="M${n2(ox)},${n2(oy)}L${n2(ex - ux * hl * 0.8)},${n2(ey - uy * hl * 0.8)}" stroke="${C.axes}" stroke-width="${n2(L.axes)}" fill="none"${a.dashed ? ` stroke-dasharray="${n2(L.axes * 3)},${n2(L.axes * 2)}"` : ''}/>`,
                `<path d="M${n2(ex)},${n2(ey)}L${n2(ex - ux * hl - uy * hl * 0.36)},${n2(ey - uy * hl + ux * hl * 0.36)}L${n2(ex - ux * hl + uy * hl * 0.36)},${n2(ey - uy * hl - ux * hl * 0.36)}Z" fill="${C.axes}"/>`
            );
            // Beyond the tip; turned aside if it would cover another one (c axis close to an edge)
            wantLabel(0, [...around(ex, ey, ux, uy, fs * 0.35, turns), ...around(ex, ey, ux, uy, fs * 1.3, turns)], a.content, C.text, a.width);
        });
        core = boundsOf(ink);
    }

    const placed = placeLabels(core);

    // ---- tight crop around everything drawn
    const frame = grow(placed.reduce(unionBox, core), fs * 0.3);
    const W = n2(frame.r - frame.l);
    const H = n2(frame.b - frame.t);
    const X = n2(frame.l);
    const Y = n2(frame.t);

    const s = state.setup;
    const desc = `${s.name || s.cutPreset} · ${s.deviceType} θ_dev=${fmtAngle(s.thetaDevDeg)} · Einbau ${s.mountMode} ψ_S=${fmtAngle(s.psiSDeg)} · φ_M=${fmtAngle(s.phiMDeg)}`;
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
        layer('guides', 'Hilfslinien', layers.guides) +
        layer('arcs', 'Winkel', layers.arcs) +
        layer('arrows', 'Pfeile', layers.arrows) +
        layer('axes', 'Kristallachsen', layers.axes) +
        layer('labels', 'Beschriftung', layers.labels) +
        `</svg>\n`
    );
}
