/**
 * Figure export: one uniform schematic of the sample and the magnetic field.
 *
 * The sample always sits in the same pose (fixed camera in the sample frame F_Sub), so figures made
 * from different setups line up 1:1 in a thesis. Drawn: slab, antenna(s) turned by theta_dev, the
 * wave vector, the field arrow with its in-plane and out-of-plane angle, and a small triad with the
 * Miller-Bravais indices of the sample edges and the surface normal.
 *
 * Geometry lives in sample units (half edge = 1, top surface at z = 0); the SVG viewBox is in pt,
 * so font sizes and line widths in the style are print sizes when the figure is placed at its width.
 */

const DEG = Math.PI / 180;
const PT_PER_MM = 72 / 25.4;

// Camera in F_Sub. Oblique: front-left corner towards the viewer; top: looking down the normal.
const VIEWS = {
    oblique: { azDeg: -58, elDeg: 36 },
    top: { azDeg: -90, elDeg: 90 }
};

const HALF = 1; // half edge of the square sample
const FIELD_LEN = 1.3; // field arrow, reaches past the sample edge
const WAVE_LEN = 0.5; // k arrow, stays in the gap between the two antennas
const ANT_IN = 0.6; // antennas sit at |k| = ANT_IN … ANT_OUT
const ANT_OUT = 0.9;
const MOMENT_LEN = 0.85;
const R_DEV = 0.27; // arc radii: antenna angle, in-plane field angle, elevation
const R_PHI = 0.46;
const R_ELEV = 0.62;
const TINY_DEG = 0.05;

// ---------------------------------------------------------------------------
// Styles: everything that should look the same in every figure of a thesis
// ---------------------------------------------------------------------------

export const BUILTIN_FIGURE_STYLES = {
    Standard: {
        widthMm: 80,
        fontFamily: 'Arial, Helvetica, sans-serif',
        fontSizePt: 9,
        vectorStyle: 'bold',
        sampleThickness: 0.1,
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
        symbols: { field: 'B', wave: 'k', moment: 'm', phi: 'φ', theta: 'θ', dev: 'α' }
    },
    'Schwarz-Weiß': {
        widthMm: 80,
        fontFamily: 'Times New Roman, Times, serif',
        fontSizePt: 10,
        vectorStyle: 'bold',
        sampleThickness: 0.1,
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
        symbols: { field: 'B', wave: 'k', moment: 'm', phi: 'φ', theta: 'θ', dev: 'α' }
    },
    Folie: {
        widthMm: 120,
        fontFamily: 'Arial, Helvetica, sans-serif',
        fontSizePt: 14,
        vectorStyle: 'bold',
        sampleThickness: 0.12,
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
        symbols: { field: 'B', wave: 'k', moment: 'm', phi: 'φ', theta: 'θ', dev: 'α' }
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

/** Miller-Bravais label "[1̅1̅20]" (combining overlines) -> tspans with real overlines, or LaTeX. */
function millerSpans(label, latex, latexLabel) {
    if (latex) return esc(`$${latexLabel}$`);
    if (/\s/.test(label)) return esc(label); // non-integer indices: plain
    let out = '';
    const chars = [...label];
    chars.forEach((c, i) => {
        if (c === '̅') return;
        if (chars[i + 1] === '̅') out += `<tspan text-decoration="overline">${esc(c)}</tspan><tspan font-size="30%"> </tspan>`;
        else out += esc(c);
    });
    return out;
}

// ---------------------------------------------------------------------------
// Projection
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

/** Fixed frame per view: sample + the whole upper half sphere the field arrow can reach. */
function frameBounds(cam, thickness) {
    const pts = [];
    for (const z of [0, -thickness]) for (const x of [-HALF, HALF]) for (const y of [-HALF, HALF]) pts.push([x, y, z]);
    const r = FIELD_LEN * 1.08;
    for (let a = 0; a < 360; a += 5) {
        for (let e = 0; e <= 90; e += 5) {
            pts.push([r * Math.cos(e * DEG) * Math.cos(a * DEG), r * Math.cos(e * DEG) * Math.sin(a * DEG), r * Math.sin(e * DEG)]);
        }
    }
    const xy = pts.map(cam.proj);
    return {
        minX: Math.min(...xy.map((p) => p[0])),
        maxX: Math.max(...xy.map((p) => p[0])),
        minY: Math.min(...xy.map((p) => p[1])),
        maxY: Math.max(...xy.map((p) => p[1]))
    };
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

    const v = state.vectors;
    const kHat = v.kHat.Sub;
    const tHat = v.tHat.Sub;
    const bHat = unit(v.bHat.Sub);
    const xEdge = [1, 0, 0];
    const nrm = [0, 0, 1];

    // ---- frame: fixed scale per view and width, so every figure lines up
    const fb = frameBounds(cam, thick);
    const pad = fs * 1.6;
    const widthPt = (Number(st.widthMm) || 80) * PT_PER_MM;
    const scale = (widthPt - 2 * pad) / (fb.maxX - fb.minX);
    let heightPt = (fb.maxY - fb.minY) * scale + 2 * pad;
    const toSvg = (p) => {
        const [x, y] = cam.proj(p);
        return [pad + (x - fb.minX) * scale, pad + (fb.maxY - y) * scale];
    };
    const pt = (p) => toSvg(p).map(n2).join(',');
    const pathOf = (pts, close = false) => 'M' + pts.map(pt).join('L') + (close ? 'Z' : '');

    const layers = { sample: [], antenna: [], guides: [], arcs: [], arrows: [], labels: [], axes: [] };
    const halo = C.background && C.background !== 'none' ? C.background : '#ffffff';

    const textAt = (x, y, content, color, anchor = 'middle', extra = '') =>
        `<text x="${n2(x)}" y="${n2(y)}" text-anchor="${anchor}" fill="${color}" font-size="${fs}"${extra}>${content}</text>`;

    // Label just beyond a screen point, pushed away along a screen direction (unit, y down)
    const labelNear = (sx, sy, dirX, dirY, content, color, gap = fs * 0.55) => {
        const x = sx + dirX * gap;
        const y = sy + dirY * gap;
        const anchor = dirX > 0.35 ? 'start' : dirX < -0.35 ? 'end' : 'middle';
        const yAdj = y + fs * 0.35 + dirY * fs * 0.35;
        return textAt(x, yAdj, content, color, anchor, ` stroke="${halo}" stroke-width="${n2(fs * 0.28)}" stroke-linejoin="round" paint-order="stroke"`);
    };

    // Box of a labelNear() text (same anchor rules), for overlap checks; widths are estimates
    const textWidth = (plain) => [...String(plain).replace(/̅/g, '')].length * fs * 0.56;
    const symbolWidth = (str) => {
        const { main, sub: s } = parseSymbol(str);
        return main.length * fs * 0.62 + s.length * fs * 0.45;
    };
    const labelBox = (sx, sy, dirX, dirY, w, gap) => {
        const x = sx + dirX * gap;
        const y = sy + dirY * gap + fs * 0.35 + dirY * fs * 0.35;
        const l = dirX > 0.35 ? x : dirX < -0.35 ? x - w : x - w / 2;
        return { l, r: l + w, t: y - fs * 0.8, b: y + fs * 0.25 };
    };

    // Labels are collected and placed at the end: each takes the first of its candidate spots that
    // stays inside the figure and clear of the labels placed before it (lower prio = placed first)
    const labelReqs = [];
    const wantLabel = (prio, cands, content, color, w) => labelReqs.push({ prio, cands, content, color, w });
    const around = (sx, sy, dx, dy, gap, turns) =>
        turns.map((t) => {
            const c = Math.cos(t * DEG);
            const s = Math.sin(t * DEG);
            return { sx, sy, dx: dx * c - dy * s, dy: dx * s + dy * c, gap };
        });
    const placeLabels = () => {
        const placed = [];
        const inside = (b) => b.l >= 1 && b.r <= widthPt - 1 && b.t >= 1 && b.b <= heightPt - 1;
        const m = fs * 0.2;
        const clash = (b) => placed.some((p) => b.l < p.r + m && p.l < b.r + m && b.t < p.b + m && p.t < b.b + m);
        labelReqs
            .sort((a, b) => a.prio - b.prio)
            .forEach((q) => {
                const boxes = q.cands.map((c) => labelBox(c.sx, c.sy, c.dx, c.dy, q.w, c.gap));
                let i = boxes.findIndex((b) => inside(b) && !clash(b));
                if (i < 0) i = Math.max(0, boxes.findIndex(inside));
                placed.push(boxes[i]);
                const c = q.cands[i];
                layers.labels.push(labelNear(c.sx, c.sy, c.dx, c.dy, q.content, q.color, c.gap));
            });
    };

    const screenDir = (from, to) => {
        const a = toSvg(from);
        const b = toSvg(to);
        const d = Math.hypot(b[0] - a[0], b[1] - a[1]) || 1;
        return [(b[0] - a[0]) / d, (b[1] - a[1]) / d];
    };

    // Arrow in 3D, head drawn in screen space so it keeps its size when foreshortened
    const arrow = (from, to, color, width, { dash = '', haloW = 0 } = {}) => {
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
        const shaft = `M${n2(a[0])},${n2(a[1])}L${n2(bx)},${n2(by)}`;
        const head = `M${n2(b[0])},${n2(b[1])}L${n2(b[0] - ux * hl - uy * hw)},${n2(b[1] - uy * hl + ux * hw)}L${n2(b[0] - ux * hl + uy * hw)},${n2(b[1] - uy * hl - ux * hw)}Z`;
        let s = '';
        if (haloW) {
            s += `<path d="${shaft}" stroke="${halo}" stroke-width="${n2(width + 2 * haloW)}" stroke-linecap="round" fill="none"/>`;
            s += `<path d="${head}" fill="${halo}" stroke="${halo}" stroke-width="${n2(2 * haloW)}" stroke-linejoin="round"/>`;
        }
        s += `<path d="${shaft}" stroke="${color}" stroke-width="${n2(width)}" fill="none"${dash ? ` stroke-dasharray="${dash}"` : ''}/>`;
        s += `<path d="${head}" fill="${color}"/>`;
        return s;
    };

    const guideLine = (from, to, color = C.angle) =>
        `<path d="${pathOf([from, to])}" stroke="${color}" stroke-width="${n2(L.guide)}" stroke-dasharray="${n2(L.guide * 1.5)},${n2(L.guide * 2.6)}" stroke-linecap="round" fill="none"/>`;

    // Arc from direction a to direction b (3D unit vectors) around the origin; returns label point too
    const arc = (a, b, r) => {
        const ang = angleDeg(a, b);
        const ua = unit(a);
        let perp = sub(unit(b), mul(ua, dot(ua, unit(b))));
        if (len(perp) < 1e-9) perp = unit([-ua[1], ua[0], 0]);
        perp = unit(perp);
        const pts = [];
        const steps = Math.max(8, Math.ceil(ang / 3));
        for (let i = 0; i <= steps; i++) {
            const t = (ang * i * DEG) / steps;
            pts.push(add(mul(ua, r * Math.cos(t)), mul(perp, r * Math.sin(t))));
        }
        const dirAt = (f) => add(mul(ua, Math.cos(ang * f * DEG)), mul(perp, Math.sin(ang * f * DEG)));
        return { ang, d: pathOf(pts), mid: dirAt(0.5), dirAt, r };
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

    // Arc with its value just outside the middle; if that spot is taken, further along the arc or further out
    let arcCount = 0;
    const drawArc = (a, b, r, symbol, color = C.angle) => {
        const A = arc(a, b, r);
        if (A.ang < TINY_DEG) return;
        layers.arcs.push(`<path d="${A.d}" stroke="${color}" stroke-width="${n2(L.arc)}" fill="none"/>`);
        const cands = [];
        for (const gap of [fs * 0.45, fs * 1.4]) {
            for (const f of [0.5, 0.3, 0.7, 0.15, 0.85]) {
                const d = A.dirAt(f);
                const [sx, sy] = toSvg(mul(d, A.r));
                const [dx, dy] = screenDir([0, 0, 0], d);
                cands.push({ sx, sy, dx, dy, gap });
            }
        }
        wantLabel(10 + arcCount++, cands, angleText(symbol, A.ang), C.text, angleWidth(symbol, A.ang));
    };

    // ---- sample slab: visible side faces (auto-shaded), then the top face
    const corners = (z) => [[-HALF, -HALF, z], [HALF, -HALF, z], [HALF, HALF, z], [-HALF, HALF, z]];
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
            layers.sample.push(`<path d="${pathOf(f.q, true)}" fill="${shade(C.sample, i % 2 ? -0.28 : -0.16)}" stroke="${sampleStroke}" stroke-width="${n2(L.outline)}" stroke-linejoin="round"/>`);
        });
    }
    layers.sample.push(`<path d="${pathOf(top, true)}" fill="${C.sample}" stroke="${sampleStroke}" stroke-width="${n2(L.outline)}" stroke-linejoin="round"/>`);

    // ---- antennas in the wave frame (k, t): emitter and receiver either side of the centre
    const wave = (kk, tt) => add(mul(kHat, kk), mul(tHat, tt));
    const rect = (k0, k1, t0, t1) => [wave(k0, t0), wave(k1, t0), wave(k1, t1), wave(k0, t1)];
    if (o.showAntenna) {
        const parts = [];
        const isIdt = state.setup.deviceType === 'IDT';
        for (const side of [-1, 1]) {
            const kIn = ANT_IN * side;
            const kOut = ANT_OUT * side;
            const [k0, k1] = [Math.min(kIn, kOut), Math.max(kIn, kOut)];
            if (isIdt) {
                const bar = 0.06;
                const tSpan = 0.4;
                parts.push(rect(k0, k1, tSpan - bar, tSpan), rect(k0, k1, -tSpan, -tSpan + bar));
                const fingers = 7;
                const pitch = (k1 - k0) / fingers;
                const fw = pitch * 0.42;
                for (let i = 0; i < fingers; i++) {
                    const kc = k0 + pitch * (i + 0.5);
                    const fromTop = i % 2 === 0;
                    parts.push(rect(kc - fw / 2, kc + fw / 2, fromTop ? -tSpan + bar + 0.07 : -tSpan + bar, fromTop ? tSpan - bar : tSpan - bar - 0.07));
                }
            } else {
                // coplanar waveguide: ground - signal - ground, all along t
                const kc = (k0 + k1) / 2;
                const tSpan = 0.45;
                const sig = 0.03;
                const gap = 0.035;
                const gnd = 0.08;
                parts.push(rect(kc - sig, kc + sig, -tSpan, tSpan));
                parts.push(rect(kc + sig + gap, kc + sig + gap + gnd, -tSpan, tSpan));
                parts.push(rect(kc - sig - gap - gnd, kc - sig - gap, -tSpan, tSpan));
            }
        }
        const stroke = shade(C.antenna, -0.35);
        layers.antenna.push(
            `<path d="${parts.map((q) => pathOf(q, true)).join('')}" fill="${C.antenna}" stroke="${stroke}" stroke-width="${n2(L.outline * 0.6)}" stroke-linejoin="round"/>`
        );
    }

    // ---- angles on the surface: antenna turn (0° edge -> k), field in plane, field out of plane
    const devDeg = angleDeg(xEdge, kHat);
    const showDev = devDeg > TINY_DEG && (o.showAntenna || o.showWave);
    const bn = dot(bHat, nrm);
    const bPlane = [bHat[0], bHat[1], 0];
    const pDir = len(bPlane) > 1e-6 ? unit(bPlane) : null;
    const elevDeg = Math.asin(Math.min(1, Math.abs(bn))) / DEG;
    const showElev = elevDeg > TINY_DEG;
    const refDir = o.angleRef === 'edge' ? xEdge : kHat;

    // Field arrow stays above the surface: pointing away from it, or (B·n < 0) coming down onto it
    const down = bn < -1e-6;
    const fTail = down ? mul(bHat, -FIELD_LEN) : [0, 0, 0];
    const fHead = down ? [0, 0, 0] : mul(bHat, FIELD_LEN);
    const fOuter = down ? fTail : fHead; // the end away from the sample centre
    const outerDir = down ? mul(bHat, -1) : bHat;
    const shadowDir = pDir && down ? mul(pDir, -1) : pDir;

    if (o.angleRef === 'edge') layers.guides.push(guideLine([0, 0, 0], mul(xEdge, 0.92)));
    if (o.angleRef === 'wave' && !o.showWave) layers.guides.push(guideLine([0, 0, 0], mul(kHat, 0.6)));

    // Antenna angle edge -> k, or the equal vertical angle (-edge -> -k) when the field crowds that side
    if (showDev) {
        const busy = pDir ? [pDir, shadowDir, arc(refDir, pDir, 1).mid] : [];
        const devMid = arc(xEdge, kHat, 1).mid;
        const clearance = (m) => Math.min(180, ...busy.map((b) => angleDeg(m, b)));
        const flip = clearance(devMid) < 40 && clearance(mul(devMid, -1)) > clearance(devMid);
        const edgeRay = flip ? mul(xEdge, -1) : xEdge;
        if (!(o.angleRef === 'edge' && !flip)) layers.guides.push(guideLine([0, 0, 0], mul(edgeRay, 0.92)));
        if (flip) layers.guides.push(guideLine([0, 0, 0], mul(kHat, -(R_DEV + 0.12))));
        drawArc(edgeRay, flip ? mul(kHat, -1) : kHat, R_DEV, st.symbols.dev);
    }

    if (pDir) {
        // in-plane angle to k (or to the 0° edge)
        drawArc(refDir, pDir, R_PHI, st.symbols.phi);
        if (showElev && !isTop) {
            const foot = mul(shadowDir, FIELD_LEN * Math.cos(elevDeg * DEG));
            layers.guides.push(guideLine([0, 0, 0], foot));
            layers.guides.push(guideLine(fOuter, foot));
            if (down) layers.guides.push(guideLine([0, 0, 0], mul(pDir, R_PHI + 0.12))); // ray for the in-plane angle
            drawArc(shadowDir, outerDir, R_ELEV, st.symbols.theta);
        }
    } else if (!isTop) {
        drawArc(refDir, outerDir, R_ELEV, st.symbols.theta); // field along the normal
    }

    // ---- arrows: k, optional m, field on top with a halo so it reads over the antenna
    if (o.showWave) {
        layers.arrows.push(arrow([0, 0, 0], mul(kHat, WAVE_LEN), C.wave, L.wave));
        const tip = toSvg(mul(kHat, WAVE_LEN));
        const [dx, dy] = screenDir([0, 0, 0], kHat);
        // beside the tip, on the side away from the field arrow
        const side = pDir && dot(tHat, pDir) > 0 ? -1 : 1;
        const [sx, sy] = screenDir([0, 0, 0], mul(tHat, side));
        const n = Math.hypot(dx + sx, dy + sy) || 1;
        wantLabel(2, around(tip[0], tip[1], (dx + sx) / n, (dy + sy) / n, fs * 0.5, [0, 40, -40, 80, -80, 130, -130]), vectorText(st.symbols.wave), C.wave, symbolWidth(st.symbols.wave));
    }
    if (o.showMoment && !state.metrics.isEasyPlaneSingular) {
        const mHat = unit(v.mHat.Sub);
        const mDown = dot(mHat, nrm) < -1e-6;
        const mFrom = mDown ? mul(mHat, -MOMENT_LEN) : [0, 0, 0];
        const mTo = mDown ? [0, 0, 0] : mul(mHat, MOMENT_LEN);
        layers.arrows.push(arrow(mFrom, mTo, C.moment, L.wave, { haloW: 0.8 }));
        const outer = mDown ? mFrom : mTo;
        const [dx, dy] = screenDir([0, 0, 0], mDown ? mul(mHat, -1) : mHat);
        const p = toSvg(outer);
        wantLabel(3, around(p[0], p[1], dx, dy, fs * 0.55, [0, 40, -40, 80, -80, 130, -130]), vectorText(st.symbols.moment), C.moment, symbolWidth(st.symbols.moment));
    }
    // Top view cannot show the tilt: the field along the normal becomes ⊙ / ⊗, a tilt goes into the label
    const fieldLabel = vectorText(st.symbols.field);
    if (isTop && !pDir) {
        const [cx, cy] = toSvg([0, 0, 0]);
        const r = fs * 0.75;
        const mark = down
            ? `<path d="M${n2(cx - r * 0.6)},${n2(cy - r * 0.6)}L${n2(cx + r * 0.6)},${n2(cy + r * 0.6)}M${n2(cx - r * 0.6)},${n2(cy + r * 0.6)}L${n2(cx + r * 0.6)},${n2(cy - r * 0.6)}" stroke="${C.field}" stroke-width="${n2(L.field * 0.8)}"/>`
            : `<circle cx="${n2(cx)}" cy="${n2(cy)}" r="${n2(L.field * 1.3)}" fill="${C.field}"/>`;
        layers.arrows.push(`<circle cx="${n2(cx)}" cy="${n2(cy)}" r="${n2(r)}" fill="${halo}" stroke="${C.field}" stroke-width="${n2(L.field * 0.8)}"/>`, mark);
        wantLabel(1, around(cx + r * 0.7, cy - r * 0.7, 0.7, -0.7, fs * 0.3, [0, 90, -90, 180]), fieldLabel, C.field, symbolWidth(st.symbols.field));
    } else {
        layers.arrows.push(arrow(fTail, fHead, C.field, L.field, { haloW: Math.max(0.8, L.field * 0.6) }));
        const p = toSvg(fOuter);
        const [dx, dy] = screenDir([0, 0, 0], outerDir);
        const tiltMode = o.angleLabels === 'value' ? 'both' : o.angleLabels;
        const tilt = isTop && showElev ? `, ${angleText(st.symbols.theta, elevDeg, tiltMode)}` : '';
        const w = symbolWidth(st.symbols.field) + (tilt ? textWidth(', ') + angleWidth(st.symbols.theta, elevDeg, tiltMode) : 0);
        wantLabel(1, around(p[0], p[1], dx, dy, fs * 0.6, [0, 40, -40, 80, -80, 130, -130]), fieldLabel + tilt, C.field, w);
    }

    // ---- triad: sample edges and normal with Miller-Bravais indices, bottom left.
    // Its place depends only on the view (fixed axis directions, typical label width), never on the setup.
    if (o.showAxes) {
        const mb = state.millerBravais;
        const axLen = fs * 2.6;
        const o0 = cam.proj([0, 0, 0]);
        const screenOf = (d) => {
            const p = cam.proj(d);
            const dx = p[0] - o0[0];
            const dy = -(p[1] - o0[1]);
            return { dx, dy, d: Math.hypot(dx, dy) };
        };
        const axes = [
            { d: nrm, label: millerSpans(mb.zsPlane.label, latex, mb.zsPlane.latex), width: textWidth(mb.zsPlane.label) },
            { d: xEdge, label: millerSpans(mb.xs.label, latex, mb.xs.latex), width: textWidth(mb.xs.label) },
            { d: [0, 1, 0], label: millerSpans(mb.ys.label, latex, mb.ys.latex), width: textWidth(mb.ys.label) }
        ];

        // Box around the triad relative to its origin: axis tips plus a typical label each
        const labelW = fs * 0.55 * 7;
        const box = { l: -fs, r: fs, t: -fs, b: fs };
        axes.forEach((a) => {
            const s = screenOf(a.d);
            if (s.d < 0.15) {
                box.l = Math.min(box.l, -fs * 0.9 - labelW);
                box.b = Math.max(box.b, fs * 1.9);
                return;
            }
            const ex = (s.dx / s.d) * axLen * Math.min(1, s.d);
            const ey = (s.dy / s.d) * axLen * Math.min(1, s.d);
            const lx = ex + (s.dx / s.d) * fs * 0.4;
            const ly = ey + (s.dy / s.d) * fs * 0.4;
            const ux = s.dx / s.d;
            const x0 = ux > 0.35 ? lx : ux < -0.35 ? lx - labelW : lx - labelW / 2;
            box.l = Math.min(box.l, x0);
            box.r = Math.max(box.r, x0 + labelW);
            box.t = Math.min(box.t, ly - fs * 1.1);
            box.b = Math.max(box.b, ly + fs * 0.6);
        });

        // Bottom left; slide down (and grow the figure) until it clears the sample
        const hull = convexHull([...corners(0), ...corners(-thick)].map(toSvg));
        const ox = pad * 0.5 - box.l;
        let oy = heightPt - pad * 0.5 - box.b;
        const rectAt = (y) => [[ox + box.l, y + box.t], [ox + box.r, y + box.t], [ox + box.r, y + box.b], [ox + box.l, y + box.b]];
        while (convexOverlap(hull, rectAt(oy))) oy += fs * 0.25;
        heightPt = Math.max(heightPt, oy + box.b + pad * 0.5);

        const cHat = unit(v.ecHat.Sub);
        if (o.showCAxis && axes.every((a) => Math.abs(dot(a.d, cHat)) < 0.999)) {
            const up = dot(cHat, nrm) >= 0;
            axes.push({ d: up ? cHat : mul(cHat, -1), label: millerSpans(up ? '[0001]' : '[0001\u0305]', latex, up ? '[0001]' : '[000\\bar{1}]'), width: textWidth('[0001]'), dashed: true });
        }
        axes.forEach((a) => {
            const s = screenOf(a.d);
            if (s.d < 0.15) {
                // axis towards the viewer (top view): dot in a circle
                layers.axes.push(`<circle cx="${n2(ox)}" cy="${n2(oy)}" r="${n2(fs * 0.38)}" fill="none" stroke="${C.axes}" stroke-width="${n2(L.axes)}"/>`);
                layers.axes.push(`<circle cx="${n2(ox)}" cy="${n2(oy)}" r="${n2(L.axes * 1.1)}" fill="${C.axes}"/>`);
                wantLabel(0, around(ox - fs * 0.2, oy + fs * 0.2, -0.75, 0.66, fs * 0.45, [0, 40, -40]), a.label, C.text, a.width);
                return;
            }
            const ux = s.dx / s.d;
            const uy = s.dy / s.d;
            const ex = ox + ux * axLen * Math.min(1, s.d);
            const ey = oy + uy * axLen * Math.min(1, s.d);
            const hl = fs * 0.55;
            layers.axes.push(
                `<path d="M${n2(ox)},${n2(oy)}L${n2(ex - ux * hl * 0.8)},${n2(ey - uy * hl * 0.8)}" stroke="${C.axes}" stroke-width="${n2(L.axes)}" fill="none"${a.dashed ? ` stroke-dasharray="${n2(L.axes * 3)},${n2(L.axes * 2)}"` : ''}/>`,
                `<path d="M${n2(ex)},${n2(ey)}L${n2(ex - ux * hl - uy * hl * 0.36)},${n2(ey - uy * hl + ux * hl * 0.36)}L${n2(ex - ux * hl + uy * hl * 0.36)},${n2(ey - uy * hl - ux * hl * 0.36)}Z" fill="${C.axes}"/>`
            );
            // Beyond the tip; turned aside if it would cover another one (c axis close to an edge)
            const turns = [0, 35, -35, 70, -70, 110, -110];
            wantLabel(0, [...around(ex, ey, ux, uy, fs * 0.35, turns), ...around(ex, ey, ux, uy, fs * 1.3, turns)], a.label, C.text, a.width);
        });
    }

    placeLabels();

    // ---- assemble
    const s = state.setup;
    const desc = `${s.name || s.cutPreset} · ${s.deviceType} θ_dev=${fmtAngle(s.thetaDevDeg)} · Einbau ${s.mountMode} ψ_S=${fmtAngle(s.psiSDeg)} · φ_M=${fmtAngle(s.phiMDeg)}`;
    const layer = (id, label, items) =>
        items.length ? `  <g id="${id}" inkscape:groupmode="layer" inkscape:label="${label}">\n    ${items.join('\n    ')}\n  </g>\n` : '';
    const W = n2(widthPt);
    const H = n2(heightPt);
    const bg = C.background && C.background !== 'none' ? `  <rect id="background" width="${W}" height="${H}" fill="${C.background}"/>\n` : '';

    return (
        `<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n` +
        `<svg xmlns="http://www.w3.org/2000/svg" xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape" ` +
        `width="${n2(widthPt / PT_PER_MM)}mm" height="${n2(heightPt / PT_PER_MM)}mm" viewBox="0 0 ${W} ${H}" ` +
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
