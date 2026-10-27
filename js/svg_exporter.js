/**
 * svg_exporter.js
 * Publication-ready layered SVG generator for Bachelor/Master/PhD theses
 * (Inkscape / Illustrator / LaTeX pdf_tex compatible).
 *
 * Generates clean vector schematics with grouped Inkscape layers:
 *   - 01_Background_Grid
 *   - 02_Crystal_EasyPlane
 *   - 03_Substrate_Chip
 *   - 04_Microstructure_IDT_CPW
 *   - 04b_Field_Sweep (only with a sweep range)
 *   - 05_Vectors_Projections
 *   - 06_Angle_Arcs
 *   - 07_Labels_Legend
 */

import {
    DEG2RAD,
    add3,
    scale3,
    dot3,
    norm3,
    normalize3,
    matVec3,
    evaluateSetup,
    FRAME_METADATA
} from './math_core.js';

function escXml(str) {
    return String(str)
        .replace(/&/g, '&amp;')
        .replace(/</g, '&lt;')
        .replace(/>/g, '&gt;')
        .replace(/"/g, '&quot;');
}

function fmt(num) {
    return Number(num).toFixed(2);
}

export function generateThesisSVG(stateA, options = {}) {
    const {
        layout = 'combined', // 'combined' | '3d' | '2d_sample' | '2d_easy'
        lockFrame = 'Sub',   // 'Lab' | 'Sub' | 'Wave' | 'Cryst'
        cameraAzimuthDeg = -38,
        cameraElevationDeg = 28,
        latexLabels = false,
        transparentBg = false,
        showEasyPlane = true,
        showMicrostructure = true,
        showCrystalAxes = true,
        showEquilibrium = true,
        showAngleArcs = true,
        showLegend = true,
        sweep = null // { fromDeg, toDeg }: draw where B and m̂ travel while φ_M runs from start to stop
    } = options;

    const width = layout === 'combined' ? 1080 : 680;
    const height = layout === 'combined' ? 540 : 540;

    // Publication color palette (optimized for crisp print & vector clarity)
    const palette = {
        bg: '#ffffff',
        grid: '#cbd5e1',
        textMain: '#0f172a',
        textMuted: '#475569',
        chipFill: '#f1f5f9',
        chipSide: '#cbd5e1',
        chipStroke: '#334155',
        easyFill: '#f3e8ff',
        easyStroke: '#7e22ce',
        goldFill: '#fef3c7',
        goldStroke: '#d97706',
        B: '#0284c7',
        B_easy: '#d97706',
        B_sample: '#0369a1',
        m: '#e11d48',
        l: '#059669',
        k: '#a21caf',
        t: '#475569',
        n: '#2563eb',
        ec: '#6d28d9',
        ea: '#7c3aed',
        em: '#4f46e5'
    };

    const lbl = (unicodeStr, texStr) => (latexLabels ? texStr : unicodeStr);

    // Layer buckets
    const layers = {
        bg: [],
        crystal: [],
        substrate: [],
        device: [],
        sweep: [],
        vectors: [],
        arcs: [],
        labels: []
    };

    if (!transparentBg) {
        layers.bg.push(`<rect x="0" y="0" width="${width}" height="${height}" fill="${palette.bg}" />`);
    }

    // Marker definitions for vector arrowheads
    const markerColors = [
        { id: 'arrow-B', color: palette.B },
        { id: 'arrow-Beasy', color: palette.B_easy },
        { id: 'arrow-Bsamp', color: palette.B_sample },
        { id: 'arrow-m', color: palette.m },
        { id: 'arrow-l', color: palette.l },
        { id: 'arrow-k', color: palette.k },
        { id: 'arrow-t', color: palette.t },
        { id: 'arrow-n', color: palette.n },
        { id: 'arrow-ec', color: palette.ec },
        { id: 'arrow-ea', color: palette.ea },
        { id: 'arrow-em', color: palette.em }
    ];

    const defsContent = markerColors
        .map(
            (m) => `
    <marker id="${m.id}" viewBox="0 0 10 10" refX="8.2" refY="5" markerWidth="6.5" markerHeight="6.5" orient="auto-start-reverse">
      <path d="M 0 1.2 L 9.5 5 L 0 8.8 z" fill="${m.color}" />
    </marker>`
        )
        .join('');

    // Helper to draw SVG arrow + label
    const addArrow2D = (x1, y1, x2, y2, color, markerId, labelText, strokeWidth = 2.4, dashed = false) => {
        const len = Math.hypot(x2 - x1, y2 - y1);
        if (len < 3) return;
        const dashAttr = dashed ? ` stroke-dasharray="5,3"` : '';
        layers.vectors.push(
            `<line x1="${fmt(x1)}" y1="${fmt(y1)}" x2="${fmt(x2)}" y2="${fmt(y2)}" stroke="${color}" stroke-width="${strokeWidth}" stroke-linecap="round"${dashAttr} marker-end="url(#${markerId})" />`
        );
        if (labelText) {
            const ang = Math.atan2(y2 - y1, x2 - x1);
            const lx = x2 + Math.cos(ang) * 14;
            const ly = y2 + Math.sin(ang) * 14 + 4;
            layers.labels.push(
                `<text x="${fmt(lx)}" y="${fmt(ly)}" fill="${color}" font-family="Inter, Arial, sans-serif" font-size="12" font-weight="700" text-anchor="middle">${escXml(labelText)}</text>`
            );
        }
    };

    // Field sweep: the setup re-evaluated along φ_M = from … to (2° steps, at most one full turn)
    const sweepStates = [];
    let sweepFrom = 0;
    let sweepTo = 0;
    if (sweep && Number.isFinite(sweep.fromDeg) && Number.isFinite(sweep.toDeg) && Math.abs(sweep.toDeg - sweep.fromDeg) >= 1) {
        sweepFrom = sweep.fromDeg;
        const span = Math.max(-360, Math.min(360, sweep.toDeg - sweep.fromDeg));
        sweepTo = sweepFrom + span;
        const n = Math.max(8, Math.ceil(Math.abs(span) / 2));
        for (let i = 0; i <= n; i++) {
            sweepStates.push(evaluateSetup({ ...stateA.setup, phiMDeg: sweepFrom + (span * i) / n }));
        }
    }
    const fmtDeg = (d) => (Number.isInteger(d) ? String(d) : d.toFixed(1));
    const degLbl = (d) => lbl(`${fmtDeg(d)}°`, `$${fmtDeg(d)}^\\circ$`);

    // Dotted trace of an arrow tip along the sweep; arrowhead = sense of rotation.
    // pointOf(state) returns the projected tip or null (e.g. m̂ undefined) – the trace breaks there.
    const addSweepTrace = (pointOf, color, markerId, center, withLabels = false) => {
        if (sweepStates.length < 2) return;
        const runs = [];
        let cur = [];
        for (const s of sweepStates) {
            const p = pointOf(s);
            if (p) cur.push(p);
            else if (cur.length) {
                runs.push(cur);
                cur = [];
            }
        }
        if (cur.length) runs.push(cur);
        const pathLen = (run) => run.reduce((acc, p, i) => (i ? acc + Math.hypot(p.x - run[i - 1].x, p.y - run[i - 1].y) : 0), 0);
        const drawn = runs.filter((run) => run.length > 1 && pathLen(run) > 4);
        drawn.forEach((run, idx) => {
            const head = idx === drawn.length - 1 ? ` marker-end="url(#${markerId})"` : '';
            layers.sweep.push(
                `<polyline points="${run.map((p) => `${fmt(p.x)},${fmt(p.y)}`).join(' ')}" fill="none" stroke="${color}" stroke-width="1.5" stroke-dasharray="0.1,3.6" stroke-linecap="round" stroke-linejoin="round"${head} />`
            );
        });
        if (!withLabels || !drawn.length) return;
        const firstRun = drawn[0];
        const lastRun = drawn[drawn.length - 1];
        const first = firstRun[0];
        const last = lastRun[lastRun.length - 1];
        // Label sideways to the trace (normal pointing away from the panel centre), so it never
        // sits on the radial label of the arrow that may end at the same spot
        const labelAt = (p, q, text) => {
            let nx = -(q.y - p.y);
            let ny = q.x - p.x;
            const nl = Math.hypot(nx, ny) || 1;
            nx /= nl;
            ny /= nl;
            if (nx * (p.x - center.x) + ny * (p.y - center.y) < 0) {
                nx = -nx;
                ny = -ny;
            }
            const anchor = Math.abs(nx) < 0.35 ? 'middle' : nx > 0 ? 'start' : 'end';
            layers.labels.push(
                `<text x="${fmt(p.x + nx * 9)}" y="${fmt(p.y + ny * 9 + 3.5 + (ny > 0.35 ? 4 : 0))}" fill="${palette.textMuted}" font-family="Inter, Arial, sans-serif" font-size="10.5" text-anchor="${anchor}">${escXml(text)}</text>`
            );
        };
        layers.sweep.push(`<circle cx="${fmt(first.x)}" cy="${fmt(first.y)}" r="2.6" fill="${color}" />`);
        if (Math.hypot(last.x - first.x, last.y - first.y) < 16) {
            labelAt(first, firstRun[1], lbl(`${fmtDeg(sweepFrom)}° → ${fmtDeg(sweepTo)}°`, `$${fmtDeg(sweepFrom)}^\\circ \\to ${fmtDeg(sweepTo)}^\\circ$`));
        } else {
            labelAt(first, firstRun[1], degLbl(sweepFrom));
            labelAt(last, lastRun[lastRun.length - 2], degLbl(sweepTo));
        }
    };
    const mDefined = (s) => showEquilibrium && !s.metrics.isEasyPlaneSingular;

    // ------------------------------------------------------------------------
    // Panel 1: 3D Isometric Schematic
    // ------------------------------------------------------------------------
    const render3DPanel = (cx, cy, scale) => {
        const F = lockFrame;
        const R_Lock_from = (src) => stateA.matrices.getRotationMatrix(F, src);
        const toLock = (v, src) => matVec3(R_Lock_from(src), v);

        const az = cameraAzimuthDeg * DEG2RAD;
        const el = cameraElevationDeg * DEG2RAD;
        const caz = Math.cos(az);
        const saz = Math.sin(az);
        const cel = Math.cos(el);
        const sel = Math.sin(el);

        const proj3 = (p) => {
            const xCam = -saz * p[0] + caz * p[1];
            const dHoriz = -caz * p[0] - saz * p[1];
            const z1 = p[2];
            const yCam = sel * dHoriz + cel * z1;
            return {
                x: cx + xCam * scale,
                y: cy - yCam * scale
            };
        };

        // Reference ellipse in locked XY plane
        const circlePts = [];
        for (let d = 0; d <= 360; d += 6) {
            const r = d * DEG2RAD;
            const pt = proj3([1.05 * Math.cos(r), 1.05 * Math.sin(r), 0]);
            circlePts.push(`${fmt(pt.x)},${fmt(pt.y)}`);
        }
        layers.bg.push(
            `<polyline points="${circlePts.join(' ')}" fill="none" stroke="${palette.grid}" stroke-width="1.2" stroke-dasharray="4,4" />`
        );

        // Basal Easy-Plane (0001) hexagon in F_Cryst
        if (showEasyPlane) {
            const hexPts = [];
            for (let i = 0; i < 6; i++) {
                const ang = i * 60 * DEG2RAD;
                const pt = proj3(toLock([0.92 * Math.cos(ang), 0.92 * Math.sin(ang), 0], 'Cryst'));
                hexPts.push(`${fmt(pt.x)},${fmt(pt.y)}`);
            }
            layers.crystal.push(
                `<polygon points="${hexPts.join(' ')}" fill="${palette.easyFill}" fill-opacity="0.55" stroke="${palette.easyStroke}" stroke-width="1.6" />`
            );
            const epLabelPt = proj3(toLock([-0.62, -0.62, 0], 'Cryst'));
            layers.labels.push(
                `<text x="${fmt(epLabelPt.x)}" y="${fmt(epLabelPt.y)}" fill="${palette.easyStroke}" font-family="Inter, Arial, sans-serif" font-size="11" font-weight="600">${escXml(lbl('(0001) Easy-Plane', '$(0001)$ Easy-Plane'))}</text>`
            );
        }

        // Substrate slab in F_Sub
        const wx = 0.62;
        const wy = 0.52;
        const topSub = [
            [-wx, -wy, 0],
            [wx, -wy, 0],
            [wx, wy, 0],
            [-wx, wy, 0]
        ];
        const topPts = topSub.map((c) => proj3(toLock(c, 'Sub')));
        layers.substrate.push(
            `<polygon points="${topPts.map((p) => `${fmt(p.x)},${fmt(p.y)}`).join(' ')}" fill="${palette.chipFill}" fill-opacity="0.82" stroke="${palette.chipStroke}" stroke-width="1.8" />`
        );

        // Line of nodes on substrate
        if (showEasyPlane && Math.abs(stateA.setup.thetaCutDeg) > 1e-3) {
            const nodeSub = stateA.vectors.nodesHat.Sub;
            const p1 = proj3(toLock(scale3(nodeSub, -0.58), 'Sub'));
            const p2 = proj3(toLock(scale3(nodeSub, 0.58), 'Sub'));
            layers.crystal.push(
                `<line x1="${fmt(p1.x)}" y1="${fmt(p1.y)}" x2="${fmt(p2.x)}" y2="${fmt(p2.y)}" stroke="${palette.easyStroke}" stroke-width="1.5" stroke-dasharray="5,3" />`
            );
        }

        // Microstructure (IDT / CPW) as true 3D polygons in F_Wave (gold width === gap width)
        if (showMicrostructure) {
            let polys = [];
            if (stateA.setup.deviceType === 'IDT') {
                const w = 0.042;
                const x0 = -5.5 * w;
                const yOuter = 0.28;
                const yBus = yOuter - 1.35 * w;
                const yTipTop = -yBus + w;
                const yTipBot = yBus - w;
                polys = [
                    [
                        [x0, yOuter], [-x0, yOuter], [-x0, yBus],
                        [x0 + 9 * w, yBus], [x0 + 9 * w, yTipTop], [x0 + 8 * w, yTipTop], [x0 + 8 * w, yBus],
                        [x0 + 5 * w, yBus], [x0 + 5 * w, yTipTop], [x0 + 4 * w, yTipTop], [x0 + 4 * w, yBus],
                        [x0 + 1 * w, yBus], [x0 + 1 * w, yTipTop], [x0 + 0 * w, yTipTop], [x0, yBus]
                    ],
                    [
                        [x0, -yOuter], [-x0, -yOuter], [-x0, -yBus],
                        [x0 + 11 * w, -yBus], [x0 + 11 * w, yTipBot], [x0 + 10 * w, yTipBot], [x0 + 10 * w, -yBus],
                        [x0 + 7 * w, -yBus], [x0 + 7 * w, yTipBot], [x0 + 6 * w, yTipBot], [x0 + 6 * w, -yBus],
                        [x0 + 3 * w, -yBus], [x0 + 3 * w, yTipBot], [x0 + 2 * w, yTipBot], [x0 + 2 * w, -yBus], [x0, -yBus]
                    ]
                ];
            } else {
                const w = 0.056;
                const yBot = -0.35;
                const yTop = 0.35;
                const yInner = yTop - w;
                polys = [[
                    [-2.5 * w, yTop], [2.5 * w, yTop], [2.5 * w, yBot],
                    [1.5 * w, yBot], [1.5 * w, yInner], [0.5 * w, yInner],
                    [0.5 * w, yBot], [-0.5 * w, yBot], [-0.5 * w, yInner],
                    [-1.5 * w, yInner], [-1.5 * w, yBot], [-2.5 * w, yBot]
                ]];
            }

            for (const poly of polys) {
                const ptsStr = poly
                    .map(([xk, yt]) => {
                        const p = proj3(toLock([xk, yt, 0], 'Wave'));
                        return `${fmt(p.x)},${fmt(p.y)}`;
                    })
                    .join(' ');
                layers.device.push(
                    `<polygon points="${ptsStr}" fill="${palette.goldStroke}" stroke="none" />`
                );
            }
        }

        const origin2 = proj3([0, 0, 0]);

        // Sweep traces: tips of m̂ (scaled with |B_easy|/|B| like its arrow) and B
        addSweepTrace(
            (s) => {
                if (!mDefined(s)) return null;
                const len = norm3(s.vectors.bEasyProjVec[F]) / (s.metrics.b0Tesla || 1);
                return len > 0.04 ? proj3(scale3(s.vectors.mHat[F], len)) : null;
            },
            palette.m,
            'arrow-m',
            origin2
        );
        addSweepTrace((s) => proj3(s.vectors.bHat[F]), palette.B, 'arrow-B', origin2, true);

        const addVec3 =(v3, scaleLen, color, markerId, text, sw = 2.5, dashed = false) => {
            const pEnd = proj3(scale3(v3, scaleLen));
            addArrow2D(origin2.x, origin2.y, pEnd.x, pEnd.y, color, markerId, text, sw, dashed);
        };

        // Crystal axes with Miller-Bravais notation
        if (showCrystalAxes) {
            addVec3(stateA.vectors.eaHat[F], 0.92, palette.ea, 'arrow-ea', lbl('a [112\u03050]', '$a\\,[11\\bar{2}0]$'), 1.6, true);
            addVec3(stateA.vectors.emHat[F], 0.88, palette.em, 'arrow-em', lbl('m [11\u030500]', '$m\\,[1\\bar{1}00]$'), 1.5, true);
            addVec3(stateA.vectors.ecHat[F], 1.00, palette.ec, 'arrow-ec', lbl('c [0001]', '$c\\,[0001]$'), 2.1, false);
        }

        // Wave vector k_hat
        addVec3(stateA.vectors.kHat[F], 1.00, palette.k, 'arrow-k', lbl('k̂', '$\\hat{k}$'), 2.8, false);

        // B projection onto Easy-Plane & drop-line from tip of B to tip of m_hat
        const bHatF = stateA.vectors.bHat[F];
        const bEasyF = scale3(stateA.vectors.bEasyProjVec[F], 1 / (stateA.metrics.b0Tesla || 1));
        const bEasyLen = norm3(bEasyF);
        const diffOut = Math.hypot(bHatF[0] - bEasyF[0], bHatF[1] - bEasyF[1], bHatF[2] - bEasyF[2]);
        if (bEasyLen > 0.04 && diffOut > 0.05) {
            const pB = proj3(bHatF);
            const pBe = proj3(bEasyF);
            layers.vectors.push(
                `<line x1="${fmt(pB.x)}" y1="${fmt(pB.y)}" x2="${fmt(pBe.x)}" y2="${fmt(pBe.y)}" stroke="${palette.B_easy}" stroke-width="1.3" stroke-dasharray="4,3" />`
            );
        }

        // Equilibrium vectors m_hat (scaled to Easy-Plane projection) and l_hat
        if (showEquilibrium && !stateA.metrics.isEasyPlaneSingular) {
            addVec3(stateA.vectors.lHat[F], 1.00, palette.l, 'arrow-l', lbl('l̂', '$\\hat{l}$'), 2.6, false);
            if (bEasyLen > 0.04) {
                addVec3(stateA.vectors.mHat[F], bEasyLen, palette.m, 'arrow-m', lbl('m̂', '$\\hat{m}$'), 2.6, false);
            }
        }

        // Applied magnetic field B
        addVec3(bHatF, 1.00, palette.B, 'arrow-B', lbl('B', '$\\mathbf{B}$'), 3.2, false);

        // 3D Angle Arc between B and k_hat
        if (showAngleArcs) {
            const u1 = normalize3(bHatF);
            const u2 = normalize3(stateA.vectors.kHat[F]);
            const omega = Math.acos(Math.max(-1, Math.min(1, dot3(u1, u2))));
            if (omega > 0.05 && omega < Math.PI - 0.05) {
                const sinO = Math.sin(omega);
                const arcPts = [];
                for (let i = 0; i <= 18; i++) {
                    const t = i / 18;
                    const a = Math.sin((1 - t) * omega) / sinO;
                    const b = Math.sin(t * omega) / sinO;
                    const pt = proj3(scale3(add3(scale3(u1, a), scale3(u2, b)), 0.42));
                    arcPts.push(`${fmt(pt.x)},${fmt(pt.y)}`);
                }
                layers.arcs.push(
                    `<polyline points="${arcPts.join(' ')}" fill="none" stroke="${palette.B}" stroke-width="1.5" stroke-dasharray="3,2" />`
                );
            }
        }

        // Panel title
        const frameMeta = FRAME_METADATA[F] || FRAME_METADATA.Sub;
        layers.labels.push(
            `<text x="${fmt(cx)}" y="28" fill="${palette.textMain}" font-family="Inter, Arial, sans-serif" font-size="14" font-weight="700" text-anchor="middle">${escXml(`3D-Geometrie (${frameMeta.name})`)}</text>`
        );
    };

    // ------------------------------------------------------------------------
    // Panel 2: 2D Sample Surface (Chip Plane F_Sub)
    // ------------------------------------------------------------------------
    const render2DSamplePanel = (cx, cy, rMax) => {
        const cw = rMax * 1.48;
        const ch = rMax * 1.25;
        layers.substrate.push(
            `<rect x="${fmt(cx - cw / 2)}" y="${fmt(cy - ch / 2)}" width="${fmt(cw)}" height="${fmt(ch)}" rx="4" fill="${palette.chipFill}" stroke="${palette.chipStroke}" stroke-width="1.8" />`
        );
        layers.bg.push(
            `<circle cx="${fmt(cx)}" cy="${fmt(cy)}" r="${fmt(rMax)}" fill="none" stroke="${palette.grid}" stroke-width="1.2" stroke-dasharray="4,4" />`
        );

        // Line of nodes
        if (showEasyPlane && Math.abs(stateA.setup.thetaCutDeg) > 1e-3) {
            const nd = stateA.vectors.nodesHat.Sub;
            layers.crystal.push(
                `<line x1="${fmt(cx - nd[0] * rMax)}" y1="${fmt(cy + nd[1] * rMax)}" x2="${fmt(cx + nd[0] * rMax)}" y2="${fmt(cy - nd[1] * rMax)}" stroke="${palette.easyStroke}" stroke-width="1.5" stroke-dasharray="5,3" />`
            );
        }

        const drawSubVec = (v, color, marker, text, sw = 2.4, dashed = false) => {
            addArrow2D(cx, cy, cx + v[0] * rMax, cy - v[1] * rMax, color, marker, text, sw, dashed);
        };

        const subTip = (v) => ({ x: cx + v[0] * rMax, y: cy - v[1] * rMax });
        addSweepTrace((s) => (mDefined(s) ? subTip(s.vectors.mHat.Sub) : null), palette.m, 'arrow-m', { x: cx, y: cy });
        addSweepTrace((s) => subTip(s.vectors.bHat.Sub), palette.B, 'arrow-B', { x: cx, y: cy }, true);

        drawSubVec([1, 0, 0], palette.textMuted, 'arrow-t', lbl('x̂_s', '$\\hat{x}_s$'), 1.3, true);
        drawSubVec([0, 1, 0], palette.textMuted, 'arrow-t', lbl('ŷ_s', '$\\hat{y}_s$'), 1.3, true);
        drawSubVec(stateA.vectors.tHat.Sub, palette.t, 'arrow-t', lbl('t̂', '$\\hat{t}$'), 1.5, true);
        drawSubVec(stateA.vectors.kHat.Sub, palette.k, 'arrow-k', lbl('k̂', '$\\hat{k}$'), 2.8, false);

        if (showEquilibrium && !stateA.metrics.isEasyPlaneSingular) {
            drawSubVec(stateA.vectors.lHat.Sub, palette.l, 'arrow-l', lbl('l̂_∥s', '$\\hat{l}_{\\parallel s}$'), 2.2, false);
            drawSubVec(stateA.vectors.mHat.Sub, palette.m, 'arrow-m', lbl('m̂_∥s', '$\\hat{m}_{\\parallel s}$'), 2.2, false);
        }
        drawSubVec(stateA.vectors.bHat.Sub, palette.B, 'arrow-B', lbl('B_∥s', '$\\mathbf{B}_{\\parallel s}$'), 3.0, false);

        layers.labels.push(
            `<text x="${fmt(cx)}" y="${fmt(cy - rMax - 22)}" fill="${palette.textMain}" font-family="Inter, Arial, sans-serif" font-size="13" font-weight="700" text-anchor="middle">${escXml('Probenebene (ℱ_Sub)')}</text>`
        );
    };

    // ------------------------------------------------------------------------
    // Panel 3: 2D Basal Easy-Plane (0001) (F_Cryst)
    // ------------------------------------------------------------------------
    const render2DEasyPanel = (cx, cy, rMax) => {
        const hexPts = [];
        for (let i = 0; i < 6; i++) {
            const ang = i * 60 * DEG2RAD;
            hexPts.push(`${fmt(cx + rMax * Math.cos(ang))},${fmt(cy - rMax * Math.sin(ang))}`);
        }
        layers.crystal.push(
            `<polygon points="${hexPts.join(' ')}" fill="${palette.easyFill}" fill-opacity="0.55" stroke="${palette.easyStroke}" stroke-width="1.8" />`
        );

        const drawCrystVec = (v, color, marker, text, sw = 2.4, dashed = false) => {
            addArrow2D(cx, cy, cx + v[0] * rMax, cy - v[1] * rMax, color, marker, text, sw, dashed);
        };

        const crystTip = (v) => ({ x: cx + v[0] * rMax, y: cy - v[1] * rMax });
        addSweepTrace((s) => (mDefined(s) ? crystTip(s.vectors.mHat.Cryst) : null), palette.m, 'arrow-m', { x: cx, y: cy });
        addSweepTrace((s) => crystTip(s.vectors.bHat.Cryst), palette.B_easy, 'arrow-Beasy', { x: cx, y: cy }, true);

        drawCrystVec([1, 0, 0], palette.ea, 'arrow-ea', lbl('ê_a [112\u03050]', '$\\hat{e}_a$'), 1.4, true);
        drawCrystVec([0, 1, 0], palette.em, 'arrow-em', lbl('ê_m [11\u030500]', '$\\hat{e}_m$'), 1.4, true);
        drawCrystVec(stateA.vectors.kHat.Cryst, palette.k, 'arrow-k', lbl('k̂_⊥c', '$\\hat{k}_{\\perp c}$'), 2.3, true);

        if (showEquilibrium && !stateA.metrics.isEasyPlaneSingular) {
            drawCrystVec(stateA.vectors.lHat.Cryst, palette.l, 'arrow-l', lbl('l̂', '$\\hat{l}$'), 2.7, false);
            drawCrystVec(stateA.vectors.mHat.Cryst, palette.m, 'arrow-m', lbl('m̂', '$\\hat{m}$'), 2.7, false);
        }
        drawCrystVec(stateA.vectors.bHat.Cryst, palette.B_easy, 'arrow-Beasy', lbl('B_⊥c', '$\\mathbf{B}_{\\perp c}$'), 3.0, false);

        layers.labels.push(
            `<text x="${fmt(cx)}" y="${fmt(cy - rMax - 22)}" fill="${palette.textMain}" font-family="Inter, Arial, sans-serif" font-size="13" font-weight="700" text-anchor="middle">${escXml('Basalebene (0001) Easy-Plane')}</text>`
        );
    };

    // Dispatch layout
    if (layout === 'combined') {
        render3DPanel(300, 265, 195);
        render2DSamplePanel(825, 140, 88);
        render2DEasyPanel(825, 370, 88);
    } else if (layout === '3d') {
        render3DPanel(width * 0.5, height * 0.5, 195);
    } else if (layout === '2d_sample') {
        render2DSamplePanel(width * 0.5, height * 0.48, 175);
    } else if (layout === '2d_easy') {
        render2DEasyPanel(width * 0.5, height * 0.48, 175);
    }

    // Bottom Summary / Legend Strip
    if (showLegend) {
        const s = stateA.setup;
        const m = stateA.metrics;
        const infoLine = `${s.name} | ${s.mountMode}-Halterung | ${s.cutPreset} (θ_cut=${s.thetaCutDeg.toFixed(2)}°) | ${s.deviceType} θ_dev=${s.thetaDevDeg.toFixed(1)}° | ψ_S=${s.psiSDeg.toFixed(1)}°, φ_M=${s.phiMDeg.toFixed(1)}° | B_⊥c/B₀=${(m.bEasyRatio * 100).toFixed(1)}%, ∠(B,k̂)=${m.angleBKDeg.toFixed(1)}°, ∠(m̂,k̂)=${m.angleMKDeg.toFixed(1)}°, ∠(l̂,k̂)=${m.angleLKDeg.toFixed(1)}°${sweepStates.length ? ` | Sweep φ_M: ${fmtDeg(sweepFrom)}°→${fmtDeg(sweepTo)}°` : ''}`;
        // Monospace 10.5px ≈ 6.4px per character: split into two rows at a " | " if it would overflow
        const parts = infoLine.split(' | ');
        let rows = [infoLine];
        if (infoLine.length * 6.4 > width - 48 && parts.length > 1) {
            let cut = 1;
            while (cut < parts.length - 1 && parts.slice(0, cut + 1).join(' | ').length <= infoLine.length / 2) cut++;
            rows = [parts.slice(0, cut).join(' | '), parts.slice(cut).join(' | ')];
        }
        const stripH = 10 + rows.length * 14;
        layers.labels.push(
            `<rect x="16" y="${height - 12 - stripH}" width="${width - 32}" height="${stripH}" rx="4" fill="#f8fafc" stroke="#cbd5e1" stroke-width="1" />`,
            ...rows.map(
                (row, i) =>
                    `<text x="${width / 2}" y="${height - 12 - stripH + 16 + i * 14}" fill="${palette.textMain}" font-family="JetBrains Mono, Consolas, monospace" font-size="10.5" text-anchor="middle">${escXml(row)}</text>`
            )
        );
    }

    return `<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<svg xmlns="http://www.w3.org/2000/svg"
     xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape"
     viewBox="0 0 ${width} ${height}"
     width="${width}" height="${height}">
  <defs>${defsContent}
  </defs>
  <g inkscape:groupmode="layer" id="layer_01_bg" inkscape:label="01_Background_Grid">
    ${layers.bg.join('\n    ')}
  </g>
  <g inkscape:groupmode="layer" id="layer_02_crystal" inkscape:label="02_Crystal_EasyPlane">
    ${layers.crystal.join('\n    ')}
  </g>
  <g inkscape:groupmode="layer" id="layer_03_substrate" inkscape:label="03_Substrate_Chip">
    ${layers.substrate.join('\n    ')}
  </g>
  <g inkscape:groupmode="layer" id="layer_04_device" inkscape:label="04_Microstructure_IDT_CPW">
    ${layers.device.join('\n    ')}
  </g>${
      layers.sweep.length
          ? `
  <g inkscape:groupmode="layer" id="layer_04b_sweep" inkscape:label="04b_Field_Sweep">
    ${layers.sweep.join('\n    ')}
  </g>`
          : ''
  }
  <g inkscape:groupmode="layer" id="layer_05_vectors" inkscape:label="05_Vectors_Projections">
    ${layers.vectors.join('\n    ')}
  </g>
  <g inkscape:groupmode="layer" id="layer_06_arcs" inkscape:label="06_Angle_Arcs">
    ${layers.arcs.join('\n    ')}
  </g>
  <g inkscape:groupmode="layer" id="layer_07_labels" inkscape:label="07_Labels_Legend">
    ${layers.labels.join('\n    ')}
  </g>
</svg>`;
}
