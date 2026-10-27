/**
 * visualizer_3d.js
 * Clean, didactic 3D visualization engine for the Hematite
 * Spin-Wave (CPW) & SAW (IDT) Coordinate Calculator.
 *
 * Designed to avoid visual overload:
 *   - Progressive layer toggles (essential vectors B & k_hat first, optional
 *     easy-plane moments m_hat & l_hat, optional crystal axes & magnet coils).
 *   - Emphasis system: hovering an item anywhere in the UI dims everything else
 *     in the 3D scene (and temporarily shows hidden helper vectors like t_hat, n_hat).
 *   - Frame-Lock support (Sub, Wave, Cryst, Lab) with smooth animated transitions,
 *     so one can literally watch how the frames are rotated against each other.
 */

import {
    DEG2RAD,
    add3,
    scale3,
    dot3,
    norm3,
    normalize3,
    matVec3,
    rotZ,
    FRAME_METADATA
} from './math_core.js';

// Which visibility flag an emphasized key forces on (so hovering e.g. "t̂" in the
// vector table shows t̂ in 3D even though it is normally hidden).
const EMPHASIS_FORCES_VISIBILITY = {
    B: 'showB',
    m: 'showM',
    l: 'showL',
    k: 'showK',
    t: 'showT',
    n: 'showN',
    xs: 'showEdges',
    ys: 'showEdges',
    easy: 'showEasyPlane',
    nodes: 'showEasyPlane',
    crystal: 'showCrystalAxes',
    ea: 'showCrystalAxes',
    em: 'showCrystalAxes',
    ec: 'showCrystalAxes',
    magnet: 'showLabMagnet',
    sample: 'showSampleSlab',
    device: 'showMicrostructure',
    custom: 'showCustom'
};

// --- Small quaternion helpers for smooth frame-lock transitions ---------------
function matToQuat(M) {
    const tr = M[0][0] + M[1][1] + M[2][2];
    let w, x, y, z;
    if (tr > 0) {
        const s = Math.sqrt(tr + 1.0) * 2;
        w = 0.25 * s;
        x = (M[2][1] - M[1][2]) / s;
        y = (M[0][2] - M[2][0]) / s;
        z = (M[1][0] - M[0][1]) / s;
    } else if (M[0][0] > M[1][1] && M[0][0] > M[2][2]) {
        const s = Math.sqrt(1.0 + M[0][0] - M[1][1] - M[2][2]) * 2;
        w = (M[2][1] - M[1][2]) / s;
        x = 0.25 * s;
        y = (M[0][1] + M[1][0]) / s;
        z = (M[0][2] + M[2][0]) / s;
    } else if (M[1][1] > M[2][2]) {
        const s = Math.sqrt(1.0 + M[1][1] - M[0][0] - M[2][2]) * 2;
        w = (M[0][2] - M[2][0]) / s;
        x = (M[0][1] + M[1][0]) / s;
        y = 0.25 * s;
        z = (M[1][2] + M[2][1]) / s;
    } else {
        const s = Math.sqrt(1.0 + M[2][2] - M[0][0] - M[1][1]) * 2;
        w = (M[1][0] - M[0][1]) / s;
        x = (M[0][2] + M[2][0]) / s;
        y = (M[1][2] + M[2][1]) / s;
        z = 0.25 * s;
    }
    const n = Math.hypot(w, x, y, z) || 1;
    return [w / n, x / n, y / n, z / n];
}

function quatToMat([w, x, y, z]) {
    return [
        [1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
        [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
        [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)]
    ];
}

function slerpQuat(a, b, t) {
    let cosHalf = a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3];
    let bb = b;
    if (cosHalf < 0) {
        bb = b.map((v) => -v);
        cosHalf = -cosHalf;
    }
    if (cosHalf > 0.9995) {
        const out = a.map((v, i) => v + t * (bb[i] - v));
        const n = Math.hypot(...out) || 1;
        return out.map((v) => v / n);
    }
    const half = Math.acos(cosHalf);
    const sinHalf = Math.sqrt(1 - cosHalf * cosHalf);
    const ra = Math.sin((1 - t) * half) / sinHalf;
    const rb = Math.sin(t * half) / sinHalf;
    return a.map((v, i) => v * ra + bb[i] * rb);
}

function easeInOutCubic(t) {
    return t < 0.5 ? 4 * t * t * t : 1 - Math.pow(-2 * t + 2, 3) / 2;
}

export const VECTOR_COLORS = {
    B: '#38bdf8',        // Sky cyan: Applied magnetic field B
    B_easy: '#fbbf24',   // Amber gold: Basal easy-plane projection B_perp_c
    B_sample: '#7dd3fc', // Light cyan: Sample-plane projection B_||s
    m: '#fb7185',        // Rose coral: Weak FM moment m_hat || B_easy
    l: '#34d399',        // Emerald green: Neel vector l_hat = m_hat x e_c
    k: '#e879f9',        // Fuchsia/magenta: Wave vector k_hat (SAW / Spin Wave)
    t: '#94a3b8',        // Slate: Transverse finger/CPW axis t_hat
    n: '#60a5fa',        // Blue: Substrate surface normal n_hat = z_s
    ec: '#a78bfa',       // Violet: Crystal c-axis [0001] (e_c)
    ea: '#c084fc',       // Light purple: Crystal a-axis [11-20]
    em: '#818cf8',       // Indigo: Crystal m-axis [1-100]
    edge: '#cbd5e1',     // Light slate: sample edges x_s / y_s
    nodes: '#c084fc',    // Purple: intersection line sample plane ∩ easy plane
    custom: '#fb923c'    // Orange: user-defined direction ("Eigene Richtung")
};

export class Scene3DController {
    constructor(canvas3d) {
        this.canvas3d = canvas3d;

        this.lockFrame = 'Sub'; // 'Sub' | 'Wave' | 'Cryst' | 'Lab'
        this._lockAnim = null;

        // null or Set of keys ('B','m','l','k','t','n','xs','ys','easy','nodes',
        // 'crystal','ea','em','ec','magnet','sample','device','custom')
        this.emphasis = null;

        // Optional user-defined direction { Lab: [x,y,z], label: '[uvtw]' }
        this.customVector = null;
        // Directions the user pinned to the list: [{ Lab, label, color, key }]
        this.pinnedVectors = [];
        this.camera = {
            azimuthDeg: -38,
            elevationDeg: 28,
            zoom: 1.0,
            perspective: true
        };
        this.onCameraChange = null;

        // Calm, uncluttered defaults: show B, k_hat, easy plane, equilibrium m/l,
        // and the magnet poles in 3D; keep crystal basis axes off until toggled.
        this.visibility = {
            showB: true,
            showBProjections: true,
            showM: true,
            showL: true,
            showK: true,
            showEasyPlane: true,
            showCrystalAxes: false,
            showLabMagnet: true,
            showSampleSlab: true,
            showMicrostructure: true,
            showCompareOverlay: true,
            showCustom: false,
            // Helper vectors that normally stay hidden and only appear when emphasized
            showT: false,
            showN: false,
            showEdges: false
        };

        this.isDragging = false;
        this.lastMouse = { x: 0, y: 0 };

        this.stateA = null;
        this.stateB = null;
        this.compareActive = false;

        this._bindEvents();
    }

    _bindEvents() {
        if (!this.canvas3d) return;

        this.canvas3d.addEventListener('pointerdown', (e) => {
            this.isDragging = true;
            this.lastMouse = { x: e.clientX, y: e.clientY };
            this.canvas3d.setPointerCapture(e.pointerId);
        });

        this.canvas3d.addEventListener('pointermove', (e) => {
            if (!this.isDragging) return;
            const dx = e.clientX - this.lastMouse.x;
            const dy = e.clientY - this.lastMouse.y;
            this.lastMouse = { x: e.clientX, y: e.clientY };

            this.camera.azimuthDeg = (this.camera.azimuthDeg - dx * 0.55) % 360;
            this.camera.elevationDeg = Math.max(-89, Math.min(89, this.camera.elevationDeg + dy * 0.5));
            if (Math.abs(dy) > 1 && this.camera.elevationDeg < 87) {
                this.camera.perspective = true;
            }
            if (typeof this.onCameraChange === 'function') {
                this.onCameraChange();
            }
            this.render();
        });

        const endDrag = (e) => {
            if (this.isDragging) {
                this.isDragging = false;
                try {
                    this.canvas3d.releasePointerCapture(e.pointerId);
                } catch (_) {
                    // ignore
                }
            }
        };
        this.canvas3d.addEventListener('pointerup', endDrag);
        this.canvas3d.addEventListener('pointercancel', endDrag);

        this.canvas3d.addEventListener(
            'wheel',
            (e) => {
                e.preventDefault();
                const factor = e.deltaY < 0 ? 1.08 : 0.92;
                this.camera.zoom = Math.max(0.55, Math.min(2.2, this.camera.zoom * factor));
                this.render();
            },
            { passive: false }
        );
    }

    isTopDown2D() {
        return !this.camera.perspective && this.camera.elevationDeg >= 88;
    }

    setTopDown2D(enable) {
        if (enable) {
            this.camera.azimuthDeg = -90;
            this.camera.elevationDeg = 90;
            this.camera.perspective = false;
        } else {
            this.camera.azimuthDeg = -38;
            this.camera.elevationDeg = 28;
            this.camera.perspective = true;
        }
        if (typeof this.onCameraChange === 'function') {
            this.onCameraChange();
        }
        this.render();
    }

    setCameraPreset(preset) {
        switch (preset) {
            case 'iso':
                this.camera.azimuthDeg = -38;
                this.camera.elevationDeg = 28;
                break;
            case 'top':
                this.camera.azimuthDeg = -90;
                this.camera.elevationDeg = 85;
                break;
            case 'front':
                this.camera.azimuthDeg = -90;
                this.camera.elevationDeg = 5;
                break;
            case 'side':
                this.camera.azimuthDeg = 0;
                this.camera.elevationDeg = 5;
                break;
            default:
                this.camera.azimuthDeg = -38;
                this.camera.elevationDeg = 28;
        }
        this.render();
    }

    updateStates(stateA, stateB = null, compareActive = false) {
        this.stateA = stateA;
        this.stateB = stateB;
        this.compareActive = compareActive;
        this.render();
    }

    /**
     * Switches the frame that stands still in the picture. With animate=true the
     * scene smoothly rotates (quaternion slerp) from the old to the new lock frame.
     */
    setLockFrame(frame, animate = true) {
        if (frame === this.lockFrame && !this._lockAnim) return;
        const fromMatrix = this._currentLockMatrix();
        this.lockFrame = frame;
        if (!animate || !fromMatrix) {
            this._lockAnim = null;
            this.render();
            return;
        }
        this._lockAnim = { qFrom: matToQuat(fromMatrix), t0: performance.now(), duration: 700 };
        const step = () => {
            if (!this._lockAnim) return;
            const done = performance.now() - this._lockAnim.t0 >= this._lockAnim.duration;
            if (done) this._lockAnim = null;
            this.render();
            if (!done) requestAnimationFrame(step);
        };
        requestAnimationFrame(step);
    }

    setEmphasis(keys) {
        this.emphasis = keys && keys.length ? new Set(keys) : null;
        this.render();
    }

    // Rotation L that maps Lab coordinates into the coordinates of the "standing still" frame.
    _targetLockMatrix(state, frame) {
        if (frame === 'Lab') {
            // Lock the rotating electromagnet: B always points along +X on screen
            return rotZ(-(state.setup.phiMDeg || 0));
        }
        return state.matrices.getRotationMatrix(frame, 'Lab');
    }

    _currentLockMatrix() {
        if (!this.stateA) return null;
        const target = this._targetLockMatrix(this.stateA, this.lockFrame);
        if (!this._lockAnim) return target;
        const t = Math.min(1, (performance.now() - this._lockAnim.t0) / this._lockAnim.duration);
        return quatToMat(slerpQuat(this._lockAnim.qFrom, matToQuat(target), easeInOutCubic(t)));
    }

    _effectiveVisibility() {
        const vis = { ...this.visibility };
        if (this.emphasis) {
            for (const key of this.emphasis) {
                const prop = EMPHASIS_FORCES_VISIBILITY[key];
                if (prop) vis[prop] = true;
            }
        }
        return vis;
    }

    _alphaFor(keys) {
        if (!this.emphasis || !keys) return 1;
        return keys.some((k) => this.emphasis.has(k)) ? 1 : 0.12;
    }

    _prepareCanvas(canvas) {
        if (!canvas) return null;
        const rect = canvas.getBoundingClientRect();
        if (rect.width < 10 || rect.height < 10) return null;
        const dpr = window.devicePixelRatio || 1;
        const w = Math.floor(rect.width);
        const h = Math.floor(rect.height);
        if (canvas.width !== Math.floor(w * dpr) || canvas.height !== Math.floor(h * dpr)) {
            canvas.width = Math.floor(w * dpr);
            canvas.height = Math.floor(h * dpr);
        }
        const ctx = canvas.getContext('2d');
        ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
        return { ctx, w, h };
    }

    _project3D(p, cx, cy, scale) {
        const az = this.camera.azimuthDeg * DEG2RAD;
        const el = this.camera.elevationDeg * DEG2RAD;
        const caz = Math.cos(az);
        const saz = Math.sin(az);
        const cel = Math.cos(el);
        const sel = Math.sin(el);

        // Right-handed turntable camera at spherical (az, el):
        // xCam = right on screen, dHoriz = horizontal distance into screen away from viewer
        const xCam = -saz * p[0] + caz * p[1];
        const dHoriz = -caz * p[0] - saz * p[1];
        const z1 = p[2];

        // When el > 0 (looking down from above), points further back (dHoriz > 0) are higher on screen (+yCam)
        const yCam = sel * dHoriz + cel * z1;
        const depth = cel * dHoriz - sel * z1;

        const perspFactor = this.camera.perspective ? 1 / (1 + depth * 0.14) : 1.0;
        return {
            x: cx + xCam * scale * perspFactor,
            y: cy - yCam * scale * perspFactor,
            depth
        };
    }

    render() {
        if (!this.stateA || !this.canvas3d) return;
        this._render3DScene();
    }

    // Offscreen render for the image export: own size, pixel scale, color scheme and background.
    // Hover emphasis is dropped so the exported picture never comes out half faded.
    renderToCanvas({ width, height, scale = 2, light = false, transparent = false }) {
        const canvas = document.createElement('canvas');
        canvas.width = Math.round(width * scale);
        canvas.height = Math.round(height * scale);
        if (!this.stateA) return canvas;
        const ctx = canvas.getContext('2d');
        ctx.setTransform(scale, 0, 0, scale, 0, 0);
        const savedEmphasis = this.emphasis;
        this.emphasis = null;
        this._exportTarget = { ctx, w: width, h: height, light, transparent };
        try {
            this._render3DScene();
        } finally {
            this._exportTarget = null;
            this.emphasis = savedEmphasis;
        }
        return canvas;
    }

    // ========================================================================
    // 3D Main Scene Renderer
    // ========================================================================
    _render3DScene() {
        const prep = this._exportTarget || this._prepareCanvas(this.canvas3d);
        if (!prep) return;
        const { ctx, w, h } = prep;
        const isLight = prep.light ?? document.body.classList.contains('theme-light');
        const cx = w * 0.5;
        const cy = h * 0.53;
        const baseScale = Math.min(w, h) * 0.33 * this.camera.zoom;

        ctx.clearRect(0, 0, w, h);
        if (!prep.transparent) {
            const bgGrad = ctx.createRadialGradient(cx, cy, 10, cx, cy, Math.max(w, h) * 0.72);
            if (isLight) {
                bgGrad.addColorStop(0, '#ffffff');
                bgGrad.addColorStop(1, '#f1f5f9');
            } else {
                bgGrad.addColorStop(0, '#131c31');
                bgGrad.addColorStop(1, '#090d16');
            }
            ctx.fillStyle = bgGrad;
            ctx.fillRect(0, 0, w, h);
        }

        const F = this.lockFrame;
        const stA = this.stateA;
        const vis = this._effectiveVisibility();
        this._vis = vis;

        // L maps Lab coordinates into the frame that currently stands still
        // (interpolated while a frame-lock transition is animating).
        const L = this._currentLockMatrix();
        const labFrom = {};
        const toLock = (vec, srcFrame) => {
            if (srcFrame === 'Lab') return matVec3(L, vec);
            if (!labFrom[srcFrame]) labFrom[srcFrame] = stA.matrices.getRotationMatrix('Lab', srcFrame);
            return matVec3(L, matVec3(labFrom[srcFrame], vec));
        };
        const vecInLock = (vecObj) => toLock(vecObj.Lab, 'Lab');
        const proj = (pt3) => this._project3D(pt3, cx, cy, baseScale);

        // 1. Reference grid in the locked frame
        this._drawLockedFrameGrid(ctx, proj, F, isLight);

        // 2. Depth-sorted 3D surfaces
        const renderQueue = [];
        if (vis.showLabMagnet) {
            this._enqueueLabMagnet(renderQueue, toLock, proj, stA, isLight);
        }
        if (vis.showEasyPlane) {
            this._enqueueEasyPlaneDisk(renderQueue, toLock, proj, isLight);
        }
        if (vis.showSampleSlab) {
            this._enqueueSubstrateAndDevice(renderQueue, toLock, proj, stA, isLight);
        }

        renderQueue.sort((a, b) => b.depth - a.depth);
        for (const item of renderQueue) {
            ctx.save();
            ctx.globalAlpha = this._alphaFor(item.keys);
            item.draw(ctx);
            ctx.restore();
        }

        // 3. Arrows (all transformed into the locked frame)
        const arrows3D = [];

        if (vis.showCrystalAxes) {
            arrows3D.push({ vec: vecInLock(stA.vectors.eaHat), baseLen: 0.92, color: VECTOR_COLORS.ea, label: 'a [112̅0]', width: 2.0, dashed: true, keys: ['ea', 'crystal'] });
            arrows3D.push({ vec: vecInLock(stA.vectors.emHat), baseLen: 0.88, color: VECTOR_COLORS.em, label: 'm [11̅00]', width: 1.8, dashed: true, keys: ['em', 'crystal'] });
            arrows3D.push({ vec: vecInLock(stA.vectors.ecHat), baseLen: 1.00, color: VECTOR_COLORS.ec, label: 'c [0001]', width: 2.3, dashed: false, keys: ['ec', 'crystal'] });
        }
        if (vis.showEdges) {
            arrows3D.push({ vec: vecInLock(stA.vectors.xsHat), baseLen: 0.8, color: VECTOR_COLORS.edge, label: 'x_s', width: 2.0, dashed: true, keys: ['xs', 'sample'] });
            arrows3D.push({ vec: vecInLock(stA.vectors.ysHat), baseLen: 0.8, color: VECTOR_COLORS.edge, label: 'y_s', width: 2.0, dashed: true, keys: ['ys', 'sample'] });
        }
        if (vis.showN) {
            arrows3D.push({ vec: vecInLock(stA.vectors.nHat), baseLen: 0.9, color: VECTOR_COLORS.n, label: 'n̂', width: 2.4, dashed: false, keys: ['n', 'sample'] });
        }
        if (vis.showT) {
            arrows3D.push({ vec: vecInLock(stA.vectors.tHat), baseLen: 0.85, color: VECTOR_COLORS.t, label: 't̂', width: 2.4, dashed: false, keys: ['t', 'device'] });
        }

        const bHat = vecInLock(stA.vectors.bHat);
        const bEasy = scale3(vecInLock(stA.vectors.bEasyProjVec), 1 / (stA.metrics.b0Tesla || 1));
        const bEasyLen = norm3(bEasy);

        // Dashed projection line from the tip of B to the tip of its easy-plane projection (= m̂ tip)
        const diffOutPlane = Math.hypot(bHat[0] - bEasy[0], bHat[1] - bEasy[1], bHat[2] - bEasy[2]);
        if (vis.showB && (vis.showM || vis.showEasyPlane) && bEasyLen > 0.04 && diffOutPlane > 0.05) {
            ctx.save();
            ctx.globalAlpha = this._alphaFor(['B', 'm', 'easy']);
            this._drawProjectionDropLine3D(ctx, proj, bHat, bEasy, VECTOR_COLORS.B_easy, 1.3);
            ctx.restore();
        }

        if (vis.showK) {
            arrows3D.push({ vec: vecInLock(stA.vectors.kHat), baseLen: 1.00, color: VECTOR_COLORS.k, label: 'k̂', width: 3.0, dashed: false, keys: ['k', 'device'] });
        }
        if (vis.showL && !stA.metrics.isEasyPlaneSingular) {
            arrows3D.push({ vec: vecInLock(stA.vectors.lHat), baseLen: 1.00, color: VECTOR_COLORS.l, label: 'l̂', width: 2.8, dashed: false, keys: ['l'] });
        }
        if (vis.showM && !stA.metrics.isEasyPlaneSingular && bEasyLen > 0.04) {
            // m̂ is scaled with |B_easy|/|B| so its tip sits where the projection line from B lands
            arrows3D.push({ vec: vecInLock(stA.vectors.mHat), baseLen: bEasyLen, color: VECTOR_COLORS.m, label: 'm̂', width: 2.9, dashed: false, keys: ['m'] });
        }
        if (vis.showB) {
            arrows3D.push({ vec: bHat, baseLen: 1.00, color: VECTOR_COLORS.B, label: 'B', width: 3.5, dashed: false, keys: ['B'] });
        }
        const pinned = vis.showCustom ? this.pinnedVectors.filter((p) => p.Lab && norm3(p.Lab) > 1e-9) : [];
        for (const p of pinned) {
            arrows3D.push({
                vec: toLock(normalize3(p.Lab), 'Lab'),
                baseLen: 0.9,
                color: p.color,
                label: p.label,
                width: 2.6,
                dashed: false,
                keys: [p.key, 'pinned']
            });
        }
        const customIsPinned = (cv) => pinned.some((p) => Math.abs(dot3(normalize3(p.Lab), normalize3(cv.Lab)) - 1) < 1e-6);
        if (vis.showCustom && this.customVector && this.customVector.Lab && norm3(this.customVector.Lab) > 1e-9 && !customIsPinned(this.customVector)) {
            arrows3D.push({
                vec: toLock(normalize3(this.customVector.Lab), 'Lab'),
                baseLen: 0.95,
                color: VECTOR_COLORS.custom,
                label: this.customVector.label || 'eigene',
                width: 3.0,
                dashed: false,
                keys: ['custom']
            });
        }

        // Compare overlay (Setup B)
        if (this.compareActive && this.stateB && vis.showCompareOverlay) {
            const stB = this.stateB;
            if (vis.showK) {
                arrows3D.push({ vec: toLock(stB.vectors.kHat.Lab, 'Lab'), baseLen: 0.92, color: '#f0abfc', label: 'k̂(B)', width: 2.0, dashed: true, keys: ['k'] });
            }
            if (vis.showB) {
                arrows3D.push({ vec: toLock(stB.vectors.bHat.Lab, 'Lab'), baseLen: 0.96, color: '#7dd3fc', label: 'B(B)', width: 2.2, dashed: true, keys: ['B'] });
            }
        }

        this._draw3DVectorsDeCollided(ctx, proj, arrows3D, isLight);
        this._draw3DHUD(ctx, w, h, stA, F, isLight);
    }

    _drawLockedFrameGrid(ctx, proj, lockFrame, isLight) {
        const meta = FRAME_METADATA[lockFrame] || FRAME_METADATA.Sub;
        const gridColor = isLight ? 'rgba(100, 116, 139, 0.16)' : 'rgba(148, 163, 184, 0.11)';
        ctx.strokeStyle = gridColor;
        ctx.lineWidth = 1;

        const r = 1.12;
        const steps = 4;
        for (let i = -steps; i <= steps; i++) {
            const t = (i / steps) * r;
            const p1 = proj([t, -r, 0]);
            const p2 = proj([t, r, 0]);
            ctx.beginPath();
            ctx.moveTo(p1.x, p1.y);
            ctx.lineTo(p2.x, p2.y);
            ctx.stroke();

            const q1 = proj([-r, t, 0]);
            const q2 = proj([r, t, 0]);
            ctx.beginPath();
            ctx.moveTo(q1.x, q1.y);
            ctx.lineTo(q2.x, q2.y);
            ctx.stroke();
        }

        // Unit circle in XY plane of locked frame
        ctx.beginPath();
        for (let deg = 0; deg <= 360; deg += 5) {
            const rad = deg * DEG2RAD;
            const pt = proj([Math.cos(rad), Math.sin(rad), 0]);
            if (deg === 0) ctx.moveTo(pt.x, pt.y);
            else ctx.lineTo(pt.x, pt.y);
        }
        ctx.strokeStyle = isLight ? 'rgba(71, 85, 105, 0.22)' : 'rgba(148, 163, 184, 0.2)';
        ctx.stroke();

        // Subtle vertical axis line without cryptic x_s/y_s/z_s text labels
        const axColor = isLight ? 'rgba(51, 65, 85, 0.3)' : 'rgba(203, 213, 225, 0.22)';
        const o = proj([0, 0, 0]);
        const pz = proj([0, 0, 1.08]);
        ctx.beginPath();
        ctx.moveTo(o.x, o.y);
        ctx.lineTo(pz.x, pz.y);
        ctx.strokeStyle = axColor;
        ctx.lineWidth = 1;
        ctx.stroke();
    }

    _enqueueLabMagnet(queue, toLock, proj, state, isLight) {
        const phiRad = state.setup.phiMDeg * DEG2RAD;
        const bDirLab = [Math.cos(phiRad), Math.sin(phiRad), 0];
        const perpLab = [-Math.sin(phiRad), Math.cos(phiRad), 0];

        for (const poleSign of [1, -1]) {
            const center = scale3(bDirLab, poleSign * 1.14);
            const wPole = 0.16;
            const hPole = 0.11;
            const cornersLab = [
                add3(center, add3(scale3(perpLab, -wPole), [0, 0, -hPole])),
                add3(center, add3(scale3(perpLab, wPole), [0, 0, -hPole])),
                add3(center, add3(scale3(perpLab, wPole), [0, 0, hPole])),
                add3(center, add3(scale3(perpLab, -wPole), [0, 0, hPole]))
            ];
            const pts = cornersLab.map((c) => proj(toLock(c, 'Lab')));
            const depth = pts.reduce((a, p) => a + p.depth, 0) / 4;
            queue.push({
                depth,
                keys: ['magnet', 'B'],
                draw: (ctx) => {
                    ctx.save();
                    ctx.beginPath();
                    pts.forEach((p, i) => (i === 0 ? ctx.moveTo(p.x, p.y) : ctx.lineTo(p.x, p.y)));
                    ctx.closePath();
                    ctx.fillStyle = poleSign > 0 ? 'rgba(239, 68, 68, 0.35)' : 'rgba(59, 130, 246, 0.35)';
                    ctx.strokeStyle = poleSign > 0 ? '#f87171' : '#60a5fa';
                    ctx.lineWidth = 1.3;
                    ctx.fill();
                    ctx.stroke();

                    const cPt = proj(toLock(center, 'Lab'));
                    ctx.fillStyle = isLight ? '#0f172a' : '#f8fafc';
                    ctx.font = '700 10px "Inter", sans-serif';
                    ctx.textAlign = 'center';
                    ctx.fillText(poleSign > 0 ? 'S' : 'N', cPt.x, cPt.y + 3);
                    ctx.restore();
                }
            });
        }
    }

    _enqueueEasyPlaneDisk(queue, toLock, proj, isLight) {
        const r = 0.9;
        const hexPts = [];
        for (let i = 0; i < 6; i++) {
            const ang = i * 60 * DEG2RAD;
            const pCryst = [r * Math.cos(ang), r * Math.sin(ang), 0];
            hexPts.push(proj(toLock(pCryst, 'Cryst')));
        }
        const depth = hexPts.reduce((a, p) => a + p.depth, 0) / 6;

        queue.push({
            depth,
            keys: ['easy', 'm', 'l', 'nodes'],
            draw: (ctx) => {
                ctx.save();
                ctx.beginPath();
                hexPts.forEach((p, i) => (i === 0 ? ctx.moveTo(p.x, p.y) : ctx.lineTo(p.x, p.y)));
                ctx.closePath();
                ctx.fillStyle = isLight ? 'rgba(168, 85, 247, 0.14)' : 'rgba(168, 85, 247, 0.16)';
                ctx.strokeStyle = isLight ? 'rgba(147, 51, 234, 0.65)' : 'rgba(192, 132, 252, 0.65)';
                ctx.lineWidth = 1.8;
                ctx.fill();
                ctx.stroke();
                ctx.restore();
            }
        });
    }

    _enqueueSubstrateAndDevice(queue, toLock, proj, state, isLight) {
        const wx = 0.62;
        const wy = 0.52;
        const zTop = 0.0;
        const zBot = -0.07;

        const topCornersSub = [
            [-wx, -wy, zTop],
            [wx, -wy, zTop],
            [wx, wy, zTop],
            [-wx, wy, zTop]
        ];
        const botCornersSub = [
            [-wx, -wy, zBot],
            [wx, -wy, zBot],
            [wx, wy, zBot],
            [-wx, wy, zBot]
        ];

        const topPts = topCornersSub.map((c) => proj(toLock(c, 'Sub')));
        const botPts = botCornersSub.map((c) => proj(toLock(c, 'Sub')));
        const topDepth = topPts.reduce((a, p) => a + p.depth, 0) / 4;

        for (let i = 0; i < 4; i++) {
            const j = (i + 1) % 4;
            const sidePoly = [topPts[i], topPts[j], botPts[j], botPts[i]];
            const d = sidePoly.reduce((a, p) => a + p.depth, 0) / 4;
            queue.push({
                depth: d + 0.02,
                keys: ['sample'],
                draw: (ctx) => {
                    ctx.save();
                    ctx.beginPath();
                    sidePoly.forEach((p, idx) => (idx === 0 ? ctx.moveTo(p.x, p.y) : ctx.lineTo(p.x, p.y)));
                    ctx.closePath();
                    ctx.fillStyle = isLight ? 'rgba(100, 116, 139, 0.35)' : 'rgba(51, 65, 85, 0.65)';
                    ctx.strokeStyle = isLight ? 'rgba(71, 85, 105, 0.6)' : 'rgba(148, 163, 184, 0.45)';
                    ctx.lineWidth = 1;
                    ctx.fill();
                    ctx.stroke();
                    ctx.restore();
                }
            });
        }

        queue.push({
            depth: topDepth,
            keys: null, // alpha is handled per part below
            draw: (ctx) => {
                const vis = this._vis || this.visibility;
                ctx.save();
                ctx.globalAlpha = this._alphaFor(['sample', 'device', 'k', 'nodes']);
                ctx.beginPath();
                topPts.forEach((p, idx) => (idx === 0 ? ctx.moveTo(p.x, p.y) : ctx.lineTo(p.x, p.y)));
                ctx.closePath();
                ctx.fillStyle = isLight ? 'rgba(226, 232, 240, 0.78)' : 'rgba(30, 41, 59, 0.78)';
                ctx.strokeStyle = isLight ? '#475569' : '#94a3b8';
                ctx.lineWidth = 1.8;
                ctx.fill();
                ctx.stroke();

                if (vis.showMicrostructure) {
                    ctx.globalAlpha = this._alphaFor(['device', 'k', 't']);
                    this._drawMicrostructure3D(ctx, toLock, proj, state.setup.deviceType);
                }

                if (vis.showEasyPlane && Math.abs(state.setup.thetaCutDeg) > 1e-3) {
                    ctx.globalAlpha = this._alphaFor(['nodes', 'easy', 'sample']);
                    const nodeSub = state.vectors.nodesHat.Sub;
                    const pN1 = proj(toLock(scale3(nodeSub, -0.58), 'Sub'));
                    const pN2 = proj(toLock(scale3(nodeSub, 0.58), 'Sub'));
                    ctx.beginPath();
                    ctx.moveTo(pN1.x, pN1.y);
                    ctx.lineTo(pN2.x, pN2.y);
                    ctx.strokeStyle = isLight ? '#9333ea' : '#c084fc';
                    ctx.setLineDash([5, 3]);
                    ctx.lineWidth = this.emphasis && this.emphasis.has('nodes') ? 2.6 : 1.6;
                    ctx.stroke();
                    ctx.setLineDash([]);
                }

                ctx.restore();
            }
        });
    }

    _getMicrostructurePolygonsWave(deviceType = 'CPW') {
        if (deviceType === 'IDT') {
            // Finger width === Gap width === w (metallization ratio 1:1)
            const w = 0.042;
            const x0 = -5.5 * w; // Total span = 11 * w across 6 fingers + 5 gaps
            const yOuter = 0.28;
            const yBus = yOuter - 1.35 * w;
            const yTipTop = -yBus + w;
            const yTipBot = yBus - w;

            const topComb = [
                [x0, yOuter],
                [-x0, yOuter],
                [-x0, yBus],
                [x0 + 9 * w, yBus],
                [x0 + 9 * w, yTipTop],
                [x0 + 8 * w, yTipTop],
                [x0 + 8 * w, yBus],
                [x0 + 5 * w, yBus],
                [x0 + 5 * w, yTipTop],
                [x0 + 4 * w, yTipTop],
                [x0 + 4 * w, yBus],
                [x0 + 1 * w, yBus],
                [x0 + 1 * w, yTipTop],
                [x0 + 0 * w, yTipTop],
                [x0, yBus]
            ];

            const botComb = [
                [x0, -yOuter],
                [-x0, -yOuter],
                [-x0, -yBus],
                [x0 + 11 * w, -yBus],
                [x0 + 11 * w, yTipBot],
                [x0 + 10 * w, yTipBot],
                [x0 + 10 * w, -yBus],
                [x0 + 7 * w, -yBus],
                [x0 + 7 * w, yTipBot],
                [x0 + 6 * w, yTipBot],
                [x0 + 6 * w, -yBus],
                [x0 + 3 * w, -yBus],
                [x0 + 3 * w, yTipBot],
                [x0 + 2 * w, yTipBot],
                [x0 + 2 * w, -yBus],
                [x0, -yBus]
            ];

            return [topComb, botComb];
        }

        // CPW (Shorted Coplanar Waveguide Antenna):
        // Conductor width === Gap width === w (metallization ratio 1:1)
        const w = 0.056;
        const yBot = -0.35;
        const yTop = 0.35;
        const yInner = yTop - w;

        const cpwPolygon = [
            [-2.5 * w, yTop],
            [2.5 * w, yTop],
            [2.5 * w, yBot],
            [1.5 * w, yBot],
            [1.5 * w, yInner],
            [0.5 * w, yInner],
            [0.5 * w, yBot],
            [-0.5 * w, yBot],
            [-0.5 * w, yInner],
            [-1.5 * w, yInner],
            [-1.5 * w, yBot],
            [-2.5 * w, yBot]
        ];
        return [cpwPolygon];
    }

    _drawMicrostructure3D(ctx, toLock, proj, deviceType = 'CPW') {
        ctx.save();
        ctx.fillStyle = '#fbbf24';
        const polys = this._getMicrostructurePolygonsWave(deviceType);

        for (const poly of polys) {
            ctx.beginPath();
            poly.forEach(([xk, yt], idx) => {
                const pt = proj(toLock([xk, yt, 0.002], 'Wave'));
                if (idx === 0) ctx.moveTo(pt.x, pt.y);
                else ctx.lineTo(pt.x, pt.y);
            });
            ctx.closePath();
            ctx.fill();
        }
        ctx.restore();
    }

    _draw3DAngleArcs(ctx, proj, state, lockFrame, isLight) {
        const bHat = state.vectors.bHat[lockFrame];
        const kHat = state.vectors.kHat[lockFrame];

        if (this.visibility.showB && this.visibility.showK) {
            this._drawSphericalArc3D(ctx, proj, bHat, kHat, 0.42, VECTOR_COLORS.B, '', isLight);
        }
    }

    _drawSphericalArc3D(ctx, proj, u1, u2, radius, color, label, isLight) {
        const v1 = normalize3(u1);
        const v2 = normalize3(u2);
        const dotVal = Math.max(-1, Math.min(1, dot3(v1, v2)));
        const omega = Math.acos(dotVal);
        if (omega < 0.04 || omega > Math.PI - 0.04) return;

        const sinOmega = Math.sin(omega);
        const N = 20;
        const pts = [];
        for (let i = 0; i <= N; i++) {
            const t = i / N;
            const a = Math.sin((1 - t) * omega) / sinOmega;
            const b = Math.sin(t * omega) / sinOmega;
            const p3 = scale3(add3(scale3(v1, a), scale3(v2, b)), radius);
            pts.push(proj(p3));
        }

        ctx.save();
        ctx.beginPath();
        pts.forEach((p, idx) => (idx === 0 ? ctx.moveTo(p.x, p.y) : ctx.lineTo(p.x, p.y)));
        ctx.strokeStyle = color;
        ctx.lineWidth = 1.6;
        ctx.setLineDash([2, 2]);
        ctx.stroke();
        ctx.setLineDash([]);
        ctx.restore();
    }

    _drawProjectionDropLine3D(ctx, proj, pStart3D, pEnd3D, color, lineWidth = 1.3) {
        const p1 = proj(pStart3D);
        const p2 = proj(pEnd3D);
        ctx.save();
        ctx.beginPath();
        ctx.moveTo(p1.x, p1.y);
        ctx.lineTo(p2.x, p2.y);
        ctx.strokeStyle = color;
        ctx.lineWidth = lineWidth;
        ctx.setLineDash([4, 3]);
        ctx.stroke();
        ctx.setLineDash([]);

        // Small 3D right-angle marker (∟) at the projection foot pEnd3D
        const uInPlane = normalize3(pEnd3D);
        const uPerp = normalize3([pStart3D[0] - pEnd3D[0], pStart3D[1] - pEnd3D[1], pStart3D[2] - pEnd3D[2]]);
        const cornerSize = Math.min(0.08, norm3(pEnd3D) * 0.25);
        if (cornerSize > 0.025) {
            const cA = proj(add3(pEnd3D, scale3(uInPlane, -cornerSize)));
            const cB = proj(add3(add3(pEnd3D, scale3(uInPlane, -cornerSize)), scale3(uPerp, cornerSize)));
            const cC = proj(add3(pEnd3D, scale3(uPerp, cornerSize)));
            ctx.beginPath();
            ctx.moveTo(cA.x, cA.y);
            ctx.lineTo(cB.x, cB.y);
            ctx.lineTo(cC.x, cC.y);
            ctx.strokeStyle = color;
            ctx.globalAlpha = 0.75;
            ctx.lineWidth = 1.1;
            ctx.stroke();
        }
        ctx.restore();
    }

    _draw3DVectorsDeCollided(ctx, proj, arrows3D, isLight) {
        const origin = proj([0, 0, 0]);
        // Group ONLY truly parallel 3D vectors (angle < 2.5°, dot > 0.999)
        const groups = [];
        for (const item of arrows3D) {
            const u = normalize3(item.vec);
            if (norm3(item.vec) < 0.03) continue;
            let foundGroup = null;
            for (const g of groups) {
                if (dot3(u, g.uRef) > 0.999) {
                    foundGroup = g;
                    break;
                }
            }
            if (!foundGroup) {
                foundGroup = { uRef: u, items: [] };
                groups.push(foundGroup);
            }
            foundGroup.items.push(item);
        }

        const labelClusters = [];
        for (const g of groups) {
            let maxTipLen = 0;
            g.items.forEach((it, idx) => {
                // Only stagger if multiple vectors truly share the exact same 3D direction
                const actualLen = g.items.length > 1
                    ? Math.max(0.18, it.baseLen - idx * 0.14)
                    : it.baseLen;
                maxTipLen = Math.max(maxTipLen, actualLen);
                const end3D = scale3(normalize3(it.vec), actualLen);
                const p2 = proj(end3D);
                const alpha = this._alphaFor(it.keys);
                const emphasized = this.emphasis && alpha === 1;
                ctx.save();
                ctx.globalAlpha = alpha;
                this._drawSingleShaftAndHead(ctx, origin.x, origin.y, p2.x, p2.y, it.color, it.width * (emphasized ? 1.35 : 1), it.dashed);
                ctx.restore();
            });
            const tipPt = proj(scale3(g.uRef, maxTipLen));
            const outerEnd = scale3(g.uRef, maxTipLen + 0.04);
            const pOuter = proj(outerEnd);
            const ang = Math.atan2(pOuter.y - origin.y, pOuter.x - origin.x);
            labelClusters.push({
                tipX: tipPt.x,
                tipY: tipPt.y,
                x: pOuter.x,
                y: pOuter.y,
                ang,
                alpha: Math.max(...g.items.map((it) => this._alphaFor(it.keys))),
                tokens: g.items.map((it) => ({ text: it.label, color: it.color }))
            });
        }

        this._drawGroupedPillBadges(ctx, labelClusters, isLight);
    }

    _drawSingleShaftAndHead(ctx, x1, y1, x2, y2, color, lineWidth = 2.5, dashed = false) {
        const dx = x2 - x1;
        const dy = y2 - y1;
        const len = Math.hypot(dx, dy);
        if (len < 3) return;

        ctx.save();
        ctx.beginPath();
        if (dashed) ctx.setLineDash([5, 3]);
        ctx.moveTo(x1, y1);
        ctx.lineTo(x2, y2);
        ctx.strokeStyle = color;
        ctx.lineWidth = lineWidth;
        ctx.lineCap = 'round';
        ctx.stroke();
        ctx.setLineDash([]);

        const headLen = Math.min(11, Math.max(7, lineWidth * 3.1));
        const angle = Math.atan2(dy, dx);
        ctx.beginPath();
        ctx.moveTo(x2, y2);
        ctx.lineTo(
            x2 - headLen * Math.cos(angle - Math.PI / 6.5),
            y2 - headLen * Math.sin(angle - Math.PI / 6.5)
        );
        ctx.lineTo(
            x2 - headLen * Math.cos(angle + Math.PI / 6.5),
            y2 - headLen * Math.sin(angle + Math.PI / 6.5)
        );
        ctx.closePath();
        ctx.fillStyle = color;
        ctx.fill();
        ctx.restore();
    }

    _parseOverbarText(rawText) {
        let clean = '';
        const overbarIndices = [];
        for (let i = 0; i < rawText.length; i++) {
            const ch = rawText[i];
            if (ch === '\u0305') continue;
            const idx = clean.length;
            clean += ch;
            if (rawText[i + 1] === '\u0305') {
                overbarIndices.push(idx);
                i++; // skip combining macron
            }
        }
        return { clean, overbarIndices };
    }

    _drawGroupedPillBadges(ctx, clusters, isLight) {
        // IMPORTANT: Do NOT merge clusters based on 2D screen proximity!
        // Two vectors that project close to each other on the 2D screen (like c [0001] and k̂)
        // are NOT parallel in 3D. Keep each cluster as its own pill badge and de-collide their boxes!
        const placedBoxes = [];
        ctx.save();
        ctx.font = '700 11.5px "Inter", sans-serif';
        const sep = ' ∥ ';
        const sepW = ctx.measureText(sep).width;

        for (const m of clusters) {
            let contentW = 0;
            m.tokens.forEach((tk, idx) => {
                const parsed = this._parseOverbarText(tk.text);
                tk.cleanText = parsed.clean;
                tk.overbarIndices = parsed.overbarIndices;
                tk.w = ctx.measureText(parsed.clean).width;
                contentW += tk.w + (idx > 0 ? sepW : 0);
            });
            const boxW = contentW + 12;
            const boxH = 20;

            const baseLx = (() => {
                let x0 = m.x + Math.cos(m.ang) * 14;
                if (Math.cos(m.ang) < -0.25) x0 -= boxW - 4;
                else if (Math.abs(Math.cos(m.ang)) <= 0.25) x0 -= boxW * 0.5;
                return x0;
            })();
            const baseLy = m.y + Math.sin(m.ang) * 14;

            // Candidate offsets (dx, dy) to avoid overlapping any already-placed pill badge
            const candidateOffsets = [
                [0, 0],
                [0, -23],
                [0, 23],
                [26, -16],
                [-26, -16],
                [26, 16],
                [-26, 16],
                [0, -45],
                [0, 45]
            ];

            let lx = baseLx;
            let ly = baseLy;
            let nudged = false;
            for (const [ox, oy] of candidateOffsets) {
                const tryX = baseLx + ox;
                const tryY = baseLy + oy;
                const collision = placedBoxes.find(
                    (b) =>
                        Math.abs(tryX + boxW / 2 - (b.x + b.w / 2)) < (boxW + b.w) * 0.5 + 4 &&
                        Math.abs(tryY - b.y) < boxH + 3
                );
                if (!collision) {
                    lx = tryX;
                    ly = tryY;
                    nudged = ox !== 0 || oy !== 0;
                    break;
                }
            }
            placedBoxes.push({ x: lx, y: ly, w: boxW, h: boxH });
            ctx.globalAlpha = m.alpha ?? 1;

            // If nudged to avoid another badge, draw a subtle leader line to the arrowhead tip
            if (nudged && m.tipX !== undefined && m.tipY !== undefined) {
                const boxCenterX = lx - 5 + boxW * 0.5;
                const boxCenterY = ly - 3;
                ctx.save();
                ctx.beginPath();
                ctx.moveTo(m.tipX, m.tipY);
                ctx.lineTo(boxCenterX, boxCenterY);
                ctx.strokeStyle = m.tokens[0]?.color || '#94a3b8';
                ctx.globalAlpha = 0.55 * (m.alpha ?? 1);
                ctx.lineWidth = 1.1;
                ctx.setLineDash([2, 2]);
                ctx.stroke();
                ctx.restore();
            }

            // Pill background
            ctx.fillStyle = isLight ? 'rgba(255, 255, 255, 0.92)' : 'rgba(15, 23, 42, 0.90)';
            ctx.strokeStyle = isLight ? 'rgba(148, 163, 184, 0.5)' : 'rgba(71, 85, 105, 0.65)';
            ctx.lineWidth = 1;
            ctx.beginPath();
            ctx.roundRect(lx - 5, ly - 13, boxW, boxH, 5);
            ctx.fill();
            ctx.stroke();

            let curX = lx + 1;
            m.tokens.forEach((tk, idx) => {
                if (idx > 0) {
                    ctx.fillStyle = isLight ? '#64748b' : '#94a3b8';
                    ctx.fillText(sep, curX, ly + 1);
                    curX += sepW;
                }
                ctx.fillStyle = tk.color;
                ctx.fillText(tk.cleanText, curX, ly + 1);
                // Draw crisp, centered overbars directly above negative Miller-Bravais digits
                if (tk.overbarIndices && tk.overbarIndices.length > 0) {
                    for (const charIdx of tk.overbarIndices) {
                        const prefixW = ctx.measureText(tk.cleanText.slice(0, charIdx)).width;
                        const charW = ctx.measureText(tk.cleanText[charIdx]).width;
                        ctx.fillRect(curX + prefixW + 0.4, ly - 9.6, Math.max(2, charW - 0.8), 1.35);
                    }
                }
                curX += tk.w;
            });
        }
        ctx.restore();
    }

    _draw3DHUD(ctx, w, h, state, lockFrame, isLight) {
        const lockNames = {
            Sub: 'Probe steht still (Magnet dreht sich mit φ_M)',
            Lab: 'Magnet & Feld B stehen still (Probe dreht sich)',
            Wave: 'Wellenrichtung k̂ steht still',
            Cryst: '(0001) Easy-Plane liegt fest horizontal'
        };
        const subtitle = lockNames[lockFrame] || lockNames.Sub;
        ctx.save();
        ctx.font = '600 11.5px "Inter", sans-serif';
        ctx.fillStyle = isLight ? '#475569' : '#94a3b8';
        // Bottom-left: the top-left corner belongs to the legend chips
        ctx.fillText(`Ansicht: ${subtitle}`, 14, h - 14);
        ctx.restore();
    }
}
