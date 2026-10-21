/**
 * visualizer_3d.js
 * Clean, didactic 3D & 2D visualization engine for the Hematite
 * Spin-Wave (CPW) & SAW (IDT) Coordinate Calculator.
 *
 * Designed to avoid visual overload:
 *   - Progressive layer toggles (essential vectors B & k_hat first, optional
 *     easy-plane moments m_hat & l_hat, optional crystal axes & magnet coils).
 *   - High-clarity 2D Compass Rose views for the Sample Surface (F_Sub / F_Wave)
 *     and the Hematite Basal Easy-Plane (0001) (F_Cryst).
 *   - Frame-Lock support (Sub, Wave, Cryst, Lab) in 3D.
 */

import {
    DEG2RAD,
    RAD2DEG,
    add3,
    scale3,
    dot3,
    norm3,
    normalize3,
    matVec3,
    FRAME_METADATA
} from './math_core.js';

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
    em: '#818cf8'        // Indigo: Crystal m-axis [1-100]
};

export class Scene3DController {
    constructor(canvas3d, canvas2dSample, canvas2dEasy) {
        this.canvas3d = canvas3d;
        this.canvas2dSample = canvas2dSample;
        this.canvas2dEasy = canvas2dEasy;

        this.lockFrame = 'Sub'; // 'Sub' | 'Wave' | 'Cryst' | 'Lab'
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
            showAngleArcs: true,
            showCompareOverlay: true
        };

        this.sample2dLockWave = false; // false -> x_s right, true -> k_hat right
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
        if (!this.stateA) return;
        if (this.canvas3d) this._render3DScene();
        if (this.canvas2dSample) this._render2DSample();
        if (this.canvas2dEasy) this._render2DEasyPlane();
    }

    // ========================================================================
    // 3D Main Scene Renderer
    // ========================================================================
    _render3DScene() {
        const prep = this._prepareCanvas(this.canvas3d);
        if (!prep) return;
        const { ctx, w, h } = prep;
        const isLight = document.body.classList.contains('theme-light');
        const cx = w * 0.5;
        const cy = h * 0.53;
        const baseScale = Math.min(w, h) * 0.33 * this.camera.zoom;

        ctx.clearRect(0, 0, w, h);
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

        const F = this.lockFrame;
        const stA = this.stateA;
        const phiMRad = (stA.setup.phiMDeg || 0) * DEG2RAD;
        const cosNegPhi = Math.cos(-phiMRad);
        const sinNegPhi = Math.sin(-phiMRad);

        // When F === 'Lab' ("Labor-Magnet"), lock the rotating Electromagnet (poles N/S & B)
        // along the fixed +X axis so the Magnet stands still and the Sample rotates when phi_M changes!
        const toLock = (vec, srcFrame) => {
            if (F === 'Lab') {
                const vLab = srcFrame === 'Lab'
                    ? vec
                    : matVec3(stA.matrices.getRotationMatrix('Lab', srcFrame), vec);
                return [
                    cosNegPhi * vLab[0] - sinNegPhi * vLab[1],
                    sinNegPhi * vLab[0] + cosNegPhi * vLab[1],
                    vLab[2]
                ];
            }
            return matVec3(stA.matrices.getRotationMatrix(F, srcFrame), vec);
        };
        const vecInLock = (vecObj) => toLock(vecObj.Lab, 'Lab');
        const proj = (pt3) => this._project3D(pt3, cx, cy, baseScale);

        // 1. Reference grid in the locked frame
        this._drawLockedFrameGrid(ctx, proj, F, isLight);

        // 2. Depth-sorted 3D surfaces
        const renderQueue = [];
        if (this.visibility.showLabMagnet) {
            this._enqueueLabMagnet(renderQueue, toLock, proj, stA, isLight);
        }
        if (this.visibility.showEasyPlane) {
            this._enqueueEasyPlaneDisk(renderQueue, toLock, proj, isLight);
        }
        if (this.visibility.showSampleSlab) {
            this._enqueueSubstrateAndDevice(renderQueue, toLock, proj, stA, isLight);
        }

        renderQueue.sort((a, b) => b.depth - a.depth);
        for (const item of renderQueue) {
            item.draw(ctx);
        }

        // 3. Collect ONLY the essential 3D arrows (transformed via vecInLock so 'Lab' locks the Magnet!)
        const arrows3D = [];

        // Optional Crystal Basis Axes with Miller-Bravais notation (when toggled on)
        if (this.visibility.showCrystalAxes) {
            const ea = vecInLock(stA.vectors.eaHat);
            const em = vecInLock(stA.vectors.emHat);
            const ec = vecInLock(stA.vectors.ecHat);
            arrows3D.push({ vec: ea, baseLen: 0.92, color: VECTOR_COLORS.ea, label: 'a [112\u03050]', width: 2.0, dashed: true });
            arrows3D.push({ vec: em, baseLen: 0.88, color: VECTOR_COLORS.em, label: 'm [11\u030500]', width: 1.8, dashed: true });
            arrows3D.push({ vec: ec, baseLen: 1.00, color: VECTOR_COLORS.ec, label: 'c [0001]', width: 2.3, dashed: false });
        }

        const bHat = vecInLock(stA.vectors.bHat);
        const bEasy = scale3(vecInLock(stA.vectors.bEasyProjVec), 1 / (stA.metrics.b0Tesla || 1));
        const bEasyLen = norm3(bEasy);

        // Dashed orthogonal projection line from the exact tip of B (1.00 * bHat)
        // to the exact tip of the Easy-Plane projection / m̂ arrow (bEasyLen * mHat)
        const diffOutPlane = Math.hypot(bHat[0] - bEasy[0], bHat[1] - bEasy[1], bHat[2] - bEasy[2]);
        if (this.visibility.showB && (this.visibility.showM || this.visibility.showEasyPlane) && bEasyLen > 0.04 && diffOutPlane > 0.05) {
            this._drawProjectionDropLine3D(ctx, proj, bHat, bEasy, VECTOR_COLORS.B_easy, 1.3);
        }

        if (this.visibility.showK) {
            arrows3D.push({ vec: vecInLock(stA.vectors.kHat), baseLen: 1.00, color: VECTOR_COLORS.k, label: 'k̂', width: 3.0, dashed: false });
        }
        if (this.visibility.showL && !stA.metrics.isEasyPlaneSingular) {
            arrows3D.push({ vec: vecInLock(stA.vectors.lHat), baseLen: 1.00, color: VECTOR_COLORS.l, label: 'l̂', width: 2.8, dashed: false });
        }
        if (this.visibility.showM && !stA.metrics.isEasyPlaneSingular && bEasyLen > 0.04) {
            // Scale m̂ arrow with the Easy-Plane projection length (bEasyLen = sin ∠(B, c))
            // so its tip always sits right where the dashed projection line from B lands!
            arrows3D.push({ vec: vecInLock(stA.vectors.mHat), baseLen: bEasyLen, color: VECTOR_COLORS.m, label: 'm̂', width: 2.9, dashed: false });
        }
        if (this.visibility.showB) {
            arrows3D.push({ vec: bHat, baseLen: 1.00, color: VECTOR_COLORS.B, label: 'B', width: 3.5, dashed: false });
        }

        // Compare Overlay (Setup B)
        if (this.compareActive && this.stateB && this.visibility.showCompareOverlay) {
            const stB = this.stateB;
            if (this.visibility.showK) {
                arrows3D.push({ vec: toLock(stB.vectors.kHat.Lab, 'Lab'), baseLen: 0.92, color: '#f0abfc', label: 'k̂(B)', width: 2.0, dashed: true });
            }
            if (this.visibility.showB) {
                arrows3D.push({ vec: toLock(stB.vectors.bHat.Lab, 'Lab'), baseLen: 0.96, color: '#7dd3fc', label: 'B(B)', width: 2.2, dashed: true });
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
            draw: (ctx) => {
                ctx.save();
                ctx.beginPath();
                topPts.forEach((p, idx) => (idx === 0 ? ctx.moveTo(p.x, p.y) : ctx.lineTo(p.x, p.y)));
                ctx.closePath();
                ctx.fillStyle = isLight ? 'rgba(226, 232, 240, 0.78)' : 'rgba(30, 41, 59, 0.78)';
                ctx.strokeStyle = isLight ? '#475569' : '#94a3b8';
                ctx.lineWidth = 1.8;
                ctx.fill();
                ctx.stroke();

                if (this.visibility.showMicrostructure) {
                    this._drawMicrostructure3D(ctx, toLock, proj, state.setup.deviceType);
                }

                if (this.visibility.showEasyPlane && Math.abs(state.setup.thetaCutDeg) > 1e-3) {
                    const nodeSub = state.vectors.nodesHat.Sub;
                    const pN1 = proj(toLock(scale3(nodeSub, -0.58), 'Sub'));
                    const pN2 = proj(toLock(scale3(nodeSub, 0.58), 'Sub'));
                    ctx.beginPath();
                    ctx.moveTo(pN1.x, pN1.y);
                    ctx.lineTo(pN2.x, pN2.y);
                    ctx.strokeStyle = isLight ? '#9333ea' : '#c084fc';
                    ctx.setLineDash([5, 3]);
                    ctx.lineWidth = 1.6;
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
                this._drawSingleShaftAndHead(ctx, origin.x, origin.y, p2.x, p2.y, it.color, it.width, it.dashed);
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

            // If nudged to avoid another badge, draw a subtle leader line to the arrowhead tip
            if (nudged && m.tipX !== undefined && m.tipY !== undefined) {
                const boxCenterX = lx - 5 + boxW * 0.5;
                const boxCenterY = ly - 3;
                ctx.save();
                ctx.beginPath();
                ctx.moveTo(m.tipX, m.tipY);
                ctx.lineTo(boxCenterX, boxCenterY);
                ctx.strokeStyle = m.tokens[0]?.color || '#94a3b8';
                ctx.globalAlpha = 0.55;
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
            Sub: 'Sample steht still (Magnet dreht sich mit φ_M)',
            Lab: 'Labor-Magnet & Feld B stehen still (Sample dreht sich)',
            Wave: 'Wellenrichtung k̂ steht still',
            Cryst: '(0001) Easy-Plane liegt fest horizontal'
        };
        const subtitle = lockNames[lockFrame] || lockNames.Sub;
        ctx.save();
        ctx.font = '600 11.5px "Inter", sans-serif';
        ctx.fillStyle = isLight ? '#475569' : '#94a3b8';
        ctx.fillText(`Ansicht: ${subtitle}`, 14, 21);
        ctx.restore();
    }

    // ========================================================================
    // 2D Sample Surface Compass View (Chip Plane)
    // ========================================================================
    _render2DSample() {
        const prep = this._prepareCanvas(this.canvas2dSample);
        if (!prep) return;
        const { ctx, w, h } = prep;
        const isLight = document.body.classList.contains('theme-light');
        const cx = w * 0.5;
        const cy = h * 0.56;
        const rMax = Math.min(w, h) * 0.30;

        ctx.clearRect(0, 0, w, h);
        ctx.fillStyle = isLight ? '#f8fafc' : '#0b1120';
        ctx.fillRect(0, 0, w, h);

        const st = this.stateA;
        const frame2D = this.sample2dLockWave ? 'Wave' : 'Sub';

        // 1. Chip outline rectangle
        ctx.save();
        ctx.translate(cx, cy);
        const chipRotRad = this.sample2dLockWave ? -st.setup.thetaDevDeg * DEG2RAD : 0;
        ctx.rotate(chipRotRad);
        const cw = rMax * 1.44;
        const ch = rMax * 1.22;
        ctx.fillStyle = isLight ? 'rgba(226, 232, 240, 0.9)' : 'rgba(30, 41, 59, 0.78)';
        ctx.strokeStyle = isLight ? '#64748b' : '#475569';
        ctx.lineWidth = 1.8;
        ctx.fillRect(-cw / 2, -ch / 2, cw, ch);
        ctx.strokeRect(-cw / 2, -ch / 2, cw, ch);
        ctx.restore();

        // 2. Compass Rose Circle & Degree Ticks with explicit axis meaning
        const lbl0 = this.sample2dLockWave ? '0° (∥ k̂)' : '0° (Sample)';
        const lbl90 = this.sample2dLockWave ? '90° (⊥ k̂)' : '90° (Sample)';
        this._drawCompassRing(ctx, cx, cy, rMax, isLight, lbl0, lbl90);

        // 3. Center Microstructure Schematic (CPW or IDT) rotated along k_hat
        const k2D = st.vectors.kHat[frame2D];
        const kAngleRad = Math.atan2(-k2D[1], k2D[0]);
        ctx.save();
        ctx.translate(cx, cy);
        ctx.rotate(kAngleRad);
        const polyScale = rMax * 0.85;
        ctx.fillStyle = 'rgba(251, 191, 36, 0.85)';
        const polys = this._getMicrostructurePolygonsWave(st.setup.deviceType);
        for (const poly of polys) {
            ctx.beginPath();
            poly.forEach(([xk, yt], idx) => {
                const px = xk * polyScale;
                const py = -yt * polyScale;
                if (idx === 0) ctx.moveTo(px, py);
                else ctx.lineTo(px, py);
            });
            ctx.closePath();
            ctx.fill();
        }
        ctx.restore();

        // 4. Line of nodes (Intersection of (0001) Easy-Plane with sample surface)
        const nodes2D = st.vectors.nodesHat[frame2D];
        if (this.visibility.showEasyPlane && Math.abs(st.setup.thetaCutDeg) > 1e-3) {
            ctx.beginPath();
            ctx.moveTo(cx - nodes2D[0] * rMax * 1.0, cy + nodes2D[1] * rMax * 1.0);
            ctx.lineTo(cx + nodes2D[0] * rMax * 1.0, cy - nodes2D[1] * rMax * 1.0);
            ctx.strokeStyle = isLight ? '#9333ea' : '#c084fc';
            ctx.setLineDash([5, 4]);
            ctx.lineWidth = 1.6;
            ctx.stroke();
            ctx.setLineDash([]);
        }

        // 5. Vectors projected onto the Sample Chip Plane
        const bEasyLen = st.metrics.bEasyRatio;
        const vecItems = [];
        if (this.visibility.showCrystalAxes) {
            vecItems.push({ v: st.vectors.eaHat[frame2D], color: VECTOR_COLORS.ea, label: 'a [112\u03050]', width: 2.0, dashed: true });
            vecItems.push({ v: st.vectors.emHat[frame2D], color: VECTOR_COLORS.em, label: 'm [11\u030500]', width: 1.8, dashed: true });
            vecItems.push({ v: st.vectors.ecHat[frame2D], color: VECTOR_COLORS.ec, label: 'c [0001]', width: 2.2, dashed: false });
        }
        if (this.visibility.showK) {
            vecItems.push({ v: st.vectors.kHat[frame2D], color: VECTOR_COLORS.k, label: 'k̂', width: 3.0, dashed: false });
        }
        if (!st.metrics.isEasyPlaneSingular) {
            if (this.visibility.showL) {
                vecItems.push({ v: st.vectors.lHat[frame2D], color: VECTOR_COLORS.l, label: 'l̂', width: 2.6, dashed: false });
            }
            if (this.visibility.showM && bEasyLen > 0.04) {
                vecItems.push({ v: scale3(st.vectors.mHat[frame2D], bEasyLen), color: VECTOR_COLORS.m, label: 'm̂', width: 2.6, dashed: false });
            }
        }
        if (this.visibility.showB) {
            vecItems.push({ v: st.vectors.bHat[frame2D], color: VECTOR_COLORS.B, label: 'B', width: 3.4, dashed: false });
        }

        this._draw2DVectorsDeCollided(ctx, cx, cy, rMax, vecItems, isLight);

        ctx.font = '700 12.5px "Inter", sans-serif';
        ctx.fillStyle = isLight ? '#0f172a' : '#f8fafc';
        ctx.fillText('1. Draufsicht auf das Sample', 14, 21);
        ctx.font = '500 11px "Inter", sans-serif';
        ctx.fillStyle = isLight ? '#64748b' : '#94a3b8';
        ctx.fillText(`Blick senkrecht auf die Sample-Oberfläche (${st.setup.cutPreset})`, 14, 37);

        ctx.fillStyle = isLight ? '#475569' : '#94a3b8';
        ctx.fillText(
            `Feld im Sample: ${(st.metrics.bInSampleRatio * 100).toFixed(0)}%  ·  Winkel zur Welle: ${st.metrics.phiBWaveDeg.toFixed(1)}°`,
            14,
            h - 12
        );
    }

    // ========================================================================
    // 2D Basal Easy-Plane (0001) Compass View
    // ========================================================================
    _render2DEasyPlane() {
        const prep = this._prepareCanvas(this.canvas2dEasy);
        if (!prep) return;
        const { ctx, w, h } = prep;
        const isLight = document.body.classList.contains('theme-light');
        const cx = w * 0.5;
        const cy = h * 0.56;
        const rMax = Math.min(w, h) * 0.30;

        ctx.clearRect(0, 0, w, h);
        ctx.fillStyle = isLight ? '#f8fafc' : '#0b1120';
        ctx.fillRect(0, 0, w, h);

        const st = this.stateA;

        // 1. Hexagonal Basal Plane (0001) & Compass Ring with Crystal Axis Labels
        this._drawCompassRing(ctx, cx, cy, rMax, isLight, '0° (a)', '90° (m)');

        ctx.save();
        ctx.beginPath();
        for (let i = 0; i < 6; i++) {
            const ang = i * 60 * DEG2RAD;
            const x = cx + rMax * Math.cos(ang);
            const y = cy - rMax * Math.sin(ang);
            if (i === 0) ctx.moveTo(x, y);
            else ctx.lineTo(x, y);
        }
        ctx.closePath();
        ctx.fillStyle = isLight ? 'rgba(168, 85, 247, 0.11)' : 'rgba(168, 85, 247, 0.15)';
        ctx.strokeStyle = isLight ? 'rgba(147, 51, 234, 0.6)' : 'rgba(192, 132, 252, 0.6)';
        ctx.lineWidth = 1.8;
        ctx.fill();
        ctx.stroke();

        // Three C2 crystal symmetry axes at 0, 60, 120 deg
        ctx.strokeStyle = isLight ? 'rgba(100, 116, 139, 0.25)' : 'rgba(148, 163, 184, 0.2)';
        ctx.setLineDash([3, 3]);
        for (const deg of [0, 60, 120]) {
            const rad = deg * DEG2RAD;
            ctx.beginPath();
            ctx.moveTo(cx - rMax * Math.cos(rad), cy + rMax * Math.sin(rad));
            ctx.lineTo(cx + rMax * Math.cos(rad), cy - rMax * Math.sin(rad));
            ctx.stroke();
        }
        ctx.setLineDash([]);
        ctx.restore();

        const bEasyLen = st.metrics.bEasyRatio;
        const vecItems = [];
        if (this.visibility.showCrystalAxes) {
            vecItems.push({ v: st.vectors.eaHat.Cryst, color: VECTOR_COLORS.ea, label: 'a [112\u03050]', width: 2.0, dashed: true });
            vecItems.push({ v: st.vectors.emHat.Cryst, color: VECTOR_COLORS.em, label: 'm [11\u030500]', width: 1.8, dashed: true });
        }
        if (this.visibility.showK) {
            vecItems.push({ v: st.vectors.kHat.Cryst, color: VECTOR_COLORS.k, label: 'k̂', width: 2.6, dashed: true });
        }
        if (!st.metrics.isEasyPlaneSingular) {
            if (this.visibility.showL) {
                vecItems.push({ v: st.vectors.lHat.Cryst, color: VECTOR_COLORS.l, label: 'l̂', width: 2.9, dashed: false });
            }
            if (this.visibility.showM && bEasyLen > 0.04) {
                vecItems.push({ v: scale3(st.vectors.mHat.Cryst, bEasyLen), color: VECTOR_COLORS.m, label: 'm̂', width: 2.9, dashed: false });
            }
        }
        if (this.visibility.showB) {
            vecItems.push({ v: st.vectors.bHat.Cryst, color: VECTOR_COLORS.B, label: 'B', width: 3.3, dashed: false });
        }

        this._draw2DVectorsDeCollided(ctx, cx, cy, rMax, vecItems, isLight);

        ctx.font = '700 12.5px "Inter", sans-serif';
        ctx.fillStyle = isLight ? '#0f172a' : '#f8fafc';
        ctx.fillText('2. Blick in die (0001) Kristall-Easy-Plane', 14, 21);

        const isSamePlaneAsSample = Math.abs(st.setup.thetaCutDeg) < 0.5;
        ctx.font = '600 11px "Inter", sans-serif';
        ctx.fillStyle = isSamePlaneAsSample
            ? (isLight ? '#7e22ce' : '#d8b4fe')
            : (isLight ? '#64748b' : '#94a3b8');
        ctx.fillText(
            isSamePlaneAsSample
                ? 'Beim c-cut identisch mit der Sample-Ebene (0° Verkippung)'
                : `Um ${st.setup.thetaCutDeg.toFixed(1)}° gegen die Sample-Ebene verkippt (${st.setup.cutPreset})`,
            14,
            37
        );

        ctx.font = '500 11px "Inter", sans-serif';
        ctx.fillStyle = isLight ? '#475569' : '#94a3b8';
        ctx.fillText(
            `Feld in Easy-Plane: ${(st.metrics.bEasyRatio * 100).toFixed(0)}%  ·  Richtung m̂: ${st.metrics.phiEasyMDeg.toFixed(1)}°`,
            14,
            h - 12
        );
    }

    _drawCompassRing(ctx, cx, cy, rMax, isLight, label0 = '0°', label90 = '90°') {
        ctx.save();
        // Subtle crosshair
        ctx.strokeStyle = isLight ? 'rgba(100, 116, 139, 0.22)' : 'rgba(148, 163, 184, 0.16)';
        ctx.lineWidth = 1;
        ctx.setLineDash([3, 3]);
        ctx.beginPath();
        ctx.moveTo(cx - rMax, cy);
        ctx.lineTo(cx + rMax, cy);
        ctx.moveTo(cx, cy - rMax);
        ctx.lineTo(cx, cy + rMax);
        ctx.stroke();
        ctx.setLineDash([]);

        ctx.beginPath();
        ctx.arc(cx, cy, rMax, 0, Math.PI * 2);
        ctx.strokeStyle = isLight ? 'rgba(100, 116, 139, 0.34)' : 'rgba(148, 163, 184, 0.26)';
        ctx.lineWidth = 1.2;
        ctx.stroke();

        ctx.font = '500 9.5px "JetBrains Mono", monospace';
        ctx.fillStyle = isLight ? '#64748b' : '#64748b';
        ctx.textAlign = 'center';
        for (let deg = 0; deg < 360; deg += 30) {
            const rad = deg * DEG2RAD;
            const x1 = cx + rMax * Math.cos(rad);
            const y1 = cy - rMax * Math.sin(rad);
            const x2 = cx + (rMax + 5) * Math.cos(rad);
            const y2 = cy - (rMax + 5) * Math.sin(rad);
            ctx.beginPath();
            ctx.moveTo(x1, y1);
            ctx.lineTo(x2, y2);
            ctx.stroke();
            if (deg === 0) {
                ctx.textAlign = 'left';
                ctx.fillText(label0, cx + rMax + 7, cy + 3);
                ctx.textAlign = 'center';
            } else if (deg === 90) {
                ctx.fillText(label90, cx, cy - rMax - 8);
            } else if (deg === 180) {
                ctx.textAlign = 'right';
                ctx.fillText('180°', cx - rMax - 7, cy + 3);
                ctx.textAlign = 'center';
            } else if (deg === 270) {
                ctx.fillText('270°', cx, cy + rMax + 15);
            }
        }
        ctx.restore();
    }

    _draw2DVectorsDeCollided(ctx, cx, cy, rMax, items, isLight) {
        // Group ONLY vectors that are truly parallel in 3D (dot3 > 0.999) into a '∥' badge!
        const groups = [];
        for (const it of items) {
            const len2 = Math.hypot(it.v[0], it.v[1]);
            if (len2 < 0.04) continue;
            const ux = it.v[0] / len2;
            const uy = it.v[1] / len2;
            const u3 = normalize3(it.v);
            let found = null;
            for (const g of groups) {
                if (dot3(u3, g.u3Ref) > 0.999) {
                    found = g;
                    break;
                }
            }
            if (!found) {
                found = { ux, uy, u3Ref: u3, maxLen: len2, items: [] };
                groups.push(found);
            }
            found.maxLen = Math.max(found.maxLen, len2);
            found.items.push({ ...it, len2 });
        }

        const labelClusters = [];
        for (const g of groups) {
            g.items.forEach((it, idx) => {
                const stagger = g.items.length > 1 ? 0.96 - idx * 0.16 : 0.94;
                const r = Math.min(1, it.len2) * stagger * rMax;
                const x2 = cx + g.ux * r;
                const y2 = cy - g.uy * r;
                this._drawSingleShaftAndHead(ctx, cx, cy, x2, y2, it.color, it.width || 2.6, it.dashed);
            });
            const rTip = Math.min(1, g.maxLen) * 0.94 * rMax;
            const rOuter = Math.min(1, g.maxLen) * 0.98 * rMax;
            const xOut = cx + g.ux * rOuter;
            const yOut = cy - g.uy * rOuter;
            const ang = Math.atan2(-g.uy, g.ux);
            labelClusters.push({
                tipX: cx + g.ux * rTip,
                tipY: cy - g.uy * rTip,
                x: xOut,
                y: yOut,
                ang,
                tokens: g.items.map((it) => ({ text: it.label, color: it.color }))
            });
        }

        this._drawGroupedPillBadges(ctx, labelClusters, isLight);
    }
}
