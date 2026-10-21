/**
 * sweep_chart.js
 * Interactive 0°..360° Magnet Rotation Sweep (phi_M) chart comparing Setup A (solid)
 * and Setup B (dashed). Supports click/drag scrubbing to set phi_M directly.
 */

import { computeMagnetSweep, solveSweepLandmarksForFixedPsi, wrapDeg360 } from './math_core.js';
import { VECTOR_COLORS } from './visualizer_3d.js';

export class SweepChartController {
    constructor(canvas, onPhiMChange) {
        this.canvas = canvas;
        this.onPhiMChange = onPhiMChange;
        this.mode = 'fractions'; // 'fractions' | 'angles' | 'projections'
        this.sweepA = null;
        this.sweepB = null;
        this.landmarksA = null;
        this.currentPhiM = 0;
        this.compareActive = false;
        this.isScrubbing = false;

        this._bindEvents();
    }

    _bindEvents() {
        if (!this.canvas) return;

        const updateFromEvent = (e) => {
            const rect = this.canvas.getBoundingClientRect();
            const padLeft = 64;
            const padRight = 24;
            const plotW = Math.max(10, rect.width - padLeft - padRight);
            const relX = Math.max(0, Math.min(1, (e.clientX - rect.left - padLeft) / plotW));
            const newPhi = wrapDeg360(Math.round(relX * 360 * 2) / 2);
            if (this.onPhiMChange) {
                this.onPhiMChange(newPhi);
            }
        };

        this.canvas.addEventListener('pointerdown', (e) => {
            this.isScrubbing = true;
            this.canvas.setPointerCapture(e.pointerId);
            updateFromEvent(e);
        });

        this.canvas.addEventListener('pointermove', (e) => {
            if (!this.isScrubbing) return;
            updateFromEvent(e);
        });

        const stopScrub = (e) => {
            if (this.isScrubbing) {
                this.isScrubbing = false;
                try {
                    this.canvas.releasePointerCapture(e.pointerId);
                } catch (_) {
                    // ignore
                }
            }
        };
        this.canvas.addEventListener('pointerup', stopScrub);
        this.canvas.addEventListener('pointercancel', stopScrub);
    }

    update(setupA, setupB = null, compareActive = false) {
        this.sweepA = computeMagnetSweep(setupA, 2);
        this.landmarksA = solveSweepLandmarksForFixedPsi(setupA);
        this.sweepB = compareActive && setupB ? computeMagnetSweep(setupB, 2) : null;
        this.currentPhiM = Number(setupA.phiMDeg ?? 0);
        this.compareActive = compareActive;
        this.render();
    }

    render() {
        if (!this.canvas || !this.sweepA) return;
        const dpr = window.devicePixelRatio || 1;
        const rect = this.canvas.getBoundingClientRect();
        if (rect.width < 10 || rect.height < 10) return;
        const w = Math.max(240, Math.floor(rect.width));
        const h = Math.max(180, Math.floor(rect.height));

        if (this.canvas.width !== Math.floor(w * dpr) || this.canvas.height !== Math.floor(h * dpr)) {
            this.canvas.width = Math.floor(w * dpr);
            this.canvas.height = Math.floor(h * dpr);
        }

        const ctx = this.canvas.getContext('2d');
        ctx.setTransform(dpr, 0, 0, dpr, 0, 0);

        const isLight = document.body.classList.contains('theme-light');
        ctx.clearRect(0, 0, w, h);
        ctx.fillStyle = isLight ? '#f8fafc' : '#0b1120';
        ctx.fillRect(0, 0, w, h);

        const padLeft = 64;
        const padRight = 24;
        const padTop = 36;
        const padBottom = 46;
        const plotW = w - padLeft - padRight;
        const plotH = h - padTop - padBottom;

        // Determine Y axis range and series based on mode
        let yMin = 0;
        let yMax = 100;
        let yTicks = [0, 25, 50, 75, 100];
        let yUnit = '%';
        let yAxisTitle = 'Feld-Anteil / Anregung (%)';
        let yScale = 1;
        let seriesDefs = [];

        if (this.mode === 'angles') {
            yMin = 0;
            yMax = 180;
            yTicks = [0, 45, 90, 135, 180];
            yUnit = '°';
            yScale = 1;
            yAxisTitle = 'Winkel zur Welle k̂ (°)';
            seriesDefs = [
                { key: 'angleBKDeg', label: '∠(B, k̂) Magnetfeld', color: VECTOR_COLORS.B },
                { key: 'angleMKDeg', label: '∠(m̂, k̂) Magnetisierung', color: VECTOR_COLORS.m },
                { key: 'angleLKDeg', label: '∠(l̂, k̂) Néel-Vektor', color: VECTOR_COLORS.l }
            ];
        } else if (this.mode === 'fractions') {
            yMin = 0;
            yMax = 100;
            yTicks = [0, 25, 50, 75, 100];
            yUnit = '%';
            yScale = 100;
            yAxisTitle = 'Feld-Anteil / Anregung (%)';
            seriesDefs = [
                { key: 'bEasyRatio', label: 'Feld in Easy-Plane (%)', color: VECTOR_COLORS.B_easy },
                { key: 'cpwTorqueM', label: 'CPW-Anregung |k̂ × m̂| (%)', color: VECTOR_COLORS.k },
                { key: 'bInSampleRatio', label: 'Feld flach im Sample (%)', color: VECTOR_COLORS.B }
            ];
        } else {
            yMin = -1.0;
            yMax = 1.0;
            yTicks = [-1.0, -0.5, 0.0, 0.5, 1.0];
            yUnit = '';
            yScale = 1;
            yAxisTitle = 'Projektion auf Welle k̂ (cos θ)';
            seriesDefs = [
                { key: 'bDotK', label: 'b̂ · k̂ (Feld längs k̂)', color: VECTOR_COLORS.B },
                { key: 'mDotK', label: 'm̂ · k̂ (Moment längs k̂)', color: VECTOR_COLORS.m },
                { key: 'lDotK', label: 'l̂ · k̂ (Néel längs k̂)', color: VECTOR_COLORS.l }
            ];
        }

        const xToPx = (phi) => padLeft + (phi / 360) * plotW;
        const yToPx = (val) => padTop + (1 - (val * yScale - yMin) / (yMax - yMin)) * plotH;

        // Grid lines & Y ticks
        ctx.font = '500 10.5px "JetBrains Mono", monospace';
        ctx.textAlign = 'right';
        for (const yt of yTicks) {
            const py = padTop + (1 - (yt - yMin) / (yMax - yMin)) * plotH;
            ctx.beginPath();
            ctx.moveTo(padLeft, py);
            ctx.lineTo(padLeft + plotW, py);
            ctx.strokeStyle = isLight ? 'rgba(100, 116, 139, 0.22)' : 'rgba(148, 163, 184, 0.14)';
            ctx.lineWidth = 1;
            ctx.stroke();

            ctx.fillStyle = isLight ? '#475569' : '#94a3b8';
            ctx.fillText(`${yt}${yUnit}`, padLeft - 8, py + 3.5);
        }

        // In 'angles' mode, annotate physical meaning of 0°, 90°, 180° on the right inside the plot
        if (this.mode === 'angles') {
            ctx.save();
            ctx.font = '600 10px "Inter", sans-serif';
            ctx.textAlign = 'right';
            ctx.fillStyle = isLight ? 'rgba(71, 85, 105, 0.75)' : 'rgba(148, 163, 184, 0.65)';
            const xHint = padLeft + plotW - 6;
            ctx.fillText('0° = längs zur Welle (∥ +k̂)', xHint, padTop + plotH - 5);
            ctx.fillText('90° = quer zur Welle (⊥ k̂)', xHint, padTop + plotH * 0.5 - 5);
            ctx.fillText('180° = entgegen der Welle (∥ −k̂)', xHint, padTop + 12);
            ctx.restore();
        }

        // X ticks
        ctx.textAlign = 'center';
        ctx.font = '500 10.5px "JetBrains Mono", monospace';
        for (const xt of [0, 45, 90, 135, 180, 225, 270, 315, 360]) {
            const px = xToPx(xt);
            ctx.beginPath();
            ctx.moveTo(px, padTop);
            ctx.lineTo(px, padTop + plotH);
            ctx.strokeStyle = isLight ? 'rgba(100, 116, 139, 0.18)' : 'rgba(148, 163, 184, 0.11)';
            ctx.stroke();

            ctx.fillStyle = isLight ? '#475569' : '#94a3b8';
            ctx.fillText(`${xt}°`, px, padTop + plotH + 16);
        }

        // X-Axis Title
        ctx.save();
        ctx.font = '600 11.5px "Inter", sans-serif';
        ctx.fillStyle = isLight ? '#334155' : '#cbd5e1';
        ctx.textAlign = 'center';
        ctx.fillText('Magnetwinkel φ_M im Labor (0° … 360°)', padLeft + plotW * 0.5, h - 8);
        ctx.restore();

        // Y-Axis Title (rotated on left)
        ctx.save();
        ctx.translate(15, padTop + plotH * 0.5);
        ctx.rotate(-Math.PI / 2);
        ctx.font = '600 11px "Inter", sans-serif';
        ctx.fillStyle = isLight ? '#334155' : '#cbd5e1';
        ctx.textAlign = 'center';
        ctx.fillText(yAxisTitle, 0, 0);
        ctx.restore();

        // Vertical Landmark Markers for B || k and B perp k
        if (this.landmarksA) {
            const marks = [
                { phi: this.landmarksA.bParallelK.phiMDeg, text: 'B ∥ k̂', color: '#38bdf8' },
                { phi: this.landmarksA.bPerpK.phiMDeg, text: 'B ⊥ k̂', color: '#94a3b8' }
            ];
            ctx.save();
            ctx.font = '600 9.5px "Inter", sans-serif';
            for (const mk of marks) {
                const mx = xToPx(wrapDeg360(mk.phi));
                ctx.beginPath();
                ctx.moveTo(mx, padTop);
                ctx.lineTo(mx, padTop + plotH);
                ctx.strokeStyle = isLight ? 'rgba(56, 189, 248, 0.35)' : 'rgba(56, 189, 248, 0.25)';
                ctx.setLineDash([2, 4]);
                ctx.lineWidth = 1.2;
                ctx.stroke();
                ctx.setLineDash([]);

                ctx.fillStyle = mk.color;
                ctx.textAlign = 'center';
                ctx.fillText(`${mk.text} (${mk.phi.toFixed(0)}°)`, mx, padTop - 4);
            }
            ctx.restore();
        }

        // Plot curves helper
        const drawCurve = (sweepData, key, color, dashed = false) => {
            const xs = sweepData.phiMArray;
            const ys = sweepData[key];
            ctx.save();
            ctx.beginPath();
            if (dashed) ctx.setLineDash([5, 4]);
            for (let i = 0; i < xs.length; i++) {
                const px = xToPx(xs[i]);
                const py = yToPx(ys[i]);
                if (i === 0) ctx.moveTo(px, py);
                else ctx.lineTo(px, py);
            }
            ctx.strokeStyle = color;
            ctx.lineWidth = dashed ? 1.8 : 2.5;
            ctx.stroke();
            ctx.restore();
        };

        // Draw Setup B first (dashed) if active, then Setup A (solid)
        if (this.compareActive && this.sweepB) {
            for (const s of seriesDefs) {
                drawCurve(this.sweepB, s.key, s.color, true);
            }
        }
        for (const s of seriesDefs) {
            drawCurve(this.sweepA, s.key, s.color, false);
        }

        // Current phi_M vertical cursor line + intersection dots
        const wrappedPhi = wrapDeg360(this.currentPhiM);
        const curPx = xToPx(wrappedPhi);
        ctx.save();
        ctx.beginPath();
        ctx.moveTo(curPx, padTop);
        ctx.lineTo(curPx, padTop + plotH);
        ctx.strokeStyle = isLight ? '#0f172a' : '#f8fafc';
        ctx.lineWidth = 1.7;
        ctx.setLineDash([3, 3]);
        ctx.stroke();
        ctx.setLineDash([]);

        // Find nearest index in sweepA for intersection dots
        const nearIdx = Math.min(
            this.sweepA.phiMArray.length - 1,
            Math.max(0, Math.round(wrappedPhi / 2))
        );
        for (const s of seriesDefs) {
            const val = this.sweepA[s.key][nearIdx];
            const py = yToPx(val);
            ctx.beginPath();
            ctx.arc(curPx, py, 4.2, 0, Math.PI * 2);
            ctx.fillStyle = s.color;
            ctx.fill();
            ctx.lineWidth = 1.5;
            ctx.strokeStyle = isLight ? '#ffffff' : '#0b1120';
            ctx.stroke();
        }

        // Cursor angle badge at top of line
        const badgeText = `φ_M = ${wrappedPhi.toFixed(1)}°`;
        ctx.font = '700 10px "JetBrains Mono", monospace';
        const bw = ctx.measureText(badgeText).width + 10;
        const bx = Math.max(padLeft, Math.min(padLeft + plotW - bw, curPx - bw * 0.5));
        ctx.fillStyle = isLight ? '#0f172a' : '#38bdf8';
        ctx.beginPath();
        ctx.roundRect(bx, padTop + 4, bw, 16, 4);
        ctx.fill();
        ctx.fillStyle = isLight ? '#ffffff' : '#090d16';
        ctx.textAlign = 'center';
        ctx.fillText(badgeText, bx + bw * 0.5, padTop + 15);
        ctx.restore();

        // Top Legend bar
        let legX = padLeft;
        ctx.font = '600 11px "Inter", sans-serif';
        ctx.textAlign = 'left';
        for (const s of seriesDefs) {
            ctx.fillStyle = s.color;
            ctx.fillRect(legX, 10, 13, 3.5);
            ctx.fillText(s.label, legX + 18, 15);
            legX += ctx.measureText(s.label).width + 34;
        }
        if (this.compareActive && this.sweepB) {
            ctx.fillStyle = isLight ? '#334155' : '#cbd5e1';
            ctx.fillText('(Durchgezogen: A | Gestrichelt: B)', legX, 15);
        }
    }
}
