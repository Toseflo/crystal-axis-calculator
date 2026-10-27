/**
 * Preview of the magnet sweep φ_M = 0…360° (other angles fixed):
 * angles to k̂ (one chart), spin-wave geometry strip, B share in the easy plane (second chart).
 * Plain SVG, hover crosshair + tooltip, click picks φ_M.
 */

import { evaluateSetup } from './math_core.js';
import { VECTOR_COLORS } from './visualizer_3d.js';

// Width follows the container (1 SVG unit = 1 px, so text keeps its size)
let W = 640;
const M = { l: 44, r: 34 };
let PW = W - M.l - M.r;
const TOP = { y: 20, h: 190 };   // angle chart
const STRIP = { y: 240, h: 18 }; // geometry strip
const BOT = { y: 284, h: 76 };   // easy-plane share chart
const H = BOT.y + BOT.h + 40;

const SERIES = [
    { key: 'B', name: '∠(B, k̂)', label: 'B', color: VECTOR_COLORS.B, get: (m) => m.angleBKDeg },
    { key: 'm', name: '∠(m̂, k̂)', label: 'm̂', color: VECTOR_COLORS.m, get: (m) => (m.isEasyPlaneSingular ? NaN : m.angleMKDeg) },
    { key: 'l', name: '∠(l̂, k̂)', label: 'l̂', color: VECTOR_COLORS.l, dash: '6 4', get: (m) => (m.isEasyPlaneSingular ? NaN : m.angleLKDeg) }
];
const EASY_COLOR = VECTOR_COLORS.B_easy;

export function computeSweepRows(setup, stepDeg = 1) {
    const rows = [];
    for (let phi = 0; phi <= 360 + 1e-9; phi += stepDeg) {
        const m = evaluateSetup({ ...setup, phiMDeg: phi }).metrics;
        rows.push({ phi, m, vals: SERIES.map((s) => s.get(m)), easy: m.bEasyRatio, geo: m.geometry });
    }
    return rows;
}

const x = (phi) => M.l + (phi / 360) * PW;
const yTop = (deg) => TOP.y + TOP.h - (deg / 180) * TOP.h;
const yBot = (r) => BOT.y + BOT.h - r * BOT.h;
const f1 = (v) => (Number.isFinite(v) ? `${v.toFixed(1)}°` : '–');

// Polyline path that breaks at NaN gaps and at jumps (wrap-arounds)
function pathOf(rows, getY) {
    let d = '';
    let prev = null;
    for (const r of rows) {
        const y = getY(r);
        if (!Number.isFinite(y)) {
            prev = null;
            continue;
        }
        const jump = prev !== null && Math.abs(y - prev) > TOP.h * 0.5;
        d += `${prev === null || jump ? 'M' : 'L'}${x(r.phi).toFixed(1)},${y.toFixed(1)}`;
        prev = y;
    }
    return d;
}

// Contiguous runs of equal geometry
function geometryRuns(rows) {
    const runs = [];
    for (const r of rows) {
        const last = runs[runs.length - 1];
        if (last && last.key === r.geo.key) last.to = r.phi;
        else runs.push({ key: r.geo.key, short: r.geo.short, label: r.geo.label, from: r.phi, to: r.phi });
    }
    return runs;
}

// Push direct labels apart so they never overlap (min. 12 px)
function spreadLabels(items) {
    const sorted = items.filter((i) => Number.isFinite(i.y)).sort((a, b) => a.y - b.y);
    for (let i = 1; i < sorted.length; i++) {
        if (sorted[i].y - sorted[i - 1].y < 12) sorted[i].y = sorted[i - 1].y + 12;
    }
    return sorted;
}

function buildSvg(rows, currentPhi) {
    const grid = [];
    for (const deg of [0, 45, 90, 135, 180]) {
        grid.push(`<line class="sw-grid" x1="${M.l}" x2="${M.l + PW}" y1="${yTop(deg)}" y2="${yTop(deg)}"/>`);
        grid.push(`<text class="sw-tick" x="${M.l - 6}" y="${yTop(deg) + 3.5}" text-anchor="end">${deg}°</text>`);
    }
    for (const r of [0, 0.5, 1]) {
        grid.push(`<line class="sw-grid" x1="${M.l}" x2="${M.l + PW}" y1="${yBot(r)}" y2="${yBot(r)}"/>`);
        grid.push(`<text class="sw-tick" x="${M.l - 6}" y="${yBot(r) + 3.5}" text-anchor="end">${r * 100}%</text>`);
    }
    for (const phi of [0, 90, 180, 270, 360]) {
        grid.push(`<line class="sw-grid sw-grid-v" x1="${x(phi)}" x2="${x(phi)}" y1="${TOP.y}" y2="${BOT.y + BOT.h}"/>`);
        grid.push(`<text class="sw-tick" x="${x(phi)}" y="${BOT.y + BOT.h + 15}" text-anchor="middle">${phi}°</text>`);
    }

    const strip = geometryRuns(rows)
        .map((run) => {
            const x0 = x(Math.max(0, run.from - 0.5)) + 1;
            const x1 = x(Math.min(360, run.to + 0.5)) - 1;
            const w = Math.max(0, x1 - x0);
            // Only the clean geometries get a label; mixed/oblique stretches stay quiet
            const strong = run.key === 'BV' || run.key === 'DE' || run.key === 'FV';
            const text = strong && w > 18 ? run.short : '';
            return `<g><title>${run.label}: φ_M ${run.from}…${run.to}°</title>`
                + `<rect class="sw-geo ${strong ? 'sw-geo-strong' : ''}" x="${x0}" y="${STRIP.y}" width="${w}" height="${STRIP.h}" rx="3"/>`
                + (text ? `<text class="sw-geo-text" x="${(x0 + x1) / 2}" y="${STRIP.y + STRIP.h / 2 + 3.5}" text-anchor="middle">${text}</text>` : '')
                + '</g>';
        })
        .join('');

    const lines = SERIES.map((s, i) =>
        `<path d="${pathOf(rows, (r) => yTop(r.vals[i]))}" fill="none" stroke="${s.color}" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"${s.dash ? ` stroke-dasharray="${s.dash}"` : ''}/>`
    ).join('');
    const easyLine = `<path d="${pathOf(rows, (r) => yBot(r.easy))}" fill="none" stroke="${EASY_COLOR}" stroke-width="2" stroke-linejoin="round"/>`;

    const last = rows[rows.length - 1];
    const ends = spreadLabels(SERIES.map((s, i) => ({ label: s.label, y: yTop(last.vals[i]) })))
        .map((l) => `<text class="sw-end" x="${M.l + PW + 6}" y="${l.y + 4}">${l.label}</text>`)
        .join('');

    const cur = Number.isFinite(currentPhi) ? ((currentPhi % 360) + 360) % 360 : null;
    const marker = cur === null
        ? ''
        : `<line class="sw-current" x1="${x(cur)}" x2="${x(cur)}" y1="${TOP.y - 4}" y2="${BOT.y + BOT.h}"/>`
          + `<text class="sw-current-text" x="${x(cur)}" y="${TOP.y - 6}" text-anchor="${cur > 330 ? 'end' : cur < 30 ? 'start' : 'middle'}">jetzt ${cur.toFixed(1)}°</text>`;

    return `<svg class="sweep-svg" viewBox="0 0 ${W} ${H}" role="img" aria-label="Winkel zu k̂ und Easy-Plane-Anteil über φ_M = 0…360°">
        ${grid.join('')}
        <text class="sw-axis-title" x="${M.l}" y="${STRIP.y - 4}">Geometrie</text>
        ${strip}
        <text class="sw-axis-title" x="${M.l}" y="${BOT.y - 5}">B in Easy-Plane</text>
        ${lines}${easyLine}${ends}${marker}
        <g class="sw-hover" visibility="hidden">
            <line class="sw-cross" y1="${TOP.y}" y2="${BOT.y + BOT.h}"/>
            ${SERIES.map((s) => `<circle class="sw-dot" r="4" fill="${s.color}"/>`).join('')}
            <circle class="sw-dot sw-dot-easy" r="4" fill="${EASY_COLOR}"/>
        </g>
        <rect class="sw-hit" x="${M.l}" y="${TOP.y}" width="${PW}" height="${BOT.y + BOT.h - TOP.y}" fill="transparent"/>
        <text class="sw-axis-title" x="${M.l + PW / 2}" y="${H - 4}" text-anchor="middle">Magnetwinkel φ_M</text>
    </svg>`;
}

/**
 * Renders the sweep into `container` and wires hover + click.
 * onPick(phiDeg) is called when the user clicks a φ_M in the chart.
 */
export function renderSweepPlot(container, setup, onPick) {
    const rows = computeSweepRows(setup);
    W = Math.max(420, Math.round(container.clientWidth || 640));
    PW = W - M.l - M.r;
    container.innerHTML = `<div class="sweep-legend">${SERIES.map(
        (s) => `<span><svg width="18" height="8" aria-hidden="true"><line x1="1" x2="17" y1="4" y2="4" stroke="${s.color}" stroke-width="2"${s.dash ? ' stroke-dasharray="4 3"' : ''}/></svg>${s.name}</span>`
    ).join('')}</div>
        <div class="sweep-plot-wrap">${buildSvg(rows, setup.phiMDeg)}<div class="sweep-tip hidden"></div></div>`;

    const svg = container.querySelector('.sweep-svg');
    const hover = svg.querySelector('.sw-hover');
    const cross = hover.querySelector('.sw-cross');
    const dots = hover.querySelectorAll('.sw-dot');
    const tip = container.querySelector('.sweep-tip');
    const hit = svg.querySelector('.sw-hit');

    const rowAt = (evt) => {
        const pt = svg.createSVGPoint();
        pt.x = evt.clientX;
        pt.y = evt.clientY;
        const p = pt.matrixTransform(svg.getScreenCTM().inverse());
        const phi = Math.round(Math.max(0, Math.min(360, ((p.x - M.l) / PW) * 360)));
        return rows[phi];
    };

    const show = (r, evt) => {
        hover.setAttribute('visibility', 'visible');
        cross.setAttribute('x1', x(r.phi));
        cross.setAttribute('x2', x(r.phi));
        SERIES.forEach((s, i) => {
            const ok = Number.isFinite(r.vals[i]);
            dots[i].setAttribute('visibility', ok ? 'visible' : 'hidden');
            dots[i].setAttribute('cx', x(r.phi));
            dots[i].setAttribute('cy', ok ? yTop(r.vals[i]) : 0);
        });
        dots[SERIES.length].setAttribute('cx', x(r.phi));
        dots[SERIES.length].setAttribute('cy', yBot(r.easy));

        tip.replaceChildren();
        const head = document.createElement('div');
        head.className = 'sweep-tip-head';
        head.textContent = `φ_M = ${r.phi}°`;
        tip.appendChild(head);
        const addRow = (color, dash, value, name) => {
            const row = document.createElement('div');
            row.className = 'sweep-tip-row';
            row.innerHTML = color
                ? `<svg width="12" height="6" aria-hidden="true"><line x1="0" x2="12" y1="3" y2="3" stroke="${color}" stroke-width="2"${dash ? ' stroke-dasharray="3 2"' : ''}/></svg>`
                : '<span></span>';
            const v = document.createElement('strong');
            v.textContent = value;
            const n = document.createElement('span');
            n.textContent = name;
            row.append(v, n);
            tip.appendChild(row);
        };
        SERIES.forEach((s, i) => addRow(s.color, s.dash, f1(r.vals[i]), s.name));
        addRow(EASY_COLOR, null, `${Math.round(r.easy * 100)} %`, 'B in Easy-Plane');
        addRow(null, null, r.geo.short, r.geo.label === r.geo.short ? 'Geometrie' : r.geo.label);

        const wrap = tip.parentElement.getBoundingClientRect();
        tip.classList.remove('hidden');
        const left = evt.clientX - wrap.left;
        tip.style.left = `${left + (left > wrap.width * 0.6 ? -tip.offsetWidth - 14 : 14)}px`;
        tip.style.top = `${Math.max(0, evt.clientY - wrap.top - tip.offsetHeight / 2)}px`;
    };
    const hide = () => {
        hover.setAttribute('visibility', 'hidden');
        tip.classList.add('hidden');
    };

    hit.addEventListener('pointermove', (e) => show(rowAt(e), e));
    hit.addEventListener('pointerleave', hide);
    hit.addEventListener('click', (e) => onPick && onPick(rowAt(e).phi));
    return rows;
}
