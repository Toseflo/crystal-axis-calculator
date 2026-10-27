/**
 * Small everyday exports for the "Teilen & Exportieren" card:
 * a plain-text lab note, a φ_M sweep as CSV and file-download helpers.
 * (The thesis SVG lives in svg_exporter.js, the Python snippet in math_core.js.)
 */

import { CUT_PRESETS, evaluateSetup } from './math_core.js';

const MOUNT_DE = { IP: 'horizontal (IP)', OOP: 'vertikal (OOP)' };

const f1 = (x) => Number(x).toFixed(1);
const pct = (x) => `${Math.round(x * 100)} %`;

// File name stem shared by all downloads, e.g. hematite_r-cut_OOP_CPW_phiM90
export function exportBaseName(setup) {
    return `hematite_${setup.cutPreset}_${setup.mountMode}_${setup.deviceType}_phiM${Math.round(setup.phiMDeg)}`;
}

// Rows (label, text) describing one evaluated setup
function setupRows(state, tag) {
    const s = state.setup;
    const mb = state.millerBravais;
    const m = state.metrics;
    const cut = CUT_PRESETS[s.cutPreset];
    const cutName = cut && s.cutPreset !== 'custom' ? cut.label : `eigener Schnitt ${mb.zsPlane.label}`;
    const spins = m.isEasyPlaneSingular ? 'm̂ unbestimmt (B ∥ c)' : `m̂ ∥ ${mb.m.label}, l̂ ∥ ${mb.l.label}`;

    const angles = [`B in Easy-Plane ${pct(m.bEasyRatio)}`];
    if (!m.isEasyPlaneSingular) {
        angles.unshift(`∠(m̂,k̂) = ${f1(m.angleMKDeg)}°`, `∠(l̂,k̂) = ${f1(m.angleLKDeg)}°`);
        if (s.deviceType === 'CPW') angles.push(`CPW-Kopplung |k̂×m̂| ${pct(m.cpwTorqueM)}`);
    }

    const rows = [
        ['Probe', `${cutName}, 0°-Kante ∥ ${mb.xs.label}${s.alphaFlatDeg ? `, Flat α = ${f1(s.alphaFlatDeg)}°` : ''}`],
        ['Halter', `${MOUNT_DE[s.mountMode] || s.mountMode}, ψ_S = ${f1(s.psiSDeg)}°`],
        ['Struktur', `${s.deviceType}, θ_dev = ${f1(s.thetaDevDeg)}°  →  k̂ ∥ ${mb.k.label}`],
        ['Magnet', `φ_M = ${f1(s.phiMDeg)}°, B₀ = ${Number(s.b0Tesla).toFixed(2)} T  →  B ∥ ${mb.B.label}`],
        ['Spins', spins],
        ['Winkel', angles.join(', ')]
    ];
    // Preset names start with the cut – only worth repeating when comparing A and B
    if (tag) rows.unshift([`Setup ${tag}`, s.name || '–']);
    else if (s.name && !(cut && s.name.startsWith(cut.shortLabel))) rows.unshift(['Setup', s.name]);
    return rows;
}

/**
 * Plain-text summary – for the lab notebook, a mail or a chat.
 * states: one evaluated setup, or A and B when comparing.
 * extras: { pinned: ['[10-11]', …], link: 'https://…' }
 */
export function buildLabNote(states, extras = {}) {
    const blocks = states.map((st, i) => setupRows(st, states.length > 1 ? 'AB'[i] : null));
    const tail = [];
    if (extras.pinned && extras.pinned.length) tail.push(['Gemerkt', extras.pinned.join(', ')]);
    if (extras.link) tail.push(['Link', extras.link]);

    const width = Math.max(...[...blocks.flat(), ...tail].map(([k]) => k.length)) + 2;
    const fmt = (rows) => rows.map(([k, v]) => k.padEnd(width) + v).join('\n');
    return [...blocks, tail]
        .filter((b) => b.length)
        .map(fmt)
        .join('\n\n');
}

/**
 * Magnet sweep φ_M = 0…360° (other angles fixed) as CSV.
 * Vector components are given in the wave frame (k̂, t̂, n̂).
 */
export function buildSweepCsv(setup, stepDeg = 1) {
    const header = [
        'phi_M_deg',
        'B_easy_ratio',
        'B_in_sample_ratio',
        'angle_B_k_deg',
        'angle_m_k_deg',
        'angle_l_k_deg',
        'cpw_coupling_kxm',
        'm_k',
        'm_t',
        'm_n',
        'l_k',
        'l_t',
        'l_n'
    ];
    const num = (x, d = 5) => (Number.isFinite(x) ? Number(x.toFixed(d)) : '');
    const lines = [header.join(',')];
    for (let phi = 0; phi <= 360 + 1e-9; phi += stepDeg) {
        const st = evaluateSetup({ ...setup, phiMDeg: phi });
        const m = st.metrics;
        const singular = m.isEasyPlaneSingular;
        const spin = (vec) => (singular ? ['', '', ''] : vec.Wave.map((x) => num(x)));
        lines.push(
            [
                num(phi, 3),
                num(m.bEasyRatio),
                num(m.bInSampleRatio),
                num(m.angleBKDeg, 3),
                singular ? '' : num(m.angleMKDeg, 3),
                singular ? '' : num(m.angleLKDeg, 3),
                singular ? '' : num(m.cpwTorqueM),
                ...spin(st.vectors.mHat),
                ...spin(st.vectors.lHat)
            ].join(',')
        );
    }
    return lines.join('\n') + '\n';
}

export function downloadBlob(blob, filename) {
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    a.remove();
    URL.revokeObjectURL(url);
}

export function downloadText(text, filename, mime = 'text/plain') {
    downloadBlob(new Blob([text], { type: `${mime};charset=utf-8` }), filename);
}
