/**
 * narrative.js
 * Builds the plain-language "Was passiert hier?" explanation for the current setup.
 * Each colored token carries data-emph so hovering it highlights the same object in 3D.
 */

import {
    RAD2DEG,
    CUT_PRESETS,
    crystalToMillerBravaisDirection,
    formatMillerBravaisLatex,
    solveSweepLandmarksForFixedPsi
} from './math_core.js';
import { VECTOR_COLORS } from './visualizer_3d.js';

const TOKEN_COLORS = {
    B: VECTOR_COLORS.B,
    m: VECTOR_COLORS.m,
    l: VECTOR_COLORS.l,
    k: VECTOR_COLORS.k,
    t: VECTOR_COLORS.t,
    n: VECTOR_COLORS.n,
    xs: VECTOR_COLORS.edge,
    ys: VECTOR_COLORS.edge,
    ea: VECTOR_COLORS.ea,
    em: VECTOR_COLORS.em,
    ec: VECTOR_COLORS.ec,
    easy: '#c084fc',
    nodes: VECTOR_COLORS.nodes,
    magnet: '#f87171',
    sample: '#94a3b8'
};

export function token(key, html) {
    const color = TOKEN_COLORS[key] || 'inherit';
    return `<span class="tok" data-emph="${key}" style="--tok:${color}">${html}</span>`;
}

const fmt = (x, d = 1) => {
    const v = Math.abs(x) < 0.5 * Math.pow(10, -d) ? 0 : x;
    return v.toFixed(d);
};

/**
 * Classifies a crystal-frame unit vector relative to the hexagonal axes.
 * Returns a short German phrase (HTML, may contain KaTeX).
 */
export function describeDirection(vCryst, tex) {
    const [x, y, z] = vCryst;
    const tiltDeg = Math.asin(Math.min(1, Math.abs(z))) * RAD2DEG; // angle out of the (0001) plane
    if (tiltDeg > 89) {
        return `zeigt entlang der ${token('ec', 'c-Achse')} – senkrecht zur Easy-Plane`;
    }
    const phi = Math.atan2(y, x) * RAD2DEG;
    // a-type axes sit at 0° + n·60°, m-type axes at 30° + n·60° (measured from e_a)
    const r = ((phi % 60) + 60) % 60;
    const distA = Math.min(r, 60 - r);
    const distM = Math.abs(r - 30);
    let inPlane;
    if (distA < 0.75) {
        inPlane = `entlang einer ${token('ea', 'a-Achse')} (Typ ${tex('\\langle 11\\bar{2}0\\rangle')})`;
    } else if (distM < 0.75) {
        inPlane = `entlang einer ${token('em', 'm-Achse')} (Typ ${tex('\\langle 1\\bar{1}00\\rangle')})`;
    } else {
        inPlane = `${fmt(distA)}° neben einer ${token('ea', 'a-Achse')}`;
    }
    if (tiltDeg < 1) {
        return `liegt in der ${token('easy', 'Easy-Plane')}, ${inPlane}`;
    }
    return `ist ${fmt(tiltDeg)}° aus der ${token('easy', 'Easy-Plane')} gekippt (Projektion ${inPlane})`;
}

// Very short classification for the summary lines ("a-Achse, in der Easy-Plane", "12.3° aus der Easy-Plane").
export function shortDirection(vCryst) {
    const [x, y, z] = vCryst;
    const tiltDeg = Math.asin(Math.min(1, Math.abs(z))) * RAD2DEG;
    const easy = token('easy', 'Easy-Plane');
    if (tiltDeg > 89) return `entlang der ${token('ec', 'c-Achse')}`;
    if (tiltDeg >= 1) return `${fmt(tiltDeg)}° aus der ${easy}`;
    const r = ((((Math.atan2(y, x) * RAD2DEG) % 60) + 60) % 60);
    if (Math.min(r, 60 - r) < 0.75) return `${token('ea', 'a-Achse')}, in der ${easy}`;
    if (Math.abs(r - 30) < 0.75) return `${token('em', 'm-Achse')}, in der ${easy}`;
    return `in der ${easy}`;
}

/**
 * Four one-liners (Probe / Welle / Feld / Spins) – the "at a glance" version of buildNarrative.
 */
export function buildSummary(state, tex) {
    const s = state.setup;
    const m = state.metrics;
    const mb = state.millerBravais;
    const lm = solveSweepLandmarksForFixedPsi(s);
    const T = (latex) => tex(latex);
    const easy = token('easy', 'Easy-Plane');

    const cutInfo = CUT_PRESETS[s.cutPreset];
    const cutName = cutInfo && s.cutPreset !== 'custom' ? cutInfo.shortLabel : T(mb.zsPlane.latex);
    const tiltNC = Math.min(m.angleNCDeg, 180 - m.angleNCDeg);
    let probe;
    if (tiltNC < 1) probe = `${easy} liegt flach in der Probe`;
    else if (Math.abs(tiltNC - 90) < 1) probe = `${easy} steht senkrecht zur Probe`;
    else probe = `${easy} um ${fmt(tiltNC)}° gekippt`;

    const welle = `${token('k', T('\\hat{k}'))} ${T('\\parallel')} ${T(mb.k.latex)} · ${shortDirection(state.vectors.kHat.Cryst)}`;

    const elev = m.elevationSampleDeg;
    const bWhere = Math.abs(elev) < 0.5 ? 'liegt in der Probe' : `${fmt(Math.abs(elev))}° aus der Probe`;
    const bk = lm.bParallelK.errorDeg < 0.5
        ? `parallel zu ${token('k', T('\\hat{k}'))} bei ${T(`\\phi_M = ${fmt(lm.bParallelK.phiMDeg)}^\\circ`)}`
        : `nie parallel zu ${token('k', T('\\hat{k}'))}`;
    const feld = `${token('B', T('\\mathbf{B}'))} ${bWhere} · ${bk}`;

    let spins;
    if (m.isEasyPlaneSingular) {
        spins = `${token('m', T('\\hat{m}'))} unbestimmt (${token('B', T('\\mathbf{B}'))} ${T('\\parallel')} ${token('ec', 'c')})`;
    } else {
        spins = `${token('m', T('\\hat{m}'))} ${T('\\parallel')} ${T(mb.m.latex)}`;
        if (s.deviceType === 'CPW') {
            const t = m.cpwTorqueM;
            const verdict = t > 0.9 ? '<span class="gist-good">stark</span>' : t < 0.15 ? '<span class="gist-bad">schwach</span>' : 'mittel';
            spins += ` · CPW-Kopplung ${verdict}`;
        } else {
            spins += ` · ${token('l', T('\\hat{l}'))} ${T('\\parallel')} ${T(mb.l.latex)}`;
        }
    }

    const line = (label, keys, body) => `<li data-emph="${keys}"><span class="gist-label">${label}</span><span>${body}</span></li>`;
    return line('Probe', 'sample,easy', `${cutName} · ${probe}`)
        + line('Welle', 'device,k', welle)
        + line('Feld', 'B,magnet', feld)
        + line('Spins', 'm,l', spins);
}

export function buildNarrative(state, tex) {
    const s = state.setup;
    const m = state.metrics;
    const v = state.vectors;
    const mb = state.millerBravais;
    const lm = solveSweepLandmarksForFixedPsi(s);
    const T = (latex) => tex(latex);

    // --- 1. Probe -----------------------------------------------------------
    const cutInfo = CUT_PRESETS[s.cutPreset];
    const cutName = cutInfo && s.cutPreset !== 'custom' ? cutInfo.shortLabel : 'eigener Schnitt';
    const angleNC = m.angleNCDeg;
    const tiltNC = Math.min(angleNC, 180 - angleNC);
    let planeRel;
    if (tiltNC < 1) planeRel = 'liegt <strong>flach</strong> in der Oberfläche';
    else if (Math.abs(tiltNC - 90) < 1) planeRel = 'steht <strong>senkrecht</strong> auf der Oberfläche';
    else planeRel = `ist um <strong>${fmt(tiltNC)}°</strong> gegen die Oberfläche gekippt`;

    let probe = `Oberfläche ${T(mb.zsPlane.latex)} (${cutName}). Die ${token('ec', 'c-Achse')} steht ${fmt(tiltNC)}° zur Normale ${token('n', T('\\hat{n}'))}, `
        + `die ${token('easy', 'Easy-Plane (0001)')} ${planeRel}.`;
    if (tiltNC >= 1) {
        const nodesMB = crystalToMillerBravaisDirection(v.nodesHat.Cryst, s.aLattice, s.cLattice);
        probe += ` Beide Ebenen schneiden sich entlang der ${token('nodes', 'Schnittlinie')} ${T(formatMillerBravaisLatex(nodesMB))}`
            + ' – nur dort kann ein Feld gleichzeitig in der Probenebene und in der Easy-Plane liegen.';
    }
    probe += ` Kanten: ${token('xs', '0°-Kante')} ${T('\\parallel')} ${T(mb.xs.latex)}, ${token('ys', '90°-Kante')} ${T('\\parallel')} ${T(mb.ys.latex)}.`;

    // --- 2. Welle -----------------------------------------------------------
    const devName = s.deviceType === 'IDT' ? 'Der IDT (SAW)' : 'Die CPW (Spinwellen)';
    const welle = `${devName} ist um ${fmt(s.thetaDevDeg)}° zur 0°-Kante gedreht. Die Welle ${token('k', T('\\hat{k}'))} ${T('\\parallel')} ${T(mb.k.latex)} `
        + `${describeDirection(v.kHat.Cryst, tex)}. Leiter bzw. Finger laufen quer dazu entlang ${token('t', T('\\hat{t}'))} ${T('\\parallel')} ${T(mb.t.latex)}.`;

    // --- 3. Einbau ----------------------------------------------------------
    let einbau;
    if (s.mountMode === 'IP') {
        einbau = `Die Probe liegt <strong>flach</strong> im ${token('magnet', 'Magneten')}: ${token('n', T('\\hat{n}'))} zeigt nach oben, das Feld dreht immer in der Probenebene. `
            + `Bei ${T('\\phi_M = 0^\\circ')} zeigt ${token('B', T('\\mathbf{B}'))} unter ${fmt(lm.atPhi0.phiBSubDeg)}° zur 0°-Kante `
            + `(${fmt(lm.atPhi0.phiBWaveDeg)}° zu ${token('k', T('\\hat{k}'))}); ${T('\\psi_S')} verschiebt diesen Nullpunkt.`;
    } else {
        einbau = `Die Probe steht <strong>aufrecht</strong> im ${token('magnet', 'Magneten')}: ${token('n', T('\\hat{n}'))} zeigt zu ${T('\\phi_M = 0^\\circ')}. `
            + `Dort steht ${token('B', T('\\mathbf{B}'))} also senkrecht auf der Probe; bei ${T('\\phi_M = 90^\\circ')} liegt es in der Probenebene, `
            + `${fmt(lm.atPhi90.phiBSubDeg)}° zur 0°-Kante (${fmt(lm.atPhi90.phiBWaveDeg)}° zu ${token('k', T('\\hat{k}'))}). `
            + `${T('\\psi_S')} dreht die Probe um ihre Normale und legt fest, welche Probenrichtung der Magnet überstreicht.`;
    }
    const ecLab = v.ecHat.Lab;
    const minEasyPct = Math.abs(ecLab[2]) * 100;
    if (minEasyPct > 99.5) {
        einbau += ` Bei einer vollen Magnetdrehung liegt das Feld <strong>immer vollständig</strong> in der ${token('easy', 'Easy-Plane')}.`;
    } else {
        einbau += ` Bei einer vollen Magnetdrehung liegt zwischen ${fmt(minEasyPct, 0)} % und 100 % des Feldes in der ${token('easy', 'Easy-Plane')} `
            + `(vollständig bei ${T('\\phi_M')} = ${fmt(lm.bInEasyPlane.phiM1)}° und ${fmt(lm.bInEasyPlane.phiM2)}°).`;
    }
    if (lm.bParallelK.errorDeg < 0.5) {
        einbau += ` ${token('B', T('\\mathbf{B}'))} ${T('\\parallel')} ${token('k', T('\\hat{k}'))} erreichst du bei ${T('\\phi_M')} = ${fmt(lm.bParallelK.phiMDeg)}°.`;
    } else {
        einbau += ` ${token('B', T('\\mathbf{B}'))} ${T('\\parallel')} ${token('k', T('\\hat{k}'))} ist in diesem Einbau <strong>nicht</strong> erreichbar `
            + `(bestenfalls ${fmt(lm.bParallelK.errorDeg)}° daneben bei ${T('\\phi_M')} = ${fmt(lm.bParallelK.phiMDeg)}°).`;
    }

    // --- 4. Feld & Spins ----------------------------------------------------
    let spins;
    const elev = m.elevationSampleDeg;
    const bPos = Math.abs(elev) < 0.5
        ? 'liegt in der Probenebene'
        : `steht ${fmt(Math.abs(elev))}° ${elev > 0 ? 'über' : 'unter'} der Probenebene`;
    spins = `Bei ${T(`\\phi_M = ${fmt(s.phiMDeg)}^\\circ`)} ${bPos.split(' ')[0]} ${token('B', T('\\mathbf{B}'))} ${bPos.split(' ').slice(1).join(' ')} `
        + `und hat ${fmt(m.angleBKDeg)}° zu ${token('k', T('\\hat{k}'))}. `;
    if (m.isEasyPlaneSingular) {
        spins += `${token('B', T('\\mathbf{B}'))} ${T('\\parallel')} ${token('ec', 'c')}: keine Komponente in der ${token('easy', 'Easy-Plane')} – im einfachen Modell ist ${token('m', T('\\hat{m}'))} dann unbestimmt.`;
    } else {
        spins += `${fmt(m.bEasyRatio * 100, 0)} % davon (${m.bEasyTesla.toFixed(2)} T) liegen in der ${token('easy', 'Easy-Plane')}. `
            + `${token('m', T('\\hat{m}'))} folgt dieser Projektion (${T('\\parallel')} ${T(mb.m.latex)}), `
            + `${token('l', T('\\hat{l}'))} steht in der Easy-Plane senkrecht dazu (${T('\\parallel')} ${T(mb.l.latex)}).`;
        if (s.deviceType === 'CPW') {
            const torque = m.cpwTorqueM;
            const verdict = torque > 0.9
                ? 'fast maximal'
                : torque < 0.15
                    ? `sehr schwach – ${token('m', T('\\hat{m}'))} liegt fast parallel zu ${token('k', T('\\hat{k}'))}`
                    : 'mittel';
            spins += ` CPW-Kopplung ${T('|\\hat{k}\\times\\hat{m}|')} = ${fmt(torque * 100, 0)} %: ${verdict}.`;
        } else {
            spins += ` ${T('\\angle(\\hat{m},\\hat{k})')} = ${fmt(m.angleMKDeg)}°, ${T('\\angle(\\hat{k},c)')} = ${fmt(m.angleKCDeg)}°.`;
        }
    }

    const section = (num, keys, title, body) => `
        <div class="story-step" data-emph-group="${keys}">
            <div class="story-step-title"><span class="story-num">${num}</span>${title}</div>
            <p>${body}</p>
        </div>`;

    return `
        ${section(1, 'sample,crystal,easy,nodes', 'Probe', probe)}
        ${section(2, 'device,k,t', 'Welle', welle)}
        ${section(3, 'sample,magnet', 'Einbau', einbau)}
        ${section(4, 'B,m,l,magnet', 'Feld &amp; Spins', spins)}
        <p class="story-model-note">
            Modell: einfache Easy-Plane-Näherung oberhalb der Morin-Temperatur (~260 K).
            ${token('m', T('\\hat{m}'))} folgt der Projektion von ${token('B', T('\\mathbf{B}'))} in die (0001)-Ebene;
            Anisotropie innerhalb der Ebene, Entmagnetisierung und Kippung durch die DMI werden vernachlässigt.
        </p>`;
}
