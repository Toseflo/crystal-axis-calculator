/**
 * app.js
 * Main UI controller for the Hematite (alpha-Fe2O3) Spin-Wave (CPW) & SAW (IDT)
 * orientation tool.
 *
 * Layout: setup steps (left) | 3D stage (center) | live overview (right).
 * Everything that names a physical object carries a data-emph key; hovering it
 * highlights the same object in the 3D view.
 */

import {
    HEMATITE_CONSTANTS,
    CUT_PRESETS,
    FRAME_METADATA,
    TARGET_PRESETS,
    RAD2DEG,
    DEG2RAD,
    createDefaultSetup,
    evaluateSetup,
    solveTargetAngles,
    solveSweepLandmarksForFixedPsi,
    computeMountingPresetPsiS,
    millerBravaisDirectionToCrystal,
    crystalToMillerBravaisDirection,
    millerBravaisPlaneToCrystalNormal,
    crystalNormalToMillerBravaisPlane,
    formatMillerBravais,
    formatMillerBravaisLatex,
    scaleMillerBravaisInts,
    cross3,
    dot3,
    normalize3,
    matVec3,
    generatePythonSnippet,
    wrapDeg360
} from './math_core.js';

import { Scene3DController, VECTOR_COLORS } from './visualizer_3d.js';
import { generateThesisSVG } from './svg_exporter.js';
import { buildNarrative, buildSummary, shortDirection, token } from './narrative.js';
import { GuidedTour } from './tour.js';

// Application State (Default: c-cut (0001) in IP orientation)
const appState = {
    activeTab: 'A', // 'A' | 'B'
    compareActive: false,
    setupA: createDefaultSetup({
        name: 'c-cut (0001) · Flach (IP)',
        mountMode: 'IP',
        cutPreset: 'c-cut',
        thetaCutDeg: 0,
        phiCutDeg: 0,
        alphaFlatDeg: 0,
        deviceType: 'CPW',
        thetaDevDeg: 0,
        psiSDeg: 0,
        phiMDeg: 0,
        b0Tesla: 0.5
    }),
    setupB: createDefaultSetup({
        name: 'r-cut (1-102) · Vertikal (OOP)',
        mountMode: 'OOP',
        cutPreset: 'r-cut',
        thetaCutDeg: HEMATITE_CONSTANTS.thetaRCutDeg,
        phiCutDeg: 90,
        alphaFlatDeg: 0,
        deviceType: 'CPW',
        thetaDevDeg: 0,
        psiSDeg: 0,
        phiMDeg: 90,
        b0Tesla: 0.5
    }),
    targetId: 'B_parallel_k',
    customTargetUVTW: [1, 1, -2, 0],
    vectorTableMode: 'angles', // 'angles' | 'xyz'
    vectorTableAll: false, // false: only B, m, l, k (+ custom)
    customDir: { mode: 'dir' } // 'dir' | 'plane' | 'cross' | 'angles'
};

const PRESET_STORAGE_KEY = 'hematite_user_presets_v1';

const DEFAULT_STARTER_PRESETS = [
    {
        id: 'preset_ccut_ip',
        name: 'Standard: c-cut (0001) · Flach (IP)',
        setup: {
            mountMode: 'IP',
            cutPreset: 'c-cut',
            thetaCutDeg: 0,
            phiCutDeg: 0,
            alphaFlatDeg: 0,
            deviceType: 'CPW',
            thetaDevDeg: 0,
            psiSDeg: 0,
            phiMDeg: 0,
            b0Tesla: 0.5
        }
    },
    {
        id: 'preset_ccut_oop',
        name: 'c-cut (0001) · Vertikal (OOP)',
        setup: {
            mountMode: 'OOP',
            cutPreset: 'c-cut',
            thetaCutDeg: 0,
            phiCutDeg: 0,
            alphaFlatDeg: 0,
            deviceType: 'CPW',
            thetaDevDeg: 0,
            psiSDeg: 0,
            phiMDeg: 90,
            b0Tesla: 0.5
        }
    },
    {
        id: 'preset_acut_ip',
        name: 'a-cut (11-20) · Flach (IP)',
        setup: {
            mountMode: 'IP',
            cutPreset: 'a-cut',
            thetaCutDeg: 90,
            phiCutDeg: 0,
            alphaFlatDeg: 0,
            deviceType: 'IDT',
            thetaDevDeg: 0,
            psiSDeg: 0,
            phiMDeg: 0,
            b0Tesla: 0.5
        }
    },
    {
        id: 'preset_rcut_oop',
        name: 'r-cut (1-102) · Vertikal (OOP)',
        setup: {
            mountMode: 'OOP',
            cutPreset: 'r-cut',
            thetaCutDeg: HEMATITE_CONSTANTS.thetaRCutDeg,
            phiCutDeg: 90,
            alphaFlatDeg: 0,
            deviceType: 'CPW',
            thetaDevDeg: 0,
            psiSDeg: 0,
            phiMDeg: 90,
            b0Tesla: 0.5
        }
    }
];

// Which 3D objects a sidebar step / slider is "about" (used for hover emphasis)
const STEP_EMPHASIS = {
    stepCard1: ['sample', 'crystal', 'easy', 'nodes'],
    stepCardDev: ['device', 'k', 't'],
    stepCard2: ['sample', 'magnet', 'n'],
    stepCard3: ['B', 'm', 'l', 'magnet']
};
const SLIDER_EMPHASIS = {
    ThetaCut: ['sample', 'crystal', 'easy', 'nodes'],
    PhiCut: ['sample', 'crystal', 'easy', 'nodes'],
    AlphaFlat: ['xs', 'ys', 'crystal'],
    ThetaDev: ['device', 'k', 't'],
    PsiS: ['sample', 'device', 'k', 'magnet'],
    PhiM: ['B', 'm', 'l', 'magnet']
};

// Compact URL keys for sharing a setup as a link
const URL_KEYS = {
    cut: 'cutPreset',
    tc: 'thetaCutDeg',
    pc: 'phiCutDeg',
    af: 'alphaFlatDeg',
    dev: 'deviceType',
    td: 'thetaDevDeg',
    mnt: 'mountMode',
    psi: 'psiSDeg',
    phi: 'phiMDeg',
    b0: 'b0Tesla'
};
const URL_STRING_PROPS = new Set(['cutPreset', 'deviceType', 'mountMode']);

let savedPresets = [];
let activePresetId = 'preset_ccut_ip';

function loadPresetsFromStorage() {
    try {
        const raw = localStorage.getItem(PRESET_STORAGE_KEY);
        if (raw) {
            const parsed = JSON.parse(raw);
            if (Array.isArray(parsed)) {
                savedPresets = parsed;
                return;
            }
        }
    } catch (_) {
        // Fallback to starter presets if localStorage is unavailable
    }
    savedPresets = DEFAULT_STARTER_PRESETS.map((p) => ({ ...p, setup: { ...p.setup } }));
}

function savePresetsToStorage() {
    try {
        localStorage.setItem(PRESET_STORAGE_KEY, JSON.stringify(savedPresets));
    } catch (_) {
        // Ignore storage quota errors
    }
}

function renderPresetDropdown() {
    const sel = document.getElementById('selectScenarioPreset');
    const delBtn = document.getElementById('btnDeletePreset');
    if (!sel) return;

    if (savedPresets.length === 0) {
        sel.innerHTML = `<option value="">— Keine Vorlagen gespeichert —</option>`;
        sel.disabled = true;
        if (delBtn) delBtn.disabled = true;
        return;
    }

    sel.disabled = false;
    if (delBtn) delBtn.disabled = false;

    if (!savedPresets.some((p) => p.id === activePresetId)) {
        activePresetId = savedPresets[0].id;
    }

    sel.innerHTML = savedPresets
        .map((p) => `<option value="${p.id}" ${p.id === activePresetId ? 'selected' : ''}>${p.name}</option>`)
        .join('');
}

function getActiveSetup() {
    return appState.activeTab === 'A' ? appState.setupA : appState.setupB;
}

// Controllers
let sceneController = null;
let sceneController2 = null;
let tour = null;
let lastActiveState = null;

function tex(latexStr) {
    if (window.katex) {
        try {
            return window.katex.renderToString(latexStr, { throwOnError: false });
        } catch (_) {
            return latexStr;
        }
    }
    return latexStr;
}

function renderStaticLatex() {
    document.querySelectorAll('[data-tex]').forEach((el) => {
        el.innerHTML = tex(el.dataset.tex);
    });
}

function showToast(msg) {
    const el = document.getElementById('toast');
    if (!el) return;
    el.textContent = msg;
    el.classList.remove('hidden');
    clearTimeout(showToast._t);
    showToast._t = setTimeout(() => el.classList.add('hidden'), 1800);
}

document.addEventListener('DOMContentLoaded', () => {
    const canvas3d = document.getElementById('canvas3D');
    const canvas3dSec = document.getElementById('canvas3DSecondary');

    sceneController = new Scene3DController(canvas3d);
    sceneController2 = new Scene3DController(canvas3dSec);
    // Share visibility flags between View 1 and View 2; default View 2 to Easy-Plane in Top-Down (2D) mode
    sceneController2.visibility = sceneController.visibility;
    sceneController2.lockFrame = 'Cryst';
    sceneController2.setTopDown2D(true);

    readStateFromUrl();
    renderStaticLatex();
    loadPresetsFromStorage();
    renderPresetDropdown();
    populateTargetSelect();
    bindUIEvents();
    syncControlsFromSetup();
    syncLegend();
    syncFramePills();
    recomputeAndRender();

    window.addEventListener('load', () => {
        renderStaticLatex();
        syncControlsFromSetup();
        recomputeAndRender();
    });

    window.addEventListener('resize', () => {
        sceneController.render();
        if (sceneController2) sceneController2.render();
    });
});

function populateTargetSelect() {
    const sel = document.getElementById('selectTargetPreset');
    if (!sel) return;
    sel.innerHTML = TARGET_PRESETS.map(
        (t) => `<option value="${t.id}">${t.label}</option>`
    ).join('');
    sel.value = appState.targetId;
}

function syncControlsFromSetup() {
    const s = getActiveSetup();

    // Tabs & Compare Header
    document.getElementById('tabSetupA').classList.toggle('active', appState.activeTab === 'A');
    document.getElementById('tabSetupB').classList.toggle('active', appState.activeTab === 'B');
    document.getElementById('lblTabA').textContent = appState.setupA.name || 'Setup A';
    document.getElementById('lblTabB').textContent = appState.setupB.name || 'Setup B';
    document.getElementById('btnToggleCompare').classList.toggle('active', appState.compareActive);
    document.getElementById('compareSetupBox').classList.toggle('hidden', !appState.compareActive);
    document.getElementById('inputSetupName').value = s.name;

    // Mount Mode
    document.getElementById('btnMountOOP').classList.toggle('active', s.mountMode === 'OOP');
    document.getElementById('btnMountIP').classList.toggle('active', s.mountMode === 'IP');

    // Sliders & Numeric Inputs
    setPair('PhiM', s.phiMDeg);
    setPair('PsiS', s.psiSDeg);
    document.getElementById('numB0').value = s.b0Tesla;

    // Cut Presets
    document.querySelectorAll('#cutPresetGrid .cut-btn').forEach((btn) => {
        btn.classList.toggle('active', btn.dataset.cut === s.cutPreset);
    });

    const cutInfo = CUT_PRESETS[s.cutPreset] || CUT_PRESETS.custom;
    document.getElementById('cutDescriptionBox').textContent = cutInfo.description;
    document.getElementById('customCutPanel').classList.toggle('hidden', s.cutPreset !== 'custom');

    setPair('ThetaCut', Number(s.thetaCutDeg).toFixed(2));
    setPair('PhiCut', Number(s.phiCutDeg).toFixed(1));
    setPair('AlphaFlat', s.alphaFlatDeg);

    // Device Type & ThetaDev
    document.getElementById('btnDevCPW').classList.toggle('active', s.deviceType === 'CPW');
    document.getElementById('btnDevIDT').classList.toggle('active', s.deviceType === 'IDT');
    setPair('ThetaDev', s.thetaDevDeg);
}

function setPair(suffix, val) {
    const sl = document.getElementById(`slider${suffix}`);
    const nm = document.getElementById(`num${suffix}`);
    if (sl) sl.value = val;
    if (nm) nm.value = val;
}

function bindSliderPair(suffix, propName, wrap360 = false) {
    const sl = document.getElementById(`slider${suffix}`);
    const nm = document.getElementById(`num${suffix}`);
    const updateVal = (raw) => {
        let v = parseFloat(raw);
        if (Number.isNaN(v)) v = 0;
        if (wrap360) v = wrapDeg360(v);
        const s = getActiveSetup();
        s[propName] = v;
        if (sl) sl.value = v;
        if (nm) nm.value = v;
        recomputeAndRender();
    };
    if (sl) sl.addEventListener('input', (e) => updateVal(e.target.value));
    if (nm) nm.addEventListener('input', (e) => updateVal(e.target.value));

    // While a slider is dragged, highlight what it moves
    if (sl && SLIDER_EMPHASIS[suffix]) {
        sl.addEventListener('pointerdown', () => {
            emphasisLock = true;
            setEmphasis(SLIDER_EMPHASIS[suffix]);
            const release = () => {
                emphasisLock = false;
                setEmphasis(null);
                window.removeEventListener('pointerup', release);
                window.removeEventListener('pointercancel', release);
            };
            window.addEventListener('pointerup', release);
            window.addEventListener('pointercancel', release);
        });
    }
}

// ============================================================================
// Emphasis (hover anywhere → highlight in 3D)
// ============================================================================

let emphasisLock = false;
let currentEmphasisSig = '';

function setEmphasis(keys) {
    const sig = keys ? keys.join(',') : '';
    if (sig === currentEmphasisSig) return;
    currentEmphasisSig = sig;
    sceneController.setEmphasis(keys);
    if (sceneController2) sceneController2.setEmphasis(keys);
}

function emphasisKeysFor(target) {
    const el = target.closest('[data-emph], .story-step-title, .step-card-header');
    if (!el) return null;
    if (el.dataset.emph) return el.dataset.emph.split(',');
    if (el.classList.contains('story-step-title')) {
        const group = el.closest('[data-emph-group]');
        return group ? group.dataset.emphGroup.split(',') : null;
    }
    const card = el.closest('.step-card');
    return card ? STEP_EMPHASIS[card.id] || null : null;
}

function bindEmphasisHover() {
    document.addEventListener('pointerover', (e) => {
        if (emphasisLock || !(e.target instanceof Element)) return;
        setEmphasis(emphasisKeysFor(e.target));
    });
    document.addEventListener('pointerleave', () => {
        if (!emphasisLock) setEmphasis(null);
    });
}

// ============================================================================
// Legend chips (= layer toggles) & frame pills
// ============================================================================

function syncLegend() {
    const vis = sceneController.visibility;
    document.querySelectorAll('#canvasLegend .legend-chip[data-vis]').forEach((chip) => {
        chip.classList.toggle('off', !vis[chip.dataset.vis]);
    });
    // The "Eigene Richtung" layer only exists while its card is open
    document.getElementById('legendChipCustom').hidden = !document.getElementById('customDirCard').open;
}

function syncFramePills() {
    document.querySelectorAll('#frameLockPills .frame-pill').forEach((p) => {
        p.classList.toggle('active', p.dataset.frame === sceneController.lockFrame);
    });
    if (lastActiveState && document.getElementById('selectVectorTableFrame').value === 'view') {
        updateVectorTable(lastActiveState);
    }
}

function rerenderScenes() {
    sceneController.render();
    if (sceneController2) sceneController2.render();
}

function bindUIEvents() {
    // Stepper Cards Accordion Toggle (Setup Cards 1, Dev, 2 auto-collapse each other for a calm sidebar)
    const setupCardIds = ['stepCard1', 'stepCardDev', 'stepCard2'];
    document.querySelectorAll('.step-card-header[data-toggle-step]').forEach((hdr) => {
        hdr.addEventListener('click', () => {
            const card = hdr.closest('.step-card');
            if (!card) return;
            const willOpen = !card.classList.contains('step-card-active');
            card.classList.toggle('step-card-active', willOpen);
            if (willOpen && setupCardIds.includes(card.id)) {
                setupCardIds.forEach((otherId) => {
                    if (otherId !== card.id) {
                        const otherCard = document.getElementById(otherId);
                        if (otherCard) otherCard.classList.remove('step-card-active');
                    }
                });
            }
        });
    });

    // Optional Dual 3D View Toggle & Top-Down (2D) Camera Toggles inside the 3D viewports
    const stageContainer = document.getElementById('stageContainer');
    const btnToggleDual3D = document.getElementById('btnToggleDual3D');
    const btnCloseDual3D = document.getElementById('btnCloseDual3D');
    const lblFrameLock1 = document.getElementById('lblFrameLock1');
    const btnCamTop1 = document.getElementById('btnCamTop1');
    const btnCamTop2 = document.getElementById('btnCamTop2');

    const setDual3DActive = (enable) => {
        if (!stageContainer) return;
        stageContainer.classList.toggle('mode-3d', !enable);
        stageContainer.classList.toggle('mode-dual-3d', enable);
        if (btnToggleDual3D) {
            btnToggleDual3D.classList.toggle('active', enable);
            const lbl = document.getElementById('lblToggleDual3D');
            if (lbl) lbl.textContent = enable ? '2. Ansicht aktiv' : '+ 2. Ansicht';
        }
        if (lblFrameLock1) {
            lblFrameLock1.textContent = enable ? 'Ansicht 1 – still steht:' : 'Was steht im Bild still?';
        }
        requestAnimationFrame(rerenderScenes);
    };

    if (btnToggleDual3D) {
        btnToggleDual3D.addEventListener('click', () => {
            setDual3DActive(!stageContainer.classList.contains('mode-dual-3d'));
        });
    }
    if (btnCloseDual3D) {
        btnCloseDual3D.addEventListener('click', () => setDual3DActive(false));
    }

    if (btnCamTop1) {
        sceneController.onCameraChange = () => {
            btnCamTop1.classList.toggle('active', sceneController.isTopDown2D());
        };
        btnCamTop1.addEventListener('click', () => {
            sceneController.setTopDown2D(!sceneController.isTopDown2D());
        });
    }
    if (btnCamTop2 && sceneController2) {
        sceneController2.onCameraChange = () => {
            btnCamTop2.classList.toggle('active', sceneController2.isTopDown2D());
        };
        btnCamTop2.addEventListener('click', () => {
            sceneController2.setTopDown2D(!sceneController2.isTopDown2D());
        });
    }

    // Saveable & Deletable Preset Manager
    const scenarioSelect = document.getElementById('selectScenarioPreset');
    const saveRow = document.getElementById('savePresetRow');
    const inputPresetName = document.getElementById('inputPresetName');
    const btnOpenSave = document.getElementById('btnOpenSavePreset');
    const btnConfirmSave = document.getElementById('btnConfirmSavePreset');
    const btnCancelSave = document.getElementById('btnCancelSavePreset');
    const btnDeletePreset = document.getElementById('btnDeletePreset');

    if (scenarioSelect) {
        scenarioSelect.addEventListener('change', (e) => {
            if (e.target.value) applySavedPreset(e.target.value);
        });
    }

    if (btnOpenSave && saveRow && inputPresetName) {
        btnOpenSave.addEventListener('click', () => {
            const s = getActiveSetup();
            inputPresetName.value = `${s.cutPreset} · ${s.mountMode} · ${s.deviceType} (${s.thetaDevDeg.toFixed(0)}°)`;
            saveRow.classList.remove('hidden');
            inputPresetName.focus();
            inputPresetName.select();
        });
    }

    const confirmSavePreset = () => {
        if (!inputPresetName) return;
        const name = inputPresetName.value.trim();
        if (!name) return;
        saveCurrentAsPreset(name);
        if (saveRow) saveRow.classList.add('hidden');
    };

    if (btnConfirmSave) btnConfirmSave.addEventListener('click', confirmSavePreset);
    if (btnCancelSave && saveRow) {
        btnCancelSave.addEventListener('click', () => saveRow.classList.add('hidden'));
    }
    if (inputPresetName && saveRow) {
        inputPresetName.addEventListener('keydown', (e) => {
            if (e.key === 'Enter') {
                e.preventDefault();
                confirmSavePreset();
            } else if (e.key === 'Escape') {
                saveRow.classList.add('hidden');
            }
        });
    }
    if (btnDeletePreset) btnDeletePreset.addEventListener('click', deleteSelectedPreset);

    // Compare mode toggle
    document.getElementById('btnToggleCompare').addEventListener('click', () => {
        appState.compareActive = !appState.compareActive;
        syncControlsFromSetup();
        recomputeAndRender();
    });

    // Setup A / B tabs
    document.getElementById('tabSetupA').addEventListener('click', () => {
        appState.activeTab = 'A';
        syncControlsFromSetup();
        recomputeAndRender();
    });
    document.getElementById('tabSetupB').addEventListener('click', () => {
        appState.activeTab = 'B';
        appState.compareActive = true; // Automatically enable compare when editing Setup B
        syncControlsFromSetup();
        recomputeAndRender();
    });

    // Copy A -> B
    document.getElementById('btnCopyAtoB').addEventListener('click', () => {
        appState.setupB = { ...appState.setupA, name: `${appState.setupA.name} (Kopie B)` };
        appState.activeTab = 'B';
        appState.compareActive = true;
        syncControlsFromSetup();
        recomputeAndRender();
    });

    // Setup Name
    document.getElementById('inputSetupName').addEventListener('input', (e) => {
        getActiveSetup().name = e.target.value;
        document.getElementById(appState.activeTab === 'A' ? 'lblTabA' : 'lblTabB').textContent =
            e.target.value || `Setup ${appState.activeTab}`;
        recomputeAndRender();
    });

    // Mount Mode (OOP / IP)
    document.getElementById('btnMountOOP').addEventListener('click', () => {
        getActiveSetup().mountMode = 'OOP';
        syncControlsFromSetup();
        recomputeAndRender();
    });
    document.getElementById('btnMountIP').addEventListener('click', () => {
        getActiveSetup().mountMode = 'IP';
        syncControlsFromSetup();
        recomputeAndRender();
    });

    // Sliders
    bindSliderPair('PhiM', 'phiMDeg', true);
    bindSliderPair('PsiS', 'psiSDeg', true);
    bindSliderPair('ThetaCut', 'thetaCutDeg', false);
    bindSliderPair('PhiCut', 'phiCutDeg', true);
    bindSliderPair('AlphaFlat', 'alphaFlatDeg', false);
    bindSliderPair('ThetaDev', 'thetaDevDeg', false);

    document.getElementById('numB0').addEventListener('input', (e) => {
        getActiveSetup().b0Tesla = Math.max(0.001, parseFloat(e.target.value) || 0.5);
        recomputeAndRender();
    });

    // Quick-Value Chips
    const chipPropMap = { PhiM: 'phiMDeg', PsiS: 'psiSDeg', ThetaDev: 'thetaDevDeg' };
    document.querySelectorAll('.quick-chips').forEach((group) => {
        const target = group.dataset.target;
        group.querySelectorAll('button').forEach((btn) => {
            btn.addEventListener('click', () => {
                getActiveSetup()[chipPropMap[target]] = parseFloat(btn.dataset.val);
                syncControlsFromSetup();
                recomputeAndRender();
            });
        });
    });

    // Cut Presets
    document.querySelectorAll('#cutPresetGrid .cut-btn').forEach((btn) => {
        btn.addEventListener('click', () => {
            const cut = btn.dataset.cut;
            const s = getActiveSetup();
            s.cutPreset = cut;
            if (CUT_PRESETS[cut] && cut !== 'custom') {
                s.thetaCutDeg = CUT_PRESETS[cut].thetaCutDeg;
                s.phiCutDeg = CUT_PRESETS[cut].phiCutDeg;
                s.alphaFlatDeg = 0; // Reset to primary 0°-edge of the newly selected cut
            }
            syncControlsFromSetup();
            recomputeAndRender();
        });
    });

    // Popular 0°-Edge Presets per Cut (alpha_flat = 0° vs. 90°)
    document.querySelectorAll('#edgePresetGrid .edge-preset-card[data-alpha]').forEach((btn) => {
        btn.addEventListener('click', () => {
            getActiveSetup().alphaFlatDeg = parseFloat(btn.dataset.alpha) || 0;
            syncControlsFromSetup();
            recomputeAndRender();
        });
    });

    // Custom Miller-Bravais (hkil) Plane Apply
    const hInput = document.getElementById('hkilH');
    const kInput = document.getElementById('hkilK');
    const iInput = document.getElementById('hkilI');
    const lInput = document.getElementById('hkilL');
    const syncI = () => {
        iInput.value = -((parseFloat(hInput.value) || 0) + (parseFloat(kInput.value) || 0));
    };
    hInput.addEventListener('input', syncI);
    kInput.addEventListener('input', syncI);

    document.getElementById('btnApplyHKIL').addEventListener('click', () => {
        const h = parseFloat(hInput.value) || 0;
        const k = parseFloat(kInput.value) || 0;
        const l = parseFloat(lInput.value) || 0;
        const nCryst = millerBravaisPlaneToCrystalNormal([h, k, -(h + k), l]);
        const s = getActiveSetup();
        s.cutPreset = 'custom';
        s.thetaCutDeg = Number((Math.acos(Math.max(-1, Math.min(1, nCryst[2]))) * RAD2DEG).toFixed(3));
        s.phiCutDeg = Number(wrapDeg360(Math.atan2(nCryst[1], nCryst[0]) * RAD2DEG).toFixed(2));
        syncControlsFromSetup();
        recomputeAndRender();
    });

    // Device Type (CPW vs IDT)
    document.getElementById('btnDevCPW').addEventListener('click', () => {
        getActiveSetup().deviceType = 'CPW';
        syncControlsFromSetup();
        recomputeAndRender();
    });
    document.getElementById('btnDevIDT').addEventListener('click', () => {
        getActiveSetup().deviceType = 'IDT';
        syncControlsFromSetup();
        recomputeAndRender();
    });

    // Inverse Solver Controls (Modular Category Tabs + Target Cards)
    document.querySelectorAll('#solverCatTabs .solver-cat-tab').forEach((tab) => {
        tab.addEventListener('click', () => selectSolverCategory(tab.dataset.solverCat, true));
    });
    document.querySelectorAll('.solver-target-card[data-target-id]').forEach((card) => {
        card.addEventListener('click', () => setSolverTarget(card.dataset.targetId));
    });
    document.getElementById('selectTargetPreset').addEventListener('change', (e) => {
        setSolverTarget(e.target.value);
    });

    const tU = document.getElementById('targetU');
    const tV = document.getElementById('targetV');
    const tT = document.getElementById('targetT');
    const tW = document.getElementById('targetW');
    const updateCustomTarget = () => {
        const u = parseFloat(tU.value) || 0;
        const v = parseFloat(tV.value) || 0;
        const w = parseFloat(tW.value) || 0;
        tT.value = -(u + v);
        appState.customTargetUVTW = [u, v, -(u + v), w];
        recomputeAndRender();
    };
    [tU, tV, tW].forEach((el) => el.addEventListener('input', updateCustomTarget));

    // Frame Lock Pills (View 1): animated transition so one can watch the frames rotate
    document.querySelectorAll('#frameLockPills .frame-pill').forEach((pill) => {
        pill.addEventListener('click', () => {
            sceneController.setLockFrame(pill.dataset.frame, true);
            const svgOpt = document.getElementById('svgOptFrame');
            if (svgOpt) svgOpt.value = pill.dataset.frame;
            syncFramePills();
        });
    });

    // Frame Lock Pills (View 2)
    document.querySelectorAll('#frameLockPills2 .frame-pill').forEach((pill) => {
        pill.addEventListener('click', () => {
            document.querySelectorAll('#frameLockPills2 .frame-pill').forEach((p) => p.classList.toggle('active', p === pill));
            if (sceneController2) sceneController2.setLockFrame(pill.dataset.frame, true);
        });
    });

    // Legend chips toggle layers (shared by both views)
    document.querySelectorAll('#canvasLegend .legend-chip[data-vis]').forEach((chip) => {
        chip.addEventListener('click', () => {
            const prop = chip.dataset.vis;
            sceneController.visibility[prop] = !sceneController.visibility[prop];
            syncLegend();
            rerenderScenes();
        });
    });

    // Vector table: frame & display mode
    document.getElementById('btnVectorTableAll').addEventListener('click', () => {
        appState.vectorTableAll = !appState.vectorTableAll;
        if (lastActiveState) updateVectorTable(lastActiveState);
    });
    document.getElementById('selectVectorTableFrame').addEventListener('change', () => {
        if (lastActiveState) updateVectorTable(lastActiveState);
    });
    document.querySelectorAll('#vectorTableMode .seg-btn').forEach((btn) => {
        btn.addEventListener('click', () => {
            appState.vectorTableMode = btn.dataset.vmode;
            document.querySelectorAll('#vectorTableMode .seg-btn').forEach((b) => b.classList.toggle('active', b === btn));
            if (lastActiveState) updateVectorTable(lastActiveState);
        });
    });

    // 1A: Sample Mounting Preset Buttons (if present)
    document.querySelectorAll('.mount-preset-btn[data-mpreset]').forEach((btn) => {
        btn.addEventListener('click', () => {
            const s = getActiveSetup();
            s.psiSDeg = Number(computeMountingPresetPsiS(s, btn.dataset.mpreset).toFixed(2));
            syncControlsFromSetup();
            recomputeAndRender();
        });
    });

    // Share link
    document.getElementById('btnCopyLink').addEventListener('click', () => {
        writeStateToUrl();
        navigator.clipboard.writeText(window.location.href).then(
            () => showToast('Link mit allen Einstellungen kopiert'),
            () => showToast('Kopieren nicht möglich – Link steht in der Adresszeile')
        );
    });

    bindCustomDirection();
    bindEmphasisHover();
    bindExportModal();
    bindTour();

    document.addEventListener('keydown', (e) => {
        if (e.key !== 'Escape') return;
        if (tour && tour.active) tour.stop();
        document.getElementById('modalExport').classList.add('hidden');
    });
}

// ============================================================================
// Target finder helpers
// ============================================================================

function setSolverTarget(targetId) {
    appState.targetId = targetId;
    const selEl = document.getElementById('selectTargetPreset');
    if (selEl) selEl.value = targetId;
    document.querySelectorAll('.solver-target-card[data-target-id]').forEach((card) => {
        card.classList.toggle('active', card.dataset.targetId === targetId);
    });
    document.getElementById('customTargetMBBox').classList.toggle('hidden', targetId !== 'custom_mb');
    recomputeAndRender();
}

function selectSolverCategory(cat, pickTarget) {
    document.querySelectorAll('#solverCatTabs .solver-cat-tab').forEach((t) => {
        t.classList.toggle('active', t.dataset.solverCat === cat);
    });
    document.querySelectorAll('.solver-target-grid[data-solver-panel]').forEach((panel) => {
        const isActivePanel = panel.dataset.solverPanel === cat;
        panel.classList.toggle('hidden', !isActivePanel);
        if (isActivePanel && pickTarget) {
            const activeCard = panel.querySelector('.solver-target-card.active') || panel.querySelector('.solver-target-card');
            if (activeCard && activeCard.dataset.targetId) setSolverTarget(activeCard.dataset.targetId);
        }
    });
}

// ============================================================================
// "Eigene Richtung": any crystal direction / plane normal / cross product / angles
// ============================================================================

function readInt(id) {
    return parseFloat(document.getElementById(id).value) || 0;
}

function bindCustomDirection() {
    const labels = {
        dir: ['u', 'v', 't', 'w'],
        plane: ['h', 'k', 'i', 'l']
    };
    document.querySelectorAll('#customDirModes .cd-mode').forEach((btn) => {
        btn.addEventListener('click', () => {
            const mode = btn.dataset.cdmode;
            appState.customDir.mode = mode;
            document.querySelectorAll('#customDirModes .cd-mode').forEach((b) => b.classList.toggle('active', b === btn));
            setIndexBrackets(document.getElementById('cdMainIndex'), mode === 'plane');
            document.querySelectorAll('.custom-dir-panel[data-cdpanel]').forEach((p) => {
                p.classList.toggle('hidden', !p.dataset.cdpanel.split(' ').includes(mode));
            });
            if (labels[mode]) {
                document.querySelectorAll('[data-cd-lbl]').forEach((lbl) => {
                    lbl.textContent = labels[mode][Number(lbl.dataset.cdLbl)];
                });
            }
            sceneController.visibility.showCustom = true;
            syncLegend();
            if (lastActiveState) updateCustomDirection(lastActiveState);
        });
    });

    [1, 2].forEach((n) => {
        document.getElementById(`crossMode${n}`).addEventListener('change', (e) => {
            setIndexBrackets(document.querySelector(`.mi-input[data-cross="${n}"]`), e.target.value === 'plane');
        });
    });

    const ids = ['cdH', 'cdK', 'cdL', 'cross1_h', 'cross1_k', 'cross1_l', 'cross2_h', 'cross2_k', 'cross2_l', 'crossMode1', 'crossMode2', 'cdAngFrame', 'cdTheta', 'cdPhi'];
    ids.forEach((id) => {
        const el = document.getElementById(id);
        if (!el) return;
        const handler = () => {
            if (lastActiveState) {
                updateCustomDirection(lastActiveState);
                rerenderScenes();
            }
        };
        el.addEventListener('input', handler);
        el.addEventListener('change', handler);
    });

    // Opening the card draws the arrow, closing it removes it again
    document.getElementById('customDirCard').addEventListener('toggle', (e) => {
        sceneController.visibility.showCustom = e.target.open;
        syncLegend();
        if (lastActiveState) updateVectorTable(lastActiveState);
        rerenderScenes();
    });

    document.getElementById('btnCustomDirDrive').addEventListener('click', driveFieldToCustomDirection);
}

// Directions are written [uvtw], planes (hkil)
function setIndexBrackets(container, isPlane) {
    const [open, close] = container.querySelectorAll('.mi-br');
    open.textContent = isPlane ? '(' : '[';
    close.textContent = isPlane ? ')' : ']';
}

// Returns the custom direction as a unit vector in the crystal frame (or null).
function computeCustomDirCryst(state) {
    const s = state.setup;
    const a = s.aLattice;
    const c = s.cLattice;
    const mode = appState.customDir.mode;

    if (mode === 'dir' || mode === 'plane') {
        const h = readInt('cdH');
        const k = readInt('cdK');
        const l = readInt('cdL');
        document.getElementById('cdI').value = -(h + k);
        const mb = [h, k, -(h + k), l];
        if (mb.every((x) => x === 0)) return null;
        return mode === 'plane' ? millerBravaisPlaneToCrystalNormal(mb, a, c) : millerBravaisDirectionToCrystal(mb, a, c);
    }

    if (mode === 'cross') {
        const read = (n) => {
            const h = readInt(`cross${n}_h`);
            const k = readInt(`cross${n}_k`);
            const l = readInt(`cross${n}_l`);
            document.getElementById(`cross${n}_i`).value = -(h + k);
            const mb = [h, k, -(h + k), l];
            return document.getElementById(`crossMode${n}`).value === 'plane'
                ? millerBravaisPlaneToCrystalNormal(mb, a, c)
                : millerBravaisDirectionToCrystal(mb, a, c);
        };
        const v = cross3(read(1), read(2));
        return Math.hypot(...v) < 1e-9 ? null : normalize3(v);
    }

    // mode === 'angles'
    const frame = document.getElementById('cdAngFrame').value;
    const th = (parseFloat(document.getElementById('cdTheta').value) || 0) * DEG2RAD;
    const ph = (parseFloat(document.getElementById('cdPhi').value) || 0) * DEG2RAD;
    const vF = [Math.sin(th) * Math.cos(ph), Math.sin(th) * Math.sin(ph), Math.cos(th)];
    return normalize3(matVec3(state.matrices.getRotationMatrix('Cryst', frame), vF));
}

// Miller-Bravais indices; integers if they reproduce the direction within 0.5°, else scaled decimals.
// (formatMillerBravais* would silently round decimals to small integers, so decimals are formatted here.)
function bestMillerBravais(vCryst, a, c, isPlane) {
    const raw = isPlane ? crystalNormalToMillerBravaisPlane(vCryst, a, c) : crystalToMillerBravaisDirection(vCryst, a, c);
    const ints = scaleMillerBravaisInts(raw);
    const back = isPlane ? millerBravaisPlaneToCrystalNormal(ints, a, c) : millerBravaisDirectionToCrystal(ints, a, c);
    const errDeg = Math.acos(Math.max(-1, Math.min(1, dot3(back, vCryst)))) * RAD2DEG;
    if (errDeg < 0.5) {
        return {
            indices: ints,
            exact: true,
            latex: formatMillerBravaisLatex(ints, isPlane),
            text: formatMillerBravais(ints, isPlane)
        };
    }
    const maxAbs = Math.max(...raw.map(Math.abs)) || 1;
    const indices = raw.map((x) => Number((x / maxAbs).toFixed(4)));
    const vals = indices.map((x) => (Math.abs(x) < 0.005 ? 0 : x));
    const [open, close] = isPlane ? ['(', ')'] : ['[', ']'];
    const dec = (x) => String(Number(x.toFixed(2)));
    return {
        indices,
        exact: false,
        latex: `${open}${vals.map((x) => (x < 0 ? `\\overline{${dec(-x)}}` : dec(x))).join('\\;')}${close}`,
        text: `${open}${vals.map(dec).join(' ')}${close}`
    };
}

const angleDeg = (u, v) => Math.acos(Math.max(-1, Math.min(1, dot3(normalize3(u), normalize3(v))))) * RAD2DEG;

function updateCustomDirection(state) {
    const out = document.getElementById('customDirResult');
    const vCryst = computeCustomDirCryst(state);
    if (!vCryst) {
        sceneController.customVector = null;
        if (sceneController2) sceneController2.customVector = null;
        out.innerHTML = '<p class="muted">Bitte eine Richtung ungleich null eingeben.</p>';
        return;
    }
    const s = state.setup;
    const v = state.vectors;
    const R = (to) => state.matrices.getRotationMatrix(to, 'Cryst');
    const vLab = matVec3(R('Lab'), vCryst);
    const vSub = matVec3(R('Sub'), vCryst);

    const dir = bestMillerBravais(vCryst, s.aLattice, s.cLattice, false);
    const plane = bestMillerBravais(vCryst, s.aLattice, s.cLattice, true);
    appState.customDir.uvtw = dir.indices;
    appState.customDir.latex = dir.latex;

    sceneController.customVector = { Lab: vLab, label: dir.text };
    if (sceneController2) sceneController2.customVector = sceneController.customVector;

    // Tilt out of the sample plane (0° = lies in the sample)
    const elevSub = Math.asin(Math.min(1, Math.abs(vSub[2]))) * RAD2DEG;
    const stat = (val, labelHtml) => `<div class="cd-stat"><strong>${val}</strong><span>${labelHtml}</span></div>`;
    const planeText = appState.customDir.mode === 'plane' ? 'Normale der Ebene' : 'steht senkrecht auf';
    out.innerHTML = `
        <div class="cd-res-head" data-emph="custom">
            <span class="cd-res-idx">${tex(dir.latex)}</span>
            <span class="cd-res-class">${shortDirection(vCryst)}</span>
        </div>
        <div class="cd-res-sub">${planeText} ${tex(plane.latex)}${dir.exact ? '' : ' · gerundet'}</div>
        <div class="cd-stats">
            ${stat(`${angleDeg(vCryst, [0, 0, 1]).toFixed(1)}°`, `zur ${token('ec', 'c-Achse')}`)}
            ${stat(`${elevSub.toFixed(1)}°`, 'aus der Probe')}
            ${stat(`${angleDeg(vSub, v.kHat.Sub).toFixed(1)}°`, `zu ${token('k', tex('\\hat{k}'))}`)}
            ${stat(`${angleDeg(vSub, v.bHat.Sub).toFixed(1)}°`, `zu ${token('B', tex('\\mathbf{B}'))}`)}
        </div>`;
}

function driveFieldToCustomDirection() {
    const uvtw = appState.customDir.uvtw;
    if (!uvtw) return;
    appState.customTargetUVTW = [...uvtw];
    document.getElementById('targetU').value = uvtw[0];
    document.getElementById('targetV').value = uvtw[1];
    document.getElementById('targetT').value = uvtw[2];
    document.getElementById('targetW').value = uvtw[3];

    selectSolverCategory('crystal', false);
    setSolverTarget('custom_mb');

    const stepCard3 = document.getElementById('stepCard3');
    stepCard3.classList.add('step-card-active');
    const details = stepCard3.querySelector('.solver-details');
    details.open = true;
    details.scrollIntoView({ behavior: 'smooth', block: 'start' });
    details.classList.add('flash');
    setTimeout(() => details.classList.remove('flash'), 1200);
}

// ============================================================================
// Guided tour
// ============================================================================

function bindTour() {
    tour = new GuidedTour({
        scene: sceneController,
        getSetup: getActiveSetup,
        evaluate: evaluateSetup,
        tex,
        refresh: recomputeAndRender,
        onVisibilityChange: () => {
            syncLegend();
            syncFramePills();
        },
        els: {
            root: document.getElementById('tourCaption'),
            title: document.getElementById('tourTitle'),
            text: document.getElementById('tourText'),
            step: document.getElementById('tourDots'),
            prev: document.getElementById('tourPrev'),
            next: document.getElementById('tourNext'),
            exit: document.getElementById('tourExit')
        }
    });
    document.getElementById('btnStartTour').addEventListener('click', () => {
        if (tour.active) tour.stop();
        else tour.start();
    });
}

// ============================================================================
// Shareable URL (#cut=r-cut&tc=…&b.cut=…&cmp=1)
// ============================================================================

function writeStateToUrl() {
    const params = new URLSearchParams();
    const put = (prefix, setup) => {
        for (const [key, prop] of Object.entries(URL_KEYS)) {
            const val = setup[prop];
            if (val === undefined) continue;
            params.set(prefix + key, URL_STRING_PROPS.has(prop) ? val : String(Number(Number(val).toFixed(3))));
        }
    };
    put('', appState.setupA);
    if (appState.compareActive) {
        put('b.', appState.setupB);
        params.set('cmp', '1');
    }
    params.set('view', sceneController ? sceneController.lockFrame : 'Sub');
    history.replaceState(null, '', `#${params.toString()}`);
}

function readStateFromUrl() {
    if (!window.location.hash || window.location.hash.length < 2) return;
    const params = new URLSearchParams(window.location.hash.slice(1));
    const apply = (prefix, setup) => {
        let found = false;
        for (const [key, prop] of Object.entries(URL_KEYS)) {
            if (!params.has(prefix + key)) continue;
            const raw = params.get(prefix + key);
            if (URL_STRING_PROPS.has(prop)) {
                setup[prop] = raw;
            } else if (Number.isFinite(parseFloat(raw))) {
                setup[prop] = parseFloat(raw);
            }
            found = true;
        }
        if (found) setup.name = prefix ? 'Setup B (aus Link)' : 'Setup A (aus Link)';
    };
    apply('', appState.setupA);
    apply('b.', appState.setupB);
    appState.compareActive = params.get('cmp') === '1';
    const view = params.get('view');
    if (view && FRAME_METADATA[view] && view !== 'Mount') sceneController.lockFrame = view;
}

// ============================================================================
// Presets
// ============================================================================

function applySavedPreset(presetId) {
    const found = savedPresets.find((p) => p.id === presetId);
    if (!found) return;
    activePresetId = found.id;

    const targetSetup = createDefaultSetup({ ...found.setup, name: found.name });
    if (appState.activeTab === 'B') {
        appState.setupB = targetSetup;
    } else {
        appState.setupA = targetSetup;
    }

    syncControlsFromSetup();
    recomputeAndRender();
}

function saveCurrentAsPreset(name) {
    const s = getActiveSetup();
    const id = `preset_${Date.now()}`;
    savedPresets.push({
        id,
        name,
        setup: {
            mountMode: s.mountMode,
            cutPreset: s.cutPreset,
            thetaCutDeg: s.thetaCutDeg,
            phiCutDeg: s.phiCutDeg,
            alphaFlatDeg: s.alphaFlatDeg,
            deviceType: s.deviceType,
            thetaDevDeg: s.thetaDevDeg,
            psiSDeg: s.psiSDeg,
            phiMDeg: s.phiMDeg,
            b0Tesla: s.b0Tesla
        }
    });
    s.name = name;
    activePresetId = id;
    savePresetsToStorage();
    renderPresetDropdown();
    syncControlsFromSetup();
    recomputeAndRender();
}

function deleteSelectedPreset() {
    const sel = document.getElementById('selectScenarioPreset');
    const idToDelete = (sel && sel.value) || activePresetId;
    if (!idToDelete) return;

    savedPresets = savedPresets.filter((p) => p.id !== idToDelete);
    activePresetId = savedPresets.length > 0 ? savedPresets[0].id : '';
    savePresetsToStorage();
    renderPresetDropdown();
    if (activePresetId) applySavedPreset(activePresetId);
}

// ============================================================================
// Core Evaluation & DOM Update Pipeline
// ============================================================================

function recomputeAndRender() {
    const stateA = evaluateSetup(appState.setupA);
    const stateB = evaluateSetup(appState.setupB);
    const activeState = appState.activeTab === 'A' ? stateA : stateB;
    const otherState = appState.activeTab === 'A' ? stateB : stateA;
    lastActiveState = activeState;

    // Sidebar
    updateStepSummaries(activeState);
    document.getElementById('badgeSubXs').innerHTML = tex(activeState.millerBravais.xs.latex);
    document.getElementById('badgeSubYs').innerHTML = tex(activeState.millerBravais.ys.latex);
    document.getElementById('badgeSubZs').innerHTML = tex(activeState.millerBravais.zsPlane.latex);
    updateDeviceMillerReadout(activeState);
    updateMountingAndSmartJumps(activeState);
    updateInverseSolver(activeState);

    // The custom direction depends on the setup (it is drawn in Lab coordinates)
    updateCustomDirection(activeState);

    // 3D views (paused while the tour animates its own temporary field sweep)
    if (!(tour && tour.active && tour._raf)) {
        sceneController.updateStates(activeState, otherState, appState.compareActive);
    }
    if (sceneController2) sceneController2.updateStates(activeState, otherState, appState.compareActive);

    // Overview column
    updateMetricsCards(activeState, stateA, stateB);
    document.getElementById('gistBody').innerHTML = buildSummary(activeState, tex);
    document.getElementById('narrativeBody').innerHTML = buildNarrative(activeState, tex);
    updateVectorTable(activeState);
    updateCompareTable(stateA, stateB);

    const exportModal = document.getElementById('modalExport');
    if (exportModal && !exportModal.classList.contains('hidden')) refreshExportPanel();

    writeStateToUrl();
}

function updateStepSummaries(activeState) {
    const s = activeState.setup;
    const mb = activeState.millerBravais;
    const sum1 = document.getElementById('summaryStep1');
    const sumDev = document.getElementById('summaryStepDev');
    const sum2 = document.getElementById('summaryStep2');
    const sum3 = document.getElementById('summaryStep3');
    if (sum1) sum1.innerHTML = `${s.cutPreset} · 0°-Kante ${tex(mb.xs.latex)}`;
    if (sumDev) sumDev.innerHTML = `${s.deviceType} · ${s.thetaDevDeg.toFixed(0)}° (${tex(mb.k.latex)})`;
    if (sum2) sum2.innerHTML = `${s.mountMode === 'OOP' ? 'Vertikal (OOP)' : 'Horizontal (IP)'} · ${tex(`\\psi_S = ${s.psiSDeg.toFixed(1)}^\\circ`)}`;
    if (sum3) sum3.innerHTML = `${tex(`\\phi_M = ${s.phiMDeg.toFixed(1)}^\\circ`)} · ${s.b0Tesla} T`;
}

function updateMountingAndSmartJumps(activeState) {
    const landmarks = solveSweepLandmarksForFixedPsi(activeState.setup);

    const jumpItems = [
        { label: `Feld ${tex('\\parallel')} Welle ${tex('\\hat{k}')}`, phi: landmarks.bParallelK.phiMDeg, emph: 'B,k' },
        { label: `Feld ${tex('\\perp')} Welle ${tex('\\hat{k}')}`, phi: landmarks.bPerpK.phiMDeg, emph: 'B,k' },
        { label: 'Feld in Easy-Plane', phi: landmarks.bInEasyPlane.phiM1, emph: 'B,easy' },
        { label: `Moment ${tex('\\hat{m} \\parallel')} Welle ${tex('\\hat{k}')}`, phi: landmarks.mAlignK.phiMDeg, emph: 'm,k' }
    ];

    const jumpsContainer = document.getElementById('smartMagnetJumps');
    jumpsContainer.innerHTML = jumpItems
        .map((j) => `
        <button type="button" class="smart-jump-btn" data-jump-phi="${j.phi}" data-emph="${j.emph}">
            <span>${j.label}</span>
            <code>${j.phi.toFixed(1)}°</code>
        </button>`)
        .join('');

    jumpsContainer.querySelectorAll('button[data-jump-phi]').forEach((btn) => {
        btn.addEventListener('click', () => {
            getActiveSetup().phiMDeg = Number(parseFloat(btn.dataset.jumpPhi).toFixed(2));
            syncControlsFromSetup();
            recomputeAndRender();
        });
    });
}

function updateDeviceMillerReadout(activeState) {
    const s = activeState.setup;
    const mb = activeState.millerBravais;

    const elK = document.getElementById('liveDevMillerK');
    const elT = document.getElementById('liveDevMillerT');
    const sumEdge0 = document.getElementById('summaryEdge0');
    if (elK) elK.innerHTML = tex(mb.k.latex);
    if (elT) elT.innerHTML = tex(mb.t.latex);
    if (sumEdge0) sumEdge0.innerHTML = tex(mb.xs.latex);

    // Update the 2 popular 0°-edge preset cards for the active cut
    const st0 = evaluateSetup({ ...s, alphaFlatDeg: 0 });
    const st90 = evaluateSetup({ ...s, alphaFlatDeg: 90 });
    const xs0 = tex(st0.millerBravais.xs.latex);
    const ys0 = tex(st0.millerBravais.ys.latex);
    const xs90 = tex(st90.millerBravais.xs.latex);
    const ys90 = tex(st90.millerBravais.ys.latex);

    const cutLabels = {
        'c-cut': {
            t0: `0° ${tex('\\parallel')} ${tex('a')}-Achse ${xs0}`,
            s0: `90°-Kante: ${tex('m')}-Achse ${ys0}`,
            t90: `0° ${tex('\\parallel')} ${tex('m')}-Achse ${xs90}`,
            s90: `90°-Kante: ${tex('a')}-Achse ${ys90}`
        },
        'a-cut': {
            t0: `0° ${tex('\\parallel')} ${tex('m')}-Achse ${xs0}`,
            s0: `In Easy-Plane (90°: ${tex('c')} ${ys0})`,
            t90: `0° ${tex('\\parallel')} ${tex('c')}-Achse ${xs90}`,
            s90: `Senkrecht zu Easy-Plane (90°: ${ys90})`
        },
        'm-cut': {
            t0: `0° ${tex('\\parallel')} ${tex('a')}-Achse ${xs0}`,
            s0: `In Easy-Plane (90°: ${tex('c')} ${ys0})`,
            t90: `0° ${tex('\\parallel')} ${tex('c')}-Achse ${xs90}`,
            s90: `Senkrecht zu Easy-Plane (90°: ${ys90})`
        },
        'r-cut': {
            t0: `0° ${tex('\\parallel')} ${tex('a')}-Achse ${xs0}`,
            s0: `In Easy-Plane (90°: ${ys0})`,
            t90: `0° ${tex('\\parallel')} ${xs90} (proj. ${tex('c')})`,
            s90: `90°-Kante in Easy-Plane: ${ys90}`
        },
        'custom': {
            t0: `0° ${tex('\\parallel')} ${xs0} (Easy-Plane)`,
            s0: `90°-Kante: ${ys0}`,
            t90: `0° ${tex('\\parallel')} ${xs90} (proj. ${tex('c')})`,
            s90: `90°-Kante: ${ys90}`
        }
    };
    const info = cutLabels[s.cutPreset] || cutLabels.custom;

    const setHtml = (id, html) => {
        const el = document.getElementById(id);
        if (el) el.innerHTML = html;
    };
    setHtml('lblEdgePreset0Title', info.t0);
    setHtml('lblEdgePreset0Sub', info.s0);
    setHtml('lblEdgePreset90Title', info.t90);
    setHtml('lblEdgePreset90Sub', info.s90);
    document.getElementById('btnEdgePreset0').classList.toggle('active', Math.abs(s.alphaFlatDeg) < 0.25);
    document.getElementById('btnEdgePreset90').classList.toggle('active', Math.abs(s.alphaFlatDeg - 90) < 0.25);

    const quickGroup = document.getElementById('deviceQuickAngleChips');
    if (quickGroup) {
        quickGroup.querySelectorAll('button[data-val]').forEach((btn) => {
            btn.classList.toggle('active', Math.abs(s.thetaDevDeg - parseFloat(btn.dataset.val)) < 0.25);
        });
    }
}

function updateMetricsCards(activeState, stateA, stateB) {
    const m = activeState.metrics;
    const isCPW = activeState.setup.deviceType === 'CPW';

    document.getElementById('valBSampleAngles').textContent = `${m.angleBKDeg.toFixed(1)}°`;
    const valBEasy = document.getElementById('valBEasy');
    valBEasy.innerHTML = `<strong>${(m.bEasyRatio * 100).toFixed(0)} %</strong> in Easy-Plane`;
    valBEasy.title = `${m.bEasyTesla.toFixed(2)} T von ${activeState.setup.b0Tesla} T liegen in der Easy-Plane`;

    document.getElementById('valAngleMK').textContent = m.isEasyPlaneSingular ? '–' : `${m.angleMKDeg.toFixed(1)}°`;
    const valMDotK = document.getElementById('valMDotK');
    valMDotK.innerHTML = isCPW ? `Kopplung <strong>${(m.cpwTorqueM * 100).toFixed(0)} %</strong>` : '';
    valMDotK.title = isCPW ? 'CPW-Kopplung |k̂ × m̂|' : '';

    document.getElementById('valAngleLK').textContent = m.isEasyPlaneSingular ? '–' : `${m.angleLKDeg.toFixed(1)}°`;

    // Compare delta badges on result tiles
    ['cmpBSample', 'cmpAngleMK', 'cmpAngleLK'].forEach((id) => {
        document.getElementById(id).classList.toggle('hidden', !appState.compareActive);
    });
    if (appState.compareActive) {
        const mA = stateA.metrics;
        const mB = stateB.metrics;
        document.getElementById('cmpBSample').textContent = `A ${mA.angleBKDeg.toFixed(1)}° · B ${mB.angleBKDeg.toFixed(1)}°`;
        document.getElementById('cmpAngleMK').textContent = `A ${mA.angleMKDeg.toFixed(1)}° · B ${mB.angleMKDeg.toFixed(1)}°`;
        document.getElementById('cmpAngleLK').textContent = `A ${mA.angleLKDeg.toFixed(1)}° · B ${mB.angleLKDeg.toFixed(1)}°`;
    }
}

function updateCompareTable(stateA, stateB) {
    const sec = document.getElementById('compareSummarySection');
    sec.classList.toggle('hidden', !appState.compareActive);
    if (!appState.compareActive) return;

    const mA = stateA.metrics;
    const mB = stateB.metrics;
    const sA = stateA.setup;
    const sB = stateB.setup;

    const rows = [
        ['Einbau / Schnitt', `${sA.mountMode} / ${sA.cutPreset}`, `${sB.mountMode} / ${sB.cutPreset}`, '—'],
        [tex('\\psi_S, \\phi_M'), `${sA.psiSDeg.toFixed(1)}°, ${sA.phiMDeg.toFixed(1)}°`, `${sB.psiSDeg.toFixed(1)}°, ${sB.phiMDeg.toFixed(1)}°`, `${(sB.phiMDeg - sA.phiMDeg).toFixed(1)}°`],
        [`Welle ${tex('\\hat{k}')}`, tex(stateA.millerBravais.k.latex), tex(stateB.millerBravais.k.latex), `${(sB.thetaDevDeg - sA.thetaDevDeg).toFixed(1)}°`],
        ['in Easy-Plane', `${(mA.bEasyRatio * 100).toFixed(1)} %`, `${(mB.bEasyRatio * 100).toFixed(1)} %`, `${((mB.bEasyRatio - mA.bEasyRatio) * 100).toFixed(1)} %`],
        ['in Probenebene', `${(mA.bInSampleRatio * 100).toFixed(1)} %`, `${(mB.bInSampleRatio * 100).toFixed(1)} %`, `${((mB.bInSampleRatio - mA.bInSampleRatio) * 100).toFixed(1)} %`],
        [tex('\\angle(\\mathbf{B}, \\hat{k})'), `${mA.angleBKDeg.toFixed(1)}°`, `${mB.angleBKDeg.toFixed(1)}°`, `${(mB.angleBKDeg - mA.angleBKDeg).toFixed(1)}°`],
        [tex('\\angle(\\hat{m}, \\hat{k})'), `${mA.angleMKDeg.toFixed(1)}°`, `${mB.angleMKDeg.toFixed(1)}°`, `${(mB.angleMKDeg - mA.angleMKDeg).toFixed(1)}°`],
        [tex('\\angle(\\hat{l}, \\hat{k})'), `${mA.angleLKDeg.toFixed(1)}°`, `${mB.angleLKDeg.toFixed(1)}°`, `${(mB.angleLKDeg - mA.angleLKDeg).toFixed(1)}°`],
        [tex('|\\hat{k}\\times\\hat{m}|'), mA.cpwTorqueM.toFixed(3), mB.cpwTorqueM.toFixed(3), (mB.cpwTorqueM - mA.cpwTorqueM).toFixed(3)]
    ];

    document.getElementById('compareTableBody').innerHTML = rows
        .map((r) => `<tr><td>${r[0]}</td><td>${r[1]}</td><td>${r[2]}</td><td>${r[3]}</td></tr>`)
        .join('');
}

// ============================================================================
// "Wo zeigt was hin?" vector table
// ============================================================================

const TABLE_AXIS_TEX = {
    Sub: ['x_s', 'y_s', '\\hat{n}'],
    Wave: ['\\hat{k}', '\\hat{t}', '\\hat{n}'],
    Cryst: ['a', 'm', 'c'],
    Lab: ['X_L', 'Y_L', 'Z_L']
};
const TABLE_FOOT = {
    Sub: 'θ = 90°: liegt in der Probe',
    Wave: 'θ = 90°: liegt in der Probe',
    Cryst: 'θ = 90°: liegt in der Easy-Plane',
    Lab: 'θ = 90°: liegt horizontal'
};
const TABLE_CORE_KEYS = ['B', 'm', 'l', 'k', 'custom'];

function updateVectorTable(activeState) {
    const sel = document.getElementById('selectVectorTableFrame').value;
    const frame = sel === 'view' ? sceneController.lockFrame : sel;
    const axes = TABLE_AXIS_TEX[frame] || TABLE_AXIS_TEX.Sub;
    const angles = appState.vectorTableMode === 'angles';

    document.getElementById('vectorTableHeaderRow').innerHTML = angles
        ? `<th>Richtung</th><th>θ zu ${tex(axes[2])}</th><th>φ ab ${tex(axes[0])}</th><th>Kristall</th>`
        : `<th>Richtung</th><th>${tex(axes[0])}</th><th>${tex(axes[1])}</th><th>${tex(axes[2])}</th><th>Kristall</th>`;
    document.getElementById('vectorTableFoot').textContent = angles ? TABLE_FOOT[frame] || '' : 'Einheitsvektoren';

    const v = activeState.vectors;
    const mb = activeState.millerBravais;
    const s = activeState.setup;
    const singular = activeState.metrics.isEasyPlaneSingular;
    const tiltNC = Math.min(activeState.metrics.angleNCDeg, 180 - activeState.metrics.angleNCDeg);

    const entries = [
        { key: 'B', name: `Feld ${tex('\\mathbf{B}')}`, vec: v.bHat, mb: mb.B.latex },
        { key: 'm', name: `Moment ${tex('\\hat{m}')}`, vec: v.mHat, mb: mb.m.latex, skip: singular },
        { key: 'l', name: `Néel ${tex('\\hat{l}')}`, vec: v.lHat, mb: mb.l.latex, skip: singular },
        { key: 'k', name: `Welle ${tex('\\hat{k}')}`, vec: v.kHat, mb: mb.k.latex },
        { key: 't', name: `Leiter ${tex('\\hat{t}')}`, vec: v.tHat, mb: mb.t.latex },
        { key: 'n', name: `Normale ${tex('\\hat{n}')}`, vec: v.nHat, mb: mb.zsDir.latex },
        { key: 'xs', name: '0°-Kante', vec: v.xsHat, mb: mb.xs.latex },
        { key: 'ys', name: '90°-Kante', vec: v.ysHat, mb: mb.ys.latex },
        { key: 'ea', name: `${tex('a')}-Achse`, vec: v.eaHat, mb: '[11\\bar{2}0]' },
        { key: 'em', name: `${tex('m')}-Achse`, vec: v.emHat, mb: '[1\\bar{1}00]' },
        { key: 'ec', name: `${tex('c')}-Achse`, vec: v.ecHat, mb: '[0001]' },
        {
            key: 'nodes',
            name: 'Schnittlinie',
            vec: v.nodesHat,
            mb: formatMillerBravaisLatex(crystalToMillerBravaisDirection(v.nodesHat.Cryst, s.aLattice, s.cLattice)),
            skip: tiltNC < 1
        }
    ];
    if (sceneController.visibility.showCustom && sceneController.customVector && appState.customDir.uvtw) {
        const lab = sceneController.customVector.Lab;
        const R = activeState.matrices.getRotationMatrix(frame, 'Lab');
        entries.push({
            key: 'custom',
            name: 'Eigene',
            vecInFrame: matVec3(R, lab),
            mb: appState.customDir.latex
        });
    }

    const fmtComp = (x) => {
        const val = Math.abs(x) < 5e-5 ? 0 : x;
        return (val >= 0 ? '+' : '−') + Math.abs(val).toFixed(3);
    };
    const color = (key) => (key === 'xs' || key === 'ys' ? VECTOR_COLORS.edge : VECTOR_COLORS[key]);

    const shown = entries.filter((e) => !e.skip);
    const visible = appState.vectorTableAll ? shown : shown.filter((e) => TABLE_CORE_KEYS.includes(e.key));
    document.getElementById('btnVectorTableAll').textContent = appState.vectorTableAll
        ? 'weniger ▴'
        : `+ ${shown.length - visible.length} weitere (Achsen, Kanten …)`;

    document.getElementById('vectorTableBody').innerHTML = visible
        .map((e) => {
            const vec = e.vecInFrame || e.vec[frame];
            const th = Math.acos(Math.max(-1, Math.min(1, vec[2]))) * RAD2DEG;
            let ph = wrapDeg360(Math.atan2(vec[1], vec[0]) * RAD2DEG);
            if (Math.abs(ph - 360) < 0.05) ph = 0;
            const phTxt = Math.hypot(vec[0], vec[1]) < 1e-6 ? '–' : `${ph.toFixed(1)}°`;
            const cells = angles
                ? `<td class="num">${th.toFixed(1)}°</td><td class="num">${phTxt}</td>`
                : `<td class="num">${fmtComp(vec[0])}</td><td class="num">${fmtComp(vec[1])}</td><td class="num">${fmtComp(vec[2])}</td>`;
            return `<tr data-emph="${e.key}"><td><span class="vt-dot" style="background:${color(e.key)}"></span>${e.name}</td>${cells}<td>${tex(e.mb)}</td></tr>`;
        })
        .join('');
}

// ============================================================================
// Target finder results
// ============================================================================

function updateInverseSolver() {
    const setup = getActiveSetup();
    const sol = solveTargetAngles(setup, appState.targetId, appState.customTargetUVTW);
    const container = document.getElementById('solverResultsContainer');

    const oop1 = sol.oopSolutions[0];
    const oop2 = sol.oopSolutions[1];
    const ip1 = sol.ipSolutions[0];
    const fixed = sol.fixedPsiSweep;

    container.innerHTML = `
        <div class="solver-card">
            <div class="solver-card-info">
                <strong>1. Vertikal-Halterung (OOP – Exakte 3D-Lösung)</strong>
                <div>${tex(`\\psi_S = ${oop1.psiSDeg.toFixed(2)}^\\circ`)} &nbsp;|&nbsp; ${tex(`\\phi_M = ${oop1.phiMDeg.toFixed(2)}^\\circ`)}</div>
                <small>Alternativ (+180°): ${tex(`\\psi_S = ${oop2.psiSDeg.toFixed(1)}^\\circ,\\; \\phi_M = ${oop2.phiMDeg.toFixed(1)}^\\circ`)}</small>
            </div>
            <button type="button" class="btn btn-xs btn-primary" data-apply-mount="OOP" data-psi="${oop1.psiSDeg}" data-phi="${oop1.phiMDeg}">
                Anwenden
            </button>
        </div>
        <div class="solver-card">
            <div class="solver-card-info">
                <strong>2. Nur Magnet drehen (aktuelles ${tex(`\\psi_S = ${fixed.psiSDeg.toFixed(1)}^\\circ`)} fixiert)</strong>
                <div>${tex(`\\phi_M = ${fixed.bestPhiMDeg.toFixed(2)}^\\circ`)} &nbsp;<small>(oder ${fixed.oppositePhiMDeg.toFixed(1)}°)</small></div>
                <small>Restabweichung: ${fixed.errorDeg.toFixed(1)}° &nbsp;|&nbsp; Easy-Plane bei ${tex(`\\phi_M = ${fixed.easyCrossingPhiM1.toFixed(1)}^\\circ`)}</small>
            </div>
            <button type="button" class="btn btn-xs btn-secondary" data-apply-mount="${fixed.mountMode}" data-psi="${fixed.psiSDeg}" data-phi="${fixed.bestPhiMDeg}">
                Magnet setzen
            </button>
        </div>
        <div class="solver-card">
            <div class="solver-card-info">
                <strong>3. Horizontal-Halterung (IP – Beste In-Plane-Näherung)</strong>
                <div>${tex(`\\psi_S = ${ip1.psiSDeg.toFixed(2)}^\\circ`)} &nbsp;|&nbsp; ${tex(`\\phi_M = ${ip1.phiMDeg.toFixed(2)}^\\circ`)}</div>
                <small>${ip1.errorDeg < 0.05 ? 'Exakt in der Sample-Ebene erreichbar' : `Ziel liegt schräg zum Sample: min. Abweichung = ${ip1.errorDeg.toFixed(1)}°`}</small>
            </div>
            <button type="button" class="btn btn-xs btn-ghost" data-apply-mount="IP" data-psi="${ip1.psiSDeg}" data-phi="${ip1.phiMDeg}">
                Anwenden
            </button>
        </div>
    `;

    container.querySelectorAll('button[data-apply-mount]').forEach((btn) => {
        btn.addEventListener('click', () => {
            const s = getActiveSetup();
            s.mountMode = btn.dataset.applyMount;
            s.psiSDeg = Number(parseFloat(btn.dataset.psi).toFixed(2));
            s.phiMDeg = Number(parseFloat(btn.dataset.phi).toFixed(2));
            syncControlsFromSetup();
            recomputeAndRender();
        });
    });
}

// ============================================================================
// Export dialog: SVG · Python · Rotation matrix
// ============================================================================

let activeExportTab = 'svg';

function refreshSvgPreview() {
    const st = evaluateSetup(getActiveSetup());
    const svgStr = generateThesisSVG(st, {
        layout: document.getElementById('svgOptLayout').value,
        lockFrame: document.getElementById('svgOptFrame').value,
        cameraAzimuthDeg: sceneController.camera.azimuthDeg,
        cameraElevationDeg: sceneController.camera.elevationDeg,
        showEasyPlane: document.getElementById('svgOptEasyPlane').checked,
        showMicrostructure: document.getElementById('svgOptMicro').checked,
        showCrystalAxes: document.getElementById('svgOptCrystalAxes').checked,
        showEquilibrium: document.getElementById('svgOptEquilibrium').checked,
        showAngleArcs: document.getElementById('svgOptArcs').checked,
        showLegend: document.getElementById('svgOptLegend').checked,
        latexLabels: document.getElementById('svgOptLatex').checked,
        transparentBg: document.getElementById('svgOptTransparent').checked
    });
    document.getElementById('svgPreviewContainer').innerHTML = svgStr;
    return svgStr;
}

const FRAME_SHORT_DE = { Sub: 'Probe', Wave: 'Welle', Cryst: 'Kristall', Lab: 'Labor', Mount: 'Halter' };
const MATRIX_AXIS_TEX = {
    Sub: ['x_s', 'y_s', 'z_s'],
    Wave: ['\\hat{k}', '\\hat{t}', '\\hat{n}'],
    Cryst: ['\\hat{e}_a', '\\hat{e}_m', '\\hat{e}_c'],
    Lab: ['X_L', 'Y_L', 'Z_L'],
    Mount: ['x_m', 'y_m', 'z_m']
};

function updateMatrixInspector() {
    const st = evaluateSetup(getActiveSetup());
    const fromF = document.getElementById('selectMatFrom').value;
    const toF = document.getElementById('selectMatTo').value;
    const R = st.matrices.getRotationMatrix(toF, fromF);

    const cells = ['<div class="matrix-cell matrix-head"></div>'];
    MATRIX_AXIS_TEX[fromF].forEach((ax) => cells.push(`<div class="matrix-cell matrix-head">${tex(ax)}</div>`));
    for (let i = 0; i < 3; i++) {
        cells.push(`<div class="matrix-cell matrix-head">${tex(MATRIX_AXIS_TEX[toF][i])}</div>`);
        for (let j = 0; j < 3; j++) {
            const val = Math.abs(R[i][j]) < 1e-9 ? 0 : R[i][j];
            cells.push(`<div class="matrix-cell">${(val >= 0 ? '+' : '') + val.toFixed(4)}</div>`);
        }
    }
    document.getElementById('matrixDisplayGrid').innerHTML = cells.join('');
}

function refreshExportPanel() {
    if (activeExportTab === 'svg') {
        refreshSvgPreview();
    } else if (activeExportTab === 'python') {
        document.getElementById('pythonSnippetCode').textContent = generatePythonSnippet(
            appState.setupA,
            appState.setupB,
            appState.compareActive
        );
    } else {
        updateMatrixInspector();
    }
}

function openExportTab(tab) {
    activeExportTab = tab;
    document.querySelectorAll('#exportTabs .export-tab').forEach((t) => t.classList.toggle('active', t.dataset.etab === tab));
    document.querySelectorAll('.export-panel[data-epanel]').forEach((p) => p.classList.toggle('hidden', p.dataset.epanel !== tab));
    refreshExportPanel();
}

function flashButton(btn, text) {
    const old = btn.innerHTML;
    btn.textContent = text;
    setTimeout(() => (btn.innerHTML = old), 1400);
}

function bindExportModal() {
    const modal = document.getElementById('modalExport');

    document.querySelectorAll('[data-open-export]').forEach((btn) => {
        btn.addEventListener('click', () => {
            document.getElementById('svgOptFrame').value = sceneController.lockFrame;
            modal.classList.remove('hidden');
            openExportTab(btn.dataset.openExport);
        });
    });
    document.querySelectorAll('#exportTabs .export-tab').forEach((t) => {
        t.addEventListener('click', () => openExportTab(t.dataset.etab));
    });
    document.getElementById('btnCloseExportModal').addEventListener('click', () => modal.classList.add('hidden'));
    modal.addEventListener('click', (e) => {
        if (e.target === modal) modal.classList.add('hidden');
    });

    [
        'svgOptLayout',
        'svgOptFrame',
        'svgOptEasyPlane',
        'svgOptMicro',
        'svgOptCrystalAxes',
        'svgOptEquilibrium',
        'svgOptArcs',
        'svgOptLegend',
        'svgOptLatex',
        'svgOptTransparent'
    ].forEach((id) => document.getElementById(id).addEventListener('change', refreshSvgPreview));

    document.getElementById('btnDownloadSvg').addEventListener('click', () => {
        const svgStr = refreshSvgPreview();
        const s = getActiveSetup();
        const blob = new Blob([svgStr], { type: 'image/svg+xml;charset=utf-8' });
        const url = URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = `hematite_${s.cutPreset}_${s.mountMode}_${s.deviceType}_phiM${Math.round(s.phiMDeg)}.svg`;
        document.body.appendChild(a);
        a.click();
        a.remove();
        URL.revokeObjectURL(url);
    });

    document.getElementById('btnCopySvgCode').addEventListener('click', (e) => {
        navigator.clipboard.writeText(refreshSvgPreview());
        flashButton(e.currentTarget, 'SVG-Code kopiert');
    });

    document.getElementById('btnCopyPythonSnippet').addEventListener('click', (e) => {
        navigator.clipboard.writeText(document.getElementById('pythonSnippetCode').textContent);
        flashButton(e.currentTarget, 'Python-Code kopiert');
    });

    document.getElementById('selectMatFrom').addEventListener('change', updateMatrixInspector);
    document.getElementById('selectMatTo').addEventListener('change', updateMatrixInspector);
    document.getElementById('btnCopyMatrixLatex').addEventListener('click', (e) => {
        const st = evaluateSetup(getActiveSetup());
        const fromF = document.getElementById('selectMatFrom').value;
        const toF = document.getElementById('selectMatTo').value;
        const R = st.matrices.getRotationMatrix(toF, fromF);
        const rows = R.map((r) => r.map((v) => v.toFixed(4)).join(' & ')).join(' \\\\\n  ');
        const texStr = `R_{\\text{${FRAME_SHORT_DE[toF]}}\\leftarrow\\text{${FRAME_SHORT_DE[fromF]}}} = \\begin{pmatrix}\n  ${rows}\n\\end{pmatrix}`;
        navigator.clipboard.writeText(texStr);
        flashButton(e.currentTarget, 'Kopiert');
    });
}
