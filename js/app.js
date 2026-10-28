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
    RAD2DEG,
    DEG2RAD,
    createDefaultSetup,
    evaluateSetup,
    solveTargetAngles,
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
import { generateFigureSVG, BUILTIN_FIGURE_STYLES, DEFAULT_FIGURE_OPTIONS, completeFigureStyle } from './svg_exporter.js';
import { buildLabNote, buildSweepCsv, downloadBlob, downloadText, exportBaseName } from './export_tools.js';
import { renderSweepPlot } from './sweep_plot.js';
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
    remountFor: null, // direction whose remount hint is open below the φ_M slider
    customTargetUVTW: [1, 1, -2, 0],
    vectorTableMode: 'angles', // 'angles' | 'xyz'
    vectorTableAll: false, // false: only B, m, l, k (+ custom)
    customDir: { mode: 'dir' }, // 'dir' | 'plane' | 'cross' | 'angles'
    pinnedDirs: [] // [{ id, cryst, label, latex, uvtw, exact, color }] – fixed in the crystal frame
};

// Arrow colours for pinned directions (distinct from the built-in vectors)
const PIN_COLORS = ['#84cc16', '#06b6d4', '#ec4899', '#ef4444', '#14b8a6', '#a8a29e'];

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

// Which 3D objects a slider is "about" (highlighted while dragging)
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
    bindUIEvents();
    syncControlsFromSetup();
    syncLegend();
    syncFramePills();
    recomputeAndRender();
    applyUrlUiState();

    // Pasting a shared link into an already open tab only changes the hash – reload to apply it
    // (our own updates use history.replaceState, which does not fire hashchange)
    window.addEventListener('hashchange', () => window.location.reload());

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
    const el = target.closest('[data-emph], .story-step-title');
    if (!el) return null;
    if (el.dataset.emph) return el.dataset.emph.split(',');
    const group = el.closest('[data-emph-group]');
    return group ? group.dataset.emphGroup.split(',') : null;
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
    // Setup tabs (Probe · Struktur · Einbau); the magnet dock below them is always visible
    document.querySelectorAll('#sideTabs .side-tab').forEach((tab) => {
        tab.addEventListener('click', () => showSideTab(tab.dataset.sideTab));
    });
    document.getElementById('sideTabs').addEventListener('keydown', (e) => {
        if (e.key !== 'ArrowLeft' && e.key !== 'ArrowRight') return;
        const tabs = [...document.querySelectorAll('#sideTabs .side-tab')];
        const i = tabs.findIndex((t) => t.classList.contains('active'));
        const next = tabs[(i + (e.key === 'ArrowRight' ? 1 : tabs.length - 1)) % tabs.length];
        showSideTab(next.dataset.sideTab);
        next.focus();
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

    // Magnet dock: a direction button (or its tick on the slider) turns the magnet there; hovering previews the new field
    const phiMarks = document.getElementById('phiMarks');
    phiMarks.addEventListener('click', (e) => {
        const mark = e.target.closest('[data-phi]');
        if (!mark) return;
        getActiveSetup().phiMDeg = Number(Number(mark.dataset.phi).toFixed(2));
        setFieldPreview([]);
        syncControlsFromSetup();
        recomputeAndRender();
    });
    phiMarks.addEventListener('pointerover', (e) => {
        const mark = e.target.closest('[data-phi]');
        if (mark) previewPhiMark(mark);
        else setFieldPreview([]);
    });
    phiMarks.addEventListener('pointerleave', () => setFieldPreview([]));
    phiMarks.addEventListener('focusin', (e) => {
        const mark = e.target.closest('[data-phi]');
        if (mark) previewPhiMark(mark);
    });
    phiMarks.addEventListener('focusout', () => setFieldPreview([]));
    // Tick angle lines depend on the dock width and on the loaded font
    window.addEventListener('resize', layoutPhiTicks);
    document.fonts?.ready.then(layoutPhiTicks);

    // Directions out of reach: a click opens (or closes) the mount that reaches them
    document.getElementById('phiOff').addEventListener('click', (e) => {
        const item = e.target.closest('[data-off-id]');
        if (item) {
            appState.remountFor = appState.remountFor === item.dataset.offId ? null : item.dataset.offId;
            recomputeAndRender();
            return;
        }
        const btn = e.target.closest('button[data-apply-mount]');
        if (!btn) return;
        const s = getActiveSetup();
        s.mountMode = btn.dataset.applyMount;
        s.psiSDeg = Number(parseFloat(btn.dataset.psi).toFixed(2));
        s.phiMDeg = Number(parseFloat(btn.dataset.phi).toFixed(2));
        syncControlsFromSetup();
        recomputeAndRender();
    });

    // Frame Lock Pills (View 1): animated transition so one can watch the frames rotate
    document.querySelectorAll('#frameLockPills .frame-pill').forEach((pill) => {
        pill.addEventListener('click', () => {
            sceneController.setLockFrame(pill.dataset.frame, true);
            syncFramePills();
            writeStateToUrl();
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
    bindExportCard();
    bindTour();

    document.addEventListener('keydown', (e) => {
        if (e.key !== 'Escape') return;
        if (tour && tour.active) tour.stop();
        document.getElementById('modalExport').classList.add('hidden');
    });
}

// ============================================================================
// Setup tabs & field-direction helpers
// ============================================================================

function showSideTab(name) {
    document.querySelectorAll('#sideTabs .side-tab').forEach((t) => {
        const on = t.dataset.sideTab === name;
        t.classList.toggle('active', on);
        t.setAttribute('aria-selected', String(on));
        t.tabIndex = on ? 0 : -1;
    });
    document.querySelectorAll('#sidePanels .side-panel').forEach((p) => {
        p.hidden = p.dataset.sidePanel !== name;
    });
    document.getElementById('sidePanels').scrollTop = 0;
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
            refreshCustomDirection();
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
        el.addEventListener('input', refreshCustomDirection);
        el.addEventListener('change', refreshCustomDirection);
    });

    // Opening the card draws the arrow, closing it removes it again
    document.getElementById('customDirCard').addEventListener('toggle', (e) => {
        sceneController.visibility.showCustom = e.target.open;
        syncLegend();
        if (lastActiveState) updateVectorTable(lastActiveState);
        rerenderScenes();
    });

    document.getElementById('btnCustomDirDrive').addEventListener('click', driveFieldToCustomDirection);

    // Pinned directions: "+ Merken" adds the current one, × removes, clicking the index loads it back
    document.getElementById('btnCustomDirPin').addEventListener('click', pinCurrentDirection);
    document.getElementById('pinnedDirList').addEventListener('click', (e) => {
        const item = e.target.closest('[data-pin-id]');
        if (!item) return;
        const id = Number(item.dataset.pinId);
        if (e.target.closest('.pin-remove')) {
            appState.pinnedDirs = appState.pinnedDirs.filter((p) => p.id !== id);
            setEmphasis(null);
        } else if (e.target.closest('.pin-load')) {
            loadPinnedDirection(appState.pinnedDirs.find((p) => p.id === id));
        }
        refreshCustomDirection();
    });
    document.getElementById('btnPinClear').addEventListener('click', () => {
        appState.pinnedDirs = [];
        refreshCustomDirection();
    });
}

// Recompute the custom + pinned directions and redraw everything that shows them
function refreshCustomDirection() {
    if (!lastActiveState) return;
    updateCustomDirection(lastActiveState);
    updatePinnedDirections(lastActiveState);
    updateVectorTable(lastActiveState);
    rerenderScenes();
    writeStateToUrl();
}

function pinCurrentDirection() {
    const cd = appState.customDir;
    if (!cd.cryst || findPinned(cd.cryst)) return;
    const used = new Set(appState.pinnedDirs.map((p) => p.color));
    const color = PIN_COLORS.find((c) => !used.has(c)) || PIN_COLORS[appState.pinnedDirs.length % PIN_COLORS.length];
    const id = (appState.pinnedDirs.reduce((m, p) => Math.max(m, p.id), 0) || 0) + 1;
    appState.pinnedDirs.push({ id, cryst: [...cd.cryst], label: cd.text, latex: cd.latex, uvtw: cd.uvtw, exact: cd.exact, color });
    refreshCustomDirection();
}

function findPinned(vCryst) {
    return appState.pinnedDirs.find((p) => dot3(p.cryst, vCryst) > 1 - 1e-6);
}

// Puts a pinned direction back into the [uvtw] input
function loadPinnedDirection(pin) {
    if (!pin) return;
    document.querySelector('#customDirModes .cd-mode[data-cdmode="dir"]').click();
    document.getElementById('cdH').value = pin.uvtw[0];
    document.getElementById('cdK').value = pin.uvtw[1];
    document.getElementById('cdL').value = pin.uvtw[3];
}

function updatePinnedDirections(state) {
    const pins = appState.pinnedDirs;
    const R = state.matrices.getRotationMatrix('Lab', 'Cryst');
    const vectors = pins.map((p) => ({ Lab: matVec3(R, p.cryst), label: p.label, color: p.color, key: `pin${p.id}` }));
    sceneController.pinnedVectors = vectors;
    if (sceneController2) sceneController2.pinnedVectors = vectors;

    const cur = appState.customDir.cryst;
    const btn = document.getElementById('btnCustomDirPin');
    const already = cur && findPinned(cur);
    btn.disabled = !cur || Boolean(already) || pins.length >= PIN_COLORS.length;
    btn.textContent = already ? '✓ gemerkt' : pins.length >= PIN_COLORS.length ? 'Liste voll' : '+ Merken';

    document.getElementById('pinnedDirs').hidden = pins.length === 0;
    document.getElementById('btnPinClear').hidden = pins.length < 2;
    document.getElementById('customDirPinCount').textContent = pins.length ? ` · ${pins.length} gemerkt` : '';

    const v = state.vectors;
    const RSub = state.matrices.getRotationMatrix('Sub', 'Cryst');
    document.getElementById('pinnedDirList').innerHTML = pins
        .map((p) => {
            const vSub = matVec3(RSub, p.cryst);
            const active = cur && dot3(p.cryst, cur) > 1 - 1e-6;
            return `
            <li class="pin-item${active ? ' active' : ''}" data-pin-id="${p.id}" data-emph="pin${p.id}">
                <span class="vt-dot" style="background:${p.color}"></span>
                <button type="button" class="pin-load" title="wieder in die Eingabe laden">${tex(p.latex)}</button>
                <span class="pin-class">${shortDirection(p.cryst)}</span>
                <span class="pin-ang" title="Winkel zu k̂ und zu B">${tex('\\hat{k}')} ${angleDeg(vSub, v.kHat.Sub).toFixed(0)}° · ${tex('\\mathbf{B}')} ${angleDeg(vSub, v.bHat.Sub).toFixed(0)}°</span>
                <button type="button" class="pin-remove" title="aus der Liste entfernen" aria-label="entfernen">×</button>
            </li>`;
        })
        .join('');
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
    appState.customDir.cryst = vCryst;
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
    appState.customDir.text = dir.text;
    appState.customDir.exact = dir.exact;

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
    const fixed = solveTargetAngles(getActiveSetup(), 'custom_mb', uvtw).fixedPsiSweep;
    if (fixed.errorDeg <= 89.9) getActiveSetup().phiMDeg = Number(fixed.bestPhiMDeg.toFixed(2));
    // Not exactly reachable: the magnet dock offers the mount that is
    appState.remountFor = fixed.errorDeg < 0.05 ? null : 'custom_mb';
    syncControlsFromSetup();
    recomputeAndRender();
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
    if (appState.compareActive && appState.activeTab === 'B') params.set('tab', 'B');
    params.set('view', sceneController ? sceneController.lockFrame : 'Sub');
    let hash = params.toString();
    // Index lists are appended by hand so they stay readable (URLSearchParams would escape , and ;)
    if (appState.pinnedDirs.length) hash += `&pins=${appState.pinnedDirs.map((p) => p.uvtw.join(',')).join(';')}`;
    history.replaceState(null, '', `#${hash}`);
}

const URL_ALLOWED = {
    cutPreset: (v) => Object.prototype.hasOwnProperty.call(CUT_PRESETS, v),
    deviceType: (v) => v === 'CPW' || v === 'IDT',
    mountMode: (v) => v === 'IP' || v === 'OOP'
};

// Parses "1,0,-1,1" into [u, v, t, w] (t recomputed), or null if malformed / all zero
function parseUvtw(raw) {
    const n = (raw || '').split(',').map(Number);
    if (n.length !== 4 || !n.every(Number.isFinite) || n.every((x) => x === 0)) return null;
    return [n[0], n[1], -(n[0] + n[1]), n[3]];
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
                if (!URL_ALLOWED[prop](raw)) continue;
                setup[prop] = raw;
            } else if (Number.isFinite(parseFloat(raw))) {
                setup[prop] = parseFloat(raw);
            } else {
                continue;
            }
            found = true;
        }
        if (found) setup.name = prefix ? 'Setup B (aus Link)' : 'Setup A (aus Link)';
    };
    apply('', appState.setupA);
    apply('b.', appState.setupB);
    appState.compareActive = params.get('cmp') === '1';
    appState.activeTab = appState.compareActive && params.get('tab') === 'B' ? 'B' : 'A';
    const view = params.get('view');
    if (view && FRAME_METADATA[view] && view !== 'Mount') sceneController.lockFrame = view;

    // Pinned directions are stored as [uvtw]; they live in the crystal frame, so A's lattice is enough
    const { aLattice: a, cLattice: c } = appState.setupA;
    appState.pinnedDirs = (params.get('pins') || '')
        .split(';')
        .map(parseUvtw)
        .filter(Boolean)
        .slice(0, PIN_COLORS.length)
        .map((uvtw, i) => {
            const cryst = normalize3(millerBravaisDirectionToCrystal(uvtw, a, c));
            const dir = bestMillerBravais(cryst, a, c, false);
            return { id: i + 1, cryst, label: dir.text, latex: dir.latex, uvtw: dir.indices, exact: dir.exact, color: PIN_COLORS[i] };
        });
}

// UI that only exists after bindUIEvents(): the pinned-directions card
function applyUrlUiState() {
    if (appState.pinnedDirs.length) document.getElementById('customDirCard').open = true;
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
    updateFieldMarks(activeState);

    // The custom direction depends on the setup (it is drawn in Lab coordinates)
    updateCustomDirection(activeState);
    updatePinnedDirections(activeState);

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
    document.getElementById('exportCardHint').textContent = appState.compareActive ? `Dateien für Setup ${appState.activeTab}` : '';

    writeStateToUrl();
}

function updateStepSummaries(activeState) {
    // One short value per tab, so the tab bar doubles as a summary of the setup
    const s = activeState.setup;
    document.getElementById('summaryStep1').textContent = s.cutPreset === 'custom' ? 'eigener Schnitt' : s.cutPreset;
    document.getElementById('summaryStepDev').textContent = `${s.deviceType} ${fmtDeg(s.thetaDevDeg)}`;
    document.getElementById('summaryStep2').textContent = `${s.mountMode} ${fmtDeg(s.psiSDeg)}`;
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
    const geo = m.geometry;
    valMDotK.innerHTML = geo.key === 'none' ? '' : `<strong>${geo.label}</strong>`;
    valMDotK.title = geo.key === 'none'
        ? ''
        : 'Spinwellen-Geometrie: BV = m̂ ∥ k̂, DE = m̂ ⊥ k̂ (beide in der Probe), FV = m̂ ⊥ Probe'
          + (isCPW
              ? `\nDrehmoment des CPW-Felds auf m̂: h_x ∥ k̂ → ${(m.cpwTorqueM * 100).toFixed(0)} %, h_z ⊥ Probe → ${(m.cpwTorqueMz * 100).toFixed(0)} %`
              : '');

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
        ['Geometrie', mA.geometry.label, mB.geometry.label, '—']
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
    if (sceneController.visibility.showCustom) {
        const RC = activeState.matrices.getRotationMatrix(frame, 'Cryst');
        const cur = appState.customDir.cryst;
        if (cur && !findPinned(cur)) {
            entries.push({ key: 'custom', name: 'Eigene', vecInFrame: matVec3(RC, cur), mb: appState.customDir.latex });
        }
        for (const p of appState.pinnedDirs) {
            entries.push({ key: `pin${p.id}`, name: 'Gemerkt', vecInFrame: matVec3(RC, p.cryst), mb: p.latex, color: p.color, pinned: true });
        }
    }

    const fmtComp = (x) => {
        const val = Math.abs(x) < 5e-5 ? 0 : x;
        return (val >= 0 ? '+' : '−') + Math.abs(val).toFixed(3);
    };
    const color = (key) => (key === 'xs' || key === 'ys' ? VECTOR_COLORS.edge : VECTOR_COLORS[key]);

    const shown = entries.filter((e) => !e.skip);
    const visible = appState.vectorTableAll ? shown : shown.filter((e) => e.pinned || TABLE_CORE_KEYS.includes(e.key));
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
            return `<tr data-emph="${e.key}"><td><span class="vt-dot" style="background:${e.color || color(e.key)}"></span>${e.name}</td>${cells}<td>${tex(e.mb)}</td></tr>`;
        })
        .join('');
}

// ============================================================================
// Magnet dock: marks on the φ_M slider
// ============================================================================

const fmtDeg = (deg) => `${Number(deg.toFixed(1))}°`;
const phiDistDeg = (a, b) => Math.abs((((a - b) % 360) + 540) % 360 - 180);
const MOUNT_WORDS = { IP: 'flach', OOP: 'stehend' };
const uvtwTex = (uvtw) => uvtw.map((v) => (v < 0 ? `\\bar{${-v}}` : `${v}`)).join('\\,');

// Field directions offered on the slider: along the wave / sample frame or along a crystal axis.
// Plain words on purpose: single symbols (t̂, n̂ …) are hard to tell apart at a glance.
const FIELD_MARKS = [
    { id: 'B_parallel_k', key: 'k', text: 'längs zur Welle', title: 'Feld in Ausbreitungsrichtung der Welle (B ∥ +k̂)' },
    { id: 'B_antiparallel_k', key: 'k', text: 'gegen die Welle', title: 'Feld entgegen der Ausbreitungsrichtung (B ∥ −k̂)' },
    { id: 'B_parallel_t', key: 't', text: 'quer zur Welle', title: 'Feld quer zur Welle, in der Probenebene (B ∥ t̂)' },
    { id: 'B_parallel_n', key: 'n', text: 'senkrecht zur Probe', title: 'Feld senkrecht zur Probenoberfläche (B ∥ n̂)' },
    { id: 'B_parallel_ea', key: 'ea', text: 'a-Achse', title: 'Feld entlang der a-Achse [112̅0]' },
    { id: 'B_parallel_em', key: 'em', text: 'm-Achse', title: 'Feld entlang der m-Achse [11̅00]' },
    { id: 'B_parallel_ec', key: 'ec', text: 'c-Achse', title: 'Feld entlang der c-Achse [0001]' }
];
// A direction from the "Eigene Richtung" card, only while its remount hint is open
const customFieldMark = () => ({ id: 'custom_mb', key: 'custom', tex: `[${uvtwTex(appState.customTargetUVTW)}]`, title: 'eigene Richtung' });

const markLabel = (m) => `<span class="phi-mark-label" style="--tok:${VECTOR_COLORS[m.key]}">${m.tex ? tex(m.tex) : m.text}</span>`;
let phiMarksSignature = '';

// Exactly reachable directions become marks at their φ_M; the others are listed below the slider,
// where a click shows the mount that reaches them.
function updateFieldMarks(activeState) {
    const setup = activeState.setup;
    const marks = appState.remountFor === 'custom_mb' ? [...FIELD_MARKS, customFieldMark()] : FIELD_MARKS;
    const reach = [];
    const off = [];
    for (const m of marks) {
        const sol = solveTargetAngles(setup, m.id, appState.customTargetUVTW);
        const exact = sol.fixedPsiSweep.errorDeg < 0.05 && m.id !== 'custom_mb';
        (exact ? reach : off).push({ ...m, sol, phi: ((sol.fixedPsiSweep.bestPhiMDeg % 360) + 360) % 360 });
    }
    renderPhiMarks(reach, setup.phiMDeg);
    renderPhiOff(off, setup);
}

function renderPhiMarks(reach, phiNow) {
    // Directions at the same φ_M share one tick and one button
    const groups = [];
    for (const m of [...reach].sort((a, b) => a.phi - b.phi)) {
        const g = groups.find((x) => phiDistDeg(x.phi, m.phi) < 1);
        if (g) g.items.push(m);
        else groups.push({ phi: m.phi, items: [m] });
    }

    const box = document.getElementById('phiMarks');
    const signature = JSON.stringify(groups.map((g) => [g.phi.toFixed(2), g.items.map((m) => m.id)]));
    if (signature !== phiMarksSignature) {
        phiMarksSignature = signature;
        const attrs = (g) => `data-phi="${g.phi}" data-emph="preview,${g.items.map((m) => m.key).join(',')}" title="${g.items.map((m) => m.title).join('\n')}"`;
        box.innerHTML = `
            <div class="phi-ticks">${groups.map((g) => `<span class="phi-tick" style="left:${(g.phi / 360) * 100}%" ${attrs(g)}><span class="phi-tick-deg">${fmtDeg(g.phi)}</span></span>`).join('')}</div>
            <div class="phi-pills">${groups.map((g) => `<button type="button" class="phi-pill" ${attrs(g)}>${g.items.map(markLabel).join('<span class="phi-mark-sep">·</span>')}<span class="phi-pill-deg">${fmtDeg(g.phi)}</span></button>`).join('')}</div>`;
        layoutPhiTicks();
    }
    box.querySelectorAll('[data-phi]').forEach((el) => el.classList.toggle('on', phiDistDeg(phiNow, Number(el.dataset.phi)) < 0.05));
}

// Tick angles sit centred under their tick; one that would touch its left neighbour moves down a line
function layoutPhiTicks() {
    const ticks = document.querySelector('#phiMarks .phi-ticks');
    if (!ticks?.clientWidth) return;
    const rowEnds = [];
    ticks.querySelectorAll('.phi-tick-deg').forEach((deg) => {
        const b = deg.getBoundingClientRect();
        let row = rowEnds.findIndex((end) => b.left >= end + 6);
        if (row < 0) row = rowEnds.length;
        rowEnds[row] = b.right;
        deg.parentElement.style.setProperty('--row', row);
    });
    ticks.style.setProperty('--rows', Math.max(rowEnds.length, 1));
}

function renderPhiOff(off, setup) {
    const box = document.getElementById('phiOff');
    if (!off.some((m) => m.id === appState.remountFor)) appState.remountFor = null;
    box.hidden = !off.length;
    if (!off.length) {
        box.innerHTML = '';
        return;
    }
    const open = off.find((m) => m.id === appState.remountFor);
    let remount = '';
    if (open) {
        // In-plane directions are reached exactly by the flat (IP) mount, everything else by the vertical (OOP) one
        const ip = open.sol.ipSolutions[0];
        const alt = setup.mountMode === 'OOP' && ip.errorDeg < 0.05 ? ip : open.sol.oopSolutions[0];
        remount = alt.errorDeg < 0.05
            ? `<div class="phi-remount"><span>${markLabel(open)} geht mit Probe <b>${MOUNT_WORDS[alt.mountMode]}</b>, auf dem Halter um <b>${fmtDeg(alt.psiSDeg)}</b> gedreht</span>
                <button type="button" class="btn btn-xs btn-secondary" data-apply-mount="${alt.mountMode}" data-psi="${alt.psiSDeg}" data-phi="${alt.phiMDeg}" title="Probe so einbauen (Tab Einbau) und das Feld dorthin drehen">So einbauen</button></div>`
            : `<div class="phi-remount"><span>${markLabel(open)} geht mit keinem Einbau genau</span></div>`;
    }
    box.innerHTML = `
        <div class="phi-off-row"><span class="phi-off-label">Nur mit anderem Einbau:</span>
        ${off.map((m) => `<button type="button" class="phi-off-item${m === open ? ' active' : ''}" data-off-id="${m.id}" data-emph="${m.key}" title="${m.title}\nKlick zeigt, wie die Probe dafür eingebaut sein muss" aria-expanded="${m === open}">${markLabel(m)}<span class="phi-off-caret">${m === open ? '▾' : '▸'}</span></button>`).join('')}</div>
        ${remount}`;
}

// Hovering a direction: dashed "B neu" arrow where the field would point, its tick lights up
function previewPhiMark(mark) {
    const phi = Number(mark.dataset.phi);
    document.querySelectorAll('#phiMarks [data-phi]').forEach((el) => el.classList.toggle('hot', el.dataset.phi === mark.dataset.phi));
    const setup = getActiveSetup();
    if (phiDistDeg(setup.phiMDeg, phi) < 0.05) return setFieldPreview([]);
    const st = evaluateSetup({ ...setup, phiMDeg: phi });
    setFieldPreview([{ Lab: st.vectors.bHat.Lab, colorKey: 'B', label: 'B neu' }]);
}

function setFieldPreview(vectors) {
    if (!vectors.length) document.querySelectorAll('#phiMarks .hot').forEach((el) => el.classList.remove('hot'));
    sceneController.setFieldPreview(vectors);
    if (sceneController2) sceneController2.setFieldPreview(vectors);
}

// ============================================================================
// Export dialog: PNG · SVG · field sweep · Python · rotation matrix
// ============================================================================

let activeExportTab = 'svg';

// Figure (SVG): options, the style being edited and the user's saved styles live in localStorage,
// so every figure of a thesis comes out the same without setting it up again
const FIGURE_STORAGE_KEY = 'hematite_figure_v1';
const figure = { options: { ...DEFAULT_FIGURE_OPTIONS }, style: completeFigureStyle(BUILTIN_FIGURE_STYLES.Standard), styleName: 'Standard', saved: {} };

function loadFigureSettings() {
    try {
        const s = JSON.parse(localStorage.getItem(FIGURE_STORAGE_KEY) || 'null');
        if (!s) return;
        figure.options = { ...DEFAULT_FIGURE_OPTIONS, ...(s.options || {}) };
        figure.saved = s.saved && typeof s.saved === 'object' ? s.saved : {};
        figure.style = completeFigureStyle(s.style);
        figure.styleName = s.styleName || '';
    } catch (_) {
        // Keep the defaults
    }
}

function saveFigureSettings() {
    try {
        localStorage.setItem(FIGURE_STORAGE_KEY, JSON.stringify(figure));
    } catch (_) {
        // Ignore storage quota errors
    }
}

const figStylePath = (path) => path.split('.');
function getFigStyleValue(path) {
    return figStylePath(path).reduce((o, k) => o?.[k], figure.style);
}
function setFigStyleValue(path, value) {
    const keys = figStylePath(path);
    const last = keys.pop();
    keys.reduce((o, k) => o[k], figure.style)[last] = value;
}

function renderFigStyleSelect() {
    const sel = document.getElementById('figStyleSelect');
    // Built with DOM nodes: saved style names are user text
    const group = (label, names, prefix) => {
        const g = document.createElement('optgroup');
        g.label = label;
        names.forEach((n) => g.append(new Option(n, `${prefix}:${n}`)));
        return g;
    };
    const own = Object.keys(figure.saved);
    sel.replaceChildren(group('Mitgeliefert', Object.keys(BUILTIN_FIGURE_STYLES), 'b'));
    if (own.length) sel.append(group('Eigene', own, 'u'));
    if (!figure.styleName) {
        const changed = new Option('(geändert)', '');
        changed.disabled = true;
        sel.append(changed);
    }
    sel.value = figure.styleName ? `${figure.saved[figure.styleName] ? 'u' : 'b'}:${figure.styleName}` : '';
    document.getElementById('btnFigStyleDelete').disabled = !figure.saved[figure.styleName];
}

function syncFigureControls() {
    document.querySelectorAll('[data-fopt]').forEach((el) => {
        if (el.type === 'checkbox') el.checked = !!figure.options[el.dataset.fopt];
        else el.value = figure.options[el.dataset.fopt];
    });
    document.querySelectorAll('[data-fs]').forEach((el) => {
        const v = getFigStyleValue(el.dataset.fs);
        if (el.type === 'color') el.value = /^#[0-9a-f]{6}$/i.test(v) ? v : '#ffffff';
        else el.value = v;
    });
    const bgNone = getFigStyleValue('colors.background') === 'none';
    document.getElementById('figBgNone').checked = bgNone;
    document.querySelector('[data-fs="colors.background"]').disabled = bgNone;
    renderFigStyleSelect();
}

function refreshSvgPreview() {
    const svgStr = generateFigureSVG(evaluateSetup(getActiveSetup()), figure.options, figure.style);
    document.getElementById('svgPreviewContainer').innerHTML = svgStr;
    return svgStr;
}

function applyFigureStyle(style, name) {
    figure.style = completeFigureStyle(JSON.parse(JSON.stringify(style)));
    figure.styleName = name;
    syncFigureControls();
    saveFigureSettings();
    refreshSvgPreview();
}

function bindFigureControls() {
    loadFigureSettings();
    syncFigureControls();

    document.querySelectorAll('[data-fopt]').forEach((el) =>
        el.addEventListener('change', () => {
            figure.options[el.dataset.fopt] = el.type === 'checkbox' ? el.checked : el.value;
            saveFigureSettings();
            refreshSvgPreview();
        })
    );

    // Any edit turns the style into an unsaved one; the dropdown then shows "(geändert)"
    const edited = () => {
        figure.styleName = '';
        renderFigStyleSelect();
        saveFigureSettings();
        refreshSvgPreview();
    };
    document.querySelectorAll('[data-fs]').forEach((el) =>
        el.addEventListener('input', () => {
            if (el.type === 'number') {
                if (el.value === '' || !Number.isFinite(Number(el.value))) return;
                setFigStyleValue(el.dataset.fs, Number(el.value));
            } else {
                setFigStyleValue(el.dataset.fs, el.value);
            }
            edited();
        })
    );
    document.getElementById('figBgNone').addEventListener('change', (e) => {
        const color = document.querySelector('[data-fs="colors.background"]');
        color.disabled = e.target.checked;
        setFigStyleValue('colors.background', e.target.checked ? 'none' : color.value);
        edited();
    });

    document.getElementById('figStyleSelect').addEventListener('change', (e) => {
        const [group, ...rest] = e.target.value.split(':');
        const name = rest.join(':');
        const style = group === 'u' ? figure.saved[name] : BUILTIN_FIGURE_STYLES[name];
        if (style) applyFigureStyle(style, name);
    });

    const nameInput = document.getElementById('figStyleName');
    const saveBtn = document.getElementById('btnFigStyleSave');
    const save = () => {
        const name = nameInput.value.trim();
        if (!name) {
            nameInput.focus();
            return;
        }
        if (BUILTIN_FIGURE_STYLES[name]) {
            flashButton(saveBtn, 'Name vergeben');
            return;
        }
        figure.saved[name] = JSON.parse(JSON.stringify(figure.style));
        figure.styleName = name;
        nameInput.value = '';
        renderFigStyleSelect();
        saveFigureSettings();
        flashButton(saveBtn, 'Gespeichert');
    };
    saveBtn.addEventListener('click', save);
    nameInput.addEventListener('keydown', (e) => e.key === 'Enter' && save());

    document.getElementById('btnFigStyleDelete').addEventListener('click', () => {
        if (!figure.saved[figure.styleName]) return;
        delete figure.saved[figure.styleName];
        applyFigureStyle(BUILTIN_FIGURE_STYLES.Standard, 'Standard');
    });

    document.getElementById('btnFigStyleCopy').addEventListener('click', (e) => {
        navigator.clipboard.writeText(JSON.stringify(figure.style, null, 2));
        flashButton(e.currentTarget, 'Kopiert');
    });
    document.getElementById('btnFigStylePaste').addEventListener('click', (e) => {
        const btn = e.currentTarget;
        navigator.clipboard.readText().then(
            (text) => {
                try {
                    const style = JSON.parse(text);
                    if (!style || typeof style !== 'object' || !style.colors) throw new Error('no style');
                    applyFigureStyle(style, '');
                    flashButton(btn, 'Übernommen');
                } catch (_) {
                    flashButton(btn, 'Kein Stil-JSON');
                }
            },
            () => flashButton(btn, 'Kein Zugriff')
        );
    });
}

// Logical size (CSS px) per PNG format; 'view' = the size of view 1 on screen
const PNG_FORMATS = { '16:9': [960, 540], '4:3': [800, 600], '1:1': [680, 680] };

function renderPngExport(scale = Number(document.getElementById('pngOptScale').value)) {
    const rect = sceneController.canvas3d.getBoundingClientRect();
    const [width, height] = PNG_FORMATS[document.getElementById('pngOptFormat').value] || [
        Math.round(rect.width),
        Math.round(rect.height)
    ];
    const colors = document.getElementById('pngOptColors').value;
    return sceneController.renderToCanvas({
        width,
        height,
        scale,
        light: colors === 'app' ? document.body.classList.contains('theme-light') : colors === 'light',
        transparent: document.getElementById('pngOptTransparent').checked
    });
}

function refreshPngPreview() {
    const scale = Number(document.getElementById('pngOptScale').value);
    // Preview at screen resolution; the download renders at the chosen scale
    const canvas = renderPngExport(window.devicePixelRatio || 1);
    const box = document.getElementById('pngPreviewContainer');
    box.replaceChildren(canvas);
    box.classList.toggle('is-transparent', document.getElementById('pngOptTransparent').checked);
    const w = Math.round((canvas.width / (window.devicePixelRatio || 1)) * scale);
    const h = Math.round((canvas.height / (window.devicePixelRatio || 1)) * scale);
    document.getElementById('pngSizeHint').textContent = `${w} × ${h} px`;
}

function pngBlob() {
    return new Promise((resolve) => renderPngExport().toBlob(resolve, 'image/png'));
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
    } else if (activeExportTab === 'png') {
        refreshPngPreview();
    } else if (activeExportTab === 'python') {
        document.getElementById('pythonSnippetCode').textContent = generatePythonSnippet(
            appState.setupA,
            appState.setupB,
            appState.compareActive
        );
    } else if (activeExportTab === 'sweep') {
        refreshSweepPreview();
    } else {
        updateMatrixInspector();
    }
}

function refreshSweepPreview() {
    document.getElementById('sweepSetupHint').textContent = appState.compareActive ? ` (Setup ${appState.activeTab})` : '';
    renderSweepPlot(document.getElementById('sweepPlot'), getActiveSetup(), (phi) => {
        getActiveSetup().phiMDeg = phi;
        syncControlsFromSetup();
        recomputeAndRender();
    });
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

    bindFigureControls();

    ['pngOptFormat', 'pngOptColors', 'pngOptScale', 'pngOptTransparent'].forEach((id) =>
        document.getElementById(id).addEventListener('change', refreshPngPreview)
    );

    document.getElementById('btnDownloadPng').addEventListener('click', () => {
        pngBlob().then((blob) => blob && downloadBlob(blob, `${exportBaseName(getActiveSetup())}.png`));
    });

    document.getElementById('btnCopyPng').addEventListener('click', (e) => {
        const btn = e.currentTarget;
        navigator.clipboard.write([new ClipboardItem({ 'image/png': pngBlob() })]).then(
            () => flashButton(btn, 'Bild kopiert'),
            () => flashButton(btn, 'Kopieren nicht möglich')
        );
    });

    document.getElementById('btnDownloadSvg').addEventListener('click', () => {
        downloadText(refreshSvgPreview(), `${exportBaseName(getActiveSetup())}.svg`, 'image/svg+xml');
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

// "Teilen & Exportieren" card: everyday exports that need no dialog
function bindExportCard() {
    document.getElementById('btnCopyNote').addEventListener('click', () => {
        writeStateToUrl();
        const states = appState.compareActive
            ? [evaluateSetup(appState.setupA), evaluateSetup(appState.setupB)]
            : [evaluateSetup(getActiveSetup())];
        const note = buildLabNote(states, {
            pinned: appState.pinnedDirs.map((p) => p.label),
            link: window.location.href
        });
        navigator.clipboard.writeText(note).then(
            () => showToast('Notiz kopiert – einfach ins Laborbuch einfügen'),
            () => showToast('Kopieren nicht möglich')
        );
    });

    document.getElementById('btnDownloadCsv').addEventListener('click', () => {
        const s = getActiveSetup();
        downloadText(buildSweepCsv(s), `${exportBaseName(s).replace(/_phiM\d+$/, '')}_sweep.csv`, 'text/csv');
    });
}
