/**
 * app.js
 * Main UI controller for the Hematite (alpha-Fe2O3) Spin-Wave (CPW) & SAW (IDT)
 * Coordinate Transformation Calculator.
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
    computeSampleAngleMillerBravais,
    solveTargetAngles,
    solveSweepLandmarksForFixedPsi,
    computeMountingPresetPsiS,
    hardwareAnglesFromWavePerspective,
    hardwareAnglesFromSubPerspective,
    hardwareAnglesFromCrystalPerspective,
    millerBravaisPlaneToCrystalNormal,
    matVec3,
    generatePythonSnippet,
    wrapDeg360
} from './math_core.js';

import { Scene3DController, VECTOR_COLORS } from './visualizer_3d.js';
import { SweepChartController } from './sweep_chart.js';
import { generateThesisSVG } from './svg_exporter.js';

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
    customTargetUVTW: [1, 1, -2, 0]
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

document.addEventListener('DOMContentLoaded', () => {
    const canvas3d = document.getElementById('canvas3D');
    const canvas3dSec = document.getElementById('canvas3DSecondary');

    sceneController = new Scene3DController(canvas3d, null, null);
    sceneController2 = new Scene3DController(canvas3dSec, null, null);
    // Share visibility flags between View 1 and View 2; default View 2 to Easy-Plane in Top-Down (2D) mode
    sceneController2.visibility = sceneController.visibility;
    sceneController2.lockFrame = 'Cryst';
    sceneController2.setTopDown2D(true);

    renderStaticLatex();
    loadPresetsFromStorage();
    renderPresetDropdown();
    populateTargetSelect();
    bindUIEvents();
    syncControlsFromSetup();
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

    const readoutBadge = document.getElementById('activeReadoutBadge');
    if (readoutBadge) {
        readoutBadge.textContent = s.name || `Setup ${appState.activeTab}`;
        readoutBadge.className = appState.activeTab === 'A' ? 'badge-a' : 'badge-b';
    }

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
            if (lbl) {
                lbl.textContent = enable ? '2. Ansicht aktiv' : '+ 2. Ansicht daneben';
            }
        }
        if (lblFrameLock1) {
            lblFrameLock1.textContent = enable
                ? 'Was soll in Ansicht 1 stillstehen?'
                : 'Was soll im Bild stillstehen?';
        }
        requestAnimationFrame(() => {
            sceneController.render();
            if (sceneController2) sceneController2.render();
        });
    };

    if (btnToggleDual3D) {
        btnToggleDual3D.addEventListener('click', () => {
            const isDual = stageContainer.classList.contains('mode-dual-3d');
            setDual3DActive(!isDual);
        });
    }
    if (btnCloseDual3D) {
        btnCloseDual3D.addEventListener('click', () => {
            setDual3DActive(false);
        });
    }

    if (btnCamTop1 && sceneController) {
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
            if (e.target.value) {
                applySavedPreset(e.target.value);
            }
        });
    }

    if (btnOpenSave && saveRow && inputPresetName) {
        btnOpenSave.addEventListener('click', () => {
            const s = getActiveSetup();
            const suggested = `${s.cutPreset} · ${s.mountMode} · ${s.deviceType} (${s.thetaDevDeg.toFixed(0)}°)`;
            inputPresetName.value = suggested;
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

    if (btnConfirmSave) {
        btnConfirmSave.addEventListener('click', confirmSavePreset);
    }
    if (btnCancelSave && saveRow) {
        btnCancelSave.addEventListener('click', () => {
            saveRow.classList.add('hidden');
        });
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
    if (btnDeletePreset) {
        btnDeletePreset.addEventListener('click', () => {
            deleteSelectedPreset();
        });
    }

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
        appState.setupB = {
            ...appState.setupA,
            name: `${appState.setupA.name} (Kopie B)`
        };
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
        const v = Math.max(0.001, parseFloat(e.target.value) || 0.5);
        getActiveSetup().b0Tesla = v;
        recomputeAndRender();
    });

    // Quick-Value Chips
    document.querySelectorAll('.quick-chips').forEach((group) => {
        const target = group.dataset.target;
        const propMap = {
            PhiM: 'phiMDeg',
            PsiS: 'psiSDeg',
            ThetaDev: 'thetaDevDeg'
        };
        group.querySelectorAll('button').forEach((btn) => {
            btn.addEventListener('click', () => {
                const val = parseFloat(btn.dataset.val);
                const s = getActiveSetup();
                s[propMap[target]] = val;
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
            const alpha = parseFloat(btn.dataset.alpha) || 0;
            const s = getActiveSetup();
            s.alphaFlatDeg = alpha;
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
        const h = parseFloat(hInput.value) || 0;
        const k = parseFloat(kInput.value) || 0;
        iInput.value = -(h + k);
    };
    hInput.addEventListener('input', syncI);
    kInput.addEventListener('input', syncI);

    document.getElementById('btnApplyHKIL').addEventListener('click', () => {
        const h = parseFloat(hInput.value) || 0;
        const k = parseFloat(kInput.value) || 0;
        const l = parseFloat(lInput.value) || 0;
        const nCryst = millerBravaisPlaneToCrystalNormal([h, k, -(h + k), l]);
        const thetaDeg = Math.acos(Math.max(-1, Math.min(1, nCryst[2]))) * RAD2DEG;
        const phiDeg = wrapDeg360(Math.atan2(nCryst[1], nCryst[0]) * RAD2DEG);
        const s = getActiveSetup();
        s.cutPreset = 'custom';
        s.thetaCutDeg = Number(thetaDeg.toFixed(3));
        s.phiCutDeg = Number(phiDeg.toFixed(2));
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
    const setSolverTarget = (targetId) => {
        appState.targetId = targetId;
        const selEl = document.getElementById('selectTargetPreset');
        if (selEl) selEl.value = targetId;
        document.querySelectorAll('.solver-target-card[data-target-id]').forEach((card) => {
            card.classList.toggle('active', card.dataset.targetId === targetId);
        });
        document
            .getElementById('customTargetMBBox')
            .classList.toggle('hidden', targetId !== 'custom_mb');
        recomputeAndRender();
    };

    document.querySelectorAll('#solverCatTabs .solver-cat-tab').forEach((tab) => {
        tab.addEventListener('click', () => {
            const cat = tab.dataset.solverCat;
            document
                .querySelectorAll('#solverCatTabs .solver-cat-tab')
                .forEach((t) => t.classList.toggle('active', t === tab));
            document.querySelectorAll('.solver-target-grid[data-solver-panel]').forEach((panel) => {
                const isActivePanel = panel.dataset.solverPanel === cat;
                panel.classList.toggle('hidden', !isActivePanel);
                if (isActivePanel) {
                    const activeCard =
                        panel.querySelector('.solver-target-card.active') ||
                        panel.querySelector('.solver-target-card');
                    if (activeCard && activeCard.dataset.targetId) {
                        setSolverTarget(activeCard.dataset.targetId);
                    }
                }
            });
        });
    });

    document.querySelectorAll('.solver-target-card[data-target-id]').forEach((card) => {
        card.addEventListener('click', () => {
            setSolverTarget(card.dataset.targetId);
        });
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

    // Frame Lock Pills (View 1)
    document.querySelectorAll('#frameLockPills .frame-pill').forEach((pill) => {
        pill.addEventListener('click', () => {
            document
                .querySelectorAll('#frameLockPills .frame-pill')
                .forEach((p) => p.classList.remove('active'));
            pill.classList.add('active');
            sceneController.lockFrame = pill.dataset.frame;
            const svgOpt = document.getElementById('svgOptFrame');
            if (svgOpt) svgOpt.value = pill.dataset.frame;
            sceneController.render();
        });
    });

    // Frame Lock Pills (View 2)
    document.querySelectorAll('#frameLockPills2 .frame-pill').forEach((pill) => {
        pill.addEventListener('click', () => {
            document
                .querySelectorAll('#frameLockPills2 .frame-pill')
                .forEach((p) => p.classList.remove('active'));
            pill.classList.add('active');
            if (sceneController2) {
                sceneController2.lockFrame = pill.dataset.frame;
                sceneController2.render();
            }
        });
    });

    // Clickable Display Toggle Buttons (Shared below 3D Plot)
    const propByChkId = {
        chkShowB: 'showB',
        chkShowM: 'showM',
        chkShowL: 'showL',
        chkShowK: 'showK',
        chkShowEasyPlane: 'showEasyPlane',
        chkShowCrystalAxes: 'showCrystalAxes',
        chkShowLabMagnet: 'showLabMagnet'
    };

    document.querySelectorAll('[data-opt-chk]').forEach((btn) => {
        const chkId = btn.dataset.optChk;
        const chk = document.getElementById(chkId);
        const prop = propByChkId[chkId];
        if (!chk || !prop) return;

        btn.addEventListener('click', () => {
            chk.checked = !chk.checked;
            sceneController.visibility[prop] = chk.checked;
            btn.classList.toggle('active', chk.checked);
            sceneController.render();
            if (sceneController2) sceneController2.render();
        });
    });

    // Right Column Selectors (Vector Table Frame & 3x3 Matrix Frame)
    document.getElementById('selectVectorTableFrame').addEventListener('change', recomputeAndRender);
    document.getElementById('selectMatFrom').addEventListener('change', recomputeAndRender);
    document.getElementById('selectMatTo').addEventListener('change', recomputeAndRender);
    ['customVecX', 'customVecY', 'customVecZ'].forEach((id) => {
        document.getElementById(id).addEventListener('input', recomputeAndRender);
    });

    // Copy 3x3 Matrix as LaTeX
    document.getElementById('btnCopyMatrixLatex').addEventListener('click', () => {
        const st = evaluateSetup(getActiveSetup());
        const fromF = document.getElementById('selectMatFrom').value;
        const toF = document.getElementById('selectMatTo').value;
        const R = st.matrices.getRotationMatrix(toF, fromF);
        const rows = R.map((r) => r.map((v) => v.toFixed(4)).join(' & ')).join(' \\\\\n  ');
        const tex = `R_{\\text{${toF}}\\leftarrow\\text{${fromF}}} = \\begin{pmatrix}\n  ${rows}\n\\end{pmatrix}`;
        navigator.clipboard.writeText(tex);
        const btn = document.getElementById('btnCopyMatrixLatex');
        const old = btn.innerHTML;
        btn.textContent = 'Kopiert';
        setTimeout(() => (btn.innerHTML = old), 1400);
    });

    // 1A: Sample Mounting Preset Buttons (if present)
    document.querySelectorAll('.mount-preset-btn[data-mpreset]').forEach((btn) => {
        btn.addEventListener('click', () => {
            const s = getActiveSetup();
            const newPsi = computeMountingPresetPsiS(s, btn.dataset.mpreset);
            s.psiSDeg = Number(newPsi.toFixed(2));
            syncControlsFromSetup();
            recomputeAndRender();
        });
    });

    // Bidirectional Angle Translator Inputs
    bindAngleTranslatorEvents();

    // Modals: SVG Generator & Python Snippet
    bindModals();
}

function bindAngleTranslatorEvents() {
    const applyTranslatedVectorSub = (uSub, hwFallback) => {
        const s = getActiveSetup();
        const lockPsi = document.getElementById('chkLockPsiInTranslator').checked;
        if (lockPsi) {
            // Keep manual sample mounting psi_S fixed; find best magnet angle phi_M!
            const st = evaluateSetup(s);
            const uLab = matVec3(st.matrices.R_Lab_Sub, uSub);
            if (Math.hypot(uLab[0], uLab[1]) > 1e-9) {
                s.phiMDeg = Number(wrapDeg360(Math.atan2(uLab[1], uLab[0]) * RAD2DEG).toFixed(2));
            }
        } else {
            s.mountMode = hwFallback.mountMode;
            s.psiSDeg = Number(hwFallback.psiSDeg.toFixed(2));
            s.phiMDeg = Number(hwFallback.phiMDeg.toFixed(2));
        }
        syncControlsFromSetup();
        recomputeAndRender();
    };

    // 1. Wave Perspective (phi_B,k and beta_OOP)
    const onWaveTransChange = () => {
        const s = getActiveSetup();
        const phiW = parseFloat(document.getElementById('transPhiWave').value) || 0;
        const elev = parseFloat(document.getElementById('transElevSub').value) || 0;
        const hw = hardwareAnglesFromWavePerspective(s, phiW, elev, true);
        const stTmp = evaluateSetup({ ...s, ...hw });
        applyTranslatedVectorSub(stTmp.vectors.bHat.Sub, hw);
    };
    document.getElementById('transPhiWave').addEventListener('change', onWaveTransChange);
    document.getElementById('transElevSub').addEventListener('change', onWaveTransChange);

    // 2. Substrate Perspective (phi_B,Sub and theta_B,n)
    const onSubTransChange = () => {
        const s = getActiveSetup();
        const phiSub = parseFloat(document.getElementById('transPhiSub').value) || 0;
        const thetaBN = parseFloat(document.getElementById('transThetaBN').value) || 90;
        const hw = hardwareAnglesFromSubPerspective(s, phiSub, thetaBN, true);
        const stTmp = evaluateSetup({ ...s, ...hw });
        applyTranslatedVectorSub(stTmp.vectors.bHat.Sub, hw);
    };
    document.getElementById('transPhiSub').addEventListener('change', onSubTransChange);
    document.getElementById('transThetaBN').addEventListener('change', onSubTransChange);

    // 3. Crystal / Easy-Plane Perspective (phi_easy(m), phi_easy(l), theta_B,c)
    const applyCrystalAzimuth = (phiEasyMDeg) => {
        const s = getActiveSetup();
        const thetaBC = parseFloat(document.getElementById('transThetaBC').value) || 90;
        const lockPsi = document.getElementById('chkLockPsiInTranslator').checked;
        if (lockPsi) {
            // Exact analytical solver for phi_M so that B_perp_c (i.e. m_hat) points along phiEasyMDeg at fixed psi_S!
            const st = evaluateSetup(s);
            const azRad = phiEasyMDeg * DEG2RAD;
            const uNorm = [Math.cos(azRad), Math.sin(azRad), 0];
            const vPerpCryst = [-Math.sin(azRad), Math.cos(azRad), 0];
            const vPerpLab = matVec3(st.matrices.R_Lab_Sub, matVec3(st.matrices.R_Sub_from.Cryst, vPerpCryst));
            const uEasyLab = matVec3(st.matrices.R_Lab_Sub, matVec3(st.matrices.R_Sub_from.Cryst, uNorm));
            if (Math.hypot(vPerpLab[0], vPerpLab[1]) > 1e-9) {
                const c1 = wrapDeg360(Math.atan2(-vPerpLab[0], vPerpLab[1]) * RAD2DEG);
                const c2 = wrapDeg360(c1 + 180);
                const r1 = c1 * DEG2RAD;
                const dot1 = Math.cos(r1) * uEasyLab[0] + Math.sin(r1) * uEasyLab[1];
                s.phiMDeg = Number((dot1 >= 0 ? c1 : c2).toFixed(2));
                syncControlsFromSetup();
                recomputeAndRender();
                return;
            }
        }
        const hw = hardwareAnglesFromCrystalPerspective(s, phiEasyMDeg, thetaBC, true);
        const stTmp = evaluateSetup({ ...s, ...hw });
        applyTranslatedVectorSub(stTmp.vectors.bHat.Sub, hw);
    };

    document.getElementById('transPhiEasyM').addEventListener('change', (e) => {
        const phiM = parseFloat(e.target.value) || 0;
        applyCrystalAzimuth(phiM);
    });
    document.getElementById('transPhiEasyL').addEventListener('change', (e) => {
        const phiL = parseFloat(e.target.value) || 0;
        // Since l_hat = m_hat x e_c = (m_m, -m_a, 0), phi_easy(m) = phi_easy(l) + 90°
        applyCrystalAzimuth(wrapDeg360(phiL + 90));
    });
    document.getElementById('transThetaBC').addEventListener('change', () => {
        const phiM = parseFloat(document.getElementById('transPhiEasyM').value) || 0;
        const s = getActiveSetup();
        const thetaBC = parseFloat(document.getElementById('transThetaBC').value) || 90;
        const hw = hardwareAnglesFromCrystalPerspective(s, phiM, thetaBC, true);
        const stTmp = evaluateSetup({ ...s, ...hw });
        applyTranslatedVectorSub(stTmp.vectors.bHat.Sub, hw);
    });
}

function applySavedPreset(presetId) {
    const found = savedPresets.find((p) => p.id === presetId);
    if (!found) return;
    activePresetId = found.id;

    const targetSetup = createDefaultSetup({
        ...found.setup,
        name: found.name
    });

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
    const newPreset = {
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
    };
    s.name = name;
    savedPresets.push(newPreset);
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
    if (activePresetId) {
        applySavedPreset(activePresetId);
    }
}

// ============================================================================
// Core Evaluation & DOM Update Pipeline
// ============================================================================

function recomputeAndRender() {
    const stateA = evaluateSetup(appState.setupA);
    const stateB = evaluateSetup(appState.setupB);
    const activeState = appState.activeTab === 'A' ? stateA : stateB;
    const otherState = appState.activeTab === 'A' ? stateB : stateA;

    // 0. Update Step Card Header Summaries
    updateStepSummaries(activeState);

    // 1. Update Substrate & Wave Miller-Bravais Badges (rendered with KaTeX overbars!)
    document.getElementById('badgeSubXs').innerHTML = tex(activeState.millerBravais.xs.latex);
    document.getElementById('badgeSubYs').innerHTML = tex(activeState.millerBravais.ys.latex);
    document.getElementById('badgeSubZs').innerHTML = tex(activeState.millerBravais.zsPlane.latex);
    document.getElementById('badgeWaveK').innerHTML = tex(activeState.millerBravais.k.latex);
    document.getElementById('badgeWaveT').innerHTML = tex(activeState.millerBravais.t.latex);
    document.getElementById('badgeAngleKC').textContent = `${activeState.metrics.angleKCDeg.toFixed(1)}°`;

    // 1a. Update Live Miller-Bravais Readout & Quick-Angle Buttons on the Device Rotation Slider
    updateDeviceMillerReadout(activeState);

    // 1b. Update Sample Mounting Summary & Smart Magnet Jumps for FIXED psi_S
    updateMountingAndSmartJumps(activeState);

    // 2. Update Primary & Optional Secondary 3D Visualizers
    // Primary rendered setup is activeState; if compareActive, pass the other setup as overlay
    sceneController.updateStates(activeState, otherState, appState.compareActive);
    if (sceneController2) {
        sceneController2.updateStates(activeState, otherState, appState.compareActive);
    }

    // 3. Update Key Physical Metric Cards & Bidirectional Angle Translator
    updateMetricsCards(activeState, stateA, stateB);
    updateAngleTranslatorInputs(activeState);

    // 4. Update Side-by-Side Compare Table
    updateCompareTable(stateA, stateB);

    // 5. Update 5-Frame Vector Table
    updateVectorTable(activeState);

    // 6. Update 3x3 Rotation Matrix & Custom Vector Transformer
    updateMatrixInspector(activeState);

    // 7. Update Inverse Target Finder
    updateInverseSolver(activeState);
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
    const s = activeState.setup;
    const landmarks = solveSweepLandmarksForFixedPsi(s);

    const jumpsContainer = document.getElementById('smartMagnetJumps');
    const jumpItems = [
        {
            label: `Feld ${tex('\\parallel')} Welle ${tex('\\hat{k}')}`,
            phi: landmarks.bParallelK.phiMDeg
        },
        {
            label: `Feld ${tex('\\perp')} Welle ${tex('\\hat{k}')}`,
            phi: landmarks.bPerpK.phiMDeg
        },
        {
            label: `Feld in Easy-Plane`,
            phi: landmarks.bInEasyPlane.phiM1
        },
        {
            label: `Moment ${tex('\\hat{m} \\parallel')} Welle ${tex('\\hat{k}')}`,
            phi: landmarks.mAlignK.phiMDeg
        }
    ];

    jumpsContainer.innerHTML = jumpItems
        .map(
            (j) => `
        <button type="button" class="smart-jump-btn" data-jump-phi="${j.phi}">
            <span>${j.label}</span>
            <code>${j.phi.toFixed(1)}°</code>
        </button>`
        )
        .join('');

    jumpsContainer.querySelectorAll('button[data-jump-phi]').forEach((btn) => {
        btn.addEventListener('click', () => {
            const targetPhi = parseFloat(btn.dataset.jumpPhi);
            getActiveSetup().phiMDeg = Number(targetPhi.toFixed(2));
            syncControlsFromSetup();
            recomputeAndRender();
        });
    });
}

function updateAngleTranslatorInputs(activeState) {
    const m = activeState.metrics;
    const setIfNotFocused = (id, val) => {
        const el = document.getElementById(id);
        if (el && document.activeElement !== el) {
            el.value = Number(val).toFixed(1);
        }
    };
    setIfNotFocused('transPhiWave', m.phiBWaveDeg);
    setIfNotFocused('transElevSub', m.elevationSampleDeg);
    setIfNotFocused('transPhiSub', m.phiBSubDeg);
    setIfNotFocused('transThetaBN', m.angleBNDeg);
    setIfNotFocused('transPhiEasyM', m.phiEasyMDeg);
    setIfNotFocused('transPhiEasyL', m.phiEasyLDeg);
    setIfNotFocused('transThetaBC', m.angleBCDeg);
}

function updateDeviceMillerReadout(activeState) {
    const s = activeState.setup;
    const mb = activeState.millerBravais;

    const elK = document.getElementById('liveDevMillerK');
    const elT = document.getElementById('liveDevMillerT');
    const elEdge0 = document.getElementById('liveChipEdge0');
    const elEdge90 = document.getElementById('liveChipEdge90');
    const sumEdge0 = document.getElementById('summaryEdge0');

    if (elK) elK.innerHTML = tex(mb.k.latex);
    if (elT) elT.innerHTML = tex(mb.t.latex);
    if (elEdge0) elEdge0.innerHTML = tex(mb.xs.latex);
    if (elEdge90) elEdge90.innerHTML = tex(mb.ys.latex);
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

    const t0El = document.getElementById('lblEdgePreset0Title');
    const s0El = document.getElementById('lblEdgePreset0Sub');
    const t90El = document.getElementById('lblEdgePreset90Title');
    const s90El = document.getElementById('lblEdgePreset90Sub');
    const btn0 = document.getElementById('btnEdgePreset0');
    const btn90 = document.getElementById('btnEdgePreset90');

    if (t0El) t0El.innerHTML = info.t0;
    if (s0El) s0El.innerHTML = info.s0;
    if (t90El) t90El.innerHTML = info.t90;
    if (s90El) s90El.innerHTML = info.s90;
    if (btn0) btn0.classList.toggle('active', Math.abs(s.alphaFlatDeg) < 0.25);
    if (btn90) btn90.classList.toggle('active', Math.abs(s.alphaFlatDeg - 90) < 0.25);

    const quickGroup = document.getElementById('deviceQuickAngleChips');
    if (quickGroup) {
        quickGroup.querySelectorAll('button[data-val]').forEach((btn) => {
            const deg = parseFloat(btn.dataset.val);
            btn.textContent = `${deg}°`;
            btn.classList.toggle('active', Math.abs(s.thetaDevDeg - deg) < 0.25);
        });
    }
}

function updateMetricsCards(activeState, stateA, stateB) {
    const m = activeState.metrics;

    // 1. Magnetfeld B zur Welle k & aktiver Anteil in der Easy-Plane
    const elBAngle = document.getElementById('valBSampleAngles');
    const elBEasy = document.getElementById('valBEasy');
    if (elBAngle) {
        elBAngle.textContent = `${m.angleBKDeg.toFixed(1)}°`;
    }
    if (elBEasy) {
        elBEasy.innerHTML = `In Easy-Plane aktiv: <strong>${m.bEasyTesla.toFixed(2)} T</strong> (${(m.bEasyRatio * 100).toFixed(0)}%)`;
    }

    // 2. Magnetisierung m zur Welle k & CPW-Anregung
    const elAngleMK = document.getElementById('valAngleMK');
    const elMDotK = document.getElementById('valMDotK');
    if (elAngleMK) {
        elAngleMK.textContent = m.isEasyPlaneSingular
            ? '–'
            : `${m.angleMKDeg.toFixed(1)}°`;
    }
    if (elMDotK) {
        elMDotK.innerHTML = `CPW-Anregung: <strong>${(m.cpwTorqueM * 100).toFixed(0)}%</strong>`;
    }

    // 3. Néel-Vektor l zur Welle k
    const elAngleLK = document.getElementById('valAngleLK');
    const elLDotK = document.getElementById('valLDotK');
    if (elAngleLK) {
        elAngleLK.textContent = m.isEasyPlaneSingular
            ? '–'
            : `${m.angleLKDeg.toFixed(1)}°`;
    }
    if (elLDotK) {
        elLDotK.innerHTML = `Projektion auf ${tex('\\hat{k}')}: <strong>${(Math.abs(m.lDotK) * 100).toFixed(0)}%</strong>`;
    }

    // Compare delta badges on result tiles
    const cmpIds = ['cmpBSample', 'cmpAngleMK', 'cmpAngleLK'];
    cmpIds.forEach((id) => {
        const el = document.getElementById(id);
        if (el) el.classList.toggle('hidden', !appState.compareActive);
    });
    if (appState.compareActive) {
        const mA = stateA.metrics;
        const mB = stateB.metrics;
        const cmpB = document.getElementById('cmpBSample');
        const cmpM = document.getElementById('cmpAngleMK');
        const cmpL = document.getElementById('cmpAngleLK');
        if (cmpB) cmpB.textContent = `A: ${mA.angleBKDeg.toFixed(1)}° vs. B: ${mB.angleBKDeg.toFixed(1)}°`;
        if (cmpM) cmpM.textContent = `A: ${mA.angleMKDeg.toFixed(1)}° vs. B: ${mB.angleMKDeg.toFixed(1)}°`;
        if (cmpL) cmpL.textContent = `A: ${mA.angleLKDeg.toFixed(1)}° vs. B: ${mB.angleLKDeg.toFixed(1)}°`;
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
        ['Einbaulage / Schnitt', `${sA.mountMode} / ${sA.cutPreset}`, `${sB.mountMode} / ${sB.cutPreset}`, '—'],
        [`Winkel (${tex('\\psi_S, \\phi_M')})`, `${sA.psiSDeg.toFixed(1)}°, ${sA.phiMDeg.toFixed(1)}°`, `${sB.psiSDeg.toFixed(1)}°, ${sB.phiMDeg.toFixed(1)}°`, `${(sB.phiMDeg - sA.phiMDeg).toFixed(1)}°`],
        [`Wellenrichtung ${tex('\\hat{k}')}`, tex(stateA.millerBravais.k.latex), tex(stateB.millerBravais.k.latex), `${(sB.thetaDevDeg - sA.thetaDevDeg).toFixed(1)}°`],
        [`Anteil in Easy-Plane`, `${(mA.bEasyRatio * 100).toFixed(1)}%`, `${(mB.bEasyRatio * 100).toFixed(1)}%`, `${((mB.bEasyRatio - mA.bEasyRatio) * 100).toFixed(1)}%`],
        [`Anteil in Sample-Ebene`, `${(mA.bInSampleRatio * 100).toFixed(1)}%`, `${(mB.bInSampleRatio * 100).toFixed(1)}%`, `${((mB.bInSampleRatio - mA.bInSampleRatio) * 100).toFixed(1)}%`],
        [`Winkel Feld zu Welle ${tex('\\angle(\\mathbf{B}, \\hat{k})')}`, `${mA.angleBKDeg.toFixed(2)}°`, `${mB.angleBKDeg.toFixed(2)}°`, `${(mB.angleBKDeg - mA.angleBKDeg).toFixed(2)}°`],
        [`Winkel Moment zu Welle ${tex('\\angle(\\hat{m}, \\hat{k})')}`, `${mA.angleMKDeg.toFixed(2)}°`, `${mB.angleMKDeg.toFixed(2)}°`, `${(mB.angleMKDeg - mA.angleMKDeg).toFixed(2)}°`],
        [`Winkel Néel zu Welle ${tex('\\angle(\\hat{l}, \\hat{k})')}`, `${mA.angleLKDeg.toFixed(2)}°`, `${mB.angleLKDeg.toFixed(2)}°`, `${(mB.angleLKDeg - mA.angleLKDeg).toFixed(2)}°`],
        [`CPW-Anregung ${tex('|\\hat{k}\\times\\hat{m}|')}`, mA.cpwTorqueM.toFixed(3), mB.cpwTorqueM.toFixed(3), (mB.cpwTorqueM - mA.cpwTorqueM).toFixed(3)]
    ];

    document.getElementById('compareTableBody').innerHTML = rows
        .map(
            (r) => `<tr><td>${r[0]}</td><td>${r[1]}</td><td>${r[2]}</td><td>${r[3]}</td></tr>`
        )
        .join('');
}

function updateVectorTable(activeState) {
    const frame = document.getElementById('selectVectorTableFrame').value;
    const meta = FRAME_METADATA[frame] || FRAME_METADATA.Wave;
    const headerRow = document.getElementById('vectorTableHeaderRow');
    headerRow.innerHTML = `
        <th>Richtung</th>
        <th>${meta.axesShort[0]}</th>
        <th>${meta.axesShort[1]}</th>
        <th>${meta.axesShort[2]}</th>
        <th>Kristall [uvtw]</th>
    `;

    const v = activeState.vectors;
    const mb = activeState.millerBravais;

    const entries = [
        { name: `Magnetfeld ${tex('\\hat{b}')}`, color: VECTOR_COLORS.B, vec: v.bHat[frame], mb: tex(mb.B.latex) },
        { name: `FM-Moment ${tex('\\hat{m}')}`, color: VECTOR_COLORS.m, vec: v.mHat[frame], mb: tex(mb.m.latex) },
        { name: `Néel-Vektor ${tex('\\hat{l}')}`, color: VECTOR_COLORS.l, vec: v.lHat[frame], mb: tex(mb.l.latex) },
        { name: `Welle ${tex('\\hat{k}')}`, color: VECTOR_COLORS.k, vec: v.kHat[frame], mb: tex(mb.k.latex) },
        { name: `Quer zur Welle ${tex('\\hat{t}')}`, color: VECTOR_COLORS.t, vec: v.tHat[frame], mb: tex(mb.t.latex) },
        { name: `Sample-Normale ${tex('\\hat{n}')}`, color: VECTOR_COLORS.n, vec: v.nHat[frame], mb: tex(mb.zsDir.latex) },
        { name: `Kristall ${tex('c')}-Achse`, color: VECTOR_COLORS.ec, vec: v.ecHat[frame], mb: tex('[0001]') }
    ];

    const fmtComp = (x) => {
        const val = Math.abs(x) < 1e-9 ? 0 : x;
        return (val >= 0 ? '+' : '') + val.toFixed(4);
    };

    document.getElementById('vectorTableBody').innerHTML = entries
        .map(
            (e) => `
        <tr>
            <td><span style="color:${e.color}">●</span> ${e.name}</td>
            <td>${fmtComp(e.vec[0])}</td>
            <td>${fmtComp(e.vec[1])}</td>
            <td>${fmtComp(e.vec[2])}</td>
            <td>${e.mb}</td>
        </tr>`
        )
        .join('');
}

function updateMatrixInspector(activeState) {
    const fromF = document.getElementById('selectMatFrom').value;
    const toF = document.getElementById('selectMatTo').value;
    const R = activeState.matrices.getRotationMatrix(toF, fromF);

    const cells = [];
    for (let i = 0; i < 3; i++) {
        for (let j = 0; j < 3; j++) {
            const val = Math.abs(R[i][j]) < 1e-9 ? 0 : R[i][j];
            cells.push(`<div class="matrix-cell">${(val >= 0 ? '+' : '') + val.toFixed(4)}</div>`);
        }
    }
    document.getElementById('matrixDisplayGrid').innerHTML = cells.join('');

    const vx = parseFloat(document.getElementById('customVecX').value) || 0;
    const vy = parseFloat(document.getElementById('customVecY').value) || 0;
    const vz = parseFloat(document.getElementById('customVecZ').value) || 0;
    const out = matVec3(R, [vx, vy, vz]);
    document.getElementById('customVecResult').textContent = `(${out.map((c) => c.toFixed(4)).join(', ')})`;
}

function updateInverseSolver(activeState) {
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
// Modals: Thesis SVG Generator & Python Snippet
// ============================================================================

function bindModals() {
    const modalSvg = document.getElementById('modalSvgExport');
    const modalPy = document.getElementById('modalPythonSnippet');

    const refreshSvgPreview = () => {
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
    };

    document.getElementById('btnOpenSvgModal').addEventListener('click', () => {
        modalSvg.classList.remove('hidden');
        refreshSvgPreview();
    });
    document.getElementById('btnCloseSvgModal').addEventListener('click', () => {
        modalSvg.classList.add('hidden');
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
    ].forEach((id) => {
        document.getElementById(id).addEventListener('change', refreshSvgPreview);
    });

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

    document.getElementById('btnCopySvgCode').addEventListener('click', () => {
        const svgStr = refreshSvgPreview();
        navigator.clipboard.writeText(svgStr);
        const btn = document.getElementById('btnCopySvgCode');
        const old = btn.innerHTML;
        btn.textContent = 'SVG-Code kopiert';
        setTimeout(() => (btn.innerHTML = old), 1500);
    });

    // Python Snippet Modal
    document.getElementById('btnOpenPythonModal').addEventListener('click', () => {
        const code = generatePythonSnippet(
            appState.setupA,
            appState.setupB,
            appState.compareActive
        );
        document.getElementById('pythonSnippetCode').textContent = code;
        modalPy.classList.remove('hidden');
    });
    document.getElementById('btnClosePythonModal').addEventListener('click', () => {
        modalPy.classList.add('hidden');
    });
    document.getElementById('btnCopyPythonSnippet').addEventListener('click', () => {
        const code = document.getElementById('pythonSnippetCode').textContent;
        navigator.clipboard.writeText(code);
        const btn = document.getElementById('btnCopyPythonSnippet');
        const old = btn.innerHTML;
        btn.textContent = 'Python-Code kopiert';
        setTimeout(() => (btn.innerHTML = old), 1500);
    });

    // Close modals on backdrop click
    [modalSvg, modalPy].forEach((m) => {
        m.addEventListener('click', (e) => {
            if (e.target === m) m.classList.add('hidden');
        });
    });
}
