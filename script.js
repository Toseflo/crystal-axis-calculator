import {
    millerBravaisToCartesian,
    cartesianToMillerBravais,
    millerBravaisPlaneToDirection,
    labToCrystal,
    crystalToLab,
    normalize,
    crossProduct,
    transposeMatrix
} from './transformation.js';

document.addEventListener('DOMContentLoaded', () => {

    // --- CONSTANTS AND VARIABLES ---
    let presets = {};

    // --- STATE VARIABLES ---
    let preciseCrystalVector = [0, 0, 0, 0]; // [h, k, i, l]
    let preciseLabVector = [0, 0, 0];     // [x, y, z]

    // --- DOM ELEMENTS ---
    const presetSelect = document.getElementById('preset-select');
    const calculateLastAxisBtn = document.getElementById('calculate-last-axis-btn');
    const transformToCrystalBtn = document.getElementById('transform-to-crystal-btn');
    const transformToLabBtn = document.getElementById('transform-to-lab-btn');

    const modeToggles = {
        x: document.getElementById('x-mode-toggle'),
        y: document.getElementById('y-mode-toggle'),
        z: document.getElementById('z-mode-toggle'),
    };

    const inputsDef = {
        x: {
            h: document.getElementById('x_h'),
            k: document.getElementById('x_k'),
            i: document.getElementById('x_i'),
            l: document.getElementById('x_l')
        },
        y: {
            h: document.getElementById('y_h'),
            k: document.getElementById('y_k'),
            i: document.getElementById('y_i'),
            l: document.getElementById('y_l')
        },
        z: {
            h: document.getElementById('z_h'),
            k: document.getElementById('z_k'),
            i: document.getElementById('z_i'),
            l: document.getElementById('z_l')
        },
    };

    const inputsLab = {
        x: document.getElementById('lab_x'),
        y: document.getElementById('lab_y'),
        z: document.getElementById('lab_z'),
        theta: document.getElementById('lab_theta'),
        phi: document.getElementById('lab_phi'),
    };

    const inputsCrystal = {
        h: document.getElementById('crystal_h'),
        k: document.getElementById('crystal_k'),
        i: document.getElementById('crystal_i'),
        l: document.getElementById('crystal_l'),
    };

    // --- UI UPDATE & STATE UPDATE FUNCTIONS ---

    // Updates the crystal display inputs from the precise state
    const updateCrystalInputs = () => {
        // Normalize only for display
        const displayVec = normalizeMillerBravais(preciseCrystalVector);
        setMbVectorToInputs(inputsCrystal, displayVec);
    };

    // Updates the lab display inputs from the precise state
    const updateLabInputs = () => {
        const displayVec = normalize(preciseLabVector);
        inputsLab.x.value = displayVec[0].toFixed(4);
        inputsLab.y.value = displayVec[1].toFixed(4);
        inputsLab.z.value = displayVec[2].toFixed(4);
        updateAnglesFromXYZ(); // also update angles
    };

    // Updates the precise crystal state from the UI fields
    const updatePreciseCrystalFromInputs = () => {
        const h = parseFloat(inputsCrystal.h.value) || 0;
        const k = parseFloat(inputsCrystal.k.value) || 0;
        const l = parseFloat(inputsCrystal.l.value) || 0;
        preciseCrystalVector = [h, k, -(h + k), l];
    };

    // Updates the precise lab state from the cartesian UI fields
    const updatePreciseLabFromCartesianInputs = () => {
        const x = parseFloat(inputsLab.x.value) || 0;
        const y = parseFloat(inputsLab.y.value) || 0;
        const z = parseFloat(inputsLab.z.value) || 0;
        preciseLabVector = [x, y, z];
        // synchronize angle display
        updateAnglesFromXYZ();
    };

    // Updates the precise lab state from the spherical angle UI fields
    const updatePreciseLabFromAngleInputs = () => {
        const thetaDeg = parseFloat(inputsLab.theta.value) || 0;
        const phiDeg = parseFloat(inputsLab.phi.value) || 0;
        const thetaRad = thetaDeg * Math.PI / 180;
        const phiRad = phiDeg * Math.PI / 180;

        const x = Math.sin(thetaRad) * Math.cos(phiRad);
        const y = Math.sin(thetaRad) * Math.sin(phiRad);
        const z = Math.cos(thetaRad);
        preciseLabVector = [x, y, z];
        // synchronize cartesian display
        inputsLab.x.value = x.toFixed(4);
        inputsLab.y.value = y.toFixed(4);
        inputsLab.z.value = z.toFixed(4);
    };

    // --- UI FUNCTIONS ---

    // Updates the i-index based on h and k
    const updateI = (hInput, kInput, iInput) => {
        const h = parseFloat(hInput.value) || 0;
        const k = parseFloat(kInput.value) || 0;
        iInput.value = (-(h + k)).toFixed(4);
    };

    // Fills definition fields for axes
    const populateDefinitionInputs = (axes) => {
        ['x', 'y', 'z'].forEach(axis => {
            if (axes[axis]) {
                const {type, value} = axes[axis];
                const [h, k, i, l] = value;
                inputsDef[axis].h.value = h;
                inputsDef[axis].k.value = k;
                inputsDef[axis].i.value = i;
                inputsDef[axis].l.value = l;

                // Update toggle button state
                const toggle = modeToggles[axis];
                toggle.dataset.mode = type;
                toggle.querySelector('span').textContent = type.charAt(0).toUpperCase() + type.slice(1);
            } else {
                // Clear inputs if axis is not defined in preset
                Object.values(inputsDef[axis]).forEach(input => input.value = '');
                const toggle = modeToggles[axis];
                toggle.dataset.mode = 'direction';
                toggle.querySelector('span').textContent = 'Direction';
            }
        });
    };

    // Loads a preset into the UI
    const loadPreset = (presetName) => {
        if (presetName === 'custom') return;
        const preset = presets[presetName];
        if (preset) {
            populateDefinitionInputs(preset.axes);
            // After loading a preset the transforms may change.
            // Reset internal vectors to avoid inconsistencies.
            preciseCrystalVector = [0, 0, 0, 0];
            preciseLabVector = [0, 0, 0];
            updateCrystalInputs();
            updateLabInputs();
            updateSpecialAxesDisplay();
        }
    };

    // Reads a Miller-Bravais vector from a group of inputs, considering the mode
    const getMbVectorFromInputs = (inputs, axis) => {
        let mb = [
            parseFloat(inputs.h.value) || 0,
            parseFloat(inputs.k.value) || 0,
            parseFloat(inputs.i.value) || 0,
            parseFloat(inputs.l.value) || 0
        ];

        const mode = modeToggles[axis].dataset.mode;
        if (mode === 'plane') {
            return millerBravaisPlaneToDirection(mb);
        }
        return mb;
    };

    // Writes a Miller-Bravais vector into a group of inputs
    const setMbVectorToInputs = (inputs, vec, precision = 4) => {
        const isInteger = vec.every(v => Math.abs(v - Math.round(v)) < 1e-9);
        const effectivePrecision = isInteger ? 0 : precision;
        inputs.h.value = vec[0].toFixed(effectivePrecision);
        inputs.k.value = vec[1].toFixed(effectivePrecision);
        inputs.i.value = vec[2].toFixed(effectivePrecision);
        inputs.l.value = vec[3].toFixed(effectivePrecision);
    };

    // Attempts to convert Miller-Bravais indices to smallest integers
    const normalizeMillerBravais = (vec, maxIndex = 10) => {
        const epsilon = 1e-6;
        if (vec.every(v => Math.abs(v) < epsilon)) return vec; // vector is [0,0,0,0]

        const smallestNonZero = Math.min(...vec.filter(v => Math.abs(v) > epsilon).map(Math.abs));
        let scaledVec = vec.map(v => v / smallestNonZero);

        for (let m = 1; m <= 100; m++) {
            const tempVec = scaledVec.map(v => v * m);
            if (tempVec.every(v => Math.abs(v - Math.round(v)) < epsilon)) {
                const finalVec = tempVec.map(Math.round);
                if (finalVec.some(v => Math.abs(v) > maxIndex)) {
                    return vec;
                }
                finalVec[2] = -(finalVec[0] + finalVec[1]);
                return finalVec;
            }
        }
        return vec;
    };


    // Updates angles based on the (x,y,z) vector
    const updateAnglesFromXYZ = () => {
        const x = parseFloat(inputsLab.x.value) || 0;
        const y = parseFloat(inputsLab.y.value) || 0;
        const z = parseFloat(inputsLab.z.value) || 0;
        const r = Math.sqrt(x * x + y * y + z * z);
        if (r === 0) {
            inputsLab.theta.value = 0;
            inputsLab.phi.value = 0;
            return;
        }
        const thetaRad = Math.acos(z / r);
        let phiRad = Math.atan2(y, x);
        if (phiRad < 0) {
            phiRad += 2 * Math.PI;
        }

        inputsLab.theta.value = (thetaRad * 180 / Math.PI).toFixed(2);
        inputsLab.phi.value = (phiRad * 180 / Math.PI).toFixed(2);
    };

    // --- NEW: Special axes display (c-axis and [11-20]) ---
    const formatVec = (v) => `(${v[0].toFixed(4)}, ${v[1].toFixed(4)}, ${v[2].toFixed(4)})`;

    const calcAnglesForVec = (v) => {
        const x = v[0], y = v[1], z = v[2];
        const r = Math.sqrt(x * x + y * y + z * z);
        if (r === 0) return {theta: '-', phi: '-'};
        const theta = Math.acos(z / r) * 180 / Math.PI;
        let phi = Math.atan2(y, x) * 180 / Math.PI;
        if (phi < 0) phi += 360;
        // Round to 2 decimals and if rounded phi equals 360 -> show 0 instead
        let thetaRounded = Math.round(theta * 100) / 100;
        let phiRounded = Math.round(phi * 100) / 100;
        if (Math.abs(phiRounded - 360) < 1e-9) phiRounded = 0;
        return {theta: thetaRounded.toFixed(2), phi: phiRounded.toFixed(2)};
    };

    const updateSpecialAxesDisplay = () => {
        // Build current transformation matrix
        const M = getTransformationMatrix();

        // c-axis in MB: [0,0,0,1]
        const c_mb = [0, 0, 0, 1];
        const c_cart_crystal = millerBravaisToCartesian(c_mb);
        const c_cart_lab = crystalToLab(c_cart_crystal, M);
        // normalize for display
        const c_cart_lab_norm = normalize(c_cart_lab);
        const c_display = formatVec(c_cart_lab_norm);
        const c_angles = calcAnglesForVec(c_cart_lab_norm);

        const caxisCartEl = document.getElementById('caxis_cart');
        const caxisThetaEl = document.getElementById('caxis_theta');
        const caxisPhiEl = document.getElementById('caxis_phi');
        if (caxisCartEl) caxisCartEl.textContent = c_display;
        if (caxisThetaEl) caxisThetaEl.textContent = c_angles.theta;
        if (caxisPhiEl) caxisPhiEl.textContent = c_angles.phi;

        // [11-20] direction MB: [1,1,-2,0]
        const a_mb = [1, 1, -2, 0];
        const a_cart_crystal = millerBravaisToCartesian(a_mb);
        const a_cart_lab = crystalToLab(a_cart_crystal, M);
        // normalize for display
        const a_cart_lab_norm = normalize(a_cart_lab);
        const a_display = formatVec(a_cart_lab_norm);
        const a_angles = calcAnglesForVec(a_cart_lab_norm);

        const aCartEl = document.getElementById('a1120_cart');
        const aThetaEl = document.getElementById('a1120_theta');
        const aPhiEl = document.getElementById('a1120_phi');
        if (aCartEl) aCartEl.textContent = a_display;
        if (aThetaEl) aThetaEl.textContent = a_angles.theta;
        if (aPhiEl) aPhiEl.textContent = a_angles.phi;
    };

    // --- PRESET / CUSTOM HANDLING ---
    // Switch preset dropdown to 'custom' when the user modifies inputs (only on real user events)
    const setPresetToCustom = () => {
        if (!presetSelect) return;
        if (presetSelect.value === 'custom') return;
        // Add a 'custom' option if it doesn't exist
        let opt = Array.from(presetSelect.options).find(o => o.value === 'custom');
        if (!opt) {
            opt = document.createElement('option');
            opt.value = 'custom';
            opt.textContent = 'Custom';
            presetSelect.appendChild(opt);
        }
        presetSelect.value = 'custom';
    };

    // Handler to attach to inputs: only mark custom on user-initiated events (isTrusted === true)
    const userChangeHandler = (e) => {
        // e may be undefined when called programmatically; ensure we only react to trusted user events
        if (e && e.isTrusted) {
            setPresetToCustom();
        }
    };

    // --- CORE LOGIC ---

    // Builds the transformation matrix M (transforms from crystal-cartesian to lab)
    const getTransformationMatrix = () => {
        const v_x_mb = getMbVectorFromInputs(inputsDef.x, 'x');
        const v_y_mb = getMbVectorFromInputs(inputsDef.y, 'y');
        const v_z_mb = getMbVectorFromInputs(inputsDef.z, 'z');

        const v_x_cart = millerBravaisToCartesian(v_x_mb);
        const v_y_cart = millerBravaisToCartesian(v_y_mb);
        const v_z_cart = millerBravaisToCartesian(v_z_mb);

        const n_x = normalize(v_x_cart);
        const n_y = normalize(v_y_cart);
        const n_z = normalize(v_z_cart);

        const M = [
            [n_x[0], n_y[0], n_z[0]],
            [n_x[1], n_y[1], n_z[1]],
            [n_x[2], n_y[2], n_z[2]]
        ];

        // Use centralized transpose helper
        return transposeMatrix(M);
    };

    const handleCalculateLastAxis = () => {
        const vectors = {
            x: {
                isSet: [inputsDef.x.h.value, inputsDef.x.k.value, inputsDef.x.l.value].some(v => v !== ''),
                mb: getMbVectorFromInputs(inputsDef.x, 'x')
            },
            y: {
                isSet: [inputsDef.y.h.value, inputsDef.y.k.value, inputsDef.y.l.value].some(v => v !== ''),
                mb: getMbVectorFromInputs(inputsDef.y, 'y')
            },
            z: {
                isSet: [inputsDef.z.h.value, inputsDef.z.k.value, inputsDef.z.l.value].some(v => v !== ''),
                mb: getMbVectorFromInputs(inputsDef.z, 'z')
            },
        };

        const definedAxes = Object.keys(vectors).filter(ax => vectors[ax].isSet);
        const undefinedAxes = Object.keys(vectors).filter(ax => !vectors[ax].isSet);

        if (definedAxes.length !== 2 || undefinedAxes.length !== 1) {
            alert("Please define exactly two axes to compute the third.");
            return;
        }

        const [v1_name, v2_name] = definedAxes;
        const v3_name = undefinedAxes[0];

        const v1_cart = millerBravaisToCartesian(vectors[v1_name].mb);
        const v2_cart = millerBravaisToCartesian(vectors[v2_name].mb);

        let v3_cart;

        if (v3_name === 'x') v3_cart = crossProduct(v1_cart, v2_cart);
        else if (v3_name === 'y') v3_cart = crossProduct(v2_cart, v1_cart);
        else v3_cart = crossProduct(v1_cart, v2_cart);

        const v3_mb = cartesianToMillerBravais(v3_cart);
        // Normalize result for display
        const normalized_v3_mb = normalizeMillerBravais(v3_mb);

        // The result is always a direction, so we set the mode accordingly
        const toggle = modeToggles[v3_name];
        toggle.dataset.mode = 'direction';
        toggle.querySelector('span').textContent = 'Direction';

        setMbVectorToInputs(inputsDef[v3_name], normalized_v3_mb);

        // Recompute the i-index and notify listeners by dispatching input events
        // so dependent UI/state updates run as if the user typed the values.
        updateI(inputsDef[v3_name].h, inputsDef[v3_name].k, inputsDef[v3_name].i);
        // Dispatch input events to trigger any attached listeners (e.g., precise state updates)
        inputsDef[v3_name].h.dispatchEvent(new Event('input', { bubbles: true }));
        inputsDef[v3_name].k.dispatchEvent(new Event('input', { bubbles: true }));
    };

    // --- EVENT LISTENERS ---

    // Initial population of presets
    fetch('presets.json')
        .then(response => response.json())
        .then(data => {
            presets = data;
            const presetNames = Object.keys(presets);
            presetNames.forEach(name => {
                const option = document.createElement('option');
                option.value = name;
                option.textContent = name;
                presetSelect.appendChild(option);
            });
            if (presetNames.length > 0) {
                presetSelect.value = presetNames[0];
                loadPreset(presetNames[0]);
            }
        });

    presetSelect.addEventListener('change', (e) => loadPreset(e.target.value));

    Object.values(modeToggles).forEach(button => {
        button.addEventListener('click', () => {
            const currentMode = button.dataset.mode;
            const newMode = currentMode === 'direction' ? 'plane' : 'direction';
            button.dataset.mode = newMode;
            button.querySelector('span').textContent = newMode.charAt(0).toUpperCase() + newMode.slice(1);
            setPresetToCustom();
            updateSpecialAxesDisplay();
        });
    });

    // Attach listeners to definition inputs
    ['x', 'y', 'z'].forEach(axis => {
        inputsDef[axis].h.addEventListener('input', (e) => {
            updateI(inputsDef[axis].h, inputsDef[axis].k, inputsDef[axis].i);
            userChangeHandler(e);
            updateSpecialAxesDisplay();
        });
        inputsDef[axis].k.addEventListener('input', (e) => {
            updateI(inputsDef[axis].h, inputsDef[axis].k, inputsDef[axis].i);
            userChangeHandler(e);
            updateSpecialAxesDisplay();
        });
        inputsDef[axis].l.addEventListener('input', userChangeHandler);
        inputsDef[axis].l.addEventListener('input', updateSpecialAxesDisplay);
    });

    // Attach listeners to crystal inputs
    inputsCrystal.h.addEventListener('input', () => updateI(inputsCrystal.h, inputsCrystal.k, inputsCrystal.i));
    inputsCrystal.k.addEventListener('input', () => updateI(inputsCrystal.h, inputsCrystal.k, inputsCrystal.i));

    // Attach listeners for transformations
    transformToLabBtn.addEventListener('click', () => {
        updatePreciseCrystalFromInputs();
        const M = getTransformationMatrix();
        const crystalCartesian = millerBravaisToCartesian(preciseCrystalVector);
        preciseLabVector = crystalToLab(crystalCartesian, M);
        updateLabInputs();
        updateSpecialAxesDisplay();
    });

    transformToCrystalBtn.addEventListener('click', () => {
        updatePreciseLabFromCartesianInputs();
        const M = getTransformationMatrix();
        const crystalCartesian = labToCrystal(preciseLabVector, M);
        preciseCrystalVector = cartesianToMillerBravais(crystalCartesian);
        updateCrystalInputs();
        updateSpecialAxesDisplay();
    });

    // Attach listeners for lab input changes
    inputsLab.x.addEventListener('input', updatePreciseLabFromCartesianInputs);
    inputsLab.y.addEventListener('input', updatePreciseLabFromCartesianInputs);
    inputsLab.z.addEventListener('input', updatePreciseLabFromCartesianInputs);
    inputsLab.theta.addEventListener('input', updatePreciseLabFromAngleInputs);
    inputsLab.phi.addEventListener('input', updatePreciseLabFromAngleInputs);

    calculateLastAxisBtn.addEventListener('click', handleCalculateLastAxis);
    // Update special axes when the calculate button completes
    calculateLastAxisBtn.addEventListener('click', updateSpecialAxesDisplay);

});
