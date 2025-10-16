import {
    millerBravaisToCartesian,
    cartesianToMillerBravais,
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
            const [h, k, i, l] = axes[axis];
            inputsDef[axis].h.value = h;
            inputsDef[axis].k.value = k;
            inputsDef[axis].i.value = i;
            inputsDef[axis].l.value = l;
        });
    };

    // Loads a preset into the UI
    const loadPreset = (presetName) => {
        const preset = presets[presetName];
        if (preset) {
            populateDefinitionInputs(preset.axes);
            // After loading a preset the transforms may change.
            // Reset internal vectors to avoid inconsistencies.
            preciseCrystalVector = [0, 0, 0, 0];
            preciseLabVector = [0, 0, 0];
            updateCrystalInputs();
            updateLabInputs();
        }
    };

    // Reads a Miller-Bravais vector from a group of inputs
    const getMbVectorFromInputs = (inputs) => {
        return [
            parseFloat(inputs.h.value) || 0,
            parseFloat(inputs.k.value) || 0,
            parseFloat(inputs.i.value) || 0,
            parseFloat(inputs.l.value) || 0
        ];
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

    // --- CORE LOGIC ---

    // Builds the transformation matrix M (transforms from crystal-cartesian to lab)
    const getTransformationMatrix = () => {
        const v_x_mb = getMbVectorFromInputs(inputsDef.x);
        const v_y_mb = getMbVectorFromInputs(inputsDef.y);
        const v_z_mb = getMbVectorFromInputs(inputsDef.z);

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
            x: {mb: getMbVectorFromInputs(inputsDef.x), isSet: getMbVectorFromInputs(inputsDef.x).some(v => v !== 0)},
            y: {mb: getMbVectorFromInputs(inputsDef.y), isSet: getMbVectorFromInputs(inputsDef.y).some(v => v !== 0)},
            z: {mb: getMbVectorFromInputs(inputsDef.z), isSet: getMbVectorFromInputs(inputsDef.z).some(v => v !== 0)},
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
        setMbVectorToInputs(inputsDef[v3_name], normalized_v3_mb);

        // Recompute the i-index and notify listeners by dispatching input events
        // so dependent UI/state updates run as if the user typed the values.
        updateI(inputsDef[v3_name].h, inputsDef[v3_name].k, inputsDef[v3_name].i);
        // Dispatch input events to trigger any attached listeners (e.g., precise state updates)
        inputsDef[v3_name].h.dispatchEvent(new Event('input', { bubbles: true }));
        inputsDef[v3_name].k.dispatchEvent(new Event('input', { bubbles: true }));
    };

    const handleTransformToCrystal = () => {
        const M = getTransformationMatrix();
        // compute using precise state
        const v_cryst_cart = labToCrystal(preciseLabVector, M);
        preciseCrystalVector = cartesianToMillerBravais(v_cryst_cart);
        // update display
        updateCrystalInputs();
    };

    const handleTransformToLab = () => {
        const M = getTransformationMatrix();
        // compute using precise state
        const v_cryst_cart = millerBravaisToCartesian(preciseCrystalVector);
        preciseLabVector = crystalToLab(v_cryst_cart, M);
        // update display
        updateLabInputs();
    };

    // --- INITIALIZATION ---
    const init = async () => {
        // Load presets from JSON
        try {
            const response = await fetch('presets.json');
            presets = await response.json();
        } catch (error) {
            console.error('Error loading presets:', error);
            alert('Error: Could not load predefined cuts/presets.');
            return;
        }

        // Load presets into dropdown
        Object.keys(presets).forEach(name => {
            const option = document.createElement('option');
            option.value = name;
            option.textContent = name;
            presetSelect.appendChild(option);
        });

        // Event Listeners
        presetSelect.addEventListener('change', (e) => loadPreset(e.target.value));
        calculateLastAxisBtn.addEventListener('click', handleCalculateLastAxis);
        transformToCrystalBtn.addEventListener('click', handleTransformToCrystal);
        transformToLabBtn.addEventListener('click', handleTransformToLab);

        // Listener for i-index
        ['x', 'y', 'z'].forEach(axis => {
            inputsDef[axis].h.addEventListener('input', () => updateI(inputsDef[axis].h, inputsDef[axis].k, inputsDef[axis].i));
            inputsDef[axis].k.addEventListener('input', () => updateI(inputsDef[axis].h, inputsDef[axis].k, inputsDef[axis].i));
        });
        inputsCrystal.h.addEventListener('input', () => updateI(inputsCrystal.h, inputsCrystal.k, inputsCrystal.i));
        inputsCrystal.k.addEventListener('input', () => updateI(inputsCrystal.h, inputsCrystal.k, inputsCrystal.i));

        // Listener that updates the precise state on input
        inputsCrystal.h.addEventListener('input', () => {
            updateI(inputsCrystal.h, inputsCrystal.k, inputsCrystal.i);
            updatePreciseCrystalFromInputs();
        });
        inputsCrystal.k.addEventListener('input', () => {
            updateI(inputsCrystal.h, inputsCrystal.k, inputsCrystal.i);
            updatePreciseCrystalFromInputs();
        });
        inputsCrystal.l.addEventListener('input', updatePreciseCrystalFromInputs);

        [inputsLab.x, inputsLab.y, inputsLab.z].forEach(el => el.addEventListener('input', updatePreciseLabFromCartesianInputs));
        [inputsLab.theta, inputsLab.phi].forEach(el => el.addEventListener('input', updatePreciseLabFromAngleInputs));


        // Load default preset if available
        if (Object.keys(presets).length > 0) {
            loadPreset(Object.keys(presets)[0]);
        }
        // Initialize precise vectors with start values from HTML
        updatePreciseCrystalFromInputs();
        updatePreciseLabFromCartesianInputs();
    };

    init();
});
