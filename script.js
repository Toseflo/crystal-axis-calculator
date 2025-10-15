document.addEventListener('DOMContentLoaded', () => {

    // --- KONSTANTEN UND VARIABLEN ---
    const A_LATTICE = 1.0;
    const C_LATTICE = 1.365 * 2; // c/a ratio for Hematite is ~2.73, often given as c/(a/sqrt(3)) = 1.365
    let presets = {};

    // --- DOM ELEMENTE ---
    const presetSelect = document.getElementById('preset-select');
    const calculateLastAxisBtn = document.getElementById('calculate-last-axis-btn');
    const transformToCrystalBtn = document.getElementById('transform-to-crystal-btn');
    const transformToLabBtn = document.getElementById('transform-to-lab-btn');

    const inputsDef = {
        x: {h: document.getElementById('x_h'), k: document.getElementById('x_k'), i: document.getElementById('x_i'), l: document.getElementById('x_l')},
        y: {h: document.getElementById('y_h'), k: document.getElementById('y_k'), i: document.getElementById('y_i'), l: document.getElementById('y_l')},
        z: {h: document.getElementById('z_h'), k: document.getElementById('z_k'), i: document.getElementById('z_i'), l: document.getElementById('z_l')},
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

    // --- HILFSFUNKTIONEN (MATHEMATIK) ---

    // Konvertiert Miller-Bravais [hkil] zu einem kartesischen Vektor [x,y,z]
    const mbToCartesian = (mb) => {
        const [h, k, _, l] = mb;
        const x = A_LATTICE * (2*h + k) / 2;
        const y = A_LATTICE * Math.sqrt(3) * k / 2;
        const z = C_LATTICE * l / 3;
        return [x, y, z];
    };

    // Konvertiert einen kartesischen Vektor [x,y,z] zu Miller-Bravais [hkil]
    const cartesianToMb = (cart) => {
        const [x, y, z] = cart;
        const k = (2 * y) / (A_LATTICE * Math.sqrt(3));
        const h = (2 * x - k * A_LATTICE) / (2 * A_LATTICE);
        const l = (3 * z) / C_LATTICE;
        const i = -(h + k);
        return [h, k, i, l];
    };

    // Vektor normalisieren
    const normalize = (v) => {
        const len = Math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        if (len === 0) return [0, 0, 0];
        return [v[0] / len, v[1] / len, v[2] / len];
    };

    // Kreuzprodukt
    const crossProduct = (a, b) => {
        return [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0]
        ];
    };

    // Matrix (3x3) mit Vektor (3x1) multiplizieren
    const multiplyMatrixVector = (m, v) => {
        return [
            m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
            m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
            m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2]
        ];
    };

    // Matrix transponieren
    const transposeMatrix = (m) => {
        return [
            [m[0][0], m[1][0], m[2][0]],
            [m[0][1], m[1][1], m[2][1]],
            [m[0][2], m[1][2], m[2][2]]
        ];
    };


    // --- UI-FUNKTIONEN ---

    // Aktualisiert den i-Index basierend auf h und k
    const updateI = (hInput, kInput, iInput) => {
        const h = parseFloat(hInput.value) || 0;
        const k = parseFloat(kInput.value) || 0;
        iInput.value = (-(h + k)).toFixed(4);
    };

    // Füllt die Felder für die Achsendefinition
    const populateDefinitionInputs = (axes) => {
        ['x', 'y', 'z'].forEach(axis => {
            const [h, k, i, l] = axes[axis];
            inputsDef[axis].h.value = h;
            inputsDef[axis].k.value = k;
            inputsDef[axis].i.value = i;
            inputsDef[axis].l.value = l;
        });
    };

    // Lädt ein Preset in die UI
    const loadPreset = (presetName) => {
        const preset = presets[presetName];
        if (preset) {
            populateDefinitionInputs(preset.axes);
        }
    };

    // Liest einen Miller-Bravais Vektor aus einer Gruppe von Inputs
    const getMbVectorFromInputs = (inputs) => {
        return [
            parseFloat(inputs.h.value) || 0,
            parseFloat(inputs.k.value) || 0,
            parseFloat(inputs.i.value) || 0,
            parseFloat(inputs.l.value) || 0
        ];
    };

    // Schreibt einen Miller-Bravais Vektor in eine Gruppe von Inputs
    const setMbVectorToInputs = (inputs, vec, precision = 4) => {
        inputs.h.value = vec[0].toFixed(precision);
        inputs.k.value = vec[1].toFixed(precision);
        inputs.i.value = vec[2].toFixed(precision);
        inputs.l.value = vec[3].toFixed(precision);
    };

    // Aktualisiert die Winkel basierend auf dem (x,y,z)-Vektor
    const updateAnglesFromXYZ = () => {
        const x = parseFloat(inputsLab.x.value) || 0;
        const y = parseFloat(inputsLab.y.value) || 0;
        const z = parseFloat(inputsLab.z.value) || 0;
        const r = Math.sqrt(x*x + y*y + z*z);
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

    // Aktualisiert den (x,y,z)-Vektor basierend auf den Winkeln
    const updateXYZFromAngles = () => {
        const thetaDeg = parseFloat(inputsLab.theta.value) || 0;
        const phiDeg = parseFloat(inputsLab.phi.value) || 0;
        const thetaRad = thetaDeg * Math.PI / 180;
        const phiRad = phiDeg * Math.PI / 180;

        inputsLab.x.value = (Math.sin(thetaRad) * Math.cos(phiRad)).toFixed(4);
        inputsLab.y.value = (Math.sin(thetaRad) * Math.sin(phiRad)).toFixed(4);
        inputsLab.z.value = (Math.cos(thetaRad)).toFixed(4);
    };

    // --- KERNLOGIK ---

    // Erstellt die Transformationsmatrix M (transformiert von Kristall-Kartesisch zu Labor)
    const getTransformationMatrix = () => {
        const v_x_mb = getMbVectorFromInputs(inputsDef.x);
        const v_y_mb = getMbVectorFromInputs(inputsDef.y);
        const v_z_mb = getMbVectorFromInputs(inputsDef.z);

        const v_x_cart = mbToCartesian(v_x_mb);
        const v_y_cart = mbToCartesian(v_y_mb);
        const v_z_cart = mbToCartesian(v_z_mb);

        const n_x = normalize(v_x_cart);
        const n_y = normalize(v_y_cart);
        const n_z = normalize(v_z_cart);

        // M hat die normalisierten kartesischen Vektoren als Spalten
        return [
            [n_x[0], n_y[0], n_z[0]],
            [n_x[1], n_y[1], n_z[1]],
            [n_x[2], n_y[2], n_z[2]]
        ];
    };

    const handleCalculateLastAxis = () => {
        const vectors = {
            x: { mb: getMbVectorFromInputs(inputsDef.x), isSet: getMbVectorFromInputs(inputsDef.x).some(v => v !== 0) },
            y: { mb: getMbVectorFromInputs(inputsDef.y), isSet: getMbVectorFromInputs(inputsDef.y).some(v => v !== 0) },
            z: { mb: getMbVectorFromInputs(inputsDef.z), isSet: getMbVectorFromInputs(inputsDef.z).some(v => v !== 0) },
        };

        const definedAxes = Object.keys(vectors).filter(ax => vectors[ax].isSet);
        const undefinedAxes = Object.keys(vectors).filter(ax => !vectors[ax].isSet);

        if (definedAxes.length !== 2 || undefinedAxes.length !== 1) {
            alert("Bitte definieren Sie genau zwei Achsen, um die dritte zu berechnen.");
            return;
        }

        const [v1_name, v2_name] = definedAxes;
        const v3_name = undefinedAxes[0];

        const v1_cart = mbToCartesian(vectors[v1_name].mb);
        const v2_cart = mbToCartesian(vectors[v2_name].mb);

        let v3_cart;

        if (v3_name === 'x') v3_cart = crossProduct(v1_cart, v2_cart);
        else if (v3_name === 'y') v3_cart = crossProduct(v2_cart, v1_cart);
        else v3_cart = crossProduct(v1_cart, v2_cart);

        const v3_mb = cartesianToMb(v3_cart);
        setMbVectorToInputs(inputsDef[v3_name], v3_mb);
    };

    const handleTransformToCrystal = () => {
        const M = getTransformationMatrix();
        const M_inv = transposeMatrix(M);

        const v_lab = [
            parseFloat(inputsLab.x.value) || 0,
            parseFloat(inputsLab.y.value) || 0,
            parseFloat(inputsLab.z.value) || 0,
        ];

        const v_cryst_cart = multiplyMatrixVector(M_inv, v_lab);
        const v_cryst_mb = cartesianToMb(v_cryst_cart);

        setMbVectorToInputs(inputsCrystal, v_cryst_mb);
    };

    const handleTransformToLab = () => {
        const M = getTransformationMatrix();
        const v_cryst_mb = getMbVectorFromInputs(inputsCrystal);
        const v_cryst_cart = mbToCartesian(v_cryst_mb);

        const v_lab = multiplyMatrixVector(M, normalize(v_cryst_cart));

        inputsLab.x.value = v_lab[0].toFixed(4);
        inputsLab.y.value = v_lab[1].toFixed(4);
        inputsLab.z.value = v_lab[2].toFixed(4);

        updateAnglesFromXYZ();
    };

    // --- INITIALISIERUNG ---
    const init = async () => {
        // Presets aus JSON laden
        try {
            const response = await fetch('presets.json');
            presets = await response.json();
        } catch (error) {
            console.error('Fehler beim Laden der Presets:', error);
            alert('Fehler: Die vordefinierten Schnittebenen konnten nicht geladen werden.');
            return;
        }

        // Presets in Dropdown laden
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

        // Listener für i-Index
        ['x', 'y', 'z'].forEach(axis => {
            inputsDef[axis].h.addEventListener('input', () => updateI(inputsDef[axis].h, inputsDef[axis].k, inputsDef[axis].i));
            inputsDef[axis].k.addEventListener('input', () => updateI(inputsDef[axis].h, inputsDef[axis].k, inputsDef[axis].i));
        });
        inputsCrystal.h.addEventListener('input', () => updateI(inputsCrystal.h, inputsCrystal.k, inputsCrystal.i));
        inputsCrystal.k.addEventListener('input', () => updateI(inputsCrystal.h, inputsCrystal.k, inputsCrystal.i));

        // Listener für Laborsystem-Inputs
        [inputsLab.x, inputsLab.y, inputsLab.z].forEach(el => el.addEventListener('input', updateAnglesFromXYZ));
        [inputsLab.theta, inputsLab.phi].forEach(el => el.addEventListener('input', updateXYZFromAngles));

        // Standard-Preset laden
        if (Object.keys(presets).length > 0) {
            loadPreset(Object.keys(presets)[0]);
        }
    };

    init();
});
