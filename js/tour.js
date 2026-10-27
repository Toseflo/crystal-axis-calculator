/**
 * tour.js
 * "Aufbau Schritt für Schritt": a short guided walk through the chain
 * Kristall → Probe → Struktur → Einbau → Feld drehen, built from the CURRENT setup.
 * Each stage switches the visible layers and the frame that stands still,
 * so the (animated) frame transitions show how the frames are rotated against each other.
 */

const ALL_OFF = {
    showB: false,
    showBProjections: false,
    showM: false,
    showL: false,
    showK: false,
    showEasyPlane: false,
    showCrystalAxes: false,
    showLabMagnet: false,
    showSampleSlab: false,
    showMicrostructure: false,
    showCustom: false,
    showT: false,
    showN: false,
    showEdges: false
};

const fmt = (x, d = 1) => (Math.abs(x) < 0.5 * Math.pow(10, -d) ? 0 : x).toFixed(d);

const STAGES = [
    {
        title: 'Kristall',
        lock: 'Cryst',
        vis: { showEasyPlane: true, showCrystalAxes: true },
        text: (st, tex) =>
            `Hämatit ist trigonal. Die drei ${tex('a')}-Achsen ${tex('\\langle 11\\bar{2}0\\rangle')} liegen alle 60° in der Basalebene (0001). `
            + `Oberhalb der Morin-Temperatur ist das die <strong>Easy-Plane</strong> (lila Sechseck): die Spins liegen darin. `
            + `Die ${tex('c')}-Achse ${tex('[0001]')} steht senkrecht darauf, die ${tex('m')}-Achsen liegen jeweils 30° zwischen den ${tex('a')}-Achsen.`
    },
    {
        title: 'Probe',
        lock: 'Sub',
        vis: { showEasyPlane: true, showCrystalAxes: true, showSampleSlab: true, showEdges: true, showN: true },
        text: (st, tex) => {
            const tilt = Math.min(st.metrics.angleNCDeg, 180 - st.metrics.angleNCDeg);
            return `Aus dem Kristall wird eine Probe mit der Oberfläche ${tex(st.millerBravais.zsPlane.latex)} geschnitten. `
                + `Das Bild hat sich gerade so gedreht, dass jetzt die <strong>Probe</strong> flach liegt. `
                + `Die Easy-Plane ist ${fmt(tilt)}° gegen die Oberfläche gekippt. `
                + `Kanten: ${tex('x_s \\parallel')} ${tex(st.millerBravais.xs.latex)}, ${tex('y_s \\parallel')} ${tex(st.millerBravais.ys.latex)}.`;
        }
    },
    {
        title: 'Struktur',
        lock: 'Sub',
        vis: { showEasyPlane: true, showSampleSlab: true, showMicrostructure: true, showK: true, showT: true },
        text: (st, tex) =>
            `Auf die Probe kommt ${st.setup.deviceType === 'IDT' ? 'der IDT' : 'die CPW'}, um ${fmt(st.setup.thetaDevDeg)}° zur 0°-Kante gedreht. `
            + `Die Welle läuft entlang ${tex('\\hat{k} \\parallel')} ${tex(st.millerBravais.k.latex)}, `
            + `Leiter bzw. Finger entlang ${tex('\\hat{t} \\parallel')} ${tex(st.millerBravais.t.latex)}. `
            + `${tex('\\hat{k}')} steht ${fmt(st.metrics.angleKCDeg)}° zur ${tex('c')}-Achse.`
    },
    {
        title: 'Einbau',
        lock: 'Lab',
        vis: { showEasyPlane: true, showSampleSlab: true, showMicrostructure: true, showK: true, showLabMagnet: true, showN: true },
        text: (st, tex) =>
            (st.setup.mountMode === 'IP'
                ? 'Jetzt steht der <strong>Magnet</strong> still. Die Probe liegt flach zwischen den Polen, ihre Normale zeigt nach oben. '
                : 'Jetzt steht der <strong>Magnet</strong> still. Die Probe steht aufrecht, ihre Normale zeigt bei '
                    + `${tex('\\phi_M = 0^\\circ')} genau auf einen Pol. `)
            + `${tex(`\\psi_S = ${fmt(st.setup.psiSDeg)}^\\circ`)} dreht die Probe auf dem Halter um ihre Normale.`
    },
    {
        title: 'Feld drehen',
        lock: 'Sub',
        sweep: true,
        vis: {
            showEasyPlane: true,
            showSampleSlab: true,
            showMicrostructure: true,
            showK: true,
            showLabMagnet: true,
            showB: true,
            showBProjections: true,
            showM: true,
            showL: true
        },
        text: (st, tex) =>
            `Der Magnet dreht sich (${tex('\\phi_M')} = <strong class="tour-phi">${fmt(st.setup.phiMDeg)}°</strong>). `
            + `${tex('\\mathbf{B}')} (cyan) dreht mit; ${tex('\\hat{m}')} (rosa) folgt der Projektion von ${tex('\\mathbf{B}')} in die Easy-Plane `
            + `(gestrichelte Lotlinie), ${tex('\\hat{l}')} (grün) steht senkrecht dazu. `
            + `Gerade liegen <strong class="tour-easy">${fmt(st.metrics.bEasyRatio * 100, 0)} %</strong> des Feldes in der Easy-Plane.`
    }
];

export class GuidedTour {
    /**
     * @param {object} opts
     *   scene      Scene3DController of view 1
     *   getSetup   () => active setup object
     *   evaluate   evaluateSetup
     *   tex        KaTeX helper
     *   refresh    () => void; re-renders the app from appState
     *   onVisibilityChange () => void
     *   els        { root, title, text, step, prev, next, exit }
     */
    constructor(opts) {
        Object.assign(this, opts);
        this.index = -1;
        this._saved = null;
        this._raf = null;
        this.els.prev.addEventListener('click', () => this.go(this.index - 1));
        this.els.next.addEventListener('click', () => {
            if (this.index >= STAGES.length - 1) this.stop();
            else this.go(this.index + 1);
        });
        this.els.exit.addEventListener('click', () => this.stop());
    }

    get active() {
        return this.index >= 0;
    }

    start() {
        if (!this.active) {
            this._saved = { visibility: { ...this.scene.visibility }, lockFrame: this.scene.lockFrame };
        }
        this.els.root.classList.remove('hidden');
        this.go(0);
    }

    stop() {
        if (!this.active) return;
        this._stopSweep();
        this.index = -1;
        this.els.root.classList.add('hidden');
        Object.assign(this.scene.visibility, this._saved.visibility);
        this.scene.setLockFrame(this._saved.lockFrame, true);
        this._saved = null;
        this.refresh();
        this.onVisibilityChange();
    }

    go(i) {
        if (i < 0 || i >= STAGES.length) return;
        this._stopSweep();
        this.index = i;
        const stage = STAGES[i];
        Object.assign(this.scene.visibility, ALL_OFF, stage.vis);
        this.scene.setLockFrame(stage.lock, true);
        this.onVisibilityChange();
        this.refresh();

        this.els.title.textContent = `${i + 1}/${STAGES.length} · ${stage.title}`;
        this.els.step.innerHTML = STAGES.map(
            (_, k) => `<span class="tour-dot${k === i ? ' active' : ''}${k < i ? ' done' : ''}"></span>`
        ).join('');
        this.els.prev.disabled = i === 0;
        this.els.next.textContent = i === STAGES.length - 1 ? 'Fertig' : 'Weiter →';
        this.els.text.innerHTML = stage.text(this.evaluate(this.getSetup()), this.tex);

        if (stage.sweep) this._startSweep(stage);
    }

    // Rotates phi_M on a temporary copy of the setup (appState stays untouched).
    _startSweep(stage) {
        const base = { ...this.getSetup() };
        const t0 = performance.now();
        const phi0 = base.phiMDeg || 0;
        let lastText = 0;
        const tick = (now) => {
            const phi = (phi0 + ((now - t0) / 1000) * 40) % 360;
            const st = this.evaluate({ ...base, phiMDeg: phi });
            this.scene.updateStates(st, null, false);
            if (now - lastText > 120) {
                this.els.text.innerHTML = stage.text(st, this.tex);
                lastText = now;
            }
            this._raf = requestAnimationFrame(tick);
        };
        this._raf = requestAnimationFrame(tick);
    }

    _stopSweep() {
        if (this._raf) cancelAnimationFrame(this._raf);
        this._raf = null;
    }
}
