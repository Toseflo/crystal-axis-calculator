# Hematite ($\alpha\text{-Fe}_2\text{O}_3$) Spin-Wave (CPW) & SAW (IDT) Geometry Calculator

Modernes, interaktives 3D/2D-Geometrie- und Koordinatentransformations-Tool für Spinwellen- (CPW) und akustische Oberflächenwellen-Experimente (SAW / IDT) auf Hämatit ($\alpha\text{-Fe}_2\text{O}_3$, Raumgruppe $R\bar{3}c$) in einem horizontal drehbaren Elektromagneten.

- **Web-App (GitHub Pages & lokal)**: [`index.html`](index.html) + modulare ES6-Module in [`js/`](js/)
- **Python-Begleitmodul für Auswerteskripte**: [`hematite_coords.py`](hematite_coords.py) (identische $3\times 3$ $SO(3)$-Rotationsmatrizen, Easy-Plane-Gleichgewicht und $0^\circ\dots 360^\circ$ Magnet-Sweeps)

---

## 1. Die 5 verknüpften Koordinatensysteme

Alle Systeme sind als **rechtshändige Orthonormalsysteme** ($R \in SO(3)$, $\det R = +1$, $R^{-1} = R^T$) definiert:

1. **Labor-System $\mathcal{F}_{\text{Lab}} = (\hat{X}_L, \hat{Y}_L, \hat{Z}_L)$**:
   - Der Elektromagnet dreht sich horizontal in der Tischebene $(\hat{X}_L, \hat{Y}_L)$ um den Winkel $\phi_M$, während $\hat{Z}_L$ vertikal nach oben zeigt:
     $$\mathbf{B}_{\text{Lab}}(\phi_M) = B_0 \begin{pmatrix} \cos\phi_M \\ \sin\phi_M \\ 0 \end{pmatrix}$$

2. **Proben-Halterung $\mathcal{F}_{\text{Mount}} \to \mathcal{F}_{\text{Sub}}$**:
   - **OOP-Konfiguration (Vertikal)**: Probe steht vertikal, sodass bei $\phi_M = 0^\circ$ das Magnetfeld exakt Out-of-Plane entlang der Probennormale zeigt ($\hat{z}_s = \hat{X}_L$, $\hat{x}_m = \hat{Y}_L$, $\hat{y}_m = \hat{Z}_L$). In dieser Halterung lässt sich die Probe um ihre Oberflächennormale $\hat{z}_s$ um $\psi_S$ drehen:
     $$R_{\text{Lab}\leftarrow\text{Sub}}^{\text{OOP}}(\psi_S) = \begin{pmatrix} 0 & 0 & 1 \\ \cos\psi_S & -\sin\psi_S & 0 \\ \sin\psi_S & \cos\psi_S & 0 \end{pmatrix} \implies \mathbf{B}_{\text{Sub}}^{\text{OOP}}(\phi_M, \psi_S) = B_0 \begin{pmatrix} \sin\phi_M\cos\psi_S \\ -\sin\phi_M\sin\psi_S \\ \cos\phi_M \end{pmatrix}$$
   - **IP-Konfiguration (Horizontal)**: Probe liegt horizontal in der Magnetebene ($\hat{z}_s = \hat{Z}_L$), drehbar um $\psi_S$:
     $$R_{\text{Lab}\leftarrow\text{Sub}}^{\text{IP}}(\psi_S) = R_z(\psi_S) \implies \mathbf{B}_{\text{Sub}}^{\text{IP}}(\phi_M, \psi_S) = B_0 \begin{pmatrix} \cos(\phi_M - \psi_S) \\ \sin(\phi_M - \psi_S) \\ 0 \end{pmatrix}$$

3. **Substrat-System $\mathcal{F}_{\text{Sub}} = (\hat{x}_s, \hat{y}_s, \hat{z}_s)$**:
   - Probenkanten $(\hat{x}_s, \hat{y}_s)$ und äußere Oberflächennormale $\hat{z}_s = \hat{n}$ (inkl. optionalem Wafer-Flat/Kanten-Offset $\alpha_{\text{flat}}$).

4. **Mikrostruktur-/Wellen-System $\mathcal{F}_{\text{Wave}} = (\hat{k}, \hat{t}, \hat{n})$**:
   - Gedreht um den Strukturwinkel $\theta_{\text{dev}}$ um $\hat{n}$ auf der Probenoberfläche ($R_{\text{Sub}\leftarrow\text{Wave}} = R_z(\theta_{\text{dev}})$):
     - $\hat{k}$: Wellen-Ausbreitungsrichtung (SAW oder Spinwelle)
     - $\hat{t} = \hat{n} \times \hat{k}$: Transversalrichtung (parallel zu IDT-Fingern bzw. CPW-Innenleiter)
     - $\hat{n} = \hat{z}_s$: Oberflächennormale

5. **Hämatit-Kristallsystem $\mathcal{F}_{\text{Cryst}} = (\hat{e}_a, \hat{e}_m, \hat{e}_c)$**:
   - $\alpha\text{-Fe}_2\text{O}_3$ ($R\bar{3}c$, $a = 5.0356\text{ \AA}$, $c = 13.7489\text{ \AA}$):
     $$\hat{e}_a \parallel [11\bar{2}0] \text{ ($C_2$-Achse)}, \quad \hat{e}_m \parallel [1\bar{1}00], \quad \hat{e}_c \parallel [0001] \text{ ($C_3$-Achse / DMI-Vektor $\mathbf{D}$)}$$
   - **Unterstützte Cut-Presets**:
     - **c-cut $(0001)$**: $\theta_{\text{cut}} = 0^\circ$ ($\hat{n} = \hat{e}_c$, Basalebene = Probenebene)
     - **a-cut $(11\bar{2}0)$**: $\theta_{\text{cut}} = 90^\circ, \phi_{\text{cut}} = 0^\circ$ ($\hat{n} = \hat{e}_a$, $\hat{x}_s^{(0)} = \hat{e}_m \parallel [1\bar{1}00]$, $\hat{y}_s^{(0)} = \hat{e}_c \parallel [0001]$)
     - **m-cut $(10\bar{1}0) / (1\bar{1}00)$**: $\theta_{\text{cut}} = 90^\circ, \phi_{\text{cut}} = 90^\circ$ ($\hat{n} = \hat{e}_m$, $\hat{x}_s^{(0)} = -\hat{e}_a \parallel [\bar{1}\bar{1}20]$, $\hat{y}_s^{(0)} = \hat{e}_c \parallel [0001]$)
     - **r-cut $(1\bar{1}02)$**: $\theta_{\text{cut}} = \arctan\left(\frac{c}{\sqrt{3}a}\right) \approx 57.61^\circ, \phi_{\text{cut}} = 90^\circ$
     - **Custom Cut**: Beliebige Winkel $(\theta_{\text{cut}}, \phi_{\text{cut}})$ oder direkte Eingabe einer Miller-Bravais-Ebene $(h\,k\,i\,l)$

---

## 2. Easy-Plane-Gleichgewicht oberhalb des Morin-Übergangs ($T > T_M$)

Durch die Dzyaloshinskii-Moriya-Wechselwirkung ($\mathbf{D} \parallel \hat{e}_c$) liegt die magnetisch leichte Ebene in der Basalebene $(0001)$ ($\perp \hat{e}_c$):
- **Aktive Feldkomponente in der Easy-Plane**:
  $$\mathbf{B}_{\perp c} = \mathbf{B} - (\mathbf{B}\cdot\hat{e}_c)\hat{e}_c, \quad B_{\perp c} = \sqrt{B_a^2 + B_m^2}$$
- **Schwaches ferromagnetisches Nettomoment $\hat{m}$** und **Néel-Vektor $\hat{l}$**:
  $$\hat{m} = \frac{\mathbf{B}_{\perp c}}{|\mathbf{B}_{\perp c}|}, \quad \hat{l} = \hat{m} \times \hat{e}_c$$

---

## 3. Lokale Nutzung & Python-Begleitmodul

### Web-App lokal starten
Da ES6-Module genutzt werden, einfach einen lokalen HTTP-Server im Projektordner starten:
```bash
python -m http.server 8000
# oder: npx serve .
```
und `http://localhost:8000` im Browser öffnen. Auf **GitHub Pages** funktioniert die Anwendung direkt statisch ohne Build-Step.

### Python-Modul (`hematite_coords.py`)
In der Web-App kann über den Button **„🐍 Python-Snippet“** jederzeit der fertige Python-Code für das aktuell eingestellte Setup kopiert werden:
```python
import numpy as np
from hematite_coords import HematiteSetup

setup = HematiteSetup(
    name="r-cut OOP",
    mount_mode="OOP",
    cut_preset="r-cut",
    theta_dev_deg=35.0,
    psi_s_deg=25.0,
    phi_m_deg=70.0,
    b0_tesla=0.5,
)

state = setup.evaluate()
print("B in F_Wave (k, t, n):", state["vectors"]["B"]["Wave"])
print("Winkel <(l, k):", state["metrics"]["angle_l_k_deg"])

# Vektorisierter 0..360 deg Magnet-Sweep:
sweep = setup.sweep_magnet(np.linspace(0, 360, 361))
```
