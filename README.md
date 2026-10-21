# Hematite Crystal Calculator

Web-Tool zur Umrechnung und 3D-Darstellung von Kristallachsen, Sample-Orientierung und Magnetfeld-Richtungen für Hämatit ($\alpha\text{-Fe}_2\text{O}_3$) bei Spinwellen- (CPW) und SAW-Messungen (IDT).

Link zur Anwendung: [Hematite Crystal Calculator](https://toseflo.github.io/crystal-axis-calculator/)

## Überblick

- **Sample & Kristallschnitt**: Vorauswahl für `c-cut (0001)`, `a-cut (11-20)`, `m-cut (10-10)`, `r-cut (1-102)` oder beliebige Ebenen `(hkil)` inkl. frei wählbarer 0°-Bezugskante des Samples.
- **Struktur & Welle**: Drehwinkel von CPW oder IDT auf dem Sample und direkte Anzeige der zugehörigen Miller-Bravais-Richtung von $\hat{k}$.
- **Einbau & Messung im Magneten**: Horizontale (In-Plane) oder vertikale (Out-of-Plane) Halterung mit Einbauwinkel $\psi_S$ und drehbarem Magnetfeld $\phi_M$, inklusive Ziel-Finder für bestimmte Feld-, Easy-Plane- oder Kristallrichtungen.
- **3D-Ansicht**: Frei drehbare 3D-Darstellung (mit einblendbarer zweiter Ansicht daneben und 2D-Draufsicht) sowie umschaltbarem festem Bezugssystem (Sample, Labor-Magnet, Welle $\hat{k}$, Easy-Plane).
- **Tabellen & Umrechner**:
  - Vektortabelle mit kartesischen Koordinaten $(x,y,z)$, Kugelwinkeln $(\theta,\phi)$ und $[uvtw]$-Richtungen.
  - Hin- und Rück-Umrechner zwischen Kristall-Indizes (`[uvtw]` / `(hkil)`) und Sample-/Labor-Koordinaten sowie Berechnung einer senkrechten 3. Kristallachse per Kreuzprodukt.
  - 3×3-Rotationsmatrizen zwischen den Systemen.
- **Export**: SVG-Export für Grafiken und [`hematite_coords.py`](hematite_coords.py) für Auswerteskripte in Python.

## Lokal starten

Da ES-Module genutzt werden, reicht ein lokaler Webserver im Projektordner:

```bash
python -m http.server 8000
```

Danach `http://localhost:8000` im Browser öffnen.
