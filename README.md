# Hematite Crystal Calculator

Web-Tool zur Umrechnung und 3D-Darstellung von Kristallachsen, Sample-Orientierung und Magnetfeld-Richtungen für Hämatit ($\alpha\text{-Fe}_2\text{O}_3$) bei Spinwellen- (CPW) und SAW-Messungen (IDT).

Link zur Anwendung: [Hematite Crystal Calculator](https://toseflo.github.io/crystal-axis-calculator/)

## Überblick

- **Sample & Kristallschnitt**: Vorauswahl für `c-cut (0001)`, `a-cut (11-20)`, `m-cut (10-10)`, `r-cut (1-102)` oder beliebige Ebenen `(hkil)` inkl. frei wählbarer 0°-Bezugskante des Samples.
- **Struktur & Welle**: Drehwinkel von CPW oder IDT auf dem Sample und direkte Anzeige der zugehörigen Miller-Bravais-Richtung von $\hat{k}$.
- **Einbau & Messung im Magneten**: Horizontale (In-Plane) oder vertikale (Out-of-Plane) Halterung mit Einbauwinkel $\psi_S$ und drehbarem Magnetfeld $\phi_M$, inklusive Ziel-Finder für bestimmte Feld-, Easy-Plane- oder Kristallrichtungen.
- **3D-Ansicht** (Mitte, bleibt beim Scrollen stehen): frei drehbar, mit 2D-Draufsicht und optionaler zweiter Ansicht. „Was steht im Bild still?“ wählt das feste Bezugssystem (Probe, Magnet, Welle $\hat{k}$, Kristall); der Wechsel ist animiert, damit man sieht, wie die Systeme gegeneinander verdreht sind. Die Legende im Bild ist zugleich der Ein/Aus-Schalter der Ebenen.
- **Aufbau erklären**: geführte Tour Kristall → Probe → Struktur → Einbau → Feld drehen, gebaut aus dem aktuellen Setup.
- **Hervorhebung**: Fährt man mit der Maus über einen Begriff, eine Tabellenzeile oder einen Setup-Schritt, wird das zugehörige Objekt im 3D-Bild hervorgehoben (alles andere wird blass).
- **Übersicht** (rechts):
  - Ergebnis-Kacheln für $\angle(\mathbf{B},\hat{k})$, $\angle(\hat{m},\hat{k})$, $\angle(\hat{l},\hat{k})$.
  - „Kurz gesagt“: vier Zeilen zu Probe, Welle, Feld und Spins; die ausführliche Erklärung lässt sich aufklappen.
  - „Wo zeigt was hin?“: $\mathbf{B}$, $\hat{m}$, $\hat{l}$, $\hat{k}$ (auf Wunsch alle Achsen und Kanten) als Winkel $(\theta,\phi)$ oder $(x,y,z)$ im gewählten System plus $[uvtw]$.
  - „Eigene Richtung“ (aufklappbar): beliebige Richtung $[uvtw]$, Ebenennormale $(hkil)$, Kreuzprodukt zweier Richtungen/Ebenen oder Winkel $(\theta,\phi)$ eingeben → wird eingezeichnet, eingeordnet, und per Klick dreht der Ziel-Finder das Feld dorthin. Mit „+ Merken“ landen mehrere Richtungen in einer Liste – jede als eigener farbiger Pfeil im 3D-Bild (fest im Kristall, folgt also jeder Änderung am Aufbau), per × wieder entfernbar.
- **Link kopieren**: Das komplette Setup (inkl. Vergleich A/B) steht in der URL und lässt sich an Kollegen weitergeben.
- **Export**: SVG-Grafik, [`hematite_coords.py`](hematite_coords.py) für Python-Auswerteskripte und 3×3-Rotationsmatrizen zwischen allen Systemen (auch als LaTeX).

## Lokal starten

Da ES-Module genutzt werden, reicht ein lokaler Webserver im Projektordner:

```bash
python -m http.server 8000
# oder: npx http-server -p 8000
```

Danach `http://localhost:8000` im Browser öffnen.
