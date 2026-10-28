# Hematite Crystal Calculator

Web-Tool zur Umrechnung und 3D-Darstellung von Kristallachsen, Sample-Orientierung und Magnetfeld-Richtungen für Hämatit ($\alpha\text{-Fe}_2\text{O}_3$) bei Spinwellen- (CPW) und SAW-Messungen (IDT).

Link zur Anwendung: [Hematite Crystal Calculator](https://toseflo.github.io/crystal-axis-calculator/)

## Überblick

Das Setup steht links in drei Tabs (Probe · Struktur · Einbau), jeder mit einer kurzen Zusammenfassung im Reiter. Darunter ist das Magnetfeld immer sichtbar.

- **Sample & Kristallschnitt**: Vorauswahl für `c-cut (0001)`, `a-cut (11-20)`, `m-cut (10-10)`, `r-cut (1-102)` oder beliebige Ebenen `(hkil)` inkl. frei wählbarer 0°-Bezugskante des Samples.
- **Struktur & Welle**: Drehwinkel von CPW oder IDT auf dem Sample und direkte Anzeige der zugehörigen Miller-Bravais-Richtung von $\hat{k}$.
- **Einbau im Magneten**: Horizontale (In-Plane) oder vertikale (Out-of-Plane) Halterung mit Einbauwinkel $\psi_S$.
- **Magnetfeld**: $\phi_M$ und $B_0$. Unter dem Regler stehen die besonderen Feldrichtungen (längs, gegen oder quer zur Welle, senkrecht zur Probe, $a$-, $m$-, $c$-Achse), die mit diesem Einbau gehen, jeweils mit Winkel und einem Strich auf dem Regler; ein Klick dreht das Feld dorthin, beim Überfahren zeigt die 3D-Ansicht gestrichelt „B neu“. Richtungen, die nur mit anderem Einbau gehen, stehen in einer Zeile darunter; ein Klick darauf nennt den Einbau, mit dem es geht (Probe flach oder stehend, um wie viel gedreht; „So einbauen“).
- **3D-Ansicht** (Mitte, bleibt beim Scrollen stehen): frei drehbar, mit 2D-Draufsicht und optionaler zweiter Ansicht. „Was steht im Bild still?“ wählt das feste Bezugssystem (Probe, Magnet, Welle $\hat{k}$, Kristall); der Wechsel ist animiert, damit man sieht, wie die Systeme gegeneinander verdreht sind. Die Legende im Bild ist zugleich der Ein/Aus-Schalter der Ebenen.
- **Aufbau erklären**: geführte Tour Kristall → Probe → Struktur → Einbau → Feld drehen, gebaut aus dem aktuellen Setup.
- **Hervorhebung**: Fährt man mit der Maus über einen Begriff, eine Tabellenzeile oder einen Setup-Tab, wird das zugehörige Objekt im 3D-Bild hervorgehoben (alles andere wird blass).
- **Übersicht** (rechts):
  - Ergebnis-Kacheln für $\angle(\mathbf{B},\hat{k})$, $\angle(\hat{m},\hat{k})$, $\angle(\hat{l},\hat{k})$ plus Spinwellen-Geometrie: Backward-Volume ($\hat{m}\parallel\hat{k}$), Damon-Eshbach ($\hat{m}\perp\hat{k}$, beide in der Probe), Forward-Volume ($\hat{m}\perp$ Probe), je ±10°. Beim CPW zeigt der Tooltip das Drehmoment beider Feldkomponenten auf $\hat{m}$: $h_x\parallel\hat{k}$ → $|\hat{k}\times\hat{m}|$, $h_z\perp$ Probe → $|\hat{n}\times\hat{m}|$.
  - „Kurz gesagt“: vier Zeilen zu Probe, Welle, Feld und Spins; die ausführliche Erklärung lässt sich aufklappen.
  - „Wo zeigt was hin?“: $\mathbf{B}$, $\hat{m}$, $\hat{l}$, $\hat{k}$ (auf Wunsch alle Achsen und Kanten) als Winkel $(\theta,\phi)$ oder $(x,y,z)$ im gewählten System plus $[uvtw]$.
  - „Eigene Richtung“ (aufklappbar): beliebige Richtung $[uvtw]$, Ebenennormale $(hkil)$, Kreuzprodukt zweier Richtungen/Ebenen oder Winkel $(\theta,\phi)$ eingeben → wird eingezeichnet, eingeordnet, und per Klick dreht sich das Feld dorthin. Mit „+ Merken“ landen mehrere Richtungen in einer Liste – jede als eigener farbiger Pfeil im 3D-Bild (fest im Kristall, folgt also jeder Änderung am Aufbau), per × wieder entfernbar.
- **Teilen & Exportieren** (Karte rechts unten):
  - **Link kopieren**: Das komplette Setup (inkl. Vergleich A/B und aktivem Tab), die Ansicht und die gemerkten Richtungen stehen in der URL (z. B. `…&pins=1,0,-1,1;0,0,0,1`) und lassen sich an Kollegen weitergeben.
  - **Notiz kopieren**: Setup als kurzer Text fürs Laborbuch (Probe, Halter, Struktur, Feld, Spins, Winkel, gemerkte Richtungen, Link); beim Vergleich A und B untereinander.
  - **Bild (PNG)**: Ansicht 1 so, wie sie gerade gedreht ist, mit Vorschau – Format (wie Ansicht, 16:9, 4:3, quadratisch), helle oder dunkle Farben, transparenter Hintergrund, Auflösung 1×/2×/4×; herunterladen oder direkt in die Zwischenablage (für PowerPoint/Teams).
  - **Grafik (SVG)**: einheitliche Abbildung für die Abschlussarbeit – Probe und Feld immer aus derselben flachen Sicht (schräg oder Draufsicht) und im selben Maßstab, eng zugeschnitten. Drauf sind die rechteckige Probe mit Sender und Empfänger in der Form der 3D-Ansicht (zwei wie im Setup, nur Sender oder gemischt IDT/Antenne), $\hat{k}$ vom Sender zum Empfänger, $\mathbf{B}$ ab der Mitte von $\hat{k}$ mit Winkel in der Ebene (ab $\hat{k}$ oder ab der 0°-Kante) und Kippwinkel, der Antennenwinkel am Sender, daneben die Miller-Bravais-Richtungen von Kanten und Normale (beim r-cut die Ebene, da ihre Normale keine einfache Richtung ist). Winkelwerte stehen direkt am Bogen. Sweeps: „Feld drehen“ zeichnet die Bahn der Feldspitze in der Drehebene des Magneten mit Drehsinn und Drehwinkel ($\phi_M$ von/bis, optional blasse Zwischenpfeile; unter der Probe gestrichelt), „Feldstärke“ schreibt den Bereich ans $\mathbf{B}$ (über null mit blassem Gegenpfeil). Schrift und Linien sind Druckgrößen (pt) bei der gewählten Probenbreite; Beschriftungen suchen sich freie Stellen statt weißer Ränder; Teile liegen in eigenen Inkscape-Ebenen, optional mit LaTeX-Labels (Inkscape → PDF+LaTeX). Stile (Farben, Linien, Schrift, Symbole) lassen sich unter eigenem Namen speichern und als JSON weitergeben.
  - **Feld-Sweep**: $\phi_M = 0\ldots360^\circ$ als Plot (Winkel zu $\hat{k}$, Geometrie-Streifen, Easy-Plane-Anteil; Klick setzt $\phi_M$) und als CSV in 1°-Schritten inkl. $\hat{m}$, $\hat{l}$ im Wellen-System – direkt für Origin/Excel/Python.
  - Für Auswerteskripte: Python-Code für [`hematite_coords.py`](hematite_coords.py) und 3×3-Rotationsmatrizen zwischen allen Systemen (auch als LaTeX).

## Lokal starten

Da ES-Module genutzt werden, reicht ein lokaler Webserver im Projektordner:

```bash
python -m http.server 8000
# oder: npx http-server -p 8000
```

Danach `http://localhost:8000` im Browser öffnen.
