# DEA Explorer per macOS

## Una tantum
1. **Python 3.11+ da python.org** (include Tk 8.6; evita il Python di sistema/Xcode).
2. **R** da https://cran.r-project.org/bin/macosx/ (build *arm64* per Apple Silicon, *x86_64* per Intel).
3. Pacchetti R:  `Rscript install_r_packages.R`

## Opzione A - eseguire dai sorgenti
```
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements_mac.txt
python main.py
```

## Opzione B - creare l'app (.app)
```
bash build_mac.sh
```
Risultato: `dist/DEA Explorer.app`. Se vuoi l'icona, metti `logo.png` (>=1024 px) in questa cartella prima del build.
L'app generata funziona solo sull'architettura del Mac che l'ha creata (arm64 o Intel).
R e i suoi pacchetti restano fuori dall'app e vanno installati su ogni Mac.
Su un altro Mac l'app non firmata va aperta con tasto destro > Apri (oppure `xattr -dr com.apple.quarantine "DEA Explorer.app"`).
