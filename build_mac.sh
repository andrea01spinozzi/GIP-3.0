#!/bin/bash
# Crea "DEA Explorer.app" (da eseguire su un Mac, dentro questa cartella):  bash build_mac.sh
set -e
cd "$(dirname "$0")"

python3 -m venv .venv_build
source .venv_build/bin/activate
pip install --upgrade pip
pip install -r requirements_mac.txt pyinstaller

ICON_ARGS=()
DATA_ARGS=(--add-data "DE_scheletro_FINALE.Rmd:." --add-data "run_pipeline.R:.")
if [ -f logo.png ]; then
  DATA_ARGS+=(--add-data "logo.png:.")
  rm -rf icon.iconset && mkdir icon.iconset
  for s in 16 32 128 256 512; do
    sips -z $s $s logo.png --out "icon.iconset/icon_${s}x${s}.png" >/dev/null
    sips -z $((s*2)) $((s*2)) logo.png --out "icon.iconset/icon_${s}x${s}@2x.png" >/dev/null
  done
  iconutil -c icns icon.iconset -o icon.icns
  ICON_ARGS=(--icon icon.icns)
fi

pyinstaller --noconfirm --clean --windowed --name "DEA Explorer" \
  --osx-bundle-identifier "org.deaexplorer.app" \
  "${ICON_ARGS[@]}" "${DATA_ARGS[@]}" \
  --collect-all tkinterdnd2 --hidden-import core_calcolo --hidden-import gui \
  main.py

echo
echo "Fatto: dist/DEA Explorer.app  (trascinala in Applicazioni)"
