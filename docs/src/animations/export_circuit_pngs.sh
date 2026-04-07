#!/usr/bin/env zsh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

if [[ -x "./venv/bin/manim" ]]; then
  MANIM_BIN="./venv/bin/manim"
elif [[ -x "./.venv/bin/manim" ]]; then
  MANIM_BIN="./.venv/bin/manim"
elif command -v manim >/dev/null 2>&1; then
  MANIM_BIN="manim"
else
  echo "Manim executable not found in ./venv/bin/manim, ./.venv/bin/manim, or PATH"
  exit 1
fi

ASSET_DIR="./assets"
mkdir -p "$ASSET_DIR"

render_png() {
  local scene_name="$1"
  local output_name="$2"

  "$MANIM_BIN" -qh -s --format=png circuit_png_render.py "$scene_name" -o "$output_name"

  local generated
  generated="$(find ./media -type f -name "$output_name" | head -n 1)"
  if [[ -z "$generated" ]]; then
    echo "Failed to locate generated PNG for $scene_name"
    exit 1
  fi

  cp "$generated" "$ASSET_DIR/$output_name"
}

render_png "QFTCircuitPNG" "qft_circuit.png"
render_png "DTCircuitPNG" "dt_circuit.png"
render_png "ZTCircuitPNG" "zt_circuit.png"

for png in "$ASSET_DIR/qft_circuit.png" "$ASSET_DIR/dt_circuit.png" "$ASSET_DIR/zt_circuit.png"; do
  if [[ ! -s "$png" ]]; then
    echo "PNG validation failed (empty file): $png"
    exit 1
  fi

  if ! file "$png" | grep -qi "PNG image data"; then
    echo "PNG validation failed (not a PNG image): $png"
    exit 1
  fi

  echo "Validated $png"
done

echo "PNG export complete:"
ls -lh "$ASSET_DIR/qft_circuit.png" "$ASSET_DIR/dt_circuit.png" "$ASSET_DIR/zt_circuit.png"
