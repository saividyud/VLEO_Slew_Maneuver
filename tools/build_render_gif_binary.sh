#!/bin/sh
set -eu

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
DIST_DIR="$SCRIPT_DIR/.dist"
BUILD_DIR="$SCRIPT_DIR/.build"
OUTPUT_BIN="$SCRIPT_DIR/bin/render_gif"

rm -rf "$DIST_DIR" "$BUILD_DIR"

uv run \
  --with pyinstaller \
  --with numpy \
  --with scipy \
  --with matplotlib \
  --with imageio \
  pyinstaller \
  --noconfirm \
  --clean \
  --onefile \
  --exclude-module tkinter \
  --exclude-module _tkinter \
  --copy-metadata imageio \
  --copy-metadata numpy \
  --copy-metadata scipy \
  --copy-metadata matplotlib \
  --name render_gif \
  --distpath "$DIST_DIR" \
  --workpath "$BUILD_DIR/pyinstaller" \
  --specpath "$BUILD_DIR/spec" \
  "$SCRIPT_DIR/render_gif.py"

install -m 755 "$DIST_DIR/render_gif" "$OUTPUT_BIN"
rm -rf "$DIST_DIR" "$BUILD_DIR"

printf 'Built native render_gif binary at: %s\n' "$OUTPUT_BIN"
