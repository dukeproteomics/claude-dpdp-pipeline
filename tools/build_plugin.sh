#!/usr/bin/env bash
# Build the shareable plugin marketplace from this project's sources:
#   pipeline/ (R scripts), .claude/skills/dpdp-* (skills), tools/plugin_src/ (manifest, hook, README)
# Output: dist/dpdp-marketplace/  (marketplace manifest + plugins/dpdp/). No data is copied.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DIST="$ROOT/dist/dpdp-marketplace"
P="$DIST/plugins/dpdp"
rm -rf "$DIST"
mkdir -p "$DIST/.claude-plugin" "$P/.claude-plugin" "$P/skills" "$P/hooks"
cp "$ROOT/tools/plugin_src/plugin.json" "$P/.claude-plugin/plugin.json"
cp "$ROOT/tools/plugin_src/README.md" "$P/README.md"
cp "$ROOT/tools/plugin_src/hooks/"* "$P/hooks/"
mkdir -p "$P/renv" && cp "$ROOT/tools/renv/renv.lock" "$ROOT/tools/renv/README.md" "$P/renv/"
rsync -a --exclude '__pycache__' --exclude '*.Rhistory' "$ROOT/pipeline/" "$P/pipeline/"
VERSION=$(python3 -c "import json;print(json.load(open('$P/.claude-plugin/plugin.json'))['version'])")
for d in "$ROOT"/.claude/skills/dpdp-*/; do
  name=$(basename "$d"); mkdir -p "$P/skills/$name"
  python3 - "$d/SKILL.md" "$P/skills/$name/SKILL.md" <<'PY'
import re, sys
src, dst = sys.argv[1], sys.argv[2]
s = open(src).read()
s = s.replace("cd /home/user/Claude_Shiny\n", "")
s = re.sub(r"Rscript pipeline/", "Rscript <pipeline>/", s)
s = re.sub(r"`pipeline/", "`<pipeline>/", s)
s = s.replace("data/<prefix>/", "<DPDP_DATA>/<prefix>/")
note = ("\n> Paths: `<pipeline>` is this plugin's `pipeline` folder, two levels above this skill's base directory "
        "(`<base>/../../pipeline`). `<DPDP_DATA>` is `$DPDP_DATA` (default `/home/user/Claude_Shiny/data`).\n")
s = re.sub(r"(\n# [^\n]+\n)", lambda m: m.group(1) + note, s, count=1)
open(dst, "w").write(s)
PY
done
cat > "$DIST/.claude-plugin/marketplace.json" <<JSON
{
  "name": "dpdp-marketplace",
  "owner": { "name": "Duke Proteomics" },
  "description": "Claude Code plugins for the Duke proteomics data processing pipeline (DPDP).",
  "plugins": [
    { "name": "dpdp", "source": "./plugins/dpdp", "version": "$VERSION",
      "description": "Claude-driven DPDP proteomics pipeline with a review page and an analyst decision at every step." }
  ]
}
JSON
echo "Built $DIST (dpdp $VERSION)"
