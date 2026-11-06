#!/bin/bash
# Runs in the web app's terminal as the logged-in user (started by web/server.py).
# First time: installs Claude Code for this user, then Claude asks them to log in with their own Claude account.
# Every time: adds/updates the shared dpdp plugin, then starts Claude Code in ~/dpdp (continuing the last conversation).
export PATH="$HOME/.local/bin:$PATH"
W="${DPDP_WEB_WORK:-$HOME/dpdp}"
mkdir -p "$W/.claude" && cd "$W" || exit 1
y=$'\e[33m'; g=$'\e[32m'; n=$'\e[0m'

if ! command -v claude >/dev/null 2>&1; then
  echo "${y}First time: installing Claude Code for $USER (about a minute)...${n}"
  if ! curl -fsSL https://claude.ai/install.sh | bash; then
    echo "${y}Claude Code could not be installed. Ask the admin. A plain shell follows.${n}"; exec bash -l
  fi
  export PATH="$HOME/.local/bin:$PATH"
fi

# The shared plugin (skills + pipeline scripts) from the server's marketplace folder
if [ -d "$DPDP_PLUGIN" ]; then
  if ! claude plugin list 2>/dev/null | grep -q "dpdp@dpdp-marketplace"; then
    echo "Adding the DPDP plugin..."
    claude plugin marketplace add "$DPDP_PLUGIN" >/dev/null 2>&1
    claude plugin install dpdp@dpdp-marketplace 2>&1 | tail -1
  else
    claude plugin marketplace update dpdp-marketplace >/dev/null 2>&1
    claude plugin update dpdp@dpdp-marketplace >/dev/null 2>&1
  fi
fi

# Folder settings: the shared data folder is part of the workspace, pipeline commands don't need a prompt each time
if [ ! -f "$W/.claude/settings.json" ]; then
  cat > "$W/.claude/settings.json" <<JSON
{
  "permissions": {
    "additionalDirectories": ["$DPDP_DATA"],
    "allow": ["Bash(Rscript:*)", "Bash(bash */pipeline/run_step.sh:*)", "Bash(python3 */pipeline/wait_decision.py:*)"]
  }
}
JSON
fi
cat > "$W/CLAUDE.md" <<MD
# DPDP web app session
This Claude Code session runs inside the DPDP web app (DPDP_WEB=1). The analyst sees the review pages in the web page's left panel.
- Run folders: \$DPDP_DATA ($DPDP_DATA), shared by the team. Uploaded inputs: \$DPDP_DATA/input/<name>/.
- Don't start review_server.py and don't xdg-open pages: the web app serves them. Tell the analyst the page is ready in the
  left panel (it opens there by itself) and run wait_decision.py in the background as usual.
- Use the dpdp-run skill for a project from start to end.
MD

echo "${g}DPDP: starting Claude Code as $USER in $W  (data: $DPDP_DATA)${n}"
proj="$HOME/.claude/projects/$(echo "$W" | sed 's#[/.]#-#g')"
if ls "$proj"/*.jsonl >/dev/null 2>&1; then claude --continue || claude; else claude; fi
echo; echo "${y}Claude exited. Type 'claude' to start it again, or click Restart Claude.${n}"
exec bash -l
