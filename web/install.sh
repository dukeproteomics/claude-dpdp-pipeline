#!/usr/bin/env bash
# Install or update the DPDP web app on this server. Run from the repo:   sudo bash web/install.sh
#   https://<server>/dpdp/  (Apache, self-signed certificate)  ->  127.0.0.1:8770 (systemd service dpdp-web, runs as root)
# Re-run after a code change: it rebuilds the plugin, copies the code, and restarts the service (open Claude sessions end).
set -euo pipefail
[ "$(id -u)" = 0 ] || { echo "Run with sudo: sudo bash web/install.sh"; exit 1; }
REPO="$(cd "$(dirname "$0")/.." && pwd)"
SRV=/srv/dpdp
GROUP="${DPDP_WEB_GROUP:-rstudio-users}"
PORT="${DPDP_WEB_PORT:-8770}"
IP="${DPDP_WEB_IP:-$(hostname -I | awk '{print $1}')}"
OWNER="${SUDO_USER:-root}"
getent group "$GROUP" >/dev/null || { echo "Group $GROUP not found"; exit 1; }
if ss -ltn | grep -q "127.0.0.1:$PORT " && ! systemctl is-active -q dpdp-web; then
  echo "Port $PORT is in use by something else (a --dev test server?). Stop it first."; exit 1
fi
echo "== DPDP web app: $SRV, https://$IP/dpdp/, group $GROUP"

# 1. Python environment (aiohttp + python-pam)
if [ ! -x "$SRV/venv/bin/python" ]; then
  python3 -m venv "$SRV/venv" 2>/dev/null || { apt-get install -y python3-venv; rm -rf "$SRV/venv"; python3 -m venv "$SRV/venv"; }
fi
"$SRV/venv/bin/pip" install -q --upgrade aiohttp python-pam six

# 2. Code: web app, plugin (rebuilt from this repo), app at the pinned submodule commit
mkdir -p "$SRV/web" "$SRV/shiny_DPDP" "$SRV/data/input"
rsync -a --delete --exclude '__pycache__' "$REPO/web/" "$SRV/web/"
sudo -u "$OWNER" bash "$REPO/tools/build_plugin.sh"
rsync -a --delete "$REPO/dist/dpdp-marketplace/" "$SRV/dpdp-marketplace/"
rm -rf "$SRV/shiny_DPDP.new" && mkdir "$SRV/shiny_DPDP.new"
git -C "$REPO/shiny_DPDP" archive HEAD | tar -x -C "$SRV/shiny_DPDP.new"
rm -rf "$SRV/shiny_DPDP" && mv "$SRV/shiny_DPDP.new" "$SRV/shiny_DPDP"
chown -R root:root "$SRV/web" "$SRV/dpdp-marketplace" "$SRV/shiny_DPDP"
chmod -R a+rX "$SRV/web" "$SRV/dpdp-marketplace" "$SRV/shiny_DPDP"

# 3. Shared run folder: group-writable, new files keep the group
chgrp "$GROUP" "$SRV/data" "$SRV/data/input"
chmod 2775 "$SRV/data" "$SRV/data/input"

# 4. Login: PAM service for the app (same passwords as RStudio)
cat > /etc/pam.d/dpdp-web <<'PAM'
@include common-auth
@include common-account
PAM

# 5. Self-signed certificate (10 years, for the server's IP)
mkdir -p /etc/ssl/dpdp
if [ ! -f /etc/ssl/dpdp/dpdp.crt ]; then
  openssl req -x509 -newkey rsa:2048 -nodes -days 3650 -keyout /etc/ssl/dpdp/dpdp.key -out /etc/ssl/dpdp/dpdp.crt \
    -subj "/CN=$IP/O=Duke Proteomics DPDP" -addext "subjectAltName=IP:$IP" 2>/dev/null
  chmod 600 /etc/ssl/dpdp/dpdp.key
fi

# 6. Apache: /dpdp/ on port 80 (redirects to https) and 443 -> the app
cat > /etc/apache2/conf-available/dpdp.conf <<CONF
# DPDP web app (installed by claude-dpdp-pipeline web/install.sh)
<Location /dpdp/>
  ProxyPreserveHost On
  ProxyPass http://127.0.0.1:$PORT/dpdp/ upgrade=websocket timeout=3600
  ProxyPassReverse http://127.0.0.1:$PORT/dpdp/
  LimitRequestBody 0
</Location>
RedirectMatch 302 ^/dpdp\$ /dpdp/
CONF
cat > /etc/apache2/sites-available/dpdp-ssl.conf <<CONF
# DPDP web app over https (installed by claude-dpdp-pipeline web/install.sh)
<VirtualHost *:443>
  ServerName $IP
  DocumentRoot /var/www/html
  SSLEngine on
  SSLCertificateFile /etc/ssl/dpdp/dpdp.crt
  SSLCertificateKeyFile /etc/ssl/dpdp/dpdp.key
  RequestHeader set X-Forwarded-Proto "https"
  ErrorLog \${APACHE_LOG_DIR}/dpdp_error.log
  CustomLog \${APACHE_LOG_DIR}/dpdp_access.log combined
</VirtualHost>
CONF
a2enmod -q ssl proxy proxy_http proxy_wstunnel headers alias >/dev/null
a2enconf -q dpdp >/dev/null
a2ensite -q dpdp-ssl >/dev/null
apache2ctl configtest
systemctl reload apache2

# 7. Service
cat > /etc/systemd/system/dpdp-web.service <<UNIT
[Unit]
Description=DPDP web app (pipeline review pages + Claude Code terminals)
After=network.target

[Service]
Environment=DPDP_DATA=$SRV/data
Environment=DPDP_APP=$SRV/shiny_DPDP
Environment=DPDP_PLUGIN=$SRV/dpdp-marketplace
Environment=DPDP_WEB_GROUP=$GROUP
Environment=DPDP_WEB_PORT=$PORT
ExecStart=$SRV/venv/bin/python $SRV/web/server.py
Restart=on-failure
KillMode=mixed

[Install]
WantedBy=multi-user.target
UNIT
systemctl daemon-reload
systemctl enable -q dpdp-web
systemctl restart dpdp-web
sleep 1
systemctl is-active dpdp-web >/dev/null && echo "== Running: https://$IP/dpdp/  (logs: journalctl -u dpdp-web)" \
  || { journalctl -u dpdp-web -n 30 --no-pager; exit 1; }
