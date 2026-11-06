# DPDP web app

The pipeline in a browser tab, like the Shiny app: `https://<server-ip>/dpdp/`.

- **Log in** with the same username and password as RStudio on this server. You must be in the `rstudio-users` group.
- **Left:** your runs and each step's review page, which opens by itself as Claude builds it (**Follow Claude**). Decide on the
  page with **Send to Claude**, or in the chat.
- **Right:** Claude Code, the official one, running on this server under *your* Linux account and *your own* Claude account. It
  keeps running if you close the tab; open the page again to continue.
- **New project:** upload the Spectronaut report and sample sheet; they go to `/srv/dpdp/data/input/<name>/`, and the request
  is typed into Claude for you (press Enter).

## First time, per user
1. The browser warns about the certificate (self-signed). Accept it once.
2. Claude Code installs itself for you (about a minute), then asks you to log in: click the link, sign in with your Duke Claude
   account, paste the code back.
3. Claude asks whether you trust the folder `~/dpdp`: choose **Yes**. That folder's settings let Claude use the shared run
   folder and run the pipeline scripts without asking each time.

## Admin
- Install or update (also after every code change): `sudo bash web/install.sh` from the repo. It rebuilds the plugin and copies
  the web app, the plugin and the app (at the pinned submodule commit) to `/srv/dpdp/`, sets up Apache (`/dpdp/` on ports 80 and
  443) and the `dpdp-web` service. Restarting the service ends open Claude sessions; users click **Restart Claude**, and
  Claude continues the last conversation.
- Shared run folder: `/srv/dpdp/data` (group `rstudio-users`, group-writable). Logs: `journalctl -u dpdp-web`.
- The service runs as root (like RStudio Server) to check passwords and start each user's terminal as that user. It listens on
  127.0.0.1:8770 only; Apache is the way in.
- Test without root: `DPDP_DATA=<data> python3 web/server.py --dev` prints a one-time login link for your own account
  (needs `aiohttp` and `python-pam`).
