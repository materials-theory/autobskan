import argparse
import os
import sys
import tempfile
import threading
import webbrowser

import autobskan


def _parse_args(argv=None):
    parser = argparse.ArgumentParser(description="Run the AutoBSKAN analysis workspace.")
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {autobskan.__version__}",
    )
    parser.add_argument("--host", default="127.0.0.1", help="Host address.")
    parser.add_argument("--port", default=8050, type=int, help="Port number.")
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Enable Dash debug mode.",
    )
    parser.add_argument(
        "--no-browser",
        action="store_true",
        help="Do not auto-open a browser tab.",
    )
    return parser.parse_args(argv)


def _open_browser(url):
    try:
        webbrowser.open_new(url)
    except Exception:
        # Browser opening can fail in restricted environments. Keep server alive.
        pass


def _prepare_matplotlib_cache():
    if os.environ.get("MPLCONFIGDIR"):
        return
    if sys.platform == "darwin":
        cache_dir = os.path.expanduser("~/Library/Caches/AutoBSKAN/matplotlib")
    else:
        cache_root = os.environ.get("XDG_CACHE_HOME", os.path.expanduser("~/.cache"))
        cache_dir = os.path.join(cache_root, "autobskan", "matplotlib")
    try:
        os.makedirs(cache_dir, exist_ok=True)
        with tempfile.NamedTemporaryFile(dir=cache_dir):
            pass
    except OSError:
        user_tag = str(os.getuid()) if hasattr(os, "getuid") else "user"
        cache_root = os.path.join(tempfile.gettempdir(), f"autobskan-{user_tag}")
        cache_dir = os.path.join(cache_root, "matplotlib")
        os.makedirs(cache_dir, exist_ok=True)
        os.environ.setdefault("XDG_CACHE_HOME", cache_root)
    os.environ["MPLCONFIGDIR"] = cache_dir


def main(argv=None):
    args = _parse_args(argv)
    _prepare_matplotlib_cache()

    from autobskan.gui.gui import build_app

    print(
        f"[autobskan-gui] AutoBSKAN {autobskan.__version__} | "
        f"http://{args.host}:{args.port}",
        flush=True,
    )
    app = build_app(
        debug_mode=args.debug,
        shutdown_on_last_client=True,
    )
    should_open_browser = (not args.no_browser) and (
        (not args.debug) or os.environ.get("WERKZEUG_RUN_MAIN") == "true"
    )
    if should_open_browser:
        url = f"http://{args.host}:{args.port}/"
        print(f"[autobskan-gui] opening browser: {url}", flush=True)
        threading.Timer(1.0, lambda: _open_browser(url)).start()
    try:
        app.run(host=args.host, port=args.port, debug=args.debug)
    except KeyboardInterrupt:
        print("[autobskan-gui] GUI closed; server stopped.", flush=True)


if __name__ == "__main__":
    main()
