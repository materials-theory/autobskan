(function () {
  "use strict";

  if (window.__autobskanClientLifecycle) {
    return;
  }
  window.__autobskanClientLifecycle = true;

  var clientId = globalThis.crypto && globalThis.crypto.randomUUID
    ? globalThis.crypto.randomUUID()
    : Date.now().toString(36) + "-" + Math.random().toString(36).slice(2);
  var closing = false;
  var lifecycleStream = null;

  function endpoint(action) {
    return "/_autobskan-client/" + action + "/" + encodeURIComponent(clientId);
  }

  function register() {
    if (typeof window.EventSource === "function") {
      lifecycleStream = new window.EventSource(endpoint("events"));
      return;
    }
    window.fetch(endpoint("open"), {
      method: "POST",
      keepalive: true,
      cache: "no-store",
    }).catch(function () {
      // A headless or embedded viewer may stop before registration completes.
    });
  }

  function closeClient() {
    if (closing) {
      return;
    }
    closing = true;
    if (lifecycleStream) {
      lifecycleStream.close();
      lifecycleStream = null;
    }
    var closeUrl = endpoint("close");
    var queued = false;
    try {
      queued = Boolean(navigator.sendBeacon(closeUrl, ""));
    } catch (_error) {
      queued = false;
    }
    if (!queued) {
      window.fetch(closeUrl, {
        method: "POST",
        keepalive: true,
        cache: "no-store",
      }).catch(function () {
        // Page teardown can finish before the fallback request resolves.
      });
    }
  }

  window.addEventListener("pagehide", closeClient, { once: true });
  window.addEventListener("beforeunload", closeClient, { once: true });
  register();
})();
