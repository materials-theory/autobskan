(function () {
  "use strict";

  if (window.__autobskanChunkedUpload) {
    return;
  }
  window.__autobskanChunkedUpload = true;

  function elements() {
    return {
      button: document.getElementById("upload-volume-button"),
      picker: document.getElementById("upload-volume-picker"),
      progress: document.getElementById("upload-volume-progress"),
      status: document.getElementById("upload-volume-status"),
      wrap: document.getElementById("upload-volume-progress-wrap"),
    };
  }

  function ensurePicker() {
    var picker = document.getElementById("upload-volume-picker");
    if (picker) {
      return picker;
    }
    picker = document.createElement("input");
    picker.id = "upload-volume-picker";
    picker.type = "file";
    picker.style.display = "none";
    document.body.appendChild(picker);
    return picker;
  }

  function updateProgress(percent, message, state) {
    var controls = elements();
    if (controls.progress) {
      controls.progress.value = Math.max(0, Math.min(100, percent));
    }
    if (controls.status) {
      controls.status.textContent = message || "";
    }
    if (controls.wrap) {
      controls.wrap.dataset.state = state || "idle";
    }
  }

  function setDashProps(componentId, props) {
    if (
      !window.dash_clientside ||
      typeof window.dash_clientside.set_props !== "function"
    ) {
      throw new Error("The GUI is not ready to receive the uploaded file.");
    }
    window.dash_clientside.set_props(componentId, props);
  }

  async function responseJson(response) {
    var payload = await response.json().catch(function () { return {}; });
    if (!response.ok) {
      throw new Error(payload.error || "Volumetric upload failed.");
    }
    return payload;
  }

  async function uploadVolume(file) {
    var controls = elements();
    var uploadId = null;
    if (controls.button) {
      controls.button.disabled = true;
    }
    updateProgress(0, "Preparing " + file.name + "...", "active");

    try {
      var startResponse = await fetch("/_autobskan-upload/start", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ filename: file.name, size: file.size }),
      });
      var start = await responseJson(startResponse);
      uploadId = start.upload_id;
      var chunkSize = Number(start.chunk_size) || 8 * 1024 * 1024;

      for (var offset = 0; offset < file.size; offset += chunkSize) {
        var chunk = file.slice(offset, Math.min(offset + chunkSize, file.size));
        var chunkResponse = await fetch(
          "/_autobskan-upload/chunk/" + encodeURIComponent(uploadId),
          {
            method: "POST",
            headers: { "X-Upload-Offset": String(offset) },
            body: chunk,
          }
        );
        var chunkResult = await responseJson(chunkResponse);
        var percent = (100 * Number(chunkResult.received)) / file.size;
        updateProgress(
          percent,
          "Uploading " + file.name + " (" + Math.floor(percent) + "%)",
          "active"
        );
      }

      var finishResponse = await fetch(
        "/_autobskan-upload/finish/" + encodeURIComponent(uploadId),
        { method: "POST" }
      );
      var finished = await responseJson(finishResponse);
      uploadId = null;

      setDashProps("volume-upload-filename-store", { data: finished.filename });
      setDashProps("volume-path", { value: "" });
      setDashProps("volume-file-store", { data: finished.path });
      setDashProps("volume-upload-ready-store", {
        data: new Date().toISOString(),
      });
      updateProgress(100, "Ready: " + finished.filename, "complete");
    } catch (error) {
      if (uploadId) {
        fetch("/_autobskan-upload/" + encodeURIComponent(uploadId), {
          method: "DELETE",
        }).catch(function () {});
      }
      updateProgress(
        0,
        error && error.message ? error.message : "Volumetric upload failed.",
        "error"
      );
    } finally {
      if (controls.button) {
        controls.button.disabled = false;
      }
      if (controls.picker) {
        controls.picker.value = "";
      }
    }
  }

  document.addEventListener("click", function (event) {
    var trigger = event.target.closest("#upload-volume-button");
    if (!trigger) {
      return;
    }
    var picker = ensurePicker();
    if (picker && !trigger.disabled) {
      picker.click();
    }
  });

  document.addEventListener("change", function (event) {
    if (event.target.id !== "upload-volume-picker") {
      return;
    }
    var file = event.target.files && event.target.files[0];
    if (file) {
      uploadVolume(file);
    }
  });
})();
