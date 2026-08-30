(function () {
  "use strict";

  var renderControlQueue = {
    timer: null,
    resolve: null,
    latest: null,
    lastSentAt: 0,
    lastSignature: null,
  };
  var renderControlIntervalMs = 120;

  function noUpdate() {
    return window.dash_clientside.no_update;
  }

  function finiteNumber(value) {
    var number = Number(value);
    return Number.isFinite(number) ? number : null;
  }

  function renderControlPayload(isosurface, gaussianSigma, atomRadius, fitRadius) {
    var payload = {
      isosurface: finiteNumber(isosurface),
      gaussian_sigma: finiteNumber(gaussianSigma),
      atom_radius_scale: finiteNumber(atomRadius),
      fit_radius: finiteNumber(fitRadius),
    };
    if (Object.keys(payload).some(function (key) { return payload[key] === null; })) {
      return null;
    }
    return payload;
  }

  function emitRenderControls(payload) {
    renderControlQueue.lastSentAt = window.performance.now();
    renderControlQueue.lastSignature = JSON.stringify(payload);
    return Object.assign({}, payload, { nonce: Date.now() });
  }

  function roleOf(trace) {
    if (!trace || !trace.meta || typeof trace.meta !== "object") {
      return null;
    }
    return trace.meta.autobskan_role || null;
  }

  function brightnessLimits(dataMin, dataMax, brightness) {
    var lower = Number(dataMin);
    var upper = Number(dataMax);
    var value = Number(brightness);
    var span = Math.max(upper - lower, 1.0e-12);
    if (value < 0) {
      upper += -value * span;
    } else {
      lower -= value * span;
    }
    if (Math.abs(upper - lower) <= 1.0e-12) {
      upper = lower + 1.0e-9;
    }
    return [lower, upper];
  }

  function colorLimits(meta, brightness, mode, manualMin, manualMax) {
    var dataMin = finiteNumber(meta && meta.data_min);
    var dataMax = finiteNumber(meta && meta.data_max);
    if (dataMin === null || dataMax === null) {
      return null;
    }
    if (String(mode).toUpperCase() === "MANUAL") {
      var lower = finiteNumber(manualMin);
      var upper = finiteNumber(manualMax);
      if (lower !== null && upper !== null && lower < upper) {
        return [lower, upper];
      }
      return [dataMin, dataMax];
    }
    return brightnessLimits(dataMin, dataMax, brightness);
  }

  function optionEnabled(options, name) {
    return Array.isArray(options) && options.indexOf(name) !== -1;
  }

  function updateHeatmap(figure, role, colorscale, limits) {
    if (!figure || !Array.isArray(figure.data)) {
      return noUpdate();
    }
    var changed = false;
    var data = figure.data.map(function (trace) {
      if (roleOf(trace) !== role) {
        return trace;
      }
      changed = true;
      return Object.assign({}, trace, {
        colorscale: colorscale,
        zmin: limits[0],
        zmax: limits[1],
      });
    });
    if (!changed) {
      return noUpdate();
    }
    return Object.assign({}, figure, { data: data });
  }

  function updateLineRange(figure, limits) {
    if (!figure || !Array.isArray(figure.data) || figure.data.length === 0) {
      return noUpdate();
    }
    var hasLine = figure.data.some(function (trace) {
      return trace && String(trace.mode || "").indexOf("lines") !== -1;
    });
    if (!hasLine) {
      return noUpdate();
    }
    var layout = Object.assign({}, figure.layout || {});
    layout.yaxis = Object.assign({}, layout.yaxis || {}, {
      range: [limits[0], limits[1]],
      autorange: false,
    });
    return Object.assign({}, figure, { layout: layout });
  }

  window.dash_clientside = window.dash_clientside || {};
  window.dash_clientside.autobskan = Object.assign(
    {},
    window.dash_clientside.autobskan,
    {
      beginRender: function () {
        return "render-overlay is-active";
      },

      surfaceErrors: function (
        renderMessage,
        renderColor,
        volumeMessage,
        volumeColor,
        structureMessage,
        structureColor
      ) {
        var candidates = [
          [renderMessage, renderColor],
          [volumeMessage, volumeColor],
          [structureMessage, structureColor],
        ];
        var error = candidates.find(function (candidate) {
          return String(candidate[1] || "").toLowerCase() === "red";
        });
        if (!error || !error[0]) {
          return ["", "user-error-banner"];
        }
        return [String(error[0]), "user-error-banner is-visible"];
      },

      scheduleRenderControls: function (
        isosurface,
        gaussianSigma,
        atomRadius,
        fitRadius,
        viewOptions,
        simulation
      ) {
        var context = window.dash_clientside.callback_context;
        var trigger = context && context.triggered_id;
        if (
          trigger === "gaussian-sigma" &&
          !optionEnabled(viewOptions, "show_blurred")
        ) {
          return noUpdate();
        }
        if (
          trigger === "atom-radius-scale" &&
          !optionEnabled(viewOptions, "show_atoms")
        ) {
          return noUpdate();
        }
        if (trigger === "fit-radius-slider" && simulation !== "PHI_APP") {
          return noUpdate();
        }
        var payload = renderControlPayload(
          isosurface,
          gaussianSigma,
          atomRadius,
          fitRadius
        );
        if (!payload) {
          return noUpdate();
        }
        var signature = JSON.stringify(payload);
        renderControlQueue.latest = payload;
        if (
          signature === renderControlQueue.lastSignature &&
          renderControlQueue.timer === null
        ) {
          return noUpdate();
        }

        if (renderControlQueue.resolve) {
          renderControlQueue.resolve(noUpdate());
          renderControlQueue.resolve = null;
        }
        if (renderControlQueue.timer !== null) {
          window.clearTimeout(renderControlQueue.timer);
          renderControlQueue.timer = null;
        }

        var elapsed = window.performance.now() - renderControlQueue.lastSentAt;
        var delay = Math.max(0, renderControlIntervalMs - elapsed);
        if (delay === 0) {
          return emitRenderControls(payload);
        }

        return new Promise(function (resolve) {
          renderControlQueue.resolve = resolve;
          renderControlQueue.timer = window.setTimeout(function () {
            var latest = renderControlQueue.latest;
            renderControlQueue.timer = null;
            renderControlQueue.resolve = null;
            resolve(emitRenderControls(latest));
          }, delay);
        });
      },

      scheduleContextRender: function (
        gamma,
        layers,
        repeatX,
        repeatY,
        atomRadiusType,
        gridLineWidth,
        gridLineColor,
        viewOptions,
        inputSource
      ) {
        var context = window.dash_clientside.callback_context;
        var trigger = context && context.triggered_id;
        var needsAtoms = trigger === "layers-input" || trigger === "atom-radius-type";
        var needsRepeat = trigger === "repeat-x" || trigger === "repeat-y";
        var needsGrid = trigger === "grid-line-width" || trigger === "grid-line-color";
        var needsGamma = trigger === "gamma-input";

        if (needsAtoms && !optionEnabled(viewOptions, "show_atoms")) {
          return noUpdate();
        }
        if (needsRepeat && !optionEnabled(viewOptions, "show_repeated")) {
          return noUpdate();
        }
        if (needsGrid && !optionEnabled(viewOptions, "show_grids")) {
          return noUpdate();
        }
        if (
          needsGamma &&
          (inputSource !== "BSKAN" || !optionEnabled(viewOptions, "show_repeated"))
        ) {
          return noUpdate();
        }
        return { trigger: trigger, nonce: Date.now() };
      },

      syncOverlayControls: function (viewOptions) {
        var repeatDisabled = !optionEnabled(viewOptions, "show_repeated");
        var atomsDisabled = !optionEnabled(viewOptions, "show_atoms");
        var blurDisabled = !optionEnabled(viewOptions, "show_blurred");
        var gridDisabled = !optionEnabled(viewOptions, "show_grids");
        return [
          repeatDisabled,
          repeatDisabled,
          atomsDisabled,
          atomsDisabled,
          atomsDisabled,
          atomsDisabled,
          blurDisabled,
          blurDisabled,
          gridDisabled,
          gridDisabled,
        ];
      },

      syncColorRangeControls: function (mode, meta, currentMin, currentMax) {
        var manual = String(mode).toUpperCase() === "MANUAL";
        var dataMin = finiteNumber(meta && meta.data_min);
        var dataMax = finiteNumber(meta && meta.data_max);
        var lower = finiteNumber(currentMin);
        var upper = finiteNumber(currentMax);

        if (!manual && dataMin !== null && dataMax !== null) {
          lower = dataMin;
          upper = dataMax;
        } else if (
          manual &&
          (lower === null || upper === null || lower >= upper) &&
          dataMin !== null &&
          dataMax !== null
        ) {
          lower = dataMin;
          upper = dataMax;
        }
        if (lower === null || upper === null || lower >= upper) {
          lower = 0.0;
          upper = 1.0;
        }
        return [lower, upper, !manual, !manual, manual, manual];
      },

      colorRangeStatus: function (mode, vmin, vmax, meta) {
        if (String(mode).toUpperCase() !== "MANUAL") {
          var dataMin = finiteNumber(meta && meta.data_min);
          var dataMax = finiteNumber(meta && meta.data_max);
          if (dataMin === null || dataMax === null) {
            return [
              "Auto range will follow the rendered data.",
              "field-hint range-status",
            ];
          }
          return [
            "Auto data range: " + dataMin.toPrecision(7) + " to " + dataMax.toPrecision(7),
            "field-hint range-status is-valid",
          ];
        }
        var lower = finiteNumber(vmin);
        var upper = finiteNumber(vmax);
        if (lower === null || upper === null || lower >= upper) {
          return [
            "Manual range requires finite values with vmin < vmax.",
            "field-hint range-status is-error",
          ];
        }
        return [
          "Manual range is applied to the map, colorbar, and line profile.",
          "field-hint range-status is-valid",
        ];
      },

      syncSliderNumber: function (sliderValue, numberValue, minimum, maximum) {
        var context = window.dash_clientside.callback_context;
        var trigger = context && context.triggered_id;
        var lower = finiteNumber(minimum);
        var upper = finiteNumber(maximum);

        function clamp(value) {
          var number = finiteNumber(value);
          if (number === null) {
            return null;
          }
          if (lower !== null) {
            number = Math.max(lower, number);
          }
          if (upper !== null) {
            number = Math.min(upper, number);
          }
          return number;
        }

        if (typeof trigger === "string" && trigger.endsWith("-input")) {
          var typed = clamp(numberValue);
          return typed === null
            ? [noUpdate(), noUpdate()]
            : [typed, noUpdate()];
        }
        var dragged = clamp(sliderValue);
        return dragged === null
          ? [noUpdate(), noUpdate()]
          : [noUpdate(), dragged];
      },

      syncSliderBounds: function (minimum, maximum, step) {
        return [minimum, maximum, step];
      },

      applySliderRange: function (
        minimum,
        maximum,
        defaultValue,
        currentValue
      ) {
        var lower = finiteNumber(minimum);
        var upper = finiteNumber(maximum);
        var current = finiteNumber(currentValue);
        if (lower === null || upper === null || lower > upper) {
          return noUpdate();
        }
        if (current !== null && current >= lower && current <= upper) {
          return noUpdate();
        }
        var fallback = finiteNumber(defaultValue);
        if (fallback === null) {
          fallback = lower;
        }
        return Math.min(upper, Math.max(lower, fallback));
      },

      updateAppearance: function (
        appearance,
        brightness,
        rangeMode,
        manualMin,
        manualMax,
        mainFigure,
        colorbarFigure,
        lineFigure,
        meta
      ) {
        if (
          !appearance ||
          !Array.isArray(appearance.colorscale) ||
          !meta ||
          !Number.isFinite(Number(meta.data_min)) ||
          !Number.isFinite(Number(meta.data_max))
        ) {
          return [noUpdate(), noUpdate(), noUpdate()];
        }
        var limits = colorLimits(
          meta,
          brightness,
          rangeMode,
          manualMin,
          manualMax
        );
        if (!limits) {
          return [noUpdate(), noUpdate(), noUpdate()];
        }
        var colorbar = updateHeatmap(
          colorbarFigure,
          "colorbar-field",
          appearance.colorscale,
          limits
        );
        return [
          updateHeatmap(
            mainFigure,
            "scalar-field",
            appearance.colorscale,
            limits
          ),
          colorbar,
          updateLineRange(lineFigure, limits),
        ];
      },

      syncExportDimensions: function (width, height, meta) {
        var currentWidth = finiteNumber(width);
        var currentHeight = finiteNumber(height);
        var xMin = finiteNumber(meta && meta.x_min);
        var xMax = finiteNumber(meta && meta.x_max);
        var yMin = finiteNumber(meta && meta.y_min);
        var yMax = finiteNumber(meta && meta.y_max);
        var ratio = null;
        if (
          xMin !== null && xMax !== null && yMin !== null && yMax !== null &&
          Math.abs(yMax - yMin) > 1.0e-12
        ) {
          ratio = Math.abs(xMax - xMin) / Math.abs(yMax - yMin);
        }
        if (!Number.isFinite(ratio) || ratio <= 0) {
          ratio = Math.max(currentWidth || 4096, 1) / Math.max(currentHeight || 2048, 1);
        }

        var trigger = window.dash_clientside.callback_context &&
          window.dash_clientside.callback_context.triggered_id;
        var nextWidth;
        var nextHeight;
        if (trigger === "export-height") {
          nextHeight = Math.round(currentHeight || 2048);
          nextWidth = Math.round(nextHeight * ratio);
        } else {
          nextWidth = Math.round(currentWidth || 4096);
          nextHeight = Math.round(nextWidth / ratio);
        }

        var minimum = 256;
        var maximum = 12000;
        if (nextWidth < minimum) {
          nextWidth = minimum;
          nextHeight = Math.round(nextWidth / ratio);
        }
        if (nextHeight < minimum) {
          nextHeight = minimum;
          nextWidth = Math.round(nextHeight * ratio);
        }
        if (nextWidth > maximum) {
          nextWidth = maximum;
          nextHeight = Math.round(nextWidth / ratio);
        }
        if (nextHeight > maximum) {
          nextHeight = maximum;
          nextWidth = Math.round(nextHeight * ratio);
        }
        var pixelLimit = 60000000;
        if (nextWidth * nextHeight > pixelLimit) {
          var reduction = Math.sqrt(pixelLimit / (nextWidth * nextHeight));
          nextWidth = Math.max(minimum, Math.floor(nextWidth * reduction));
          nextHeight = Math.max(minimum, Math.floor(nextHeight * reduction));
        }

        return [
          currentWidth === nextWidth ? noUpdate() : nextWidth,
          currentHeight === nextHeight ? noUpdate() : nextHeight,
          nextWidth.toLocaleString() + " x " + nextHeight.toLocaleString() + " px",
        ];
      },

      toggleControls: function (
        toggleClicks,
        closeClicks,
        backdropClicks,
        renderClicks,
        keyboardRequest,
        currentlyOpen
      ) {
        var context = window.dash_clientside.callback_context;
        var trigger = context && context.triggered_id;
        var mobile = window.matchMedia("(max-width: 1260px)").matches;
        var rail = document.getElementById("control-rail");
        var currentlyVisible = mobile
          ? Boolean(rail && rail.classList.contains("controls-open"))
          : Boolean(rail && !rail.classList.contains("controls-collapsed"));
        var opened = currentlyVisible;
        if (trigger === "controls-toggle-button") {
          opened = !currentlyVisible;
        } else if (trigger === "controls-keyboard-store") {
          var action = keyboardRequest && keyboardRequest.action;
          opened = action === "open" ? true : !currentlyVisible;
        } else if (
          trigger === "controls-close-button" ||
          trigger === "controls-backdrop"
        ) {
          opened = false;
        } else if (trigger === "render-button") {
          opened = mobile ? false : currentlyVisible;
        } else {
          return [noUpdate(), noUpdate(), noUpdate(), noUpdate()];
        }
        var railClass = mobile
          ? opened
            ? "control-rail controls-open"
            : "control-rail"
          : opened
            ? "control-rail"
            : "control-rail controls-collapsed";
        return [
          opened,
          railClass,
          mobile && opened
            ? "controls-backdrop controls-open"
            : "controls-backdrop",
          !mobile && !opened
            ? "workspace-grid controls-collapsed"
            : "workspace-grid",
        ];
      },
    }
  );

  if (!window.__autobskanControlKeys) {
    window.__autobskanControlKeys = true;
    window.addEventListener("keydown", function (event) {
      var rail = document.getElementById("control-rail");
      if (!rail) {
        return;
      }
      if (document.querySelector('[role="dialog"]')) {
        return;
      }
      var target = event.target;
      var tag = target && target.tagName ? target.tagName.toLowerCase() : "";
      var isFormControl =
        ["input", "textarea", "select", "button"].indexOf(tag) !== -1 ||
        Boolean(target && target.isContentEditable);
      var mobile = window.matchMedia("(max-width: 1260px)").matches;
      var collapsed = mobile
        ? !rail.classList.contains("controls-open")
        : rail.classList.contains("controls-collapsed");

      if (event.key === "Escape") {
        event.preventDefault();
        if (isFormControl && typeof target.blur === "function") {
          target.blur();
        }
        var closeButton = document.getElementById("controls-close-button");
        if (closeButton) {
          closeButton.click();
        } else if (window.dash_clientside.set_props) {
          window.dash_clientside.set_props("controls-keyboard-store", {
            data: { action: "toggle", nonce: Date.now() },
          });
        }
        return;
      }
      if (event.key === "Tab" && collapsed && !isFormControl) {
        event.preventDefault();
        var toggleButton = document.getElementById("controls-toggle-button");
        if (toggleButton) {
          toggleButton.click();
        } else if (window.dash_clientside.set_props) {
          window.dash_clientside.set_props("controls-keyboard-store", {
            data: { action: "open", nonce: Date.now() },
          });
        }
      }
    });
  }
})();
