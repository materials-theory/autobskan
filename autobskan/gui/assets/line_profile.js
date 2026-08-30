(function () {
  "use strict";

  function isLineSelection(trace) {
    return Boolean(
      trace &&
        trace.meta &&
        typeof trace.meta === "object" &&
        trace.meta.autobskan_role === "line-selection"
    );
  }

  function normalisePoints(points) {
    if (!Array.isArray(points)) {
      return [];
    }
    return points.slice(0, 2).reduce(function (result, point) {
      if (!Array.isArray(point) || point.length < 2) {
        return result;
      }
      var x = Number(point[0]);
      var y = Number(point[1]);
      if (Number.isFinite(x) && Number.isFinite(y)) {
        result.push([x, y]);
      }
      return result;
    }, []);
  }

  function nextPoints(currentPoints, clickData) {
    var clicked = clickData && clickData.points && clickData.points[0];
    if (!clicked) {
      return null;
    }

    var x = Number(clicked.x);
    var y = Number(clicked.y);
    if (!Number.isFinite(x) || !Number.isFinite(y)) {
      return null;
    }

    var point = [x, y];
    var current = normalisePoints(currentPoints);
    if (current.length !== 1) {
      return [point];
    }

    var dx = point[0] - current[0][0];
    var dy = point[1] - current[0][1];
    if (Math.hypot(dx, dy) <= 1.0e-12) {
      return [point];
    }
    return [current[0], point];
  }

  function selectionTraces(points) {
    var traces = [];
    if (points.length === 2) {
      traces.push({
        type: "scatter",
        x: [points[0][0], points[1][0]],
        y: [points[0][1], points[1][1]],
        mode: "lines",
        line: { color: "#0f766e", width: 3 },
        hoverinfo: "skip",
        cliponaxis: true,
        showlegend: false,
        meta: { autobskan_role: "line-selection" },
      });
    }
    if (points.length > 0) {
      traces.push({
        type: "scatter",
        x: points.map(function (point) {
          return point[0];
        }),
        y: points.map(function (point) {
          return point[1];
        }),
        mode: "markers+text",
        text: ["P1", "P2"].slice(0, points.length),
        textposition: "top center",
        marker: {
          size: 13,
          color: ["#0f766e", "#d97706"].slice(0, points.length),
          line: { color: "#ffffff", width: 2 },
        },
        textfont: { color: "#134e4a", size: 13 },
        hoverinfo: "skip",
        cliponaxis: true,
        showlegend: false,
        meta: { autobskan_role: "line-selection" },
      });
    }
    return traces;
  }

  function updateFigure(figure, points) {
    if (!figure || !Array.isArray(figure.data)) {
      return window.dash_clientside.no_update;
    }
    var data = figure.data.filter(function (trace) {
      return !isLineSelection(trace);
    });
    data = data.concat(selectionTraces(points));
    return Object.assign({}, figure, { data: data });
  }

  window.dash_clientside = window.dash_clientside || {};
  window.dash_clientside.autobskan = Object.assign(
    {},
    window.dash_clientside.autobskan,
    {
      updateLineSelection: function (
        clickData,
        clearClicks,
        lineEnabled,
        currentPoints,
        figure
      ) {
        var context = window.dash_clientside.callback_context;
        var trigger = context && context.triggered_id;
        var points;

        if (trigger === "clear-line-button" || !lineEnabled) {
          points = [];
        } else if (trigger === "main-image") {
          points = nextPoints(currentPoints, clickData);
          if (points === null) {
            return [
              window.dash_clientside.no_update,
              window.dash_clientside.no_update,
            ];
          }
        } else {
          return [
            window.dash_clientside.no_update,
            window.dash_clientside.no_update,
          ];
        }

        return [points, updateFigure(figure, points)];
      },

      applyLineDisplay: function (points, figure) {
        return updateFigure(figure, normalisePoints(points));
      },
    }
  );
})();
