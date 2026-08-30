from __future__ import annotations

from collections.abc import Callable

import dash_mantine_components as dmc
from dash import dcc, html

GRAPH_CONFIG = {
    "displaylogo": False,
    "responsive": True,
    "displayModeBar": "hover",
    "modeBarButtonsToRemove": ["lasso2d", "select2d"],
}


def _field_label(text, component_id=None):
    props = {"className": "field-label"}
    if component_id is not None:
        props["id"] = component_id
    return dmc.Text(text, **props)


def _precision_slider(
    *,
    slider_id,
    input_id,
    value,
    minimum,
    maximum,
    step,
    decimal_scale,
    color="teal",
):
    return html.Div(
        [
            dmc.Slider(
                id=slider_id,
                min=minimum,
                max=maximum,
                step=step,
                value=value,
                updatemode="drag",
                labelAlwaysOn=False,
                color=color,
                marks=None,
                className="precision-slider",
            ),
            dmc.NumberInput(
                id=input_id,
                value=value,
                min=minimum,
                max=maximum,
                step=step,
                decimalScale=decimal_scale,
                debounce=250,
                clampBehavior="strict",
                size="xs",
                className="precision-number-input",
            ),
        ],
        className="precision-slider-row",
    )


def _upload_button(label, upload_id):
    return dcc.Upload(
        id=upload_id,
        children=dmc.Button(
            label,
            variant="light",
            color="teal",
            size="sm",
            radius="sm",
            fullWidth=True,
        ),
        multiple=False,
        className="upload-control",
    )


def _chunked_volume_upload():
    return html.Div(
        [
            dmc.Button(
                "Upload volumetric file",
                id="upload-volume-button",
                variant="light",
                color="teal",
                size="sm",
                radius="sm",
                fullWidth=True,
            ),
            html.Div(
                [
                    html.Progress(
                        id="upload-volume-progress",
                        value=0,
                        max=100,
                    ),
                    html.Span(
                        "",
                        id="upload-volume-status",
                        role="status",
                        **{"aria-live": "polite"},
                    ),
                ],
                id="upload-volume-progress-wrap",
                className="chunk-upload-progress",
            ),
        ],
        className="upload-control",
    )


def _debug_console(debug_mode, status_text, status_wrap):
    alert_style = {"whiteSpace": "pre-line", "fontSize": "0.78rem"}

    def _alert(wrapper_id, alert_id, title, message):
        return html.Div(
            id=wrapper_id,
            style=status_wrap(True),
            children=dmc.Alert(
                id=alert_id,
                title=title,
                children=status_text("WAIT", message),
                color="gray",
                variant="light",
                radius="sm",
                style=alert_style,
            ),
        )

    return html.Div(
        [
            dmc.Text("Debug status", className="field-label"),
            _alert(
                "render-status-wrap",
                "render-status",
                "Render",
                "No rendered image yet.",
            ),
            _alert(
                "volume-status-wrap",
                "volume-status",
                "Volumetric parse",
                "No volumetric file parsed yet.",
            ),
            _alert(
                "structure-status-wrap",
                "structure-status",
                "Structure parse",
                "No structure file parsed yet.",
            ),
        ],
        className="debug-console",
        style={"display": "grid" if debug_mode else "none"},
    )


def _data_controls(defaults):
    return dmc.Stack(
        [
            _field_label("Data source"),
            dmc.SegmentedControl(
                id="input-source",
                value=defaults["input_source"],
                data=[
                    {"label": "bSKAN", "value": "BSKAN"},
                    {"label": "VASP", "value": "VASP"},
                ],
                fullWidth=True,
                color="teal",
                radius="sm",
                className="mode-control",
            ),
            dmc.Text(
                "CURRENT for bSKAN; PARCHG, LOCPOT, or CHGCAR for VASP",
                id="volume-source-hint",
                className="field-hint",
            ),
            _field_label("CURRENT file", "volume-file-label"),
            dmc.TextInput(
                id="volume-path",
                placeholder="Path to volumetric data",
                debounce=450,
                size="sm",
            ),
            dmc.Text(
                "No volumetric file selected.",
                id="volume-path-feedback",
                className="path-feedback",
            ),
            _chunked_volume_upload(),
            html.Div(
                [
                    dmc.Text("Embedded structure", className="embedded-structure-title"),
                    dmc.Text(
                        "VASP structure metadata will appear here after loading.",
                        id="embedded-structure-message",
                        className="embedded-structure-message",
                    ),
                ],
                id="embedded-structure-banner",
                className="embedded-structure-banner",
                style={"display": "none"},
            ),
            html.Div(className="control-divider"),
            _field_label("Structure file"),
            dmc.TextInput(
                id="structure-path",
                placeholder="Path to POSCAR or CONTCAR",
                debounce=450,
                size="sm",
            ),
            dmc.Text(
                "No structure file selected.",
                id="structure-path-feedback",
                className="path-feedback",
            ),
            _upload_button("Upload structure file", "upload-structure"),
            dmc.NumberInput(
                id="gamma-input",
                label="Lattice gamma (deg)",
                value=90.0,
                step=0.001,
                decimalScale=6,
                debounce=350,
                size="sm",
            ),
        ],
        gap="xs",
    )


def _simulation_controls(defaults):
    return dmc.Stack(
        [
            _field_label("Simulation"),
            dmc.SegmentedControl(
                id="simulation-type",
                value=defaults["simulation"],
                data=[
                    {"label": "STM", "value": "STM"},
                    {
                        "label": dcc.Markdown(
                            r"$\Phi_{\rm app}$",
                            mathjax=True,
                            className="simulation-math-label",
                        ),
                        "value": "PHI_APP",
                    },
                    {
                        "label": dcc.Markdown(
                            r"$\Phi_{\rm loc}$",
                            mathjax=True,
                            className="simulation-math-label",
                        ),
                        "value": "LWF",
                    },
                ],
                fullWidth=True,
                color="teal",
                radius="sm",
                className="mode-control simulation-mode-control",
            ),
            dcc.Markdown(
                "CURRENT or PARCHG surface analysis.",
                id="simulation-definition",
                mathjax=True,
                className="field-hint simulation-definition",
            ),
            html.Div(
                [
                    _field_label("STM acquisition mode"),
                    dmc.SegmentedControl(
                        id="stm-mode",
                        value=defaults["mode"],
                        data=[
                            {"label": "Constant current", "value": "CONSTANT CURRENT"},
                            {"label": "Constant height", "value": "CONSTANT HEIGHT"},
                        ],
                        fullWidth=True,
                        color="teal",
                        radius="sm",
                        className="mode-control mode-control-wide",
                    ),
                ],
                id="mode-container",
                className="control-stack",
            ),
            html.Div(
                [
                    dmc.NumberInput(
                        id="fermi-level-input",
                        label="Fermi level E_F (eV)",
                        value=defaults.get("fermi_level"),
                        placeholder="Auto from sibling OUTCAR",
                        step=0.01,
                        decimalScale=8,
                        debounce=350,
                        size="sm",
                    ),
                    dmc.Text(
                        "Enter E_F or keep this empty to read a sibling OUTCAR.",
                        id="fermi-level-help",
                        className="field-hint",
                    ),
                ],
                id="fermi-level-container",
                className="control-stack",
                style={"display": "none"},
            ),
            _field_label("log10(Isosurface)", "isosurface-label"),
            _precision_slider(
                slider_id="isosurface-slider",
                input_id="isosurface-input",
                value=-2.5,
                minimum=-5.0,
                maximum=0.0,
                step=0.1,
                decimal_scale=6,
            ),
            dmc.Text(
                "Range set after file load.",
                id="isosurface-help",
                className="field-hint",
            ),
            html.Div(
                [
                    _field_label("Fit half-width (A)"),
                    _precision_slider(
                        slider_id="fit-radius-slider",
                        input_id="fit-radius-input",
                        value=defaults["fit_radius"],
                        minimum=0.05,
                        maximum=3.0,
                        step=0.05,
                        decimal_scale=4,
                        color="orange",
                    ),
                    dmc.Text(
                        "Regression uses z = fit height +/- this value.",
                        className="field-hint",
                    ),
                ],
                id="fit-radius-container",
                className="control-stack",
                style={"display": "none"},
            ),
        ],
        gap="sm",
    )


def _appearance_controls(defaults, colormaps):
    return dmc.Stack(
        [
            dmc.Select(
                id="colormap-dropdown",
                label="Colormap",
                data=[{"label": cmap, "value": cmap} for cmap in colormaps],
                value=defaults["colormap_stm"],
                allowDeselect=False,
                searchable=True,
                size="sm",
            ),
            html.Div(
                [
                    html.Div(
                        [
                            _field_label("Color range"),
                            dmc.Button(
                                "Auto range",
                                id="color-range-auto-button",
                                variant="subtle",
                                color="teal",
                                size="compact-xs",
                                radius="sm",
                            ),
                        ],
                        className="field-heading-row",
                    ),
                    dmc.SegmentedControl(
                        id="color-range-mode",
                        value=defaults["color_range_mode"],
                        data=[
                            {"label": "Auto", "value": "AUTO"},
                            {"label": "Manual", "value": "MANUAL"},
                        ],
                        fullWidth=True,
                        color="teal",
                        radius="sm",
                        size="xs",
                        className="mode-control",
                    ),
                    html.Div(
                        [
                            dmc.NumberInput(
                                id="color-vmin-input",
                                label="vmin",
                                value=0.0,
                                step=0.01,
                                decimalScale=8,
                                debounce=200,
                                disabled=True,
                                size="sm",
                            ),
                            dmc.NumberInput(
                                id="color-vmax-input",
                                label="vmax",
                                value=1.0,
                                step=0.01,
                                decimalScale=8,
                                debounce=200,
                                disabled=True,
                                size="sm",
                            ),
                        ],
                        className="two-column-fields",
                    ),
                    dmc.Text(
                        "Auto range follows the rendered data.",
                        id="color-range-help",
                        className="field-hint range-status",
                    ),
                ],
                className="control-stack color-range-controls",
            ),
            _field_label("Interactive display"),
            dmc.SegmentedControl(
                id="display-quality",
                value=defaults["display_quality"],
                data=[
                    {"label": "Fast", "value": "FAST"},
                    {"label": "Balanced", "value": "BALANCED"},
                    {"label": "Full", "value": "FULL"},
                ],
                fullWidth=True,
                color="teal",
                radius="sm",
                size="xs",
            ),
            _field_label("Brightness"),
            _precision_slider(
                slider_id="brightness-slider",
                input_id="brightness-input",
                value=defaults["brightness"],
                minimum=-1.0,
                maximum=1.0,
                step=0.01,
                decimal_scale=4,
            ),
            _field_label("Contrast"),
            _precision_slider(
                slider_id="contrast-slider",
                input_id="contrast-input",
                value=defaults["contrast"],
                minimum=-0.95,
                maximum=0.95,
                step=0.01,
                decimal_scale=4,
            ),
            dmc.NumberInput(
                id="xy-upsample",
                label="XY interpolation factor",
                min=1,
                max=8,
                step=1,
                allowDecimal=False,
                value=defaults["xy_upsample"],
                debounce=300,
                size="sm",
            ),
        ],
        gap="sm",
    )


def _overlay_controls(defaults):
    return dmc.Stack(
        [
            dmc.CheckboxGroup(
                id="view-options",
                value=defaults["view_options"],
                children=dmc.Stack(
                    [
                        html.Div(
                            [
                                dmc.Checkbox(
                                    value="show_blurred",
                                    label="Gaussian blur",
                                    color="orange",
                                ),
                                _field_label("Gaussian sigma"),
                                _precision_slider(
                                    slider_id="gaussian-sigma",
                                    input_id="gaussian-sigma-input",
                                    value=defaults["gaussian_sigma"],
                                    minimum=0.0,
                                    maximum=30.0,
                                    step=0.1,
                                    decimal_scale=3,
                                    color="orange",
                                ),
                            ],
                            className="blur-control-block",
                        ),
                        html.Div(
                            [
                                dmc.Checkbox(
                                    value="show_atoms", label="Atoms", color="teal"
                                ),
                                dmc.Checkbox(
                                    value="show_repeated",
                                    label="Repeated cells",
                                    color="teal",
                                ),
                                dmc.Checkbox(
                                    value="show_grids",
                                    label="Unit-cell grid",
                                    color="orange",
                                ),
                            ],
                            className="checkbox-grid",
                        ),
                    ],
                    gap="sm",
                ),
            ),
            html.Div(
                [
                    dmc.NumberInput(
                        id="repeat-x",
                        label="Repeat x",
                        min=1,
                        step=1,
                        allowDecimal=False,
                        value=defaults["repeat_x"],
                        debounce=300,
                        size="sm",
                    ),
                    dmc.NumberInput(
                        id="repeat-y",
                        label="Repeat y",
                        min=1,
                        step=1,
                        allowDecimal=False,
                        value=defaults["repeat_y"],
                        debounce=300,
                        size="sm",
                    ),
                ],
                className="two-column-fields",
            ),
            html.Div(
                [
                    dmc.NumberInput(
                        id="layers-input",
                        label="Outmost layers",
                        min=1,
                        max=50,
                        step=1,
                        value=defaults["layers"],
                        allowDecimal=False,
                        debounce=300,
                        size="sm",
                    ),
                    dmc.NumberInput(
                        id="scalebar-length",
                        label="Scalebar (A)",
                        min=0.0,
                        step=0.1,
                        value=defaults["scalebar"],
                        debounce=300,
                        size="sm",
                    ),
                ],
                className="two-column-fields",
            ),
            _field_label("Atom radius basis"),
            dmc.SegmentedControl(
                id="atom-radius-type",
                value=defaults["atom_radius_type"],
                data=[
                    {"label": "Atomic", "value": "ATOMIC"},
                    {"label": "vdW", "value": "VDW"},
                    {"label": "Ionic", "value": "IONIC"},
                ],
                fullWidth=True,
                color="teal",
                radius="sm",
            ),
            _field_label("Atom radius scale"),
            _precision_slider(
                slider_id="atom-radius-scale",
                input_id="atom-radius-scale-input",
                value=defaults["atom_radius_scale"],
                minimum=1.0,
                maximum=30.0,
                step=0.1,
                decimal_scale=3,
            ),
            html.Div(
                [
                    dmc.NumberInput(
                        id="grid-line-width",
                        label="Grid width",
                        min=0.1,
                        max=10.0,
                        step=0.1,
                        value=defaults["grid_line_width"],
                        debounce=300,
                        size="sm",
                    ),
                    dmc.ColorInput(
                        id="grid-line-color",
                        label="Grid color",
                        value=defaults["grid_line_color"],
                        format="hex",
                        size="sm",
                    ),
                ],
                className="two-column-fields",
            ),
        ],
        gap="sm",
    )


def _export_controls():
    return dmc.Stack(
        [
            _field_label("Surface image"),
            dmc.Select(
                id="export-format",
                label="Output format",
                value="png",
                data=[
                    {"label": "PNG", "value": "png"},
                    {"label": "JPEG", "value": "jpeg"},
                    {"label": "WebP", "value": "webp"},
                    {"label": "SVG", "value": "svg"},
                    {"label": "PDF", "value": "pdf"},
                ],
                allowDeselect=False,
                size="sm",
            ),
            html.Div(
                [
                    dmc.NumberInput(
                        id="export-width",
                        label="Width (px)",
                        value=4096,
                        min=256,
                        max=12000,
                        step=128,
                        allowDecimal=False,
                        debounce=300,
                        size="sm",
                    ),
                    dmc.NumberInput(
                        id="export-height",
                        label="Height (px)",
                        value=2048,
                        min=256,
                        max=12000,
                        step=128,
                        allowDecimal=False,
                        debounce=300,
                        size="sm",
                    ),
                ],
                className="two-column-fields",
            ),
            dmc.Text(
                "4096 x 2048 px",
                id="export-size-status",
                className="field-hint export-size-status",
            ),
            _field_label("Colorbar strip"),
            html.Div(
                [
                    dmc.NumberInput(
                        id="colorbar-width",
                        label="Width (px)",
                        value=1600,
                        min=256,
                        max=12000,
                        step=128,
                        allowDecimal=False,
                        debounce=300,
                        size="sm",
                    ),
                    dmc.NumberInput(
                        id="colorbar-height",
                        label="Height (px)",
                        value=80,
                        min=16,
                        max=512,
                        step=8,
                        allowDecimal=False,
                        debounce=300,
                        size="sm",
                    ),
                ],
                className="two-column-fields",
            ),
            html.Div(
                [
                    dmc.Button(
                        "Save image",
                        id="save-image-button",
                        variant="light",
                        color="teal",
                        size="sm",
                        radius="sm",
                    ),
                    dmc.Button(
                        "Save colorbar",
                        id="save-cbar-button",
                        variant="light",
                        color="teal",
                        size="sm",
                        radius="sm",
                    ),
                    dmc.Button(
                        "Export bskan.in",
                        id="export-bskanin-button",
                        variant="light",
                        color="orange",
                        size="sm",
                        radius="sm",
                    ),
                    dmc.Button(
                        "Reset defaults",
                        id="set-default-button",
                        variant="subtle",
                        color="gray",
                        size="sm",
                        radius="sm",
                    ),
                ],
                className="action-grid",
            ),
        ],
        gap="sm",
    )


def _control_rail(defaults, colormaps, debug_console):
    return html.Aside(
        html.Div(
            [
                html.Div(
                    [
                        html.Div(
                            [
                                dmc.Text("Controls", className="panel-title"),
                                dmc.Text("Render settings", className="panel-kicker"),
                            ]
                        ),
                        dmc.Button(
                            "Update image",
                            id="render-button",
                            color="teal",
                            size="sm",
                            radius="sm",
                            className="render-button",
                        ),
                        dmc.Button(
                            "Close",
                            id="controls-close-button",
                            variant="subtle",
                            color="gray",
                            size="xs",
                            className="controls-close-button",
                        ),
                    ],
                    className="control-rail-header",
                ),
                dmc.Tabs(
                    [
                        dmc.TabsList(
                            [
                                dmc.TabsTab("Data", value="data"),
                                dmc.TabsTab("Model", value="simulation"),
                                dmc.TabsTab("Image", value="appearance"),
                                dmc.TabsTab("Overlay", value="overlays"),
                                dmc.TabsTab("Export", value="export"),
                            ],
                            grow=True,
                            className="control-tabs-list",
                        ),
                        dmc.TabsPanel(
                            _data_controls(defaults),
                            value="data",
                            className="control-tab-panel",
                        ),
                        dmc.TabsPanel(
                            _simulation_controls(defaults),
                            value="simulation",
                            className="control-tab-panel",
                        ),
                        dmc.TabsPanel(
                            _appearance_controls(defaults, colormaps),
                            value="appearance",
                            className="control-tab-panel",
                        ),
                        dmc.TabsPanel(
                            _overlay_controls(defaults),
                            value="overlays",
                            className="control-tab-panel",
                        ),
                        dmc.TabsPanel(
                            _export_controls(),
                            value="export",
                            className="control-tab-panel",
                        ),
                    ],
                    id="control-tabs",
                    value="data",
                    keepMounted=True,
                    className="control-tabs",
                ),
                debug_console,
            ],
            className="control-panel",
        ),
        id="control-rail",
        className="control-rail",
        **{"aria-label": "Analysis controls"},
    )


def _analysis_panel(empty_figure, defaults):
    return html.Section(
        dmc.Tabs(
            [
                dmc.TabsList(
                    [
                        dmc.TabsTab("Line profile", value="profile"),
                        dmc.TabsTab("FFT", value="fft"),
                    ],
                    grow=False,
                    className="analysis-tabs-list",
                ),
                dmc.TabsPanel(
                    [
                        html.Div(
                            [
                                dmc.Group(
                                    [
                                        dmc.Switch(
                                            id="line-profile-enabled",
                                            label="Enabled",
                                            checked=defaults["line_profile_enabled"],
                                            color="teal",
                                            size="sm",
                                        ),
                                        dmc.Switch(
                                            id="line-magnet-enabled",
                                            label="Snap to extrema",
                                            checked=defaults["line_magnet_enabled"],
                                            color="orange",
                                            size="sm",
                                        ),
                                        dmc.Button(
                                            "Clear",
                                            id="clear-line-button",
                                            variant="subtle",
                                            color="gray",
                                            size="xs",
                                            radius="sm",
                                        ),
                                        dmc.Button(
                                            "Export CSV",
                                            id="save-line-profile-button",
                                            variant="light",
                                            color="teal",
                                            size="xs",
                                            radius="sm",
                                        ),
                                    ],
                                    gap="md",
                                    className="analysis-actions",
                                ),
                                html.Div(
                                    dcc.Graph(
                                        id="line-profile-graph",
                                        figure=empty_figure(
                                            "Line profile",
                                            "Select two points on the surface map.",
                                            height=280,
                                        ),
                                        config=GRAPH_CONFIG,
                                        responsive=True,
                                        className="analysis-graph",
                                    ),
                                    className="analysis-graph-frame",
                                    role="region",
                                    **{"aria-label": "Line profile plot"},
                                ),
                            ],
                            className="analysis-tab-content",
                        )
                    ],
                    value="profile",
                ),
                dmc.TabsPanel(
                    [
                        html.Div(
                            [
                                dmc.Group(
                                    [
                                        dmc.Switch(
                                            id="fft-enabled",
                                            label="Enabled",
                                            checked=defaults["fft_enabled"],
                                            color="teal",
                                            size="sm",
                                        ),
                                        dmc.Button(
                                            "Save FFT",
                                            id="save-fft-button",
                                            variant="subtle",
                                            color="gray",
                                            size="xs",
                                            radius="sm",
                                        ),
                                    ],
                                    gap="md",
                                    className="analysis-actions",
                                ),
                                html.Div(
                                    dcc.Graph(
                                        id="fft-graph",
                                        figure=empty_figure(
                                            "2D FFT",
                                            "Render a surface map first.",
                                            height=280,
                                        ),
                                        config=GRAPH_CONFIG,
                                        responsive=True,
                                        className="analysis-graph",
                                    ),
                                    className="analysis-graph-frame",
                                    role="region",
                                    **{"aria-label": "Two-dimensional Fourier transform plot"},
                                ),
                            ],
                            className="analysis-tab-content",
                        )
                    ],
                    value="fft",
                ),
            ],
            value="profile",
            id="analysis-tabs",
            keepMounted=True,
            color="teal",
            className="analysis-tabs",
        ),
        className="panel analysis-panel",
    )


def _visual_workspace(empty_figure, defaults):
    return html.Div(
        [
            html.Section(
                [
                    html.Div(
                        [
                            html.Div(
                                [
                                    dmc.Text("Surface map", className="panel-title"),
                                    dmc.Text(
                                        "Interactive image", className="panel-kicker"
                                    ),
                                ]
                            ),
                            html.Div(
                                dmc.Group(
                                    [
                                        dmc.Loader(size="xs", color="teal"),
                                        dmc.Text("Rendering", size="xs", fw=700),
                                    ],
                                    gap="xs",
                                ),
                                id="processing-indicator",
                                className="processing-indicator",
                                style={"display": "none"},
                                role="status",
                                **{"aria-live": "polite"},
                            ),
                        ],
                        className="panel-header",
                    ),
                    html.Div(
                        [
                            dcc.Graph(
                                id="main-image",
                                figure=empty_figure(
                                    "Surface map",
                                    "Select CURRENT, PARCHG, LOCPOT, or CHGCAR data.",
                                    height=580,
                                ),
                                config=GRAPH_CONFIG,
                                mathjax=True,
                                responsive=True,
                                className="main-graph",
                            ),
                            html.Div(
                                [
                                    dmc.Text(
                                        "Preparing data and updating surface",
                                        fw=700,
                                        size="sm",
                                    ),
                                    dmc.Progress(
                                        value=100,
                                        animated=True,
                                        striped=True,
                                        color="teal",
                                        size="sm",
                                        radius="sm",
                                        className="render-progress",
                                    ),
                                ],
                                id="render-overlay",
                                className="render-overlay",
                                role="status",
                                **{"aria-live": "polite", "aria-label": "Rendering progress"},
                            ),
                        ],
                        className="main-graph-loading",
                        role="region",
                        **{"aria-label": "Interactive surface map"},
                    ),
                    html.Div(
                        dcc.Graph(
                            id="colorbar-graph",
                            figure=empty_figure(
                                "Color scale",
                                "Available after rendering.",
                                height=145,
                            ),
                            config={**GRAPH_CONFIG, "displayModeBar": False},
                            mathjax=True,
                            responsive=True,
                            className="colorbar-graph",
                        ),
                        className="colorbar-strip",
                        role="region",
                        **{"aria-label": "Surface map color scale"},
                    ),
                ],
                className="panel canvas-panel",
            ),
            _analysis_panel(empty_figure, defaults),
        ],
        className="visual-workspace",
    )


def build_layout(
    *,
    defaults: dict,
    colormaps: list[str],
    version_label: str,
    debug_mode: bool,
    empty_figure: Callable,
    status_text: Callable,
    status_wrap: Callable,
):
    stores = [
        dcc.Store(id="line-points-store", data=[]),
        dcc.Store(id="line-display-points-store", data=[]),
        dcc.Store(id="plot-meta-store", data={}),
        dcc.Store(id="surface-key-store", data=""),
        dcc.Store(id="appearance-store", data={}),
        dcc.Store(
            id="isosurface-render-store",
            data={
                "isosurface": -2.5,
                "gaussian_sigma": defaults["gaussian_sigma"],
                "atom_radius_scale": defaults["atom_radius_scale"],
                "fit_radius": defaults["fit_radius"],
            },
        ),
        dcc.Store(id="isosurface-range-default-store", data=-2.5),
        dcc.Store(id="render-request-store", data={}),
        dcc.Store(id="controls-open-store", data=False),
        dcc.Store(id="controls-keyboard-store", data={}),
        dcc.Store(id="volume-file-store", data=""),
        dcc.Store(id="structure-file-store", data=""),
        dcc.Store(id="volume-upload-ready-store", data=""),
        dcc.Store(id="structure-upload-ready-store", data=""),
        dcc.Store(id="volume-upload-filename-store", data=""),
        dcc.Store(id="structure-upload-filename-store", data=""),
        dcc.Store(id="embedded-structure-store", data={}),
        dcc.Store(id="structure-override-pending-store", data={}),
        dcc.Store(id="structure-override-decision-store", data={}),
        dcc.Store(id="startup-log-store", data={}),
        dcc.Interval(id="startup-ping", interval=400, max_intervals=1, n_intervals=0),
        dcc.Download(id="download-image"),
        dcc.Download(id="download-colorbar"),
        dcc.Download(id="download-line-profile"),
        dcc.Download(id="download-fft-image"),
        dcc.Download(id="download-bskanin"),
    ]
    debug_console = _debug_console(debug_mode, status_text, status_wrap)

    return dmc.MantineProvider(
        theme={
            "fontFamily": "Inter, ui-sans-serif, system-ui, -apple-system, Segoe UI, sans-serif",
            "primaryColor": "teal",
            "defaultRadius": "sm",
            "cursorType": "pointer",
        },
        children=html.Div(
            [
                *stores,
                html.Header(
                    [
                        html.Div(
                            [
                                dmc.Title("AutoBSKAN", order=1, className="brand-name"),
                                dmc.Badge(
                                    version_label,
                                    color="gray",
                                    variant="light",
                                    size="sm",
                                    radius="sm",
                                    className="version-badge",
                                ),
                            ],
                            className="brand-row",
                        ),
                        dmc.Text(
                            "STM, apparent barrier, and local work function analysis",
                            className="brand-subtitle",
                        ),
                        dmc.Group(
                            [
                                html.A(
                                    "New window",
                                    id="new-window-button",
                                    href="/",
                                    target="_blank",
                                    rel="noopener noreferrer",
                                    className="new-window-link",
                                    **{"aria-label": "Open an independent analysis window"},
                                ),
                                dmc.Button(
                                    "Controls",
                                    id="controls-toggle-button",
                                    variant="light",
                                    color="teal",
                                    size="xs",
                                    radius="sm",
                                    className="mobile-controls-trigger",
                                ),
                            ],
                            gap="xs",
                            className="header-actions",
                        ),
                    ],
                    className="app-header",
                ),
                html.Div(
                    "",
                    id="user-error-banner",
                    className="user-error-banner",
                    role="alert",
                    **{"aria-live": "assertive", "aria-atomic": "true"},
                ),
                html.Main(
                    [
                        _visual_workspace(empty_figure, defaults),
                        _control_rail(defaults, colormaps, debug_console),
                        html.Button(
                            "",
                            id="controls-backdrop",
                            type="button",
                            **{"aria-label": "Close controls"},
                            className="controls-backdrop",
                        ),
                    ],
                    id="app-interactive-wrap",
                    className="workspace-grid",
                    style={"pointerEvents": "auto"},
                ),
                html.Footer(
                    [
                        dmc.Text("AutoBSKAN", fw=700),
                        html.Div(
                            [
                                dmc.Text("Materials Theory Group@Yonsei University"),
                                dmc.Text(
                                    "Condensed matter theory@The University of Sydney"
                                ),
                            ],
                            className="footer-affiliations",
                        ),
                    ],
                    className="app-footer",
                ),
                dmc.Modal(
                    [
                        dmc.Text(
                            "The VASP volumetric file already contains an embedded "
                            "structure. Proceeding will replace it with the selected "
                            "structure for overlays and cell-dependent processing.",
                            id="structure-override-message",
                            size="sm",
                        ),
                        dmc.Group(
                            [
                                dmc.Button(
                                    "Use embedded structure",
                                    id="structure-override-cancel-button",
                                    variant="subtle",
                                    color="gray",
                                ),
                                dmc.Button(
                                    "Proceed",
                                    id="structure-override-proceed-button",
                                    color="teal",
                                ),
                            ],
                            justify="flex-end",
                            mt="md",
                        ),
                    ],
                    id="structure-override-modal",
                    title="Replace embedded VASP structure?",
                    opened=False,
                    centered=True,
                    withCloseButton=False,
                    closeOnClickOutside=False,
                    closeOnEscape=False,
                    radius="sm",
                ),
            ],
            className="app-shell",
        ),
    )
