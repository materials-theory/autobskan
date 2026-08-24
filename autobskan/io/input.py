import datetime

import numpy as np


def TorF(boolean):
    boolean = boolean.upper()
    if boolean == ".TRUE." or boolean.startswith("T"):
        return True
    elif boolean == ".FALSE." or boolean.startswith("F"):
        return False
    else:
        raise ValueError("Boolean value must be either True or False.")

def str_to_range(X, remove_overlap = True, include_endpoint = True) -> list:
    '''
    make array from string.
    When the data separated by "_", it regards as range.
    When the data separated by ",", it adds single data.
    "a_b_c,d" -> [a~b(step=c), d]
    ex1) "-0.01_0.01_0.01,0.02" -> [-0.01, 0.00, 0.01, 0.02]
    ex2) "1e3,1e4,1e5"      -> [1000., 10000., 100000.]
    '''
    assert isinstance(X, str), f"Input must be string, not {type(X)}: {X}" # TODO: 이거 왜 둔거야? 수정했는디 에러가 나는 경우가 있는지 확인하자

    if "," not in X:
        X_split_ = X.split("_")
        if len(X_split_) == 1:
            # If there is no "_", it regards as single data.
            return [] if X_split_[0] == "" else [float(X_split_[0])]

        X_float = list(map(lambda y: float(y), X_split_))
        if len(X_float) == 2:
            st, ed = X_float
            step = 1
        else:
            assert len(X_float) == 3, "If it is divided by '_' symbol, start and end value must be given"
            assert X_float[2] + 1 != 1, f"Step size must not be zero. Current io: {X}"
            st, ed, step = X_float

        X_range = list(np.arange(st, ed, step))

        if include_endpoint:
            if X_float[1] not in X_range:
                X_range.append(X_float[1])

        return sorted(np.unique(X_range)) if remove_overlap else sorted(X_range)

    else:
        X_split_and = X.split(",")
        return_list = []

        for fragment in X_split_and:
            return_list += str_to_range(fragment, remove_overlap, include_endpoint)

        return sorted(np.unique(return_list)) if remove_overlap else sorted(return_list)

# Import io options from file
def get_kv_pair_from_file(filename:str, read_option:dict):
    fileobj = open(filename, 'r')
    option_dict_raw = {}
    for line_number, line in enumerate(fileobj, start=1):
        if "#" in line:
            if line.strip().startswith("#"):
                continue
            new_line = ""
            for letter in line:
                if letter != "#":
                    new_line += letter
                else:
                    break
            if new_line != "":
                line = new_line

        if line.strip() == "":
            continue

        try:
            key, value = line.split("=", 1)
        except ValueError as exc:
            raise ValueError(
                f"Invalid input line {line_number} in {filename}: expected KEY = VALUE."
            ) from exc
        key = "".join(key.upper().split())
        r = read_option[key]

        value = value.strip()
        if r[0]:
            # Delete spaces
            value = value.replace(" ", "")

        if r[1]:
            # Equate upper/lower cases
            value = value.upper()

        option_dict_raw[key] = value
    fileobj.close()
    return option_dict_raw


def parse_input(option_dict_raw:dict, data_preset:dict):
    raw_options = option_dict_raw.copy()
    option_dict = {}

    for key in raw_options:
        data_value = data_preset[key]

        if data_value is None:
            raw_value = raw_options[key].strip()
            option_dict[key] = (
                None
                if raw_value.upper() in {"", "AUTO", "NONE"}
                else float(raw_value)
            )

        elif isinstance(data_value, list)|isinstance(data_value, tuple):
            dtype = type(data_value[0])
            if dtype in [float, int]:
                option_dict[key] = str_to_range(raw_options[key], remove_overlap = False, include_endpoint = True)
                if dtype is int:
                    option_dict[key] = list(map(int, option_dict[key]))
            else:
                option_dict[key] = [dtype(x) for x in raw_options[key].split(",")]

        elif isinstance(data_value, bool):
            option_dict[key] = TorF(raw_options[key])

        else:
            # str, float, int
            dtype = type(data_value)
            if dtype is int:
                try:
                    option_dict[key] = int(raw_options[key])
                except ValueError:
                    _ = float(raw_options[key])
                    if int(_) != _:
                        print(f"{key} value should be integer, not float. Switched from {_} to {int(_)}.")
                    option_dict[key] = int(float(raw_options[key]))
            else:
                option_dict[key] = dtype(raw_options[key])
    return option_dict

class Options:
    # Dictionary which contains information of all Keywords
    # Keyword : [Read option, Default value]

    # Read options
    #  |___ TT: Delete spaces (T) & equate upper/lower cases (T)
    #  |___ TF: Delete spaces (T), distinguish upper/lower case (F)
    #  |___ FT: Do not delete spaces (F) & equate upper/lower cases (T)
    #  |___ FF: Do not delete spaces (F), distinguish upper/lower case (F)

    TT = True, True
    TF = True, False
    FT = False, True
    FF = False, False

    keys_calc  = ["VASP_COMMAND", "BSKAN_COMMAND", "SCF_PATH", "TIP_PATH", "METHOD", "BIAS"]
    keys_image = ["SIMULATION", "INPUT_SOURCE", "CURRENT", "VOLUME", "FERMI_LEVEL", "ISO_AUTO", "ISO",
                  "FIT_RADIUS", "CMAP", "CONTRAST", "BRIGHTNESS", "EXT",
                  "IMAGE_MODE", "CONTOUR_RESOLUTION"]
    keys_post = ["POSCAR", "ITERATION", "BLUR_METHOD", "BLUR_SIGMA", "GAMMA", "DISPLAY_CBAR", "DISPLAY_CBAR_ON_SEPARATE_PLOT",
                 "DISPLAY_ATOMS", "LAYERS", "DIAGONAL_TRANSFORM_FOR_HEXAGONAL", "ATOMS_INFO", "ATOM_ADDINFO",
                 "RADIUS_TYPE", "SIZE_RATIO", "DISPLAY_CELL"]

    keywords_info = dict(
        TASK = [TT, "IMAGE"],

        # for VASP & bSKAN calculation
        VASP_COMMAND = [FF, "Needs_to_be_specified when TASK = CALCULATION"],
        # GVEC_PARAM = [],
        BSKAN_COMMAND = [FF, "Needs_to_be_specified when TASK = CALCULATION"],
        SCF_PATH = [FF, "Needs_to_be_specified when TASK = CALCULATION"],
        TIP_PATH = [FF, "Needs_to_be_specified when TASK = CALCULATION"],
        METHOD = [TT, "CHEN"],
        BIAS = [FF, [0.0]],

        # for Image making
        SIMULATION = [TT, "STM"],
        INPUT_SOURCE = [TT, "BSKAN"],
        IMAGE_MODE=[FT, "CONSTANT CURRENT"],
        CURRENT = [FF, "CURRENT"],
        VOLUME = [FF, ""],
        FERMI_LEVEL = [TT, None],
        ISO_AUTO = [TT, "LOGSCALE"], # LOGSCALE function added (2020.08.29)
        ISO = [TT, [5.0]],
        FIT_RADIUS = [TT, 0.5],
        CMAP = [FF, "afmhot"],
        CONTRAST = [TT, 0.0],
        BRIGHTNESS = [TT, 0.0],
        EXT = [TT, "PNG"],
        CONTOUR_RESOLUTION=[TT, 200],

        # for Postprocessing
        POSCAR = [FF, "POSCAR"],
        ITERATION = [FF, [1, 1]],
        GAMMA = [TT, 90.0],
        DIAGONAL_TRANSFORM_FOR_HEXAGONAL = [TT, True],

        DISPLAY_ATOMS = [TT, False],
        ATOMS_INFO = [TT, "VESTA"], # TODO: 이거 제대로 안 해뒀음. VESTA, ASE 등 뭐 쓰는겨?
        ATOM_ADDINFO = [FF, "elements.ini"],
        LAYERS = [TT, 1],
        RADIUS_TYPE = [TT, "ATOMIC"],
        SIZE_RATIO = [TT, 30],

        DISPLAY_CELL  = [TT, False],
        DISPLAY_CBAR  = [TT, False],
        DISPLAY_CBAR_ON_SEPARATE_PLOT = [TT, False],
        BLUR_METHOD = [TT, "GAUSSIAN"],
        BLUR_SIGMA = [TT, 0.0],
    )

    atoms = None # ase.Atoms object
    current = None # autobskan.image.stmplot.Current object
    elements_info = None # list of elements info: colors, radius, and so on

    def __init__(self, filename = None):
        self._filename = filename
        ro, dt = {}, {}
        for key in self.keywords_info:
            ro[key] = self.keywords_info[key][0]
            dt[key] = self.keywords_info[key][1]

        if filename is not None:
            self.option_dict_raw = get_kv_pair_from_file(filename, ro)
            self.option_dict = parse_input(self.option_dict_raw, dt)
        else:
            self.option_dict_raw = None
            self.option_dict = {}
            for key in self.keywords_info:
                self.option_dict[key] = self.keywords_info[key][1]

        # For default values
        for keyword in self.keywords_info:
            if keyword not in self.option_dict:
                self.option_dict[keyword] = self.keywords_info[keyword][1]

        # More detailed modification in options

        ## Autobskan task selection
        task = self.option_dict["TASK"].upper()
        if task.startswith("C"):
            self.option_dict["TASK"] = "CALCULATION"
        elif task.startswith("P"):
            self.option_dict["TASK"] = "POST_PROCESSING_ONLY"
        # elif self.options["TASK"].upper().startswith("TE"):
        #     self.task = "TEST"
        else:
            self.option_dict["TASK"] = "IMAGE"

        ## (2) Method; STM simulation methodolgy selection
        bskan_method = self.option_dict["METHOD"].upper()
        if bskan_method.startswith("T"):
            self.option_dict["METHOD"] = "TH"
        elif bskan_method.startswith("C"):
            self.option_dict["METHOD"] = "CHEN"
        elif bskan_method.startswith("B") or bskan_method.startswith("N"):
            self.option_dict["METHOD"] = "BARDEEN"
        else:
            raise IOError("Input Method is not supported. Choose among TH/CHEN/BARDEEN")

        ## BIAS
        self.option_dict["BIAS"] = np.unique(self.option_dict["BIAS"])

        ## ISO_AUTO and ISO
        current_iso_auto = str(self.option_dict["ISO_AUTO"]).strip().upper()
        if (
            current_iso_auto.startswith("T")
            or current_iso_auto == ".TRUE."
            or current_iso_auto.startswith("LOG")
            or current_iso_auto == "L"
        ):
            self.option_dict["ISO_AUTO"] = "LOGSCALE"
            self.option_dict["ISO"] = self.keywords_info["ISO"][1]
        elif current_iso_auto.startswith("LIN") or current_iso_auto.startswith("EQUAL"):
            self.option_dict["ISO_AUTO"] = "LINEAR"
            if len(self.option_dict["ISO"]) != 1:
                raise ValueError("ISO_AUTO=LINEAR requires one integer image count in ISO.")
            iso_count = float(self.option_dict["ISO"][0])
            if not iso_count.is_integer() or iso_count < 1:
                raise ValueError("ISO_AUTO=LINEAR requires ISO to be a positive integer.")
            self.option_dict["ISO"] = int(iso_count)
        elif current_iso_auto.startswith("F") or current_iso_auto == ".FALSE.":
            self.option_dict["ISO_AUTO"] = False
            self.option_dict["ISO"] = np.unique(self.option_dict["ISO"])
        else:
            raise ValueError(
                "ISO_AUTO must be TRUE/LOGSCALE, LINEAR, or FALSE. "
                f"Received: {self.option_dict['ISO_AUTO']}"
            )

        ## CMAP
        self.option_dict["CMAP"] = str(self.option_dict["CMAP"])

        ## Simulation type
        simulation = self.option_dict["SIMULATION"].upper()
        if simulation.startswith("P") or simulation.startswith("APP"):
            self.option_dict["SIMULATION"] = "PHI_APP"
        elif simulation.startswith("L"):
            self.option_dict["SIMULATION"] = "LWF"
        else:
            self.option_dict["SIMULATION"] = "STM"

        ## Input source
        input_source = self.option_dict["INPUT_SOURCE"].upper()
        if input_source.startswith("V"):
            self.option_dict["INPUT_SOURCE"] = "VASP"
        else:
            self.option_dict["INPUT_SOURCE"] = "BSKAN"

        if not self.option_dict["VOLUME"]:
            self.option_dict["VOLUME"] = self.option_dict["CURRENT"]

        ## IMAGE MODE | CONSTANT CURRENT or CONSTANT HEIGHT ?
        image_mode = self.option_dict["IMAGE_MODE"].upper()
        if len(image_mode.split())==2:
            if image_mode.split()[-1].startswith("CU"):
                self.option_dict["IMAGE_MODE"] = "CONSTANT CURRENT"
            elif image_mode.split()[-1].startswith("H"):
                self.option_dict["IMAGE_MODE"] = "CONSTANT HEIGHT"
            else:
                raise IOError(f"Invalid IMAGE_MODE: {image_mode}."
                              f"Choose one between CONSTANT CURRENT or CONSTANT HEIGHT")
        else:
            assert len(image_mode.split())==1, (f"Invalid IMAGE_MODE: {image_mode}."
                                                  f"Choose one of 'CONSTANT CURRENT' or 'CONSTANT HEIGHT' mode")
            if image_mode.startswith("CU") or image_mode.startswith("CC"):
                self.option_dict["IMAGE_MODE"] = "CONSTANT CURRENT"
            elif image_mode.startswith("H"):
                self.option_dict["IMAGE_MODE"] = "CONSTANT HEIGHT"
            else:
                raise IOError(f"Invalid IMAGE_MODE: {image_mode}."
                              f"Choose one between CONSTANT CURRENT or CONSTANT HEIGHT")

        ## RADIUS TYPE
        current_radius_type = self.option_dict["RADIUS_TYPE"].upper()
        if current_radius_type.startswith("A"):
            self.option_dict["RADIUS_TYPE"] = "ATOMIC"
        elif current_radius_type.startswith("V"):
            self.option_dict["RADIUS_TYPE"] = "VDW"
        elif current_radius_type.startswith("I"):
            self.option_dict["RADIUS_TYPE"] = "IONIC"
        else:
            raise IOError("Wrong input of atomic radius type")


        ## ITERATION
        it_len = len(self.option_dict["ITERATION"])
        assert it_len <= 2, "Iteration only along a and b axis is available"
        if it_len == 1:
            self.option_dict["ITERATION"] = [self.option_dict["ITERATION"][0]] * 2

        ## BLUR_METHOD
        assert self.option_dict["BLUR_METHOD"].upper().startswith("GAU"), "Only Gaussian blurring is supported for now"

    def __repr__(self):
        _repr = f"Options parsed from: {self._filename}\n\n"
        _repr += "Input options: \n"
        for key, value in sorted(self.option_dict.items()):
            _repr += f"* {key}: {repr(value)}\n"
        return _repr

    def export_to_file(self, filename = "bskan.in_exported", openmode="w", indent = 0):
        assert isinstance(indent, int), "Indent must be integer"
        with open(filename, openmode) as fileobj:
            task = self.option_dict["TASK"]
            print(f"### Exported bskan.in at {datetime.datetime.now()}", file=fileobj)
            print(f"TASK = {task}", file=fileobj)
            _header_calc_written = False
            _header_image_written = False
            _header_post_processing_written = False

            for keyword in self.keywords_info:
                if keyword == "TASK":
                    continue

                k_info = self.keywords_info[keyword]
                value = self.option_dict[keyword] if keyword in self.option_dict else k_info[1]

                _write = False
                if keyword in self.keys_calc:
                    if not _header_calc_written:
                        _header_calc_written = True
                        print("\n### Calculation", file=fileobj)
                    _write = True if task=="CALCULATION" else False
                elif keyword in self.keys_image:
                    if not _header_image_written:
                        _header_image_written = True
                        print("\n### Image", file=fileobj)
                    _write = True if task=="IMAGE" else False
                elif keyword in self.keys_post:
                    if not _header_post_processing_written:
                        _header_post_processing_written = True
                        print("\n### Post processing", file=fileobj)
                    _write = True if task in ["IMAGE", "POST_PROCESSING_ONLY"] else False
                else:
                    _write = False
                    raise KeyError(f"Unrecognized keyword {keyword}")

                if _write:
                    ## Should I skip writing the default settings?
                    # if value == self.keywords_info[keyword][1]:
                    #     continue
                    if isinstance(value, list) or isinstance(value, np.ndarray):
                        value = [str(v) for v in value]
                        value = ", ".join(value)
                    print(f"{' '*indent}{keyword} = {value}", file=fileobj)

def main(input_file='bskan.in'):
    return Options(input_file)
