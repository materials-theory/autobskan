import numpy as np
import os, copy

# Defining bskan_input class.
def raw_input(bskan_in='bskan.in'):
    bskan_in = open(bskan_in, 'r')
    option_dict_raw = {}
    for line in bskan_in:
        if "#" in line:
            new_line = ""
            for letter in line:
                if letter != "#":
                    new_line += letter
                else:
                    break
            if new_line != "":
                line = new_line
        try:
            a, b = line.split("=")
            a = "".join(a.upper().split())
            if a in ["BSKAN", "VASP", "CURRENT", "POSCAR", "ATOM_ADDINFO", "CMAP", "EXT"]:
                option_dict_raw[a] = b.strip()
            else:
                option_dict_raw[a] = "".join(b.upper().split())
        except:
            pass
    return option_dict_raw

def parse_input(raw_options): # 파일 인풋 형식에 따라 구분해둠
    '''raw_options : type_dictionary'''
    option_dict = {}
    for i in ["MODE", "METHOD", "CMAP", "BSKAN", "VASP", "BLUR_METHOD", "CURRENT", "EXT", "POSCAR", "RADIUS_TYPE", "ATOM_ADDINFO"]:
        # 처리 필요 없는 경우. 모든 공백이 제거되어서 나오고, BSKAN과 VASP과 같은 대소문자 구분이 필요한 경우는 raw_input에서 앞뒤 공백만 잘려서 나옴
        if i in raw_options:
            option_dict[i] = raw_options[i]
    for i in ["BIAS", "ISO"]: # _와 &를 사용해 여러 값을 넣는 경우
        if i in raw_options:
            option_dict[i] = str_to_range(raw_options[i])
    for i in ["ISO_AUTO", "POST_PROCESSING"]: # True or False를 받는 경우
        if i == "ISO_AUTO" and raw_options[i].startswith("L"):
            option_dict[i] = "LOGSCALE"
            # LOGSCALE 추가됨. (2020.08.29) LOGSCALE이 아닌 경우 TorF 실행
        elif i in raw_options:
            option_dict[i] = TorF(raw_options[i])
    for i in ["ITERATION"]:
        if i in raw_options:
            option_dict[i] = list(map(lambda x:int(x), raw_options[i].split(",")))
    for i in ["BRIGHTNESS", "GAMMA", "CONTRAST"]: # float값만 받는 경우
        if i in raw_options:
            option_dict[i] = float(raw_options[i])
    for i in ["BLUR_SIGMA", "CONTOUR_RESOLUTION", "LAYERS", "SIZE_RATIO"]: # int값만 받는 경우
        if i in raw_options:
            option_dict[i] = int(raw_options[i])
    # for i in ["CURRENT"]:  # 정규표현식 사용 가능케하면서 & 구분은 아래에서 수행하게 됨.
    #     if i in raw_options:
    #         option_dict[i] = raw_options[i].split("&")
    return option_dict

def str_to_range(X):
    '''
    make array from string.
    When the data separated by "_", it regards as range.
    When the data separated by "&", it adds single data.
    "a_b_c&d" -> [a~b(step=c), d]
    ex1) "-0.01_0.01_0.01&0.02" -> [-0.01, 0.00, 0.01, 0.02]
    ex2) "1e3&1e4&1e5"      -> [1000., 10000., 100000.]
    '''
    from numpy import array, arange, sort, ndarray, append
    if type(X) != ndarray:
        if len(X.split("&")) == 1:
            if len(X.split("_")) == 3:
                X = list(map(lambda y: float(y), X.split("_")))
                if X[2] + 1 == 1:
                    X = array([X[0]])
                else:
                    X_temp = arange(X[0], X[1], X[2])
                    if X[1] not in X_temp:
                        X = append(X_temp, X[1])
                    else:
                        X = X_temp.copy()
            else:
                X = [float(X)]
        else:
            X_split = X.split("&")
            X = array([])
            for y in X_split:
                if len(y.split("_")) == 3:
                    y_range = list(map(lambda y: float(y), y.split("_")))
                    y_range_temp = arange(y_range[0], y_range[1], y_range[2])
                    if y_range[1] not in y_range_temp:
                        y_range = append(y_range_temp, y_range[1])
                    else:
                        y_range = y_range_temp.copy()
                    X = append(X, y_range)

                else:
                    X = append(X, float(y))
    X = sort(X)
    return X

def TorF(boolean):
    boolean = boolean.upper()
    if boolean == ".TRUE." or boolean.startswith("T"):
        return True
    else:
        return False

class Bskan_input:
    def __init__(self, filename="bskan.in"):
        '''
        filename can be either filename or dictionary type.
        '''
        # Default Setting
        self.mode = "IMAGE"
        #for VASP & bSKAN calculation
        self.vasp = None
        self.bskan = None
        self.method = None
        self.bias = None

        # for Image making
        self.current = "CURRENT"
        self.iso_auto = "LOGSCALE"
        self.iso = [5]
        self.cmap = "afmhot"
        self.contrast = 0
        self.brightness = 0
        self.ext = "png"
        self.poscar = None        
        self.atom_addinfo = None
        self.layers = 1
        self.radius_type = "ATOMIC"
        self.size_ratio = 60
        self.contour_resolution = 200

        # for Postprocessing
        self.post_processing = False
        self.iteration = [4, 4]
        self.blur_method = "GAUSSIAN"
        self.blur_sigma = 10
        self.gamma = 90


        if type(filename) in [dict, str]:
            if type(filename) == str:
                if os.path.exists(filename):
                    self.options = parse_input(raw_input(filename))
                else:
                    raise IOError("Wrong input filename. Default = bskan.in")
            else:
                self.options = parse_input(filename)
            if "MODE" in self.options:
                if self.options["MODE"].upper().startswith("CA"):
                    self.mode = "CALCULATION"
                # elif self.options["MODE"].upper().startswith("TE"):
                #     self.mode = "TEST"
                elif self.options["MODE"].upper().startswith("PO"):
                    self.mode = "POST_PROCESSING"
                else:
                    self.mode = "IMAGE"

            # (1) Calculation Parts
            if "VASP" in self.options:
                self.vasp = self.options["VASP"]
            if "BSKAN" in self.options:
                self.bskan = self.options["BSKAN"]
            if "METHOD" in self.options:
                tmp = self.options["METHOD"].upper()
                if tmp.startswith("TE") or tmp.startswith("TH"):
                    self.method = "TH"
                elif tmp.startswith("CH"):
                    self.method = "CHEN"
                elif tmp.startswith("BA") or tmp.startswith("NU"):
                    self.method = "BARDEEN"
                else:
                    raise IOError("Input Method is not supported. Choose among TH/CHEN/BARDEEN")
            if "BIAS" in self.options:
                self.bias = self.options["BIAS"]

            # (2) Image Parts
            if "CURRENT" in self.options:
                if "&" in self.options["CURRENT"]:
                    self.current = self.options["CURRENT"].split("&")
                elif self.options["CURRENT"].upper().startswith("AL"):
                    self.current = "ALL"
                else:
                    self.current = self.options["CURRENT"]
            if "ISO_AUTO" in self.options:
                self.iso_auto = self.options["ISO_AUTO"]
            if "ISO" in self.options:
                self.iso = self.options["ISO"]
            if "CMAP" in self.options:
                self.cmap = self.options["CMAP"].lower()
            if "CONTRAST" in self.options:
                self.contrast = self.options["CONTRAST"]
            if "CONTOUR_RESOLUTION" in self.options:
                self.contour_resolution = self.options["CONTOUR_RESOLUTION"]
            if "BRIGHTNESS" in self.options:
                self.brightness = self.options["BRIGHTNESS"]                
            if "EXT" in self.options:
                self.ext = self.options["EXT"]
            if "POSCAR" in self.options:
                self.poscar = self.options["POSCAR"]
            if "ATOM_ADDINFO" in self.options:
                self.atom_addinfo = self.options["ATOM_ADDINFO"]                
            if "LAYERS" in self.options:
                self.layers = self.options["LAYERS"]
            if "RADIUS_TYPE" in self.options:
                tmp = self.options["RADIUS_TYPE"].upper()
                if tmp.startswith("A"):
                    self.radius_type = "ATOMIC"
                elif tmp.startswith("V"):
                    self.radius_type = "VDW"
                elif tmp.startswith("I"):
                    self.radius_type = "IONIC"
                else:
                    raise IOError("wrong input of atomic radius type")

            # (3) Post Processing
            if "POST_PROCESSING" in self.options:
                self.post_processing = self.options["POST_PROCESSING"]
            if "ITERATION" in self.options:
                self.iteration = self.options["ITERATION"]
            if "GAMMMA" in self.options:
                self.gamma = self.options["GAMMA"]
            if "BLUR_METHOD" in self.options:
                if self.options["BLUR_METHOD"].upper().startswith("GAU"):
                    self.blur_method = "GAUSSIAN"
                else:
                    raise IOError("Only Gaussian blurring is supported for now.")
            if "BLUR_SIGMA" in self.options:
                self.blur_sigma = self.options["BLUR_SIGMA"]

    def export(self):
        # TODO: for GUI to CLI
        pass


def main(input_file='bskan.in'):
    return Bskan_input(input_file)

# if __name__ == "__main__":
#     import argparse
#     import numpy as np
#     pars = argparse.ArgumentParser()
#     pars.add_argument('-r', type=str, default='bskan.in', help='input file\'s name. Default = bskan.in')
#     args = pars.parse_args()
#     Bskan_input(args.r)
