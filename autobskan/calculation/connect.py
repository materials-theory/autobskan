import json, time, os, subprocess, math, re, shutil
from glob import glob
import numpy as np

from autobskan.image import AR
from ase.io.vasp import read_vasp, write_vasp
from ase import Atoms, Atom
from ase.geometry import get_layers
from ase.calculators.vasp.create_input import GenerateVaspInput as GVI

from autobskan.io.input import Options


##_functions_##

def gvector(model):
    AUTOA = 0.529177249
    RYTOEV = 13.605826
    HSQDTM = RYTOEV * AUTOA**2
    TPI = 2.0 * math.pi

    def fft_indices(n):
        half = n // 2
        return [i if i <= half else i - n for i in range(n)]

    def read_outcar_info(outcar_path="OUTCAR"):
        encut = None
        ngx = ngy = ngz = None
        kpoints = []

        with open(outcar_path, "r", errors="ignore") as fr:
            lines = fr.readlines()

        for line in lines:
            if encut is None:
                m = re.search(r"ENCUT\s*=\s*([0-9.]+)", line)
                if m is not None:
                    encut = float(m.group(1))

            if ngx is None:
                m = re.search(
                    r"NGX\s*=\s*(\d+)\s+NGY\s*=\s*(\d+)\s+NGZ\s*=\s*(\d+)",
                    line
                )
                if m is not None:
                    ngx, ngy, ngz = [int(x) for x in m.groups()]

        for i, line in enumerate(lines):
            if re.search(r"irreducible k-points:", line):
                for trial in lines[i + 1:]:
                    vals = trial.split()
                    if len(vals) == 4:
                        try:
                            kpoints.append([float(vals[0]), float(vals[1]), float(vals[2])])
                        except:
                            pass
                    elif len(kpoints) > 0:
                        break
                break

        if encut is None:
            raise RuntimeError("gvector: Cannot find ENCUT at OUTCAR")
        if ngx is None:
            raise RuntimeError("gvector: Cannot find NGX/NGY/NGZ at OUTCAR")
        if len(kpoints) == 0:
            kpoints = [[0.0, 0.0, 0.0]]

        return encut, ngx, ngy, ngz, kpoints

    encut, ngx, ngy, ngz, kpoints = read_outcar_info("OUTCAR")

    cell = np.array(model.get_cell(), dtype=float)

    reciprocal_no_2pi = np.linalg.inv(cell).T

    max_valid_energy = 0.0

    for kpoint in kpoints:
        kx, ky, kz = kpoint

        for ix in fft_indices(ngx):
            for iy in fft_indices(ngy):
                g_frac = np.array([ix + kx, iy + ky, 0.0 + kz], dtype=float)
                g_cart = TPI * np.dot(g_frac, reciprocal_no_2pi)
                energy = HSQDTM * np.dot(g_cart, g_cart)

                if energy < encut and energy > max_valid_energy:
                    max_valid_energy = energy

    if max_valid_energy <= 0:
        raise RuntimeError("gvector calculation failed")

    return -round(max_valid_energy, 6)

def stm(POSCAR):
    ans = [0,1,0.052918,3,2.50,2.50,6]
    model = read_vasp(POSCAR)

    ##_surface_height_##

    a, b = get_layers(model, (0,0,1))
    if len(b) == 1:
        ans[0] = b[0] + 0.5
    else:
        for i in range(0,len(b)-1):
            layer_dis = b[1] - b[0]
            if (b[i+1] - b[i]) <=1.1* layer_dis and (b[i+1]-b[i]) >= 0.9*layer_dis:
                surf_atom = i
                continue
            else:
                surf_atom = i
                break
        ans[0] = b[i] + 0.5

    ##_[0]+5.2918_##

    ans[1] = ans[0] + 5.2918

    ##_g_vector_##

    ans[3] = gvector(model)

    ##_fermi_energy_level_##

    Outcar = open("OUTCAR", 'r')
    outcar = Outcar.readlines()
    Outcar.close()
    for line in outcar:
        fermi_energy = re.search('E-fermi', line)
        if fermi_energy != None:
            fermi_line = line
    fenergy = float(fermi_line.split()[2])
    ans[6] = fenergy

    return ans


def inscan_info(model, STm, bias, inp):

    ##_finding_longer_axis_&_grid_calculation_##

    cell = model.get_cell()
    bskan_cell = AR.bskan_cell_transform(cell)
    axis_1 = bskan_cell[0][0]
    axis_2 = bskan_cell[1][1]
    if axis_1 >= axis_2:
        longaxis = axis_1
    else:
        longaxis = axis_2
    grid = math.floor(longaxis * 10 + 1)
    gridpoint = np.min([grid,201])

    ##_INSCAN_data_##

    if inp.option_dict['METHOD']  == "TH":
        method = "TERsoff Hamann"
    elif inp.option_dict['METHOD'] == "CHEN":
        method = "CHEN"
    elif inp.option_dict['METHOD'] == "BARDEEN":
        method = "BARDEEN"

    inscan_list=[
    "{0}\n".format(method),
    "BIAS = {0}\n".format(bias),
    "LIMITS = -0.00 0.00\n",
    "# ANGLES = 0 0 0\n",
    "# ORBITALS = 1 0 0 0 0 0 0 0 0\n",
    "NKELDYSH = -1\n",
    "GRIDPOINTS = {0}\n".format(gridpoint),
    "PIVOT = 0.000 0.000\n",
    "CELL = 1.0 1.0\n",
    "ZVACUUM = {0}\n".format(STm)
    ]

    return inscan_list


def th_calc(inp):
    subprocess.call(inp.option_dict['BSKAN_COMMAND'], shell=True)
    this_pwd = os.getcwd()
    os.mkdir("current")
    os.symlink(f"{this_pwd}/ASAMPLE",f"{this_pwd}/current/ASAMPLE")
    shutil.copy(f"{this_pwd}/INSCAN",f"{this_pwd}/current/INSCAN")
    os.symlink(f"{this_pwd}/WAVSAMPLE",f"{this_pwd}/current/WAVSAMPLE")
    os.symlink(f"{this_pwd}/CURMAT",f"{this_pwd}/current/CURMAT")
    os.chdir("current")
    fr = open("./INSCAN",'r')
    lines = fr.readlines()
    fr.close()
    fw = open("./INSCAN",'w')
    for line in lines:
        fw.write(line)
    fw.write("CURRENT = 0.0\n")
    fw.close()
    subprocess.call(inp.option_dict['BSKAN_COMMAND'], shell=True)


def traffic(inpjson_path, calculation_title, job, dict_data):
    original_path = os.getcwd()
    os.chdir(inpjson_path)
    title = f".traffic_light_{time.time()}"
    fw = open(title,'w')
    fw.close()
    light = glob(".traffic_light*")
    light.sort()
    while len(light) > 1:
        light = glob(".traffic_light*")
        light.sort()
        if light.index(title) == 0:
            break
        time.sleep(1)
    dict_data[calculation_title]=job
    save_json(inpjson_path,dict_data)
    if os.path.exists(f"{inpjson_path}/{title}"):
        os.remove(f"{inpjson_path}/{title}")
    os.chdir(original_path)
    return title

def read_json(json_path):
    with open(f"{json_path}/status.json", 'r') as solid_file:
        data = json.load(solid_file)
        return data

def save_json(json_path, dict_data):
    if os.path.exists(f"{json_path}/status.json"):
        data = read_json(json_path)
    else:
        data={}
    for keys in dict_data:
        # if value of key already exists, data will be overlapped.
        data[keys] = dict_data[keys]
    with open(f"{json_path}/status.json", 'w') as solid_file:
        json.dump(data, solid_file, indent=4, sort_keys=True)


def may_I_calculate(calculation_title, json_path):
    original_path = os.getcwd()
    os.chdir(json_path)
    light = glob(".traffic_light*")
    while len(light) > 0:
        light = glob(".traffic_light*")
        time.sleep(1)
    if os.path.exists(f"{json_path}/.{calculation_title}_traffic_light"):
        return False
    else:
        fw = open(f".{calculation_title}_traffic_light",'w')
        fw.close()
    os.chdir(original_path)
    return True

def clean_json():
    current_pwd = os.getcwd()
    try:
        os.mkdir("__status__")
    except:
        pass
    inpjson_path = f"{current_pwd}/__status__"
    os.chdir(inpjson_path)
    data=dict()
    data['start']='calculation'
    if os.path.exists(f"{inpjson_path}/status.json"):
        pass
    else:
        fw = open("status.json",'w')
        json.dump(data,fw,indent=4,sort_keys=True)
        fw.close()
    os.chdir(current_pwd)
    return inpjson_path


def calculation(json_path, calculation_title, inp, method, dict_data):
    traffic(json_path, calculation_title, 'calculating', dict_data)
    if method == 'NONSCF':
        subprocess.call(inp.option_dict['VASP_COMMAND'], shell=True)
    elif method == 'TH':
        th_calc(inp)
    elif method == 'CHEN':
        subprocess.call(inp.option_dict['BSKAN_COMMAND'], shell=True)
    traffic(json_path, calculation_title, 'done', dict_data)


def nonscf_incar(INCAR,POSCAR):
    incar = GVI()
    incar.read_incar(INCAR)
    if incar.int_params['kpar'] != None:
        incar.int_params['kpar'] = None
    if incar.bool_params['lwave'] != False:
        incar.bool_params['lwave'] = False
    if incar.bool_params['lcharg'] != False:
        incar.bool_params['lcharg'] = False
    STM = stm(POSCAR)
    incar.list_float_params['stm'] = STM
    incar.write_incar(atoms=None,directory='./')


def outcar_to_kpoints(OUTCAR):
    fr = open(OUTCAR,'r')
    fw = open("KPOINTS",'w')
    k = None
    kpoint = 0
    i = 0
    for line in fr:
        if k == None:
            k = re.search("irreducible k-points:",line)
        if k != None:
            if kpoint == 0:
                kpoint_line = line
                kpoint = int(line.split()[1])
                fw.write("Explicit k-point list by auto-bskan\n")
                fw.write(f"{kpoint}\n")
                fw.write("reciprocal\n")
                continue
            if len(line.split()) != 4 and i == 1:
                break
            elif len(line.split()) != 4 and i == 0:
                continue
            else:
                fw.write(line)
                i = 1
    fr.close()
    fw.close()


def contcar_to_asample(CONTCAR):
    poscar_fr = open(CONTCAR, 'r')
    poscar_line = poscar_fr.readlines()
    poscar_fr.close()
    poscar_fw = open('ASAMPLE','w')
    line_count = 1
    for line in poscar_line:
        if len(line) > 1:
            if line_count <= 5:
                poscar_fw.write(line)
                line_count += 1
            elif line_count == 6:
                line_count += 1
            elif line_count == 7:
                poscar_fw.write(line)
                line_count += 1
            elif line_count == 8:
                if line.split(maxsplit=len(line))[0] == 'S':
                    poscar_fw.write(line)
                    line_count += 1
                else :
                    poscar_fw.write("Selective dynamics\n")
                    poscar_fw.write(line)
                    line_count += 1
            elif line_count == 9:
                poscar_fw.write(line)
    poscar_fw.close()


def bias_directory(inp, model, STM, until_tip, th_for_chen=False):
    bias_dir = inp.option_dict['BIAS']
    if th_for_chen == True:
        bias_dir = [0]
    for bias in bias_dir:
        try:
            os.chdir(f"{until_tip}")
            bias_title = f"bias_{bias}"
            os.mkdir(bias_title)
            os.symlink(f"{until_tip}/ASAMPLE",f"{until_tip}/{bias_title}/ASAMPLE")
            os.symlink(f"{until_tip}/WAVSAMPLE",f"{until_tip}/{bias_title}/WAVSAMPLE")
            os.chdir(bias_title)
            STm = STM[0]
            inscan_fw = open('./INSCAN', 'w')
            for line in inscan_info(model, STm, bias, inp):
                inscan_fw.write(line)
            inscan_fw.close()
            os.chdir("..")
        except:
            continue

def ctoc():
    fr = open("CURSAVE", 'r')
    cursave_line = fr.readlines()
    fr.close()
    cur_pos = list()
    fr = open('./current', 'r')
    for i in range(0, 6):
        a = fr.readline()
        cur_pos.append(a)
    atom_number = int(a)
    while True:
        b = fr.readline()
        if len(b.split()) >= 5:
            break
        cur_pos.append(b)
        c = b
    fr.close()
    n_x, n_y, n_z = c.rstrip().split()
    nx = int(n_x)
    ny = int(n_y)
    nz = int(n_z)
    cursave_info = list()
    i = 0
    a = 0
    while i < len(cursave_line)-1:
        if len(cursave_line[i]) == 21:
            x = float(cursave_line[i].split()[0])
            y = float(cursave_line[i].split()[1])
            cursave_info.append([x,y,[]])
            a += 1
        i += 1
        while len(cursave_line[i]) == 13:
            cursave_info[a-1][2].append(float(cursave_line[i].strip()))
            i += 1
            if i == len(cursave_line):
                break
    current_info = list()
    for z in range(0, nz):
        for i in range(0, ny):
            b = i
            for a in range(0, nx):
                current_info.append(cursave_info[b][2][z])
                b += ny
    current_fw = open("CURRENT", 'w')
    for line in cur_pos:
        current_fw.write(line)
    i=0
    while i < len(current_info):
        info = "    " + "%.4e"%current_info[i]
        current_fw.write(info)
        i += 1
        if i%5 == 0:
            current_fw.write("\n")
    current_fw.close()

###############################################################################

def main(bskan_input = "bskan.in"):

    inp = Options(bskan_input)

    ##_current_pwd_##

    current_pwd = os.getcwd()

    ##_status.json_##

    inpjson_path = clean_json()

    ##_NONSCF_directory_&_input_files_##

    dict_data=dict()
    model = read_vasp(f"{inp.option_dict['SCF_PATH']}/POSCAR")

    os.chdir(f"{inp.option_dict['SCF_PATH']}")
    STM = stm("POSCAR")
    os.chdir(f"{current_pwd}")

    try:
        os.mkdir("./1_nonscf")
        os.chdir("./1_nonscf")
        os.symlink(f"{inp.option_dict['SCF_PATH']}/POSCAR",f"{current_pwd}/1_nonscf/POSCAR")
        os.symlink(f"{inp.option_dict['SCF_PATH']}/POTCAR",f"{current_pwd}/1_nonscf/POTCAR")
        shutil.copy(f"{inp.option_dict['SCF_PATH']}/OUTCAR",f"{current_pwd}/1_nonscf/OUTCAR")
        os.symlink(f"{inp.option_dict['SCF_PATH']}/INCAR",f"{current_pwd}/1_nonscf/_INCAR")
        os.symlink(f"{inp.option_dict['SCF_PATH']}/WAVECAR",f"{current_pwd}/1_nonscf/WAVECAR")
        os.symlink(f"{inp.option_dict['SCF_PATH']}/CHGCAR",f"{current_pwd}/1_nonscf/CHGCAR")

        if os.path.exists(f"{inp.option_dict['SCF_PATH']}/KPOINTS"):
            os.symlink(f"{inp.option_dict['SCF_PATH']}/KPOINTS",f"{current_pwd}/1_nonscf/KPOINTS")
        else:
            pass

        ##_NONSCF_INCAR_##

        nonscf_incar("./_INCAR","./POSCAR")

        ##_NONSCF_KPOINTS_##

        try:
            fr = open("KPOINTS",'r')
            fr.close()
            pass
        except:
            outcar_to_kpoints("OUTCAR")

    except:
        os.chdir(f'{current_pwd}/1_nonscf')
        while True:
            files_list = glob("*")
            all_in_path = True
            for f in ["POSCAR", "POTCAR", "INCAR", "WAVECAR", "CHGCAR", "KPOINTS", "OUTCAR"]:
                if f not in files_list:
                    all_in_path = False
            if all_in_path:
                os.chdir(f"{current_pwd}")
                break
            else:
                time.sleep(1)

    ##_NONSCF_calculation_##

    os.chdir(f"{current_pwd}/1_nonscf")
    calculation_title = "NONSCF"
    answer = may_I_calculate(calculation_title,inpjson_path)

    if answer == True:
        calculation(inpjson_path, calculation_title, inp,'NONSCF', dict_data)

    else:
        while True:
            if os.path.exists(f"{current_pwd}/2_bskan"):
                break
            else:
                time.sleep(1)

    ##_bskan_directory_##

    os.chdir(f"{current_pwd}")

    stm_title = f"STM_{inp.option_dict['METHOD']}"
    if inp.option_dict['METHOD'] == "TH":
        tip_title = "NO_TIP"
    elif inp.option_dict['METHOD'] == "CHEN" or inp.option_dict['METHOD'] == "BARDEEN":
        tip_title = inp.option_dict['TIP_PATH'].split('/')[-1]

    until_tip = f"{current_pwd}/2_bskan/{stm_title}/{tip_title}"

    try:
        os.mkdir("2_bskan")
        os.mkdir("3_result")
        os.chdir("2_bskan")

        ##_STM_directory_##

        os.mkdir(stm_title)
        os.chdir(stm_title)

        ##_tip_directory_##

        os.mkdir(tip_title)
        os.chdir(tip_title)

        ##_STM_file_preparation_without_INSCAN_(ASAMPLE,WAVSAMPLE)_##

        write_vasp("CONTCAR_to_read", model, direct = True)
        contcar_to_asample("CONTCAR_to_read")
        os.remove("CONTCAR_to_read")

        os.symlink(f"{current_pwd}/1_nonscf/STM",f"{until_tip}/WAVSAMPLE")
    except:
        pass

    try:
        os.chdir(f"{until_tip}")
        bias_directory(inp,model,STM,until_tip)
    except:
        while True:
            try:
                os.chdir(until_tip)
            except:
                continue
            a = 0
            while a < len(inp.option_dict['BIAS']):
                time.sleep(1)
                a = 0
                bias_list = glob("bias*")
                for b in inp.option_dict['BIAS']:
                    if b in bias_list:
                        a += 1
            for bias in bias_list:
                os.chdir(f'./{bias}')
                while True:
                    file_list = glob("*")
                    if 'WAVSAMPLE' in file_list and 'ASAMPLE' in file_list and 'INSCAN' in file_list:
                        break
                    else:
                        time.sleep(1)
                os.chdir("..")
            break
            os.chdir(f"{current_pwd}")


    ##_calculation_for_TH_##

    if inp.option_dict['METHOD'] == "TH":
        cal_list = glob("bias*")
        for cal in cal_list:
            os.chdir(f"{until_tip}/{cal}")
            calculation_title = "STM" + f"_{cal}"
            answer = may_I_calculate(calculation_title,inpjson_path)
            if answer == True:
                calculation(inpjson_path, calculation_title, inp, 'TH', dict_data)
                os.symlink(f"{until_tip}/{cal}/current/CURRENT",f"{current_pwd}/3_result/{inp.option_dict['METHOD']}_{tip_title}_{cal}_CURRENT")
            else:
                time.sleep(1)
                continue

    ##_calculation_for_CHEN_##

    if inp.option_dict['METHOD'] == "CHEN":

        ##_TH_with_bias_0_##
        if os.path.exists(f"{current_pwd}/2_bskan/STM_TH/NO_TIP"):
            pass
        else:
            shutil.copytree(f"{current_pwd}/2_bskan/STM_CHEN/{tip_title}", f"{current_pwd}/2_bskan/STM_TH/NO_TIP", symlinks=True)
            os.chdir(f"{current_pwd}/2_bskan/STM_TH/NO_TIP")
            bias_directory(inp,model,STM,f"{current_pwd}/2_bskan/STM_TH/NO_TIP",th_for_chen=True)
            os.chdir(f"{current_pwd}/2_bskan/STM_TH/NO_TIP/bias_0")
            fr = open("./INSCAN",'r')
            lines = fr.readlines()
            fr.close()
            fw = open("./INSCAN",'w')
            for line in lines:
                if line == "CHEN\n":
                    fw.write("TERsoff Hamann\n")
                else:
                    fw.write(line)
            fw.close()
            calculation_title = "STM_TH_for_CHEN"
            answer = may_I_calculate(calculation_title,inpjson_path)
            if answer == True:
                calculation(inpjson_path, calculation_title, inp, 'TH', dict_data)
            else:
                pass

        ##_ready_for_CHEN_##

        os.chdir(f"{until_tip}")
        cal_list = glob("bias*")
        for cal in cal_list:
            os.chdir(f"{until_tip}/{cal}")
            calculation_title = "STM" + f"_{cal}"
            answer = may_I_calculate(calculation_title,inpjson_path)
            if answer == True:
                os.symlink(f"{inp.option_dict['TIP_PATH']}/ATIP", f"{until_tip}/{cal}/ATIP")
                os.symlink(f"{inp.option_dict['TIP_PATH']}/WAVTIP", f"{until_tip}/{cal}/WAVTIP")
                os.symlink(f"{inp.option_dict['TIP_PATH']}/PROCARtip", f"{until_tip}/{cal}/PROCARtip")
                calculation(inpjson_path, calculation_title, inp, 'CHEN', dict_data)
            else:
                time.sleep(1)
                continue
            if os.path.exists(f"{until_tip}/{cal}/CURRENT"):
                pass
            else:
                while True:
                    if os.path.exists(f"{current_pwd}/2_bskan/STM_TH/NO_TIP/bias_0/current/CURRENT"):
                        os.symlink(f"{current_pwd}/2_bskan/STM_TH/NO_TIP/bias_0/current/CURRENT",f"{until_tip}/{cal}/current")
                        ctoc()
                        break
                    else:
                        time.sleep(1)
                        continue
            os.symlink(f"{until_tip}/{cal}/CURRENT",f"{current_pwd}/3_result/{inp.option_dict['METHOD']}_{tip_title}_{cal}_CURRENT")
