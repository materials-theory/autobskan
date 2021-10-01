# %matplotlib tk
import tkinter as tk
import tkinter.ttk as ttk
from tkinter import filedialog
from tkinter import PhotoImage
import matplotlib.pyplot as plt
from matplotlib.widgets import Cursor, Button
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from ase.io.vasp import read_vasp

from autobskan.input.input import Bskan_input
from autobskan.image import stmplot, post_processing, AR

def main():
    # TODO: bskan_input 인가 쨌든 이 클래스에서 가져와야함
    bskan_input = Bskan_input(None)

    iso_min, iso_max = -1, 1
    iso_suitable = 0 # logarithm 에서 제일 작은 것 중 2번째로 정하자
    max_nlayer = 5

    def open_cur_file():
        curfile.delete(0, tk.END)
        curfile.insert(0, filedialog.askopenfilename())
        current = stmplot.Current(curfile.get())
        bskan_input.current = curfile.get()
        pass

    def open_str_file():
        strfile.delete(0, tk.END)
        strfile.insert(0, filedialog.askopenfilename())
        link_str_file()

    def link_str_file(event=None):
        filename = strfile.get()
        bskan_input.poscar = read_vasp(filename)
        bskan_input.gamma = bskan_input_poscar.cell.cellpar()[-1]

    def export_bskanin():
        pass

    main_window = tk.Tk()
    main_window.title("autobskan GUI")
    main_window.geometry("1280x720")
    main_window.resizable(True, True)


    ######## 1. Ouput images
    image_frame = ttk.LabelFrame(main_window,
                                 text="STM images",
                                 #                             relief="solid",
                                 )
    # image_frame.pack(side="left", fill="y", expand=True)
    image_frame.grid(row=0, column=0, padx=10, pady=10, sticky="news")

    img1_lbl = ttk.Label(image_frame, text="Unit cell")
    img1_lbl.pack(side="top", fill="x", expand=True)
    # img1
    img1 = ttk.Frame(image_frame)
    img1.pack(side="top", fill="x", expand=True)
    temp = ttk.Button(img1)
    temp.pack()

    sep1 = ttk.Separator(image_frame)
    sep1.pack(side="top", fill="x", expand=True)

    img2_lbl = ttk.Label(image_frame, text="Repeated image (nx, ny)")
    img2_lbl.pack(side="top", fill="x", expand=True)
    # img2

    sep2 = ttk.Separator(image_frame)
    sep2.pack(side="top", fill="x", expand=True)

    img3_lbl = ttk.Label(image_frame, text="Blurred image")
    img3_lbl.pack(side="top", fill="x", expand=True)
    # img3


    ####### 2. FILE inputs
    file_frame = ttk.Frame(main_window)
    file_frame.grid(row=0, column=1, padx=10, pady=10, sticky="news")

    filelbl_1 = ttk.LabelFrame(file_frame, text="CURRENT file")
    filelbl_1.grid(row=0, column=0, sticky="news")
    curfile = ttk.Entry(filelbl_1)
    curfile.pack(side="left", fill="both", expand=True)
    curfile_btn = tk.Button(filelbl_1, text="Open",
                            width=2, height=1, command=open_cur_file)
    curfile_btn.pack(side="right", fill="both", expand=True)

    filelbl_2 = ttk.LabelFrame(file_frame, text="Structure file")
    filelbl_2.grid(row=1, column=0, sticky="news")
    strfile = ttk.Entry(filelbl_2)
    strfile.bind("<Return>", link_str_file)
    strfile.pack(side="left", fill="both", expand=True)
    strfile_btn = tk.Button(filelbl_2, text="Open",
                            width=2, height=1, command=open_str_file)
    strfile_btn.pack(side="right", fill="both", expand=True)


    cmap_frame = ttk.LabelFrame(file_frame, text="Colormap")
    cmap_frame.grid(row=2, column=0, sticky="news")
    cmap = ttk.Combobox(cmap_frame, state="readonly", height=5, values=plt.colormaps())
    cmap.pack(side="top", fill="both", expand=True)
    cmap.set("afmhot")

    toggle_frame = ttk.LabelFrame(file_frame, text="View options")
    toggle_frame.grid(row=3, column=0, sticky="news", pady=20)
    showatoms = tk.BooleanVar()
    showunitcell = tk.BooleanVar()
    showcursor = tk.BooleanVar()
    check_showatoms = tk.Checkbutton(toggle_frame, text="Show atoms",
                                     variable=showatoms, anchor="w")
    check_showatoms.pack(side="top", fill="both", expand=True)
    check_showunitcell = tk.Checkbutton(toggle_frame, text="Show unit cell",
                                        variable=showunitcell, anchor="w")
    check_showunitcell.pack(side="top", fill="both", expand=True)
    check_showcursor = tk.Checkbutton(toggle_frame, text="Show cursor",
                                      variable = showcursor, anchor="w")
    check_showcursor.pack(side="top", fill="both", expand=True)

    check_showatoms.select()
    check_showunitcell.select()
    check_showcursor.select()
    # check_showcursor.state(["selected"]) # ttk 썼을 때는 어떻게 select?

    def rp_cmd(event):
        img2_lbl.config(text=f"Repeated image ({repeat.get()})")
    repeat_lbl = ttk.LabelFrame(file_frame, text="Repeat (x, y)")
    repeat_lbl.grid(row=4, column=0, sticky="news", pady=20)
    repeat = ttk.Entry(repeat_lbl)
    repeat.pack()
    repeat.delete(0, tk.END)
    repeat.insert(0, '2, 2')
    img2_lbl.config(text=f"Repeated image ({repeat.get()})")
    repeat.bind("<Return>", rp_cmd)

    def return_to_default():
        # TODO cmd2를 사용해서 값을 원래대로 돌리는 경우는 오케이지만, 다른 경우엔 이미지를 업데이트하는 함수를 만들 필요가 있을 수 있다.
        cmap.set("afmhot")
        repeat.delete(0, tk.END)
        repeat.insert(0, '2, 2')
        rp_cmd(0)
        check_showatoms.select()
        check_showunitcell.select()
        check_showcursor.select()
        brightness.delete(0, tk.END)
        brightness.insert(0, 0.0)
        br_cmd2()
        contrast.delete(0, tk.END)
        contrast.insert(0, 0.0)
        cont_cmd2()
        isosurface.delete(0, tk.END)
        isosurface.insert(0, iso_suitable)
        iso_cmd2()
        gaus_sigma.delete(0, tk.END)
        gaus_sigma.insert(0, 10)
        gaus_cmd2()
        nlayer.delete(0, tk.END)
        nlayer.insert(0, 1)
        nlayer_cmd2()
        ratom.delete(0, tk.END)
        ratom.insert(0, 60)
        ratom_cmd2()

    set_as_default = ttk.Button(file_frame, text="SET AS DEFAULT",
                                command = return_to_default)
    set_as_default.grid(row=5, column=0, sticky="news", pady=20)


    # variable 확인해 봄. 잘 되는구만
    # def tempprt():
    #     print(showatoms.get(), showunitcell.get(), showcursor.get())
    # temp_prtbtn = tk.Button(toggle_frame, text="print", command=tempprt)
    # temp_prtbtn.pack()

    ####### 3. Detailed option manipulation
    opt_frame = ttk.Frame(main_window)
    opt_frame.grid(row=0, column=2, padx=10, pady=10, sticky="news")


    def br_cmd(value):
        # TODO : 여기에 image update 시키는 것도 넣어야 겠네
        brightness.delete(0, tk.END)
        brightness.insert(0, value)

    def br_cmd2():
        # TODO : 여기에 image update 시키는 것도 넣어야 겠네
        br_scale.set(brightness.get())

    def br_cmd3(event):
        # TODO : 여기에 image update 시키는 것도 넣어야 겠네
        br_cmd2()

    br_label = ttk.LabelFrame(opt_frame, text="Brightness")
    br_label.pack(side="top", fill="x", expand=True)
    br_var = tk.DoubleVar()
    br_scale = tk.Scale(br_label, variable = br_var, orient=tk.HORIZONTAL,
                        from_=-1, to=1, resolution=0.1, command = br_cmd,
                        showvalue = False)
    br_scale.pack(side="bottom", fill="x", expand=True)
    # brightness = ttk.Entry(br_label)
    brightness = ttk.Spinbox(br_label, from_=-1, to=1, increment=0.1,
                             command = br_cmd2)
    brightness.bind("<Return>", br_cmd3)
    brightness.delete(0, tk.END)
    brightness.insert(0, 0.0)
    br_cmd2()
    brightness.pack(fill="x", expand=True)


    def cont_cmd(value):
        # TODO : 여기에 image update 시키는 것도 넣어야 겠네
        contrast.delete(0, tk.END)
        contrast.insert(0, value)

    def cont_cmd2():
        # TODO : 여기에 image update 시키는 것도 넣어야 겠네
        cont_scale.set(contrast.get())

    def cont_cmd3(event):
        # TODO : 여기에 image update 시키는 것도 넣어야 겠네
        cont_cmd2()
    cont_label = ttk.LabelFrame(opt_frame, text="Contrast")
    cont_label.pack(side="top", fill="x", expand=True)
    cont_var = tk.DoubleVar()
    cont_scale = tk.Scale(cont_label, variable=cont_var, orient=tk.HORIZONTAL,
                          from_=-1, to=1, resolution=0.1, command = cont_cmd,
                          showvalue=False)
    cont_scale.pack(side="bottom", fill="x", expand=True)
    contrast = ttk.Spinbox(cont_label, from_=-1, to=1, increment=0.1,
                           command = cont_cmd2)
    contrast.bind("<Return>", cont_cmd3)
    contrast.delete(0, tk.END)
    contrast.insert(0, 0.0)
    cont_cmd2
    contrast.pack(fill="x", expand=True)


    def iso_cmd(value):
        # TODO : 여기에 image update 시키는 것도 넣어야 겠네
        isosurface.delete(0, tk.END)
        isosurface.insert(0, value)

    def iso_cmd2():
        # TODO : 여기에 image update 시키는 것도 넣어야 겠네
        iso_scale.set(isosurface.get())

    def iso_cmd3(event):
        # TODO : 여기에 image update 시키는 것도 넣어야 겠네
        iso_cmd2()
    iso_label = ttk.LabelFrame(opt_frame, text="Isosurface")
    iso_label.pack(side="top", fill="x", expand=True)
    iso_var = tk.DoubleVar()
    iso_scale = tk.Scale(iso_label, variable=iso_var, orient=tk.HORIZONTAL,
                         from_=iso_min, to=iso_max, resolution=(iso_max-iso_min)/10, command = iso_cmd,
                         showvalue=False)
    iso_scale.pack(side="bottom", fill="x", expand=True)
    isosurface = ttk.Spinbox(iso_label, from_=iso_min, to=iso_max, increment=(iso_max-iso_min)/10,
                             command = iso_cmd2)
    isosurface.bind("<Return>", iso_cmd3)
    isosurface.delete(0, tk.END)
    isosurface.insert(0, iso_suitable)
    iso_cmd2()
    isosurface.pack(fill="x", expand=True)


    def gaus_cmd(value):
        # TODO : 여기에 image update 시키는 것도 넣어야 겠네
        gaus_sigma.delete(0, tk.END)
        gaus_sigma.insert(0, value)

    def gaus_cmd2():
        # TODO : 여기에 image update 시키는 것도 넣어야 겠네
        gaus_scale.set(gaus_sigma.get())

    def gaus_cmd3(event):
        # TODO : 여기에 image update 시키는 것도 넣어야 겠네
        gaus_cmd2()
    gaus_label = ttk.LabelFrame(opt_frame, text="Gaussian Blurring sigma")
    gaus_label.pack(side="top", fill="x", expand=True)
    gaus_var = tk.DoubleVar()
    gaus_scale = tk.Scale(gaus_label, variable=gaus_var, orient=tk.HORIZONTAL,
                          from_=0, to=1000, resolution=10, command = gaus_cmd,
                          showvalue=False)
    gaus_scale.pack(side="bottom", fill="x", expand=True)
    gaus_sigma = ttk.Spinbox(gaus_label, from_=0, to=1000, increment=5,
                             command = gaus_cmd2)
    gaus_sigma.bind("<Return>", gaus_cmd3)
    gaus_sigma.delete(0, tk.END)
    gaus_sigma.insert(0, 10)
    gaus_cmd2()
    gaus_sigma.pack(fill="x", expand=True)


    def nlayer_cmd(value):
        # TODO : 여기에 image update 시키는 것도 넣어야 겠네
        nlayer.delete(0, tk.END)
        nlayer.insert(0, value)

    def nlayer_cmd2():
        # TODO : 여기에 image update 시키는 것도 넣어야 겠네
        nlayer_scale.set(nlayer.get())

    def nlayer_cmd3(event):
        # TODO : 여기에 image update 시키는 것도 넣어야 겠네
        nlayer_cmd2()
    nlayer_label = ttk.LabelFrame(opt_frame, text="Number of outmost layers")
    nlayer_label.pack(side="top", fill="x", expand=True)
    nlayer_var = tk.DoubleVar()
    nlayer_scale = tk.Scale(nlayer_label, variable=nlayer_var, orient=tk.HORIZONTAL,
                            from_=1, to=max_nlayer, resolution=1, command = nlayer_cmd,
                            showvalue=False)
    nlayer_scale.pack(side="bottom", fill="x", expand=True)
    nlayer = ttk.Spinbox(nlayer_label, from_=1, to=max_nlayer, increment=1,
                         command = nlayer_cmd2)
    nlayer.bind("<Return>", nlayer_cmd3)
    nlayer.delete(0, tk.END)
    nlayer.insert(0, 1)
    nlayer_cmd2()
    nlayer.pack(fill="x", expand=True)

    def ratom_cmd(value):
        # TODO : 여기에 image update 시키는 것도 넣어야 겠네
        ratom.delete(0, tk.END)
        ratom.insert(0, value)

    def ratom_cmd2():
        # TODO : 여기에 image update 시키는 것도 넣어야 겠네
        ratom_scale.set(ratom.get())

    def ratom_cmd3(event):
        # TODO : 여기에 image update 시키는 것도 넣어야 겠네
        ratom_cmd2()
    ratom_label = ttk.LabelFrame(opt_frame, text="Atom radius scale")
    ratom_label.pack(side="top", fill="x", expand=True)
    ratom_var = tk.DoubleVar()
    ratom_scale = tk.Scale(ratom_label, variable=ratom_var, orient=tk.HORIZONTAL,
                           from_=0, to=500, resolution=10, command = ratom_cmd,
                           showvalue=False)
    ratom_scale.pack(side="bottom", fill="x", expand=True)
    ratom = ttk.Spinbox(ratom_label, from_=0, to=500, increment=10,
                        command = ratom_cmd2)
    ratom.bind("<Return>", ratom_cmd3)
    ratom.delete(0, tk.END)
    ratom.insert(0, 60)
    ratom_cmd2()
    ratom.pack(fill="x", expand=True)


    # Menu
    # menu = tk.Menu(main_window)

    # ## Menu1: File
    # menu_file = tk.Menu(menu, tearoff=0)
    # menu_file.add_command(label="Open Current File", command=open_cur_file)
    # menu_file.add_command(label="Open Structure File", command=open_struc_file)
    # menu_file.add_separator()
    # menu_file.add_command(label="Save Image as")
    # menu_file.add_separator()
    # menu_file.add_command(label="Export bskan.in", command=export_bskanin)
    # menu_file.add_separator()
    # menu_file.add_command(label="Exit", command=main_window.quit)

    # menu.add_cascade(label="File", menu=menu_file)

    # main_window.config(menu = menu)
    main_window.mainloop()

if __name__ == "__main__":
	main()