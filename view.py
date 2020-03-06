# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 15:41:16 2020

@author: Mark
"""
import tkinter as tk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backend_bases import MouseEvent, LocationEvent
from tkinter import filedialog
from tkinter.ttk import Combobox, Scrollbar
from tkinter import ttk
import functools
from matplotlib.widgets import RectangleSelector

class View:
    def __init__(self, controller, root):
        self.controller = controller
        self.container = root
        self.container.title('Scatter Simulator')
        self.selected_spectrum = ''
        self.setup()

    def setup(self):
        
        self.createWidgets()
        self.setupLayout()
        
    def createWidgets(self):
        self.bcolor = 'grey'
        self.bthickness = 0
        self.axis_label_fontsize = 12
        self.axis_tick_font = 8
        
        # controls frame
        self.left_frame = tk.Frame(self.container, borderwidth=2, width=300, height=600, highlightbackground=self.bcolor, highlightcolor=self.bcolor, highlightthickness=self.bthickness)

        # load spectra frame
        self.load_spectra_frame = tk.Frame(self.left_frame, borderwidth=10, width=300, height=500, highlightbackground=self.bcolor, highlightcolor=self.bcolor, highlightthickness=self.bthickness)
        self.load_spectra_frame = tk.Frame(self.left_frame, borderwidth=2,width=400,height=200, highlightbackground=self.bcolor, highlightcolor=self.bcolor, highlightthickness=self.bthickness)
        self.step1_label = tk.Label(self.load_spectra_frame, text='1. Load XPS spectra',font=("Helvetica", 12))
        # Load spectrum buttons
        self.btn1 = tk.Button(self.load_spectra_frame, text = "Load spectrum", width = 15, command = self.loadSpectrum,borderwidth=2)
        self.btn2 = tk.Button(self.load_spectra_frame, text = "Build spectrum",width = 15, borderwidth=2, command = self.addSynthSpec)
        
        # Loss buttons frame
        self.loss_buttons_frame = tk.Frame(self.left_frame, borderwidth=2,width=400,height=600, highlightbackground=self.bcolor, highlightcolor=self.bcolor, highlightthickness=self.bthickness)
        self.step2_label = tk.Label(self.loss_buttons_frame, text='2. Choose a loss function',font=("Helvetica", 12))

        # Experimental parameter inputs
        self.exp_param_frame = tk.Frame(self.left_frame, borderwidth=10, width=300, height=500, highlightbackground=self.bcolor, highlightcolor=self.bcolor, highlightthickness=self.bthickness)
        self.exp_param_subframe = tk.Frame(self.exp_param_frame, borderwidth=10, width=300, height=500, highlightbackground=self.bcolor, highlightcolor=self.bcolor, highlightthickness=self.bthickness)
        self.step3_label = tk.Label(self.exp_param_frame, text='3. Set parameters',font=("Helvetica", 12))
        self.pressure = tk.DoubleVar()
        self.pressure.set(1.0)
        self.pressure_label = tk.Label(self.exp_param_subframe, text="Pressure [mbar]",borderwidth=2)
        self.pressure_entry = tk.Entry(self.exp_param_subframe, width=15,borderwidth=2, textvariable = self.pressure)
        self.distance = tk.DoubleVar()
        self.distance.set(0.8)
        self.distance_label = tk.Label(self.exp_param_subframe, text="Distance [mm]",borderwidth=2)
        self.distance_entry = tk.Entry(self.exp_param_subframe, width=15,borderwidth=2, textvariable = self.distance)
        
        self.n_iter = tk.IntVar()
        self.n_iter_label = tk.Label(self.exp_param_subframe, text="Nr. iterations",borderwidth=2)
        self.n_iter_entry = tk.Entry(self.exp_param_subframe, width=15,borderwidth=2, textvariable = self.n_iter)
        
        self.prob = tk.DoubleVar()
        self.prob_label = tk.Label(self.exp_param_subframe, text="Probability",borderwidth=2)
        self.prob_entry = tk.Entry(self.exp_param_subframe, width=15,borderwidth=2, textvariable = self.prob)

        self.advanced = tk.IntVar()
        self.advanced_chk = tk.Checkbutton(self.exp_param_frame, text="Advanced", variable=self.advanced, command = self.showAdvanced)
        
        # Execute simulation
        self.simulate_frame = tk.Frame(self.left_frame,borderwidth=10, width=300, height=500, highlightbackground=self.bcolor, highlightcolor=self.bcolor, highlightthickness=self.bthickness)
        self.step4_label = tk.Label(self.simulate_frame, text='4. Run simulation',font=("Helvetica", 12))
        self.simulate_label = tk.Label(self.simulate_frame, text="Simulation")
        self.scatter_btn = tk.Button(self.simulate_frame, text = "Scatter", command = self.controller.scatterSpectrum, borderwidth=2, width=15, pady=2)
        self.unscatter_btn = tk.Button(self.simulate_frame, text = "Un-scatter", borderwidth=2, width=15, pady=2)    
        self.bulk = tk.IntVar()
        self.bulk_chk = tk.Checkbutton(self.simulate_frame, text="Bulk spectrum", variable=self.bulk)
        
        # XPS spectra Frame
        self.mid_frame = tk.Frame(self.container, borderwidth=2,width=400,height=600, highlightbackground=self.bcolor, highlightcolor=self.bcolor, highlightthickness=self.bthickness)

        # Figure for XPS spectra
        x = []
        y = []
        self.fig1 = Figure(x,y)
        self.fig1.ax.tick_params(direction='out', length=4, width=1, colors='black',
               grid_color='black', labelsize=self.axis_tick_font, grid_alpha=0.5)
        self.fig1.ax.set_xlabel('Energy [eV]', fontsize=self.axis_label_fontsize)
        self.fig1.ax.set_ylabel('Intensity [cts./sec.]', fontsize=self.axis_label_fontsize)
        self.fig1.ax.set_title('Spectra')
        self.fig1.fig.tight_layout()
        self.chart1 = FigureCanvasTkAgg(self.fig1.fig, self.mid_frame)
        self.chart1.mpl_connect('button_press_event', functools.partial(self.controller.doubleClkChart, ax=self.fig1.ax, chart = self.chart1))
        rectprops = dict(facecolor='white', edgecolor = 'black',
                           alpha=0.2, fill=True)
        self.RS1 = RectangleSelector(self.fig1.ax, functools.partial(self.controller.selector, ax=self.fig1.ax, chart = self.chart1),
                                       drawtype='box', useblit=True,
                                       button=[1, 3],  # don't use middle button
                                       minspanx=5, minspany=5,
                                       spancoords='pixels',
                                       interactive=True,
                                       rectprops = rectprops) 
        self.RS1.set_active(True)
        # XPS spectra table
        self.bottom_middle_frame = tk.Frame(self.mid_frame, borderwidth=2,width=100,height=1, highlightbackground=self.bcolor, highlightcolor=self.bcolor, highlightthickness=self.bthickness)
        columns = ('Nr.','Type', 'Visibility')
        self.spectra_table = ttk.Treeview(self.mid_frame, height=4,show='headings',columns=columns, selectmode='browse')
        self.spectra_table.name = 'spectra'
        self.spectra_table.column('Nr.',width=50,anchor=tk.W)
        self.spectra_table.heading('Nr.', text='Nr.', anchor=tk.W)
        self.spectra_table.column('Type',width=150,anchor=tk.W)
        self.spectra_table.heading('Type', text='Type', anchor=tk.W)        
        self.spectra_table.column('Visibility',width=200,anchor=tk.W)
        self.spectra_table.heading('Visibility', text='Visibility', anchor=tk.W)
        self.spectra_table.entry_choices = {'0':['none','Peak','VacuumExcitation'],'1':['visible','hidden']}
        self.spectra_table.bind('<Button-3>',functools.partial(self.controller.tablePopup, table=self.spectra_table, table_choices = self.controller.spectra_table_choices))
        self.spectra_table.bind('<Double-1>',functools.partial(self.controller.doubleClkTable, table=self.spectra_table))
        self.spectra_table.bind('<Delete>', self.removeSpectrum)
        
        # normalize spectra box
        self.normalize = tk.IntVar()
        self.normalize_chk = tk.Checkbutton(self.mid_frame, text="Normalize", variable=self.normalize, command = self.controller.rePlotFig1)

        # Loss function Frame
        self.right_frame = tk.Frame(self.container,borderwidth=2,width=400,height=600, highlightbackground=self.bcolor, highlightcolor=self.bcolor, highlightthickness=self.bthickness)

        # Loss function figure frame
        self.bottom_right_frame = tk.Frame(self.right_frame, borderwidth=2,width=400,height=100, highlightbackground=self.bcolor, highlightcolor=self.bcolor, highlightthickness=self.bthickness)

        # Load file of loss functions
        self.btn3 = tk.Button(self.loss_buttons_frame, text="Load loss functions",
                              borderwidth=2, width=15, command=self.loadScatterers)

        # Select loss function (from loaded file)
        self.load_loss_frame = tk.Frame(self.loss_buttons_frame, borderwidth=2,width=400,height=600, highlightbackground=self.bcolor, highlightcolor=self.bcolor, highlightthickness=self.bthickness)
        self.load_loss_label = tk.Label(self.load_loss_frame, text='Select loss function')
        self.cbox = Combobox(self.load_loss_frame, width=15, textvariable = self.controller.selected_scatterer,
                             values=self.controller.scatterer_choices)
        #self.cbox['values'] = self.controller.scatterer_choices
        self.cbox.bind("<<ComboboxSelected>>", self.controller.setCurrentScatterer)

        # Build loss function
        self.btn4 = tk.Button(self.loss_buttons_frame, text = "New loss function", borderwidth=2, width=15)
        # Save loss function
        self.btn5 = tk.Button(self.loss_buttons_frame, text = "Save loss functions", borderwidth=2, width=15, command = self.saveScatterers)

        # cross sections frame
        self.cross_section_frame = tk.Frame(self.bottom_right_frame)
        self.cross_section_label = tk.Label(self.cross_section_frame, text='Inelastic probaility:')
        self.cross_section = tk.StringVar()
        self.cross_section_entry = tk.Entry(self.cross_section_frame, width = 10,borderwidth = 2, textvariable = self.cross_section)
        self.cross_section.trace('w',self.controller.updateCrossSection)
        self.angle_factor = tk.StringVar()
        self.angle_factor_label = tk.Label(self.cross_section_frame, text='Angle Factor:')
        self.angle_factor_entry = tk.Entry(self.cross_section_frame, width = 10,borderwidth = 2, textvariable = self.angle_factor)
        self.angle_factor.trace('w',self.controller.updateAngle)
        self.gas_diameter_label = tk.Label(self.cross_section_frame, text='Gas diameter (nm):')
        self.gas_diameter = tk.StringVar()
        self.gas_diameter_entry = tk.Entry(self.cross_section_frame,width = 10,borderwidth = 2, textvariable = self.gas_diameter)
        self.gas_diameter.trace('w',self.controller.updateDiameter)
        
        # Loss function figure (Figure2)
        x = []
        y = []
        self.fig2 = Figure(x,y)
        self.fig2.ax.tick_params(direction='out', length=4, width=1, colors='black',
               grid_color='black', labelsize=self.axis_tick_font, grid_alpha=0.5)
        self.fig2.ax.set_xlabel('Energy Loss [eV]', fontsize=self.axis_label_fontsize)
        self.fig2.ax.set_ylabel('Probability', fontsize=self.axis_label_fontsize)
        self.fig2.ax.set_title('Loss Function')
        self.fig2.fig.tight_layout()
        self.chart2 = FigureCanvasTkAgg(self.fig2.fig, self.right_frame)
        self.chart2.mpl_connect('button_press_event', functools.partial(self.controller.doubleClkChart, ax=self.fig2.ax, chart = self.chart2))
        rectprops = dict(facecolor='white', edgecolor = 'black',
                           alpha=0.2, fill=True)
        self.RS2 = RectangleSelector(self.fig2.ax, functools.partial(self.controller.selector, ax=self.fig2.ax, chart = self.chart2),
                                       drawtype='box', useblit=True,
                                       button=[1, 3],  # don't use middle button
                                       minspanx=5, minspany=5,
                                       spancoords='pixels',
                                       interactive=True,
                                       rectprops = rectprops) 
        self.RS2.set_active(True)

        # Scatterers Table
        columns = ('Nr.','Type')
        self.scatterers_table = ttk.Treeview(self.bottom_right_frame, height=4,show='headings',columns=columns, selectmode='browse')
        self.scatterers_table.column('Nr.',width=50,anchor=tk.W)
        self.scatterers_table.heading('Nr.', text='Nr.', anchor=tk.W)
        self.scatterers_table.column('Type',width=350,anchor=tk.W)
        self.scatterers_table.heading('Type', text='Type', anchor=tk.W)
        self.scatterers_table.bind('<Double-1>',self.controller.callLossEditor)
        self.scatterers_table.bind('<Delete>', self.removeComponent)
        self.add_comp_btn = tk.Label(self.bottom_right_frame, text = "+ Add component", borderwidth=2, padx = 2, pady=0, relief ='raised')
        self.add_comp_btn.bind('<Button-1>', self.addComponent)
        
    def setupLayout(self):
        self.left_frame.pack(side=tk.LEFT, fill=None, anchor="nw")
        self.load_spectra_frame.pack(side=tk.TOP, expand=False, fill=tk.Y, anchor='center', padx=15, pady=5)
        self.step1_label.pack(side=tk.TOP, pady=5)
        self.loss_buttons_frame.pack(side=tk.TOP, expand=False, fill=tk.Y, anchor='center', padx=15, pady=5)
        self.step2_label.pack(side=tk.TOP, pady=5)
        
        self.exp_param_frame.pack(side=tk.TOP, fill = None, anchor='center',pady=5)
        self.step3_label.pack(side=tk.TOP, pady=5)
        self.exp_param_subframe.pack(side=tk.TOP)
        self.distance_label.pack(side=tk.TOP)
        self.distance_entry.pack(side=tk.TOP)
        self.pressure_label.pack(side=tk.TOP)
        self.pressure_entry.pack(side=tk.TOP)
        
        self.simulate_frame.pack(side=tk.TOP, fill = None, anchor='center', pady=5)
        self.step4_label.pack(side=tk.TOP, pady=5)
        self.advanced_chk.pack(side=tk.TOP)
        self.simulate_label.pack(side=tk.TOP)
        self.scatter_btn.pack(side=tk.TOP)
        self.unscatter_btn.pack(side=tk.TOP)
        self.bulk_chk.pack(side=tk.TOP)
        
        self.mid_frame.pack(side=tk.LEFT, fill=tk.Y, anchor='n')
        self.btn1.pack(side=tk.TOP, fill = None, pady=2)
        self.btn2.pack(side=tk.TOP, fill = None, pady=2)
        
        # Figure XPS
        self.chart1.get_tk_widget().pack(side=tk.TOP)
        self.bottom_middle_frame.pack(side=tk.TOP, pady=10) # This is just for some extra space
        self.spectra_table.pack(side=tk.TOP, anchor="ne", pady=10, padx=15)
        self.normalize_chk.pack(side=tk.TOP, anchor="ne", pady=10, padx=15)
        
        self.right_frame.pack(side=tk.LEFT, fill = None, anchor='n')
        self.btn3.pack(side=tk.TOP)
        self.load_loss_frame.pack(side=tk.TOP)
        self.load_loss_label.pack(side=tk.TOP)
        self.cbox.pack(side=tk.TOP)
        self.btn4.pack(side=tk.TOP)
        self.btn5.pack(side=tk.TOP)
        self.chart2.get_tk_widget().pack(side=tk.TOP)
        
        self.bottom_right_frame.pack(side=tk.TOP, anchor="ne")
        self.cross_section_frame.pack(side=tk.TOP, padx=15)
        self.cross_section_label.pack(side=tk.LEFT)
        self.cross_section_entry.pack(side=tk.LEFT)
        self.angle_factor_label.pack(side=tk.LEFT)
        self.angle_factor_entry.pack(side=tk.LEFT)
        self.gas_diameter_label.pack(side=tk.LEFT)
        self.gas_diameter_entry.pack(side=tk.LEFT)
        self.scatterers_table.pack(side=tk.TOP, anchor="ne", pady=10, padx=15)
        self.add_comp_btn.pack(side=tk.TOP, anchor="ne", padx=15)
        
    def showAdvanced(self):
        if self.advanced.get() == 1:
            #self.exp_param_subframe.pack_forget()
            self.distance_label.pack_forget()
            self.distance_entry.pack_forget()
            self.pressure_label.pack_forget()
            self.pressure_entry.pack_forget()
            
            self.n_iter_label.pack(side=tk.TOP)
            self.n_iter_entry.pack(side=tk.TOP)
            self.prob_label.pack(side=tk.TOP)
            self.prob_entry.pack(side=tk.TOP)
            
            self.cross_section_entry.configure(state='disable')
            self.angle_factor_entry.configure(state='disable')
            self.gas_diameter_entry.configure(state='disable')
            
        elif self.advanced.get() == 0:
            self.n_iter_label.pack_forget()
            self.n_iter_entry.pack_forget()
            self.prob_label.pack_forget()
            self.prob_entry.pack_forget()
            
            self.distance_label.pack(side=tk.TOP)
            self.distance_entry.pack(side=tk.TOP)
            self.pressure_label.pack(side=tk.TOP)
            self.pressure_entry.pack(side=tk.TOP)
            
            self.cross_section_entry.configure(state='normal')
            self.angle_factor_entry.configure(state='normal')
            self.gas_diameter_entry.configure(state='normal')
        
    def loadSpectrum(self):
        file = filedialog.askopenfilename(initialdir = self.controller.datapath)
        if file:
            self.controller.loadSpectrum(file)
            self.controller.datapath = (file.rsplit('/',maxsplit = 1)[0])
    
    def loadScattered(self):
        file = filedialog.askopenfilename(initialdir = self.controller.datapath)
        self.controller.loadSpectrumToFit(file)

    def loadScatterers(self):
        file = filedialog.askopenfilename(initialdir=self.controller.datapath,
                                          title='Select loss function file',
                                          filetypes=[('json', '*.json')])
        self.controller.loadScatterers(file)

    def saveScatterers(self):
        file = filedialog.asksaveasfilename(initialdir=self.controller.datapath,
                                            title='Save as',
                                            filetypes=[('json', '*.json')])
        self.controller.saveScatterers(file)

    def updateScattererChoice(self):
        self.cbox.set('')
        self.cbox.config(values=self.controller.scatterer_choices)

    def noScatterer(self):
        tk.messagebox.showerror('Error','Please select an Unscattered spectrum and a Loss Function')        
        
    def createSpectrum(self):
        self.controller.createSynthetic()
        
    def removeSpectrum(self, event):
        selected_item = self.spectra_table.selection()
        cur_item = self.spectra_table.item(selected_item[0])['values'][0]
        self.spectra_table.delete(selected_item[0])
        self.controller.removeSpectrum(cur_item)
        
    def removeComponent(self, event):
        selected_item = self.scatterers_table.selection()
        cur_item = self.scatterers_table.item(selected_item[0])['values'][0]
        self.scatterers_table.delete(selected_item[0])
        self.controller.removeComponent(cur_item)
        
    def addComponent(self, event):
        def setChoice(choice):
            self.controller.addComponent(choice)
        popup = tk.Menu(self.container, tearoff=0)
        choices = self.controller.component_choices
        for i,j in enumerate(choices):
                # This line below is a bit tricky. Needs to be done this way because the i in the loop is only scoped for the loop, and does not persist
                popup.add_command(command = lambda choice = choices[i]: setChoice(choice), label=j)
        popup.post(event.x_root, event.y_root)
        
    def addSynthSpec(self):
        self.controller.addSynthSpec()
        self.controller.fillTable1()
        self.selected_spectrum = len(self.spectra_table.get_children())-1
        self.spec_builder = SpecBuilder(self.controller,self.selected_spectrum)
                   
class Figure:

    def __init__(self, x, y, dpi = 100):
        self.x = x
        self.y = y
        self.xlabel = ''
        self.ylabel = ''
        self.title = ''
        self.fig, self.ax = plt.subplots(figsize=(5,4), dpi=100)
        self.fig.patch.set_facecolor('0.9411')
        self.ax.plot(x, y)
        self.ax.set_xlabel(self.xlabel)
        self.ax.set_ylabel(self.ylabel)
        self.ax.set_title(self.title)
        self.fig.tight_layout()
        
class LossEditor:
    
    def __init__(self, controller, params, comp_nr):
        self.controller = controller
        self.bcolor = 'grey'
        self.bthickness = 0
        self.params = params
        self.comp_nr = comp_nr
        window = tk.Toplevel()
        title = 'Component: ' + str(comp_nr)
        window.wm_title(title)
        window.attributes("-topmost", True)
        header = tk.Label(window, text = 'Component: ' + str(comp_nr))
        header.pack(side=tk.TOP)
        self.labels = []
        self.entries = []
        self.stringvars = []
        i=0
        for key in params:
            subsubframe = tk.Frame(window, highlightbackground=self.bcolor, highlightcolor=self.bcolor, highlightthickness=self.bthickness)
            subsubframe.pack(side=tk.TOP, padx=10, pady=10)
            self.stringvars += [tk.StringVar()]
            self.stringvars[i].set(params[key])
            self.stringvars[i].trace('w', self.callBack)
            self.labels += [tk.Label(subsubframe, text = key, borderwidth = 2)]
            self.labels[i].pack(side=tk.TOP)
            self.entries += [tk.Entry(subsubframe, width = 20,borderwidth = 2, textvariable = self.stringvars[i])]
            self.entries[i].pack(side=tk.TOP)
            i+=1

    def callBack(self,event, *args):
        if self.buildDict() == 1:
            self.controller.changeLossFunction()
            
    def buildDict(self):
        ready = 0
        for i, key in enumerate(self.params):
            if len(self.stringvars[i].get()) == 0:
                break
            else:
                self.params[key] = float(self.stringvars[i].get())
                ready = 1
        return ready
    
class SpecBuilder:
    def __init__(self, controller, spec_idx):
        self.controller = controller
        self.spec_idx = spec_idx
        self.selected_peak = None
        self.bcolor = 'grey'
        self.bthickness = 0
        self.window = tk.Toplevel()
        title = 'Component: '
        self.window.wm_title(title)
        self.window.attributes("-topmost", True)
        header = tk.Label(self.window, text = 'Component: ')
        header.pack(side=tk.TOP)
        self.labels = []
        self.entries = []
        self.stringvars = []
        self.createWidgets()
        self.setupLayout()
        
    def createWidgets(self):
        self.done = tk.Button(self.window, text='Done', command=self.Done)
        columns = ('Nr.','Type', 'Position', 'Width', 'Intensity')
        self.peak_table = ttk.Treeview(self.window, height=8,show='headings',columns=columns, selectmode='browse')
        self.peak_table.name = 'spectra'
        self.peak_table.column('Nr.',width=50,anchor=tk.W)
        self.peak_table.heading('Nr.', text='Nr.', anchor=tk.W)
        self.peak_table.column('Type',width=100,anchor=tk.W)
        self.peak_table.heading('Type', text='Type', anchor=tk.W)        
        self.peak_table.column('Position',width=50,anchor=tk.W)
        self.peak_table.heading('Position', text='Position', anchor=tk.W)
        self.peak_table.column('Width',width=50,anchor=tk.W)
        self.peak_table.heading('Width', text='Width', anchor=tk.W)
        self.peak_table.column('Intensity',width=60,anchor=tk.W)
        self.peak_table.heading('Intensity', text='Intensity', anchor=tk.W)
        self.peak_table.bind('<Delete>', self.removeComponent)
        self.add_comp_btn = tk.Label(self.window, text = "+ Add component", borderwidth=2, padx = 2, pady=0, relief ='raised')
        self.add_comp_btn.bind('<Button-1>', self.addComponent)
        self.peak_table.bind('<ButtonRelease-1>',self.selectComponent)
        
        self.entries_frame = tk.Frame(self.window)
        self.position_frame = tk.Frame(self.entries_frame)
        self.position = tk.DoubleVar()
        self.position.trace('w', self.modPeak)
        self.position_label = tk.Label(self.position_frame, text = 'Position')
        self.position_entry = tk.Entry(self.position_frame, width = 10, textvariable = self.position)
        #self.position_entry.config(validate='key', validatecommand = self.modPeak)
        
        self.width_frame = tk.Frame(self.entries_frame)
        self.width = tk.DoubleVar()
        self.width.trace('w', self.modPeak)
        self.width_label = tk.Label(self.width_frame, text = 'Width')
        self.width_entry = tk.Entry(self.width_frame, width = 10, textvariable = self.width)
        #self.width_entry.config(validate='key', validatecommand = self.modPeak)
        
        self.intensity_frame = tk.Frame(self.entries_frame)
        self.intensity = tk.DoubleVar()
        self.intensity.trace('w', self.modPeak)
        self.intensity_label = tk.Label(self.intensity_frame, text = 'Intensity')
        self.intensity_entry = tk.Entry(self.intensity_frame, width = 10, textvariable = self.intensity)
        #self.intensity_entry.config(validate='key', validatecommand = self.modPeak)
        
    def setupLayout(self):
        padding = 5
        self.entries_frame.pack(side=tk.TOP, padx=padding, pady=padding)
        self.position_frame.pack(side=tk.LEFT, padx=padding, pady=padding)
        self.width_frame.pack(side=tk.LEFT, padx=padding, pady=padding)
        self.intensity_frame.pack(side=tk.LEFT, padx=padding, pady=padding)
        
        self.position_label.pack(side=tk.TOP, padx=padding, pady=padding)
        self.position_entry.pack(side=tk.TOP, padx=padding, pady=padding)
        
        self.width_label.pack(side=tk.TOP, padx=padding, pady=padding)
        self.width_entry.pack(side=tk.TOP, padx=padding, pady=padding)

        self.intensity_label.pack(side=tk.TOP, padx=padding, pady=padding)
        self.intensity_entry.pack(side=tk.TOP, padx=padding, pady=padding)

        self.peak_table.pack(side=tk.TOP, padx=15, pady=0)
        self.add_comp_btn.pack(side=tk.TOP,anchor='ne', padx=15)
        self.done.pack(side=tk.TOP, pady=15)
        
    def Done(self):
        self.window.destroy()
              
    def removeComponent(self):
        return
    
    def addComponent(self, event):
        def setChoice(choice):
            self.controller.addPeak(self.spec_idx, choice)
            self.refreshTable()
            self.controller.rePlotFig1()
        popup = tk.Menu(self.window, tearoff=0)
        choices = self.controller.peak_choices
        for i,j in enumerate(choices):
                # This line below is a bit tricky. Needs to be done this way because the i in the loop is only scoped for the loop, and does not persist
                popup.add_command(command = lambda choice = choices[i]: setChoice(choice), label=j)
        popup.post(event.x_root, event.y_root)
        
    def refreshTable(self):
        for row in self.peak_table.get_children():
            self.peak_table.delete(row)
        values = self.controller.getComps(self.spec_idx)
        for value in values:
            self.peak_table.insert('',value,values=values[value], iid=str(value))
            
    def selectComponent(self, event):
        sel = self.peak_table.selection()
        self.selected_peak = sel[0]
        cur_item = self.peak_table.item(sel[0])['values']
        #peak_type = cur_item[1]
        pos = cur_item[2]
        width = cur_item[3]
        intensity = cur_item[4]
        self.position.set(pos)
        self.width.set(width)
        self.intensity.set(intensity)

    def modPeak(self, *args):
        if self.selected_peak == '':
            return
        else:
            position = self.position.get()
            width = self.width.get()
            intensity = self.intensity.get()
            new_values = {'spec_idx':self.spec_idx,'peak_idx':self.selected_peak,'position': position, 'width':width, 'intensity': intensity}
            self.controller.updatePeak(new_values)
            self.updateEntry(new_values)
            
    def updateEntry(self, new_values):
        peak_idx = new_values['peak_idx']
        peak_type = self.peak_table.item(self.selected_peak)['values'][1]
        values = (new_values['peak_idx'],peak_type, new_values['position'],new_values['width'],new_values['intensity'])
        self.peak_table.item(peak_idx, values = values)
            

if __name__ == "__main__":
    mainwin = tk.Tk()
    app = View(mainwin)
    mainwin.mainloop()
