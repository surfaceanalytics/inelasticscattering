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
        self.default_start = 0
        self.default_stop = 100
        self.default_step = 0.1

    def setup(self):
        
        self.createWidgets()
        self.setupLayout()
        
    def createWidgets(self):
        self.bcolor = 'grey'
        self.bthickness = 0
        self.axis_label_fontsize = 12
        self.axis_tick_font = 8
        # frames are leveld in a heirarchy
        # a frame embedded into another frame will have a name like f1_1_2
        
        # f1 : left highest-level frame
        # controls frame
        # contains the user controls for loading, saving files, for changing parameters of simulation, and for running simulation
        self.f1 = tk.Frame(self.container, borderwidth=2, width=300, height=600, highlightbackground=self.bcolor, highlightcolor=self.bcolor, highlightthickness=self.bthickness)

        # f2 : middle highest-level frame
        # XPS spectra Frame
        # contains the plot of the XPS spectra, and the list of spectra
        self.f2 = tk.Frame(self.container, borderwidth=2,width=400,height=600, highlightbackground=self.bcolor, highlightcolor=self.bcolor, highlightthickness=self.bthickness)

        # f3 : right highest-level frame
        # Loss function Frame
        # contains plot of loss function, and loss functions table
        self.f3 = tk.Frame(self.container,borderwidth=2,width=400,height=600, highlightbackground=self.bcolor, highlightcolor=self.bcolor, highlightthickness=self.bthickness)

        # f3_1
        # Loss function figure frame
        self.f3_1 = tk.Frame(self.f3, borderwidth=2,width=400,height=100, highlightbackground=self.bcolor, highlightcolor=self.bcolor, highlightthickness=self.bthickness)

        # f1_1 
        # load spectra frame
        self.f1_1 = tk.Frame(self.f1, borderwidth=10, width=300, height=500, highlightbackground=self.bcolor, highlightcolor=self.bcolor, highlightthickness=self.bthickness)
        self.f1_1 = tk.Frame(self.f1, borderwidth=2,width=400,height=200, highlightbackground=self.bcolor, highlightcolor=self.bcolor, highlightthickness=self.bthickness)
        self.step1_label = tk.Label(self.f1_1, text='1. Load XPS spectra',font=("Helvetica", 12))
        self.btn1 = tk.Button(self.f1_1, text = "Load spectrum", width = 15, command = self.loadSpectrum,borderwidth=2)
        self.btn2 = tk.Button(self.f1_1, text = "Build spectrum",width = 15, borderwidth=2, command = self.addSynthSpec)
        
        # f1_2
        # Loss buttons frame
        self.f1_2 = tk.Frame(self.f1, borderwidth=2,width=400,height=600, highlightbackground=self.bcolor, highlightcolor=self.bcolor, highlightthickness=self.bthickness)
        self.step2_label = tk.Label(self.f1_2, text='2. Choose a scatterer',font=("Helvetica", 12))

        # f1_3
        # Parameter inputs frame
        self.f1_3 = tk.Frame(self.f1, borderwidth=2, width=300, height=500, highlightbackground=self.bcolor, highlightcolor=self.bcolor, highlightthickness=self.bthickness)
        self.step3_label = tk.Label(self.f1_3, text='3. Set parameters',font=("Helvetica", 12))

        # f1_3_1
        # Parameters subframe
        self.f1_3_1 = tk.Frame(self.f1_3, borderwidth=2, width=300, height=500, highlightbackground=self.bcolor, highlightcolor=self.bcolor, highlightthickness=self.bthickness)
        
        # f1_3_1_1
        # parameters sub-sub-frame
        self.f1_3_1_1 = tk.Frame(self.f1_3_1)
        
        # f1_3_1_2
        self.f1_3_1_2 = tk.Frame(self.f1_3_1)
        
        # f1_3_1_3
        self.f1_3_1_3 = tk.Frame(self.f1_3_1)
        
        # f1_3_2
        # Variants subframe
        self.f1_3_2 = tk.Frame(self.f1_3, borderwidth=2, width=300, height=30, highlightbackground=self.bcolor, highlightcolor=self.bcolor, highlightthickness=self.bthickness)
        
        
        # f1_4
        # run simulation frame
        self.f1_4 = tk.Frame(self.f1,borderwidth=10, width=300, height=500, highlightbackground=self.bcolor, highlightcolor=self.bcolor, highlightthickness=self.bthickness)
        
        # f2_1
        # spectra table frame
        self.f2_1 = tk.Frame(self.f2)

        
        # parameter variant 0
        self.entry_width = 7
        self.bord_width = 2
        self.pressure = tk.DoubleVar()
        self.pressure.set(1.0)
        self.pressure_label = tk.Label(self.f1_3_1_1, text="P [mbar]", borderwidth = self.bord_width)
        self.pressure_entry = tk.Entry(self.f1_3_1_1, width = self.entry_width, borderwidth = self.bord_width, textvariable = self.pressure)
        
        self.distance = tk.DoubleVar()
        self.distance.set(0.8)
        self.distance_label = tk.Label(self.f1_3_1_1, text="D [mm]", borderwidth = self.bord_width)
        self.distance_entry = tk.Entry(self.f1_3_1_1, width = self.entry_width, borderwidth = self.bord_width, textvariable = self.distance)
        
        self.inelastic_xsect_label = tk.Label(self.f1_3_1_2, text='Inelastic \n X-sect.:')
        self.inelastic_xsect = tk.StringVar()
        self.inelastic_xsect_entry = tk.Entry(self.f1_3_1_2, width = self.entry_width, borderwidth = self.bord_width, textvariable = self.inelastic_xsect)
        self.inelastic_xsect.trace('w',self.controller.updateInelasticXSect)
        
        self.inel_angle_factor = tk.StringVar()
        self.inel_angle_factor_label = tk.Label(self.f1_3_1_2, text='f(Angle):')
        self.inel_angle_factor_entry = tk.Entry(self.f1_3_1_2, width = self.entry_width, borderwidth = self.bord_width, textvariable = self.inel_angle_factor)
        self.inel_angle_factor.trace('w',self.controller.updateAngle)
        
        self.elastic_xsect_label = tk.Label(self.f1_3_1_3, text='Elastic \n X-sect.:')
        self.elastic_xsect = tk.StringVar()
        self.elastic_xsect_entry = tk.Entry(self.f1_3_1_3,width = self.entry_width,borderwidth = self.bord_width, textvariable = self.elastic_xsect)
        self.elastic_xsect.trace('w',self.controller.updateElasticXSect)

        self.el_angle_factor = tk.StringVar()
        self.el_angle_factor_label = tk.Label(self.f1_3_1_3, text = 'f(Angle):')
        self.el_angle_factor_entry = tk.Entry(self.f1_3_1_3, width = self.entry_width,borderwidth=self.bord_width, textvariable = self.el_angle_factor)
        
        # parameter variant 1
        self.n_iter = tk.IntVar()
        self.n_iter_label = tk.Label(self.f1_3_1_1, text="Nr. iter.",borderwidth=2)
        self.n_iter_entry = tk.Entry(self.f1_3_1_1, width = self.entry_width, borderwidth = self.bord_width, textvariable = self.n_iter)
  
        self.inel_prob = tk.DoubleVar()
        self.inel_prob_label = tk.Label(self.f1_3_1_2, text="Inelastic \n Prob.",borderwidth=2)
        self.inel_prob_entry = tk.Entry(self.f1_3_1_2, width = self.entry_width, borderwidth = self.bord_width, textvariable = self.inel_prob)
        
        self.el_prob = tk.DoubleVar()
        self.el_prob_label = tk.Label(self.f1_3_1_3, text="Elastic \n Prob.",borderwidth=2)
        self.el_prob_entry = tk.Entry(self.f1_3_1_3, width = self.entry_width, borderwidth = self.bord_width, textvariable = self.el_prob)

        self.variant = tk.IntVar()
        self.variant.set(0)
        # variant 2 has the same parameters as variant 0
        
        #self.variant_chk = tk.Checkbutton(self.f1_3, text="Advanced", variable=self.variant, command = self.toggleParamVariants)
        # dictionary to hold all of the parameters frame variants
        self.variants_dict = {
                            0:[
                            self.distance_label, self.distance_entry, 
                            self.pressure_label, self.pressure_entry,
                            self.inelastic_xsect_label, self.inelastic_xsect_entry,
                            self.inel_angle_factor_label, self.inel_angle_factor_entry,
                            self.elastic_xsect_label, self.elastic_xsect_entry,
                            self.el_angle_factor_label, self.el_angle_factor_entry
                            ],
                            1:[
                            self.n_iter_label,
                            self.n_iter_entry,
                            self.inel_prob_label,
                            self.inel_prob_entry,
                            self.inel_angle_factor_label,
                            self.inel_angle_factor_entry,
                            self.el_prob_label,
                            self.el_prob_entry,
                            self.el_angle_factor_label,
                            self.el_angle_factor_entry
                            ],
                            2:[
                            self.distance_label, self.distance_entry, 
                            self.pressure_label, self.pressure_entry,
                            self.inelastic_xsect_label, self.inelastic_xsect_entry,
                            self.inel_angle_factor_label, self.inel_angle_factor_entry,
                            self.elastic_xsect_label, self.elastic_xsect_entry,
                            self.el_angle_factor_label, self.el_angle_factor_entry
                            ]
                            }
        # make radio buttons for variants
        for k in self.variants_dict.keys():
            tk.Radiobutton(self.f1_3_2, text=str(k), variable=self.variant, indicatoron=0, command = self.toggleParamVariants, value=k).pack(side=tk.LEFT)
        
        # Run simulation
        self.step4_label = tk.Label(self.f1_4, text='4. Run simulation',font=("Helvetica", 12))
        self.simulate_label = tk.Label(self.f1_4, text="Simulation")
        self.scatter_btn = tk.Button(self.f1_4, text = "Scatter", command = self.controller.scatterSpectrum, borderwidth=2, width=15, pady=2)
        self.unscatter_btn = tk.Button(self.f1_4, text = "Un-scatter", borderwidth=2, width=15, pady=2)    
        self.bulk = tk.IntVar()
        self.bulk_chk = tk.Checkbutton(self.f1_4, text="Bulk spectrum", variable=self.bulk)
        
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
        self.chart1 = FigureCanvasTkAgg(self.fig1.fig, self.f2)
        self.chart1.mpl_connect('button_press_event', functools.partial(self.doubleClkChart, ax=self.fig1.ax, chart = self.chart1))
        rectprops = dict(facecolor='white', edgecolor = 'black',
                           alpha=0.2, fill=True)
        self.RS1 = RectangleSelector(self.fig1.ax, functools.partial(self.selector, ax=self.fig1.ax, chart = self.chart1),
                                       drawtype='box', useblit=True,
                                       button=[1, 3],  # don't use middle button
                                       minspanx=5, minspany=5,
                                       spancoords='pixels',
                                       interactive=True,
                                       rectprops = rectprops) 
        self.RS1.set_active(True)
        # XPS spectra table
        columns = ('Nr.','Type', 'Visibility')
        #self.spectra_table = ttk.Treeview(self.f2, height=4,show='headings',columns=columns, selectmode='browse')
        self.spectra_table = ttk.Treeview(self.f2_1, height=4, columns=columns, selectmode='browse')
        self.spectra_table.name = 'spectra'
        self.spectra_table.column('#0', width=60, anchor=tk.W)
        self.spectra_table.column('Nr.',width=50,anchor=tk.W)
        self.spectra_table.heading('Nr.', text='Nr.', anchor=tk.W)
        self.spectra_table.column('Type',width=145,anchor=tk.W)
        self.spectra_table.heading('Type', text='Type', anchor=tk.W)        
        self.spectra_table.column('Visibility',width=145,anchor=tk.W)
        self.spectra_table.heading('Visibility', text='Visibility', anchor=tk.W)

        self.spectra_table.entry_choices = {'0':['none','Peak','VacuumExcitation'],'1':['visible','hidden']}
        self.spectra_table.bind('<Button-3>',functools.partial(self.controller.tablePopup, table=self.spectra_table, table_choices = self.controller.spectra_table_choices))
        self.spectra_table.bind('<Double-1>',functools.partial(self.controller.doubleClkTable, table=self.spectra_table))
        self.spectra_table.bind('<Delete>', self.removeSpectrum)

        # normalize spectra box
        self.normalize = tk.IntVar()
        self.normalize_chk = tk.Checkbutton(self.f2_1, text="Normalize", variable=self.normalize, command = self.controller.rePlotFig1)

        # Load file of loss functions
        self.btn3 = tk.Button(self.f1_2, text="Load scatterer",
                              borderwidth=2, width=15, command=self.loadScatterers)

        # Select loss function (from loaded file)
        self.load_loss_frame = tk.Frame(self.f1_2, borderwidth=2,width=400,height=600, highlightbackground=self.bcolor, highlightcolor=self.bcolor, highlightthickness=self.bthickness)
        self.load_loss_label = tk.Label(self.load_loss_frame, text='Select loss function')
        self.selected_scatterer = tk.StringVar()
        self.scatterer_choices = []
        self.cbox = Combobox(self.load_loss_frame, width=15, textvariable = self.selected_scatterer,
                             values=self.scatterer_choices)
        #self.cbox['values'] = self.controller.scatterer_choices
        self.cbox.bind("<<ComboboxSelected>>", self.setCurrentScatterer)

        # Build loss function
        self.btn4 = tk.Button(self.f1_2, text = "New loss function", borderwidth=2, width=15, command = self.newScatterer)
        # Save loss function
        self.btn5 = tk.Button(self.f1_2, text = "Save loss functions", borderwidth=2, width=15, command = self.saveScatterers)
        
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
        self.chart2 = FigureCanvasTkAgg(self.fig2.fig, self.f3)
        self.chart2.mpl_connect('button_press_event', functools.partial(self.doubleClkChart, ax=self.fig2.ax, chart = self.chart2))
        rectprops = dict(facecolor='white', edgecolor = 'black',
                           alpha=0.2, fill=True)
        self.RS2 = RectangleSelector(self.fig2.ax, functools.partial(self.selector, ax=self.fig2.ax, chart = self.chart2),
                                       drawtype='box', useblit=True,
                                       button=[1, 3],  # don't use middle button
                                       minspanx=5, minspany=5,
                                       spancoords='pixels',
                                       interactive=True,
                                       rectprops = rectprops) 
        self.RS2.set_active(True)

        # Scatterers Table
        columns = ('Nr.','Type')
        self.scatterers_table = ttk.Treeview(self.f3_1, height=4,show='headings',columns=columns, selectmode='browse')
        self.scatterers_table.column('Nr.',width=50,anchor=tk.W)
        self.scatterers_table.heading('Nr.', text='Nr.', anchor=tk.W)
        self.scatterers_table.column('Type',width=350,anchor=tk.W)
        self.scatterers_table.heading('Type', text='Type', anchor=tk.W)
        self.scatterers_table.bind('<Double-1>',self.controller.callLossEditor)
        self.scatterers_table.bind('<Delete>', self.removeComponent)
        self.add_comp_btn = tk.Label(self.f3_1, text = "+ Add component", borderwidth=2, padx = 2, pady=0, relief ='raised')
        self.add_comp_btn.bind('<Button-1>', self.addComponent)
        
    def setupLayout(self):
        self.f1.pack(side=tk.LEFT, fill=None, anchor="nw")
        self.f1_1.pack(side=tk.TOP, expand=False, fill=tk.Y, anchor='center', padx=15, pady=5)
        self.step1_label.pack(side=tk.TOP, pady=5)
        self.f1_2.pack(side=tk.TOP, expand=False, fill=tk.Y, anchor='center', padx=15, pady=5)
        self.step2_label.pack(side=tk.TOP, pady=1)
        
        self.f1_3.pack(side=tk.TOP, fill = None, anchor='center',pady=1)
        self.step3_label.pack(side=tk.TOP, pady=1)
        self.f1_3_1.pack(side=tk.TOP)
        self.f1_3_1_1.pack(side=tk.LEFT, anchor = tk.N)
        self.f1_3_1_2.pack(side=tk.LEFT)
        self.f1_3_1_3.pack(side=tk.LEFT)
        self.f1_3_2.pack(side=tk.TOP)
        
        self.f1_4.pack(side=tk.TOP, fill = None, anchor='center', pady=1)
        self.step4_label.pack(side=tk.TOP, pady=1)
        #self.variant_chk.pack(side=tk.TOP)
        self.simulate_label.pack(side=tk.TOP)
        self.scatter_btn.pack(side=tk.TOP)
        self.unscatter_btn.pack(side=tk.TOP)
        self.bulk_chk.pack(side=tk.TOP)
        
        self.f2.pack(side=tk.LEFT, fill=tk.Y, anchor='n')
        self.btn1.pack(side=tk.TOP, fill = None, pady=2)
        self.btn2.pack(side=tk.TOP, fill = None, pady=2)
        
        # Figure XPS
        self.chart1.get_tk_widget().pack(side=tk.TOP)
        self.f2_1.pack(side=tk.TOP, fill=tk.Y, anchor='ne')
        self.spectra_table.pack(side=tk.TOP, anchor="ne", pady=10, padx=15)
        self.normalize_chk.pack(side=tk.TOP, anchor="ne", pady=10, padx=15)
        
        self.f3.pack(side=tk.LEFT, fill = None, anchor='n')
        self.btn3.pack(side=tk.TOP)
        self.load_loss_frame.pack(side=tk.TOP)
        self.load_loss_label.pack(side=tk.TOP)
        self.cbox.pack(side=tk.TOP)
        self.btn4.pack(side=tk.TOP)
        self.btn5.pack(side=tk.TOP)
        self.chart2.get_tk_widget().pack(side=tk.TOP)
        
        self.f3_1.pack(side=tk.TOP, anchor="ne")
        self.scatterers_table.pack(side=tk.TOP, anchor="ne", pady=10, padx=15)
        self.add_comp_btn.pack(side=tk.TOP, anchor="ne", padx=15)
        
        self.toggleParamVariants()
        
    def toggleParamVariants(self):
        variant = self.variant.get()
        widgets = self.variants_dict[variant]
        
        for child in self.f1_3_1.children.values():
            for next_child in child.children.values():
                next_child.pack_forget()
        
        x_pad = 2
        self.f1_3_1_1.pack(side=tk.LEFT, anchor = tk.S, padx=x_pad)
        self.f1_3_1_2.pack(side=tk.LEFT, anchor = tk.S, padx=x_pad)
        self.f1_3_1_3.pack(side=tk.LEFT, anchor = tk.S, padx=x_pad)
                
        for widget in widgets:
            widget.pack(side=tk.TOP, anchor = tk.S)   
      
    def loadSpectrum(self):
        file = filedialog.askopenfilename(initialdir = self.controller.datapath)
        if file:
            self.controller.loadSpectrum(file)
            self.controller.datapath = (file.rsplit('/',maxsplit = 1)[0])

    def loadScatterers(self):
        file = filedialog.askopenfilename(initialdir=self.controller.datapath,
                                          title='Select loss function file',
                                          filetypes=[('json', '*.json')])
        self.controller.loadScatterers(file)
        try:
            self.cbox.set('default')
            self.controller.setCurrentScatterer()
        except:
            pass
        
    def setCurrentScatterer(self, event):
        self.controller.setCurrentScatterer()
        
    def saveScatterers(self):
        file = filedialog.asksaveasfilename(initialdir=self.controller.datapath,
                                            title='Save as',
                                            filetypes=[('json', '*.json')])
        self.controller.saveScatterers(file)

    def updateScattererChoices(self, choices):
        self.cbox.set('')
        self.scatterer_choices = choices
        self.cbox.config(values=self.scatterer_choices)

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
        if len(self.spectra_table.get_children()) == 0:
            start = self.default_start
            stop = self.default_stop
        else:
            start, stop = self.fig1.ax.get_xlim()
        step = self.default_step
        self.controller.addSynthSpec(start, stop, step)
        self.controller.fillTable1()
        self.selected_spectrum = len(self.spectra_table.get_children())-1
        self.spec_builder = SpecBuilder(self.controller,self.selected_spectrum)
        self.spec_builder.start = start
        self.spec_builder.stop = stop
        self.spec_builder.step = step
        
    def editSynthSpec(self, idx):
        self.spec_builder = SpecBuilder(self.controller,idx)
        self.spec_builder.setSelection(idx)
        self.default_step = self.spec_builder.step

    def selector(self, eclick, erelease, ax, chart):
        
        if eclick.dblclick: # in the case that the user double clicks, we dont want to use the selector. We use zoomOut
            pass
        else:
            xvals = [eclick.xdata,erelease.xdata]
            yvals = [eclick.ydata,erelease.ydata]
            ax.set_xlim(min(xvals),max(xvals))
            ax.set_ylim(min(yvals),max(yvals))
            chart.draw()
        
    def doubleClkChart(self, event, ax, chart): 
        if event.dblclick:
            self.zoomOut(ax=ax, chart=chart)
            
    def zoomOut(self, ax, chart):
        self.press = None
        ax.relim()
        ax.autoscale()
        chart.draw()
        
    def newScatterer(self):
        self.new_scatterer = newScatterer(self, self.controller)
            
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
        self.refreshTable()
        self.start = 0
        self.stop = 0
        self.step = 0
        
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
        self.position = tk.StringVar()
        self.position.trace('w', self.modPeak)
        self.position_label = tk.Label(self.position_frame, text = 'Position')
        self.position_entry = tk.Entry(self.position_frame, width = 10, textvariable = self.position)
        #self.position_entry.config(validate='key', validatecommand = self.modPeak)
        
        self.width_frame = tk.Frame(self.entries_frame)
        self.width = tk.StringVar()
        self.width.trace('w', self.modPeak)
        self.width_label = tk.Label(self.width_frame, text = 'Width')
        self.width_entry = tk.Entry(self.width_frame, width = 10, textvariable = self.width)
        #self.width_entry.config(validate='key', validatecommand = self.modPeak)
        
        self.intensity_frame = tk.Frame(self.entries_frame)
        self.intensity = tk.StringVar()
        self.intensity.trace('w', self.modPeak)
        self.intensity_label = tk.Label(self.intensity_frame, text = 'Intensity')
        self.intensity_entry = tk.Entry(self.intensity_frame, width = 10, textvariable = self.intensity)
        #self.intensity_entry.config(validate='key', validatecommand = self.modPeak)
        
        self.edit_range_frame = tk.Frame(self.entries_frame)
        self.edit_range_label = tk.Label(self.edit_range_frame, text = '  ')
        self.edit_range_button = tk.Button(self.edit_range_frame, text='Edit range',command=self.editRange)
        
    def setupLayout(self):
        padding = 5
        self.entries_frame.pack(side=tk.TOP, padx=padding, pady=padding)
        self.position_frame.pack(side=tk.LEFT, padx=padding, pady=padding)
        self.width_frame.pack(side=tk.LEFT, padx=padding, pady=padding)
        self.intensity_frame.pack(side=tk.LEFT, padx=padding, pady=padding)
        self.edit_range_frame.pack(side=tk.LEFT, padx=padding, pady=padding)
        
        self.position_label.pack(side=tk.TOP, padx=padding, pady=padding)
        self.position_entry.pack(side=tk.TOP, padx=padding, pady=padding)
        
        self.width_label.pack(side=tk.TOP, padx=padding, pady=padding)
        self.width_entry.pack(side=tk.TOP, padx=padding, pady=padding)

        self.intensity_label.pack(side=tk.TOP, padx=padding, pady=padding)
        self.intensity_entry.pack(side=tk.TOP, padx=padding, pady=padding)
        
        self.edit_range_label.pack(side=tk.TOP, padx=padding, pady=padding)
        self.edit_range_button.pack(side=tk.TOP, padx=padding, pady=padding)

        self.peak_table.pack(side=tk.TOP, padx=15, pady=0)
        self.add_comp_btn.pack(side=tk.TOP,anchor='ne', padx=15)
        self.done.pack(side=tk.TOP, pady=15)
        self.edit_range_button.pack(side=tk.TOP)
        
    def Done(self):
        self.window.destroy()
              
    def removeComponent(self):
        return
    
    def addComponent(self, event):
        def setChoice(choice):
            self.controller.addPeak(self.spec_idx, choice)
            self.refreshTable()
            self.controller.rePlotFig1()
            self.setSelection(len(self.peak_table.get_children())-1)
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
        pos = str(cur_item[2])
        width = str(cur_item[3])
        intensity = str(cur_item[4])
        self.position.set(pos)
        self.width.set(width)
        self.intensity.set(intensity)

    def modPeak(self, *args):
        if self.window.focus_get().__class__.__name__ == 'Entry': 
            position = self.position.get()
            width = self.width.get()
            intensity = self.intensity.get()
            if (position != '') & (width != '') & (intensity != ''):  
                position = float(position)
                width = float(width)
                intensity = float(intensity)
                new_values = {'spec_idx':self.spec_idx,'peak_idx':self.selected_peak,'position': position, 'width':width, 'intensity': intensity}
                self.controller.updatePeak(new_values)
                self.updateEntry(new_values)
        else:
            return
            
    def updateEntry(self, new_values):
        peak_idx = new_values['peak_idx']
        peak_type = self.peak_table.item(self.selected_peak)['values'][1]
        values = (new_values['peak_idx'],peak_type, new_values['position'],new_values['width'],new_values['intensity'])
        self.peak_table.item(peak_idx, values = values)
        
    def editRange(self):
        self.range_editor = RangeEditor(self.controller, self.spec_idx)
        self.range_editor.setParams(self.start, self.stop, self.step)
        
    def setSelection(self, idx):  
        self.selected_peak = idx
        self.peak_table.selection_set(str(idx))
        cur_item = self.peak_table.item(idx)['values']
        #peak_type = cur_item[1]
        pos = str(cur_item[2])
        width = str(cur_item[3])
        intensity = str(cur_item[4])
        self.position.set(pos)
        self.width.set(width)
        self.intensity.set(intensity)
    
class RangeEditor:
    def __init__(self, controller, spec_idx):
        self.controller = controller
        self.spec_idx = spec_idx
        self.bcolor = 'grey'
        self.bthickness = 0
        self.window = tk.Toplevel()
        title = 'Range: '
        self.window.wm_title(title)
        self.window.attributes("-topmost", True)
        header = tk.Label(self.window, text = 'Enter range for spectrum')
        header.pack(side=tk.TOP)
        self.createWidgets()
        self.setupLayout()
        
    def createWidgets(self):
        self.container = tk.Frame(self.window)
        self.f1 = tk.Frame(self.container)
        self.f2 = tk.Frame(self.container)
        self.f3 = tk.Frame(self.container)
        
        self.start = tk.DoubleVar()
        self.start_label = tk.Label(self.f1, text = 'Start')
        self.start_entry = tk.Entry(self.f1, width = 10, textvariable = self.start)
        
        self.stop = tk.DoubleVar()
        self.stop_label = tk.Label(self.f2, text = 'Stop')
        self.stop_entry = tk.Entry(self.f2, width = 10, textvariable = self.stop)
  
        self.step = tk.DoubleVar()
        self.step_label = tk.Label(self.f3, text = 'Step')
        self.step_entry = tk.Entry(self.f3, width = 10, textvariable = self.step)
        
        self.done_button = tk.Button(self.window, text = 'Done', command = self.Done)
        
    def setupLayout(self):
        self.container.pack(side=tk.TOP, pady=5)
        self.f1.pack(side=tk.LEFT, padx=5)
        self.f2.pack(side=tk.LEFT, padx=5)  
        self.f3.pack(side=tk.LEFT, padx=5)
        self.start_label.pack(side=tk.TOP, pady=5)
        self.start_entry.pack(side=tk.TOP, pady=5)
        self.stop_label.pack(side=tk.TOP, pady=5)
        self.stop_entry.pack(side=tk.TOP, pady=5)
        self.step_label.pack(side=tk.TOP, pady=5)
        self.step_entry.pack(side=tk.TOP, pady=5)
        self.done_button.pack(side=tk.TOP,pady=5)
        
    def Done(self):
        start = self.start.get()
        stop = self.stop.get()
        step = self.step.get()
        idx = self.spec_idx
        self.controller.editRange(idx, start, stop, step)
        self.window.destroy()
        
    def setParams(self, start, stop, step):
        self.start.set(start)
        self.stop.set(stop)
        self.step.set(step)
        
        
class newScatterer:
    def __init__(self, parent, controller):
        self.parent = parent
        self.controller = controller
        self.bcolor = 'grey'
        self.bthickness = 0
        self.window = tk.Toplevel()
        title = 'New Scatterer'
        self.window.wm_title(title)
        self.window.attributes("-topmost", True)
        self.header = tk.Label(self.window, text = 'Enter a name for the scatterer.')
        self.header.pack(side=tk.TOP, padx = 15, pady = 15)
        self.name = tk.StringVar()
        self.name_entry = tk.Entry(self.window, width = 20, textvariable = self.name)
        self.name_entry.pack(side=tk.TOP,pady=5)
        self.done = tk.Button(self.window, text='Done', command = self.Done)
        self.done.pack(side=tk.TOP, pady=5)
        self.name_entry.focus()
        
    def Done(self):
        self.controller.newScatterer(self.name.get())
        self.parent.cbox.set(self.name.get())
        self.controller.setCurrentScatterer()
        self.window.destroy()
        

        

if __name__ == "__main__":
    mainwin = tk.Tk()
    app = View(mainwin)
    mainwin.mainloop()
