# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 15:41:16 2020

@author: Mark
"""
import tkinter
import tkinter.ttk as tk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backend_bases import MouseEvent, LocationEvent
from tkinter import filedialog
from tkinter.ttk import Combobox, Scrollbar
from tkinter import DoubleVar, StringVar, IntVar, LEFT, TOP, W, Y, N, S, Toplevel, Menu, CENTER
import functools
from matplotlib.widgets import RectangleSelector
from inputs_frame import InputsFrame

from specbuilder import SpecBuilder

import threading
import queue

import time

class View:
    def __init__(self, controller, root):
        self.controller = controller
        self.container = root
        self.queue = queue.Queue()
        self.container.title('Scatter Simulator')
        self.selected_spectrum = ''
        self.setup()
        self.default_start = 0
        self.default_stop = 100
        self.default_step = 0.1
        self.s = tk.Style()
        self.s.theme_use('vista')
        self.s.configure("vista.TFrame", padding=100)
        
        #print(self.s.lookup('TFrame', 'font'))
        #print(tk.Style().lookup("TButton", "font"))

    def setup(self):
        self._createWidgets()
        self._setupLayout()
        self._addTooltips()
           
    def _createWidgets(self):
        self.bcolor = 'grey'
        self.bthickness = 0
        self.axis_label_fontsize = 12
        self.axis_tick_font = 8
        ''' Frames are leveled in a heirarchy
         a frame embedded into another frame will have a name like f1_1_2
        '''
        ''' f1 : left highest-level frame
        controls frame
        contains the user controls for loading, saving files, for changing 
        parameters of simulation, and for running simulation'''
        frame1_padding = 10
        self.f1 = tk.Frame(self.container, 
                           width=300, height=600, 
                           padding=frame1_padding)
        # f2 : middle highest-level frame
        # XPS spectra Frame
        # contains the plot of the XPS spectra, and the list of spectra
        self.f2 = tk.Frame(self.container, 
                           width=400,height=600, 
                           padding=frame1_padding)
        # f3 : right highest-level frame
        # Loss function Frame
        # contains plot of loss function, and loss functions table
        self.f3 = tk.Frame(self.container,
                           width=400,height=600, 
                           padding=frame1_padding)
        # f3_1
        # Loss function figure frame
        self.f3_1 = tk.Frame(self.f3, 
                             width=400,height=100, 
                             padding=frame1_padding)
        # f1_1
        # Loss buttons frame
        self.f1_1 = tk.Frame(self.f1, 
                             width=400,height=600, 
                             padding=frame1_padding)
        self.step2_label = tk.Label(self.f1_1, 
                                    text='1. Scatterers',
                                    font=("Helvetica", 12))
        # f1_1 
        # load spectra frame
        self.f1_2 = tk.Frame(self.f1,
                             width=400,height=200, 
                             padding=frame1_padding)
        self.step1_label = tk.Label(self.f1_2, 
                                    text='2. Spectra',
                                    font=("Helvetica", 12))
        self.btn1 = tk.Button(self.f1_2, 
                              text = "Load spectrum", 
                              width = 15, 
                              command = self.loadSpectrum)
        self.btn2 = tk.Button(self.f1_2, 
                              text = "Build spectrum",
                              width = 15, 
                              command = self.addSynthSpec)
        self.export = tk.Button(self.f1_2,
                                text = 'Export',
                                width = 15,
                                command = self.export)
        
        # f1_3
        # Parameter inputs frame
        self.f1_3 = tk.Frame(self.f1, 
                             width=300, height=500, 
                             padding=frame1_padding)
        self.step3_label = tk.Label(self.f1_3, 
                                    text='3. Parameters',
                                    font=("Helvetica", 12))

        # f1_3_1
        # Parameters subframe
        self.f1_3_1 = tk.Frame(self.f1_3, width=300, height=500)
        self.inputs_frame = InputsFrame(self.f1_3_1)

        # f1_3_2
        # Variants subframe
        self.f1_3_2 = tk.Frame(self.f1_3, width=300, height=30)
                
        # f1_4
        # run simulation frame
        self.f1_4 = tk.Frame(self.f1, width=300, height=500)
        
        # f2_1
        # spectra table frame
        self.f2_1 = tk.Frame(self.f2)
        
        # Run simulation
        self.step4_label = tk.Label(self.f1_4, 
                                    text='4. Simulation',
                                    font=("Helvetica", 12))
        self.scatter_btn = tk.Button(self.f1_4, text = "Scatter", 
                                     command = self.controller.scatterSpectrum, 
                                     width=15)
        self.unscatter_btn = tk.Button(self.f1_4, 
                                       text = "Un-scatter", 
                                       width=15)    
        self.bulk = IntVar()
        self.bulk_chk = tk.Checkbutton(self.f1_4, 
                                       text="Bulk spectrum", 
                                       variable=self.bulk)
        
        # Figure for XPS spectra
        x = []
        y = []
        self.fig1 = Figure(x,y)
        self.fig1.ax.tick_params(direction='out', length=4, width=1, 
                                 colors='black', grid_color='black', 
                                 labelsize=self.axis_tick_font, grid_alpha=0.5)
        self.fig1.ax.set_xlabel('Energy [eV]', 
                                fontsize=self.axis_label_fontsize)
        self.fig1.ax.set_ylabel('Intensity [cts./sec.]', 
                                fontsize=self.axis_label_fontsize)
        self.fig1.ax.set_title('Spectra')
        self.fig1.fig.tight_layout()
        self.chart1 = FigureCanvasTkAgg(self.fig1.fig, self.f2)
        self.chart1.mpl_connect('button_press_event', 
                                functools.partial(self.doubleClkChart, 
                                                  ax=self.fig1.ax, 
                                                  chart = self.chart1))
        rectprops = dict(facecolor='white', edgecolor = 'black',
                           alpha=0.2, fill=True)
        self.RS1 = RectangleSelector(self.fig1.ax, 
                                     functools.partial(self.selector, 
                                                       ax=self.fig1.ax, 
                                                       chart = self.chart1),
                                       drawtype='box', useblit=True,
                                       button=[1, 3],  # don't use middle button
                                       minspanx=5, minspany=5,
                                       spancoords='pixels',
                                       interactive=True,
                                       rectprops = rectprops) 
        self.RS1.set_active(True)
        # XPS spectra table
        columns = ('Nr.','Type', 'Visibility')
        self.spectra_table = tk.Treeview(self.f2_1, 
                                         height=4, 
                                         columns=columns, 
                                         selectmode='browse')
        self.spectra_table.name = 'spectra'
        self.spectra_table.column('#0', width=60, anchor=W)
        self.spectra_table.column('Nr.',width=50,anchor=W)
        self.spectra_table.heading('Nr.', text='Nr.', anchor=W)
        self.spectra_table.column('Type',width=145,anchor=W)
        self.spectra_table.heading('Type', text='Type', anchor=W)        
        self.spectra_table.column('Visibility',width=145,anchor=W)
        self.spectra_table.heading('Visibility', text='Visibility', anchor=W)

        self.spectra_table.entry_choices = {
                '0':['none','Peak','VacuumExcitation'],
                '1':['visible','hidden']}
        self.spectra_table.bind('<Button-3>',
                                functools.partial(self.controller.tablePopup, 
                                                  table=self.spectra_table, 
                                                  table_choices = 
                                                  self.controller.spectra_table_choices))
        self.spectra_table.bind('<Double-1>',
                                functools.partial(self.controller.doubleClkTable, 
                                                  table=self.spectra_table))
        self.spectra_table.bind('<Delete>', self.removeSpectrum)

        # Normalize spectra box
        self.normalize = IntVar()
        self.normalize_chk = tk.Checkbutton(self.f2_1, 
                                            text="Normalize", 
                                            variable=self.normalize, 
                                            command = self.controller.rePlotFig1)

        # Load file of loss functions
        self.btn3 = tk.Button(self.f1_1, text="Load scatterers",
                              width=15, command=self.loadScatterers)

        # Select loss function (from loaded file)
        self.load_loss_frame = tk.Frame(self.f1_1, width=400,height=600)
        self.load_loss_label = tk.Label(self.load_loss_frame, 
                                        text='Select scatterer')
        self.selected_scatterer = StringVar()
        self.scatterer_choices = []
        self.cbox = Combobox(self.load_loss_frame, width=13, 
                             textvariable = self.selected_scatterer,
                             values=self.scatterer_choices)
        self.cbox.bind("<<ComboboxSelected>>", self.setCurrentScatterer)

        # Build loss function
        self.btn4 = tk.Button(self.f1_1, 
                              text = "New scatterer", 
                              width=15, 
                              command = self.newScatterer)
        # Save loss function
        self.btn5 = tk.Button(self.f1_1, 
                              text = "Save scatterer", 
                              width=15, 
                              command = self.saveScatterers)
        
        # Loss function figure (Figure2)
        x = []
        y = []
        self.fig2 = Figure(x,y)
        self.fig2.ax.tick_params(direction='out', 
                                 length=4, 
                                 width=1, 
                                 colors='black',
                                 grid_color='black', 
                                 labelsize=self.axis_tick_font, 
                                 grid_alpha=0.5)
        self.fig2.ax.set_xlabel('Energy Loss [eV]', 
                                fontsize=self.axis_label_fontsize)
        self.fig2.ax.set_ylabel('Probability', 
                                fontsize=self.axis_label_fontsize)
        self.fig2.ax.set_title('Loss Function')
        self.fig2.fig.tight_layout()
        self.chart2 = FigureCanvasTkAgg(self.fig2.fig, self.f3)
        self.chart2.mpl_connect('button_press_event', 
                                functools.partial(self.doubleClkChart, 
                                                  ax=self.fig2.ax, 
                                                  chart = self.chart2))
        rectprops = dict(facecolor='white', edgecolor = 'black',
                           alpha=0.2, fill=True)
        self.RS2 = RectangleSelector(self.fig2.ax, 
                                     functools.partial(self.selector, 
                                                       ax=self.fig2.ax, 
                                                       chart = self.chart2),
                                       drawtype='box', 
                                       useblit=True,button=[1, 3],  # don't use middle button
                                       minspanx=5, minspany=5,
                                       spancoords='pixels',
                                       interactive=True,
                                       rectprops = rectprops) 
        self.RS2.set_active(True)

        # Scatterers Table
        columns = ('Nr.','Type')
        self.scatterers_table = tk.Treeview(self.f3_1,
                                            height=4,show='headings',
                                            columns=columns, 
                                            selectmode='browse')
        self.scatterers_table.column('Nr.',width=50,anchor=W)
        self.scatterers_table.heading('Nr.', text='Nr.', anchor=W)
        self.scatterers_table.column('Type',width=350,anchor=W)
        self.scatterers_table.heading('Type', text='Type', anchor=W)
        self.scatterers_table.bind('<Double-1>',self.controller.callLossEditor)
        self.scatterers_table.bind('<Delete>', self.removeLossComponent)
        self.add_comp_btn = tk.Label(self.f3_1, 
                                     text = "+ Add component", 
                                     relief ='raised')
        self.add_comp_btn.bind('<Button-1>', 
                               self.addComponent)
        
    def _setupLayout(self):
        self.f1.pack(side=LEFT, fill=None, anchor="nw")
        self.f1_1.pack(side=TOP, expand=False, fill=Y, anchor='center')
        self.step1_label.pack(side=TOP)
        self.f1_2.pack(side=TOP, expand=False, fill=Y, anchor='center')
        self.step2_label.pack(side=TOP)
        
        self.f1_3.pack(side=TOP, fill = None, anchor='center')
        self.step3_label.pack(side=TOP)
        self.f1_3_1.pack(side=TOP)
        self.f1_3_2.pack(side=TOP)
        
        self.f1_4.pack(side=TOP, fill = None, anchor='center')
        self.step4_label.pack(side=TOP)
        self.scatter_btn.pack(side=TOP)
        self.unscatter_btn.pack(side=TOP)
        self.bulk_chk.pack(side=TOP)
        
        self.f2.pack(side=LEFT, fill=Y, anchor='n')
        self.btn1.pack(side=TOP, fill = None)
        self.btn2.pack(side=TOP, fill = None)
        
        self.export.pack(side=TOP, fill = None)
        
        # Figure XPS
        self.chart1.get_tk_widget().pack(side=TOP)
        self.f2_1.pack(side=TOP, fill=Y, anchor='ne')
        self.spectra_table.pack(side=TOP, anchor="ne", pady=10, padx=15)
        self.normalize_chk.pack(side=TOP, anchor="ne", padx=15)
        
        self.f3.pack(side=LEFT, fill = None, anchor='n')
        self.btn3.pack(side=TOP)
        self.load_loss_frame.pack(side=TOP)
        self.load_loss_label.pack(side=TOP)
        self.cbox.pack(side=TOP)
        self.btn4.pack(side=TOP)
        self.btn5.pack(side=TOP)
        self.chart2.get_tk_widget().pack(side=TOP)
        
        self.f3_1.pack(side=TOP, anchor="ne")
        self.scatterers_table.pack(side=TOP, anchor="ne", padx=5)
        self.add_comp_btn.pack(side=TOP, anchor="se", pady=10, padx=5)
        
    def _addTooltips(self):
        tool_tips = [[self.btn3, "This will load scatterers from a file."],
                     [self.btn4, "Create a new scatterer. \nGive it a name. \nThen build your own Loss Function"],
                     [self.btn5, "Save the scatterer you built"],
                     [self.btn1, "Load a measured spectrum"],
                     [self.btn2, "Create your own spectrum"],
                     [self.export, "Save the simulated spectra to Vamas or Excel."],
                     [self.step3_label, '''Here you can enter the parameters used in the simulation algorithm.
P is the pressure in mbar.
D is the distance between the sample and the detector.
Inel.X is the inelastic scattering cross section in nm^2
El.X is the eleastic scattering cross section.
f(Decay) are decay factors. A value of 1 means no decay.'''],
                      [self.scatter_btn, "Run the calculation to simulate inelastic scattering."],
                      [self.bulk_chk, '''Checking this box will configure the calculation to simulate the bulk signal.
In this case, P and D have no effect.'''],
                      [self.normalize_chk, '''Checking this box will normalize all spectra.''']
                     ]
        
        for t in tool_tips:
            CreateToolTip(t[0], self.controller.thread_pool_executer, t[1])

    def addVariantChoices(self, params):
        # make radio buttons for variants
        self.variant = StringVar()
        for p in params:
            tk.Radiobutton(self.f1_3_2, 
                           text=str(p), 
                           variable=self.variant, 
                           command = self.toggleVariant, 
                           value=p).pack(side=LEFT)
        
    def toggleVariant(self):
        new_var = self.variant.get()
        self.controller.toggleVariant(new_var)
        
    def buildAlgorithmFrame(self, params):
        self.inputs_frame.buildFrame(params)
    
    def loadSpectrum(self):
        file = filedialog.askopenfilename(initialdir = self.controller.datapath)
        if file:
            self.controller.loadSpectrum(file)
            self.controller.datapath = (file.rsplit('/',maxsplit = 1)[0])

    def loadScatterers(self):
        def getFile():
            self.file = filedialog.askopenfilename(initialdir=self.controller.datapath,
                                              title='Select loss function file',
                                              filetypes=[('json', '*.json')])
        
        threading.Thread(target=getFile()).start()
        self.controller.loadScatterers(self.file)
        
    def setCurrentScatterer(self, event):
        label = self.selected_scatterer.get()
        self.controller.setCurrentScatterer(label)
        
    def export(self):
        file = filedialog.asksaveasfilename(initialdir=self.controller.datapath,
                                    title='Save as',
                                    filetypes=[('Vamas','*.vms'),
                                               ('Excel', '*.xlsx')])
        self.controller.export(file)
        
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
        tkinter.messagebox.showerror('Error',
                                     'Please select an Unscattered spectrum and a Loss Function')        
        
    def createSpectrum(self):
        self.controller.createSynthetic()
        
    def removeSpectrum(self, event):
        selected_item = self.spectra_table.selection()
        cur_item = self.spectra_table.item(selected_item[0])['values'][0]
        self.spectra_table.delete(selected_item[0])
        self.controller.removeSpectrum(cur_item)
        
    def removeLossComponent(self, event):
        selected_item = self.scatterers_table.selection()
        cur_item = self.scatterers_table.item(selected_item[0])['values'][0]
        self.scatterers_table.delete(selected_item[0])
        self.controller.removeLossComponent(cur_item)
        
    def addComponent(self, event):
        def setChoice(choice):
            self.controller.addComponent(choice)
        popup = Menu(self.container, tearoff=0)
        choices = self.controller.component_choices
        for i,j in enumerate(choices):
                '''This line below is a bit tricky. Needs to be done this way 
                because the i in the loop is only scoped for the loop, and 
                does not persist'''
                popup.add_command(command = lambda choice = choices[i]: 
                    setChoice(choice), label=j)
        popup.post(event.x_root, event.y_root)
        
    def addSynthSpec(self):
        """ This function calls the controller to tell it a new synthetic
        spectrum should be built.
        """
        if len(self.spectra_table.get_children()) == 0:
            start = self.default_start
            stop = self.default_stop
        else:
            start, stop = self.fig1.ax.get_xlim()
        step = self.default_step

        self.selected_spectrum = len(self.spectra_table.get_children())-1
        
        params = self.controller.addSynthSpec(start, stop, step)
        self.spec_builder = SpecBuilder(self.controller, params)
        
        self.controller.fillTable1()
        
        self.spec_builder.start = start
        self.spec_builder.stop = stop
        self.spec_builder.step = step
        
    def editSynthSpec(self, params):
        self.spec_builder = SpecBuilder(self.controller,params)

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
    ''' This is a pop-up window that is called when the user adds or edits a  
    component of the loss function
    '''
    def __init__(self, controller, params, comp_nr, comp_type):
        self.controller = controller
        self.bcolor = 'grey'
        self.bthickness = 0
        self.params = params
        self.comp_nr = comp_nr
        padding = 15
        window = Toplevel(padx=padding, pady=padding)
        
      
        
        title = str(comp_nr) + ": " + str(comp_type)
        window.wm_title(title)
        window.attributes("-topmost", True)
        header = tk.Label(window, text = 'Component Nr.: ' 
                          + str(comp_nr) + '\n Type: ' 
                          + str(comp_type), anchor="center", justify="center")

        header.pack(side=TOP)
        
        self.labels = []
        self.entries = []
        self.stringvars = []
        i=0
        for key in params:
            subsubframe = tk.Frame(window)
            subsubframe.pack(side=TOP)
            self.stringvars += [StringVar()]
            self.stringvars[i].set(params[key])
            self.stringvars[i].trace('w', self.callBack)
            self.labels += [tk.Label(subsubframe, text = key)]
            self.labels[i].pack(side=TOP)
            self.entries += [tk.Entry(subsubframe, width = 20, 
                                      textvariable = self.stringvars[i])]
            self.entries[i].pack(side=TOP)
            i+=1
            
        x = self.controller.root.winfo_x()
        y = self.controller.root.winfo_y()
                
        w = self.controller.root.winfo_width()
        h = self.controller.root.winfo_height()
        dx = 600
        dy = 200
        window.geometry("+%d+%d" % (x+w-dx,y+h-dy))

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
            
class newScatterer:
    ''' This is a pop-up window that is called when the user wants to create a
    new scatterer.
    '''
    def __init__(self, parent, controller):
        self.parent = parent
        self.controller = controller
        self.bcolor = 'grey'
        self.bthickness = 0
        padding = 15
        self.window = Toplevel(padx=padding, pady=padding)
        title = 'New Scatterer'
        self.window.wm_title(title)
        self.window.attributes("-topmost", True)
        self.header = tk.Label(self.window, 
                               text = 'Enter a name for the scatterer.')
        self.header.pack(side=TOP)
        self.name = StringVar()
        self.name_entry = tk.Entry(self.window, 
                                   width = 20, 
                                   textvariable = self.name)
        self.name_entry.pack(side=TOP)
        self.done = tk.Button(self.window, 
                              text='Done', 
                              command = self.Done)
        self.done.pack(side=TOP)
        self.name_entry.focus()
        
    def Done(self):
        self.controller.newScatterer(self.name.get())
        self.parent.cbox.set(self.name.get())
        self.controller.setCurrentScatterer()
        self.window.destroy()
        
        
class CreateToolTip(object):
    '''
    create a tooltip for a given widget
    '''
    def __init__(self, widget, thread, text='widget info'):
        self.widget = widget
        self.text = text
        self.widget.bind("<Enter>", self.timer)
        self.widget.bind("<Leave>", self.close)
        self.thread = thread
        
    def timer(self, event=None):
        self.inside = True
        def go():
            time.sleep(0.1)
            if self.inside == True:
                self.enter()
                
        self.thread.submit(go())
        #t.start()
         
    def enter(self): 
        x = y = 0
        x, y, cx, cy = self.widget.bbox("insert")
        x += self.widget.winfo_rootx() + 25
        y += self.widget.winfo_rooty() + 20
        # creates a toplevel window
        self.tw = Toplevel(self.widget)
        # Leaves only the label and removes the app window
        self.tw.wm_overrideredirect(True)
        self.tw.wm_attributes("-alpha",0.9)
        self.tw.wm_geometry("+%d+%d" % (x, y))
        label = tk.Label(self.tw, text=self.text, justify='left',
                       background='white', relief='solid', borderwidth=1,
                       font=("TkDefautlFont", "8", "normal"))
        
        label.pack(padx=0, pady=0)

    def close(self, event=None):
        self.inside = False
        if hasattr(self, 'tw'): #self.tw:
            self.tw.destroy()

if __name__ == "__main__":
    mainwin = tk.Tk()
    app = View(mainwin)
    mainwin.mainloop()
