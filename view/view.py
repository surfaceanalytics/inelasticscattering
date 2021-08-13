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
import threading
import queue

from view.inputs_frame import InputsFrame
from view.specbuilder import SpecBuilder
from view.figure import Figure
from view.table import Table
from view.specselector import SpecSelector
from view.tooltip import CreateToolTip
from view.popupmenu import PopUpMenu

class View:
    def __init__(self, controller, root):
        self.controller = controller
        self.container = root
        self.colours = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
                        '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
        self.queue = queue.Queue()
        self.container.title('Scatter Simulator')
        self.selected_spectrum = ''
        self._setup()
        
        self.s = tk.Style()
        self.s.theme_use('vista')
        self.s.configure("vista.TFrame", padding=100)
           
        #print(self.s.lookup('TFrame', 'font'))
        #print(tk.Style().lookup("TButton", "font"))

    def _setup(self):
        self._createWidgets()
        self._setupLayout()
        self._addTooltips()
           
    def _createWidgets(self):
        """ This function instantiates all of the widgets.
        
        Parameters
        ----------
        None
        
        Returns
        -------
        None
        """
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
                                width = 15)#,
                                #command = self.exportFile)
        self.export.bind("<Button-1>", self.exportPopUp)
        
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
                                       command = self.controller.unScatterSpectrum,
                                       width=15)    
        self.bulk = IntVar()
        self.bulk_chk = tk.Checkbutton(self.f1_4, 
                                       text="Bulk spectrum", 
                                       variable=self.bulk)
        
        # Figure for XPS spectra
        x = []
        y = []
        params = {'xlabel':'Energy [eV]', 'ylabel': 'Intensity [cts./sec.]',
               'title':'Spectra', 'colours':self.colours}
        self.fig1 = Figure(self.f2,x,y,params)
        
        # Loss function figure (Figure2)
        params = {'xlabel':'Energy [eV]', 'ylabel':'Probability', 
               'title':'Loss Function', 'colours':self.colours}
        self.fig2 = Figure(self.f3, x,y,params)

        self.figures = [self.fig1, self.fig2]

        # XPS spectra table
        
        spec_table_params = {'table_name':'spectra', 
                             'on_double_click':'visibility',
                  'columns':[{'name':'Nr.', 
                              'width':50},
                             {'name':'Type',
                              'width':125,
                              'popup_choices':['none',
                                               'Scattered',
                                               'Unscattered']
                              },
                             {'name':'Visibility',
                              'width':120,
                              'popup_choices':['visible', 'hidden']},
                             {'name':'Edit',
                              'width':45,
                              'link':'edit'}
                             ]
                  }
        
        loss_table_params = {'table_name':'loss_function', 
                             'on_double_click':'edit',
                             'colour_keys':False,
                             'columns':[
                                     {'name':'Nr.', 'width':50},
                                     {'name':'Type','width':350}]
                             }
             
        self.spectra_table = Table(self.container,
                                   self.controller, 
                                   self.f2_1, spec_table_params)
        
        self.loss_function_table = Table(self.container,
                                         self.controller,
                                         self.f3_1, loss_table_params)
        
        self.tables = [self.spectra_table, self.loss_function_table]

        # Normalize spectra box
        self.normalize = IntVar()
        self.normalize.set(0)
        
        self.normalize_chk = tk.Checkbutton(self.f2_1, 
                                            text="Normalize", 
                                            variable=self.normalize, 
                                            command = lambda: 
                                                self.controller.changeNormalization(self.normalize.get()))
        

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
        


        ''' This is the button for adding a new component to the loss function.
        '''
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
        
        self.f2_1.pack(side=TOP, fill=Y, anchor='ne')
        self.spectra_table.table.pack(side=TOP, anchor="ne", pady=10, padx=15)
        self.normalize_chk.pack(side=TOP, anchor="ne", padx=15)
        
        self.f3.pack(side=LEFT, fill = None, anchor='n')
        self.btn3.pack(side=TOP)
        self.load_loss_frame.pack(side=TOP)
        self.load_loss_label.pack(side=TOP)
        self.cbox.pack(side=TOP)
        self.btn4.pack(side=TOP)
        self.btn5.pack(side=TOP)
        
        self.f3_1.pack(side=TOP, anchor="ne")
        self.loss_function_table.table.pack(side=TOP, anchor="ne", padx=5)
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
                      [self.normalize_chk, '''Checking this box will normalize all spectra.'''],
                      [self.unscatter_btn, '''Run the calculation to remove inelastic scattering signal.''']
                     ]
        
        for t in tool_tips:
            CreateToolTip(t[0], self.controller.thread_pool_executer, t[1])

    def addVariantChoices(self, params):
        """ This function make radio buttons for each of the algorithm
        variants."""
        self.variant = StringVar()
        for p in params:
            tk.Radiobutton(self.f1_3_2, 
                           text=str(p), 
                           variable=self.variant, 
                           command = self.toggleVariant, 
                           value=p).pack(side=LEFT)
        
    def toggleVariant(self):
        """ This function tells the controller when the choice of Algorithm
        has changed."""
        new_var = self.variant.get()
        self.controller.toggleVariant(new_var)
        
    def buildAlgorithmFrame(self, params):
        """ This method is called by the controller.
        """
        self.inputs_frame.buildFrame(params)
    
    def loadSpectrum(self):
        file = filedialog.askopenfilename(initialdir = self.controller.datapath)
        if file:
            self.controller.loadFile(file)
            self.controller.datapath = (file.rsplit('/',maxsplit = 1)[0])
            
    def callSpecSelector(self, table_params, fig_params):
        self.specselector = SpecSelector(self.controller,
                                         table_params, 
                                         fig_params)

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
        
    def exportPopUp(self, event):
        """
        This creates a PopUp menu on the export button.

        Parameters
        ----------
        event : event

        Returns
        -------
        None.

        """
        choices = self.controller.getExportFormats()  
        popup = PopUpMenu(self.container, choices, self.exportFile)
        popup.pop(event)
        
    def exportFile(self, file_format):
        """
        This is the callback function that is used in the file export PopUp
        menu.

        Parameters
        ----------
        file_format : STRING
            A key that identifies the desired export file format.

        Returns
        -------
        None.

        """
        file = filedialog.asksaveasfilename(initialdir=self.controller.datapath,
                                    title='Save as')
        self.controller.export(file, file_format)
        
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
        
    def showWarning(self, warning_msg):
        tkinter.messagebox.showerror('Warning',
                                     warning_msg)
        
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
        """ This function adds a new loss function component."""
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
        self.selected_spectrum = len(self.spectra_table.table.get_children())-1
        self.controller.addSynthSpec()
        self.spectra_table.fillTable()
        
    def editSynthSpec(self, params):
        self.spec_builder = SpecBuilder(self.controller,params)
  
    def newScatterer(self):
        self.new_scatterer = newScatterer(self, self.controller)

class LossEditor:
    """ This is a pop-up window that is called when the user adds or edits a  
    component of the loss function
    """
    def __init__(self, controller, params, comp_nr, comp_type):
        self.controller = controller
        self.bcolor = 'grey'
        self.bthickness = 0
        self.params = params
        self.comp_nr = comp_nr
        padding = 15
        self.window = Toplevel(padx=padding, pady=padding)
        
        title = str(comp_nr) + ": " + str(comp_type)
        self.window.wm_title(title)
        self.window.attributes("-topmost", True)
        header = tk.Label(self.window, text = 'Component Nr.: ' 
                          + str(comp_nr) + '\n Type: ' 
                          + str(comp_type), anchor="center", justify="center")

        header.pack(side=TOP)
        
        self._buildEntryFields(params)
            
        x = self.controller.root.winfo_x()
        y = self.controller.root.winfo_y()
                
        w = self.controller.root.winfo_width()
        h = self.controller.root.winfo_height()
        dx = 600
        dy = 200
        self.window.geometry("+%d+%d" % (x+w-dx,y+h-dy))
        
    def _buildEntryFields(self, params):
        """
        This function constructs the user entry input fields.

        Parameters
        ----------
        params : A dictionary that contains keys (STRING) that represent the
        attribute names of the component, and values. 

        Returns
        -------
        None.

        """
        self.labels = []
        self.entries = []
        self.stringvars = []
        i=0
        for key in params:
            subsubframe = tk.Frame(self.window)
            subsubframe.pack(side=TOP)
            self.stringvars += [StringVar()]
            self.stringvars[i].set(params[key])
            self.stringvars[i].trace('w', self._callBack)
            self.labels += [tk.Label(subsubframe, text = key)]
            self.labels[i].pack(side=TOP)
            self.entries += [tk.Entry(subsubframe, width = 20, 
                                      textvariable = self.stringvars[i])]
            self.entries[i].pack(side=TOP)
            i+=1

    def _callBack(self,event, *args):
        if self.buildDict() == 1:
            self.controller.modifyLossLineshape()
            
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
    """ This is a pop-up window that is called when the user wants to create a
    new scatterer.

    """
    def __init__(self, parent, controller):
        """
        Parameters
        ----------
        parent : An instance of a View object that contains a combobox, called
        cbox.
        controller : An instance of a Controller class.

        Returns
        -------
        None.

        """           
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
        name = self.name.get() 
        self.controller.newScatterer(name)
        self.parent.cbox.set(name)
        self.controller.setCurrentScatterer(name)
        self.window.destroy()

if __name__ == "__main__":
    mainwin = tk.Tk()
    app = View(mainwin)
    mainwin.mainloop()
