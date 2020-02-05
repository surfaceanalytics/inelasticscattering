# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 22:31:45 2020

@author: Mark
"""
import tkinter as tk
import tkinter.ttk as ttk
import functools

class Controller():
    def __init__(self):
        self.root = tk.Tk()
        self.table_values = [(0,'A','1'),(1,'A','2'),(2,'B','3')]
        self.table1_values = [(0,'D','4'),(1,'D','5'),(2,'E','6')]
        self.table_choices = {1:['A','B','C'],2:['1','2','3']}
        self.table1_choices = {1:['D','E','F'],2:['4','5','6']}
        self.view = View(self, self.root)
        self.fillTable(self.view.table,self.table_values)
        self.fillTable(self.view.table1,self.table1_values)
        
    def fillTable(self, table, values):
        for i,j in enumerate(values):
            table.insert('',i,values=j)
        
    def tablePopup(self, event, table, table_choices):
        row = table.identify_row(event.y)
        col = table.identify_column(event.x)
        col = int(col.lstrip('#'))-1
        def setChoice(choice):
            table.set(row,column=col,value=choice)
        popup = tk.Menu(self.view.root, tearoff=0)
        if col > 0:
            choices = table_choices[col]
            for i,j in enumerate(choices):
                # This line below is crazy. Needs to be done this way because the i in the loop is only scoped for the loop, and does not persist
                popup.add_command(command = lambda choice = choices[i]: setChoice(choice), label=j)
        try:
            popup.tk_popup(event.x_root, event.y_root, 0)
        finally:
            popup.grab_release()

class View():
    def __init__(self, controller, root):
        self.root = root
        self.controller = controller
        
        columns = ('ID','Col1', 'Col2')
        self.table = ttk.Treeview(self.root, height=4,show='headings',columns=columns, selectmode='browse')
        self.table.column('ID',width=50,anchor=tk.W)
        self.table.heading('ID', text='ID', anchor=tk.W)
        self.table.column('Col1',width=150,anchor=tk.W)
        self.table.heading('Col1', text='Col1', anchor=tk.W)        
        #self.table.bind('<<TreeviewSelect>>',self.controller.tableSelection)
        self.table.column('Col2',width=200,anchor=tk.W)
        self.table.heading('Col2', text='Col2', anchor=tk.W)
        self.table.pack(side=tk.TOP, anchor="center")
        self.table.bind('<Button-3>',functools.partial(self.controller.tablePopup, table=self.table, table_choices = self.controller.table_choices))

        columns = ('ID','Col3', 'Col4')
        self.table1 = ttk.Treeview(self.root, height=4,show='headings',columns=columns, selectmode='browse')
        self.table1.column('ID',width=50,anchor=tk.W)
        self.table1.heading('ID', text='ID', anchor=tk.W)
        self.table1.column('Col3',width=150,anchor=tk.W)
        self.table1.heading('Col3', text='Col3', anchor=tk.W)        
        #self.table.bind('<<TreeviewSelect>>',self.controller.tableSelection)
        self.table1.column('Col4',width=200,anchor=tk.W)
        self.table1.heading('Col4', text='Col4', anchor=tk.W)
        self.table1.pack(side=tk.TOP, anchor="center")
        self.table1.bind('<Button-3>',functools.partial(self.controller.tablePopup, table=self.table1, table_choices = self.controller.table1_choices))





if __name__ == "__main__":
    app = Controller()
    app.root.mainloop()


