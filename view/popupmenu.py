# -*- coding: utf-8 -*-
"""
Created on Tue Aug 25 15:27:52 2020

@author: Mark
"""

from tkinter import Menu

class PopUpMenu():
    def __init__(self, root, choices, callback):
        """
        Parameters
        ----------
        root: TKINTER parent
        choices :LIST
            A list of choices that should be displayed in the pop-up menu

        Returns
        -------
        None.

        """
        self.root = root
        self.choices = choices
        self.choice = ''
        self.callback = callback
        
    def pop(self, event):

        popup = Menu(self.root, tearoff=0)
        
        def setChoice(choice):
            self.callback(choice)
        
        for c in self.choices:
            '''This line below is a bit tricky. Needs to be done this way 
            because the i in the loop is only scoped for the loop, and 
            does not persist.'''
            popup.add_command(command = lambda choice = c: 
                              setChoice(choice), label=c)
        try:
            popup.tk_popup(event.x_root, event.y_root, 0)
        finally:
            popup.grab_release()
            return self.choice
        
            
