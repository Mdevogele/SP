#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 17 21:42:00 2019

@author: maximedevogele
"""
from tkinter import *
import tkinter as tk

LARGE_FONT = ("Verdana",12)

class SPapp(tk.Tk):
    
    def __init__(self, *args, **kwargs):
        
        tk.Tk.__init__(self, *args, **kwargs)
        container = tk.Frame(self)
        
        container.pack(side ="top", fill = "both", expand = True)
        
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)
        
        self.frames = {}
        
        frame = StartPage(container, self)
        
        self.frames[StartPage] = frame
        
        frame.grid(row=0, column=0, sticky="nsew")
        
        self.show_frame(StartPage)
        
    def show_frame(self, cont):
        
        frame = self.frames[cont]
        frame.tkraise()
        
class StartPage(tk.Frame):

    def Select_GUI():
        print('bla')
#        print(GUI_SELECT.get())  
        
#        APP = GUI_SELECT.get()
#        
#        if APP == 'GOODMAN':
#            os.system('SP_SOAR.py')
#    #        app = SP.simpleapp_tk(None)
#    #        app.title('my application')
#    #        app.mainloop()    
#    
#        if APP == 'GMOSS' or APP == 'GMOSN':
#            os.system('SP_GMOS.py')
#    #        app = SP.simpleapp_tk(None)
#    #        app.title('my application')
#    #        app.mainloop()            
    
    
    
        return None
    
    def __init__(self, parent, controller):
        
        tk.Frame.__init__(self,parent)
#        self.title("Spectroscopy pipeline")
        label = tk.Label(self, text="Select your instrument",font=LARGE_FONT)
        label.grid(row=0,column=0)
        
#        self.geometry('350x200')
        
        
        GUI_SELECT = StringVar()

        rad1 = Radiobutton(self,text='DeVeny@DCT', value='Deveny',justify=tk.LEFT,variable = GUI_SELECT)
     
        rad2 = Radiobutton(self,text='GOODMAN@SOAR', value='GOODMAN',justify=tk.LEFT,variable = GUI_SELECT)
     
        rad3 = Radiobutton(self,text='GMOS@Gemini North', value='GMOSS',justify=tk.LEFT,variable = GUI_SELECT)
    
        rad4 = Radiobutton(self,text='GMOS@Gemini South', value='GMOSN',justify=tk.LEFT,variable = GUI_SELECT)
    
     
        rad1.grid(column=0, row=1,sticky=tk.W+tk.S)
     
        rad2.grid(column=0, row=2,sticky=tk.W+tk.S)
    
        rad3.grid(column=0, row=3,sticky=tk.W+tk.S)
        
        rad4.grid(column=0, row=4,sticky=tk.W+tk.S)
    
    
        btn = Button(self, text="Select",command=self.Select_GUI)
     
        btn.grid(column=0, row=6)





if __name__ == "__main__":
    
    # create Tk object instance
    app = SPapp(None)
    app.mainloop()

        