import tkinter
from tkinter import HORIZONTAL, VERTICAL, ttk
from os import path

windowwidth = 500
windowheight = 400

window = tkinter.Tk()
window.geometry("%dx%d" % (windowwidth,windowheight))
window.title("PhoreMic version 0.1")

def callback():
    return None

menu = tkinter.Menu(window)
new_item = tkinter.Menu(menu,tearoff=0)
new_item.add_command(label='Absorption spectrum')
new_item.add_command(label='Emission spectrum')
new_item.add_command(label='Filter profile')
menu.add_cascade(label='Import new', menu=new_item)
menu.add_command(label='Compute PDF', command=callback)
menu.add_command(label='About', command=callback)
window.config(menu=menu)

tab_parent = ttk.Notebook(window,height=windowheight,width=windowwidth)
def conf(event):
    tab_parent.config(height=window.winfo_height(),width=window.winfo_width())
window.bind("<Configure>",conf)

tab1 = ttk.Frame(tab_parent)
tab2 = ttk.Frame(tab_parent)
tab3 = ttk.Frame(tab_parent)
tab4 = ttk.Frame(tab_parent)
tab_parent.add(tab1, text="Parameters")
tab_parent.add(tab2, text="Fluorophores")
tab_parent.add(tab3, text="Spectra")
tab_parent.add(tab4, text="Results")
tab_parent.pack()

#tab 1 widgets
setuprow = 1
setupentrywidth = 8

ttk.Separator(tab1, orient=VERTICAL).grid(column=0, row=0, rowspan=8, sticky='N,S')
ttk.Separator(tab1, orient=VERTICAL).grid(column=3, row=0, rowspan=setuprow+16, sticky='N,S')

ttk.Separator(tab1, orient=HORIZONTAL).grid(row=setuprow-1, columnspan=3, sticky='W,E')
labelsetup = ttk.Label(tab1,text ="Detector setup",font="Brown-gothic 11")
labelsetup.grid(column=1, row=setuprow)
ttk.Separator(tab1, orient=HORIZONTAL).grid(row=setuprow+1, columnspan=3, sticky='W,E')
NAlabel = ttk.Label(tab1,text="Numerical aperture",justify="left",anchor="w")
NAlabel.grid(column=1,row=setuprow+2)
NAentry = ttk.Entry(tab1,width=setupentrywidth,justify="right")
NAentry.grid(column=2,row=setuprow+2)
nlabel = ttk.Label(tab1,text="Refractive index",justify="left",anchor="w")
nlabel.grid(column=1,row=setuprow+3)
nentry = ttk.Entry(tab1,width=setupentrywidth,justify="right")
nentry.grid(column=2,row=setuprow+3)
qefflabel = ttk.Label(tab1,text="Detector qu. eff.",justify="left",anchor="w")
qefflabel.grid(column=1,row=setuprow+4)
qeffentry = ttk.Entry(tab1,width=setupentrywidth,justify="right")
qeffentry.grid(column=2,row=setuprow+4)
pixellabel = ttk.Label(tab1,text="Pixel size [\u03bcm]",justify="left",anchor="w")
pixellabel.grid(column=1,row=setuprow+5)
pixelentry = ttk.Entry(tab1,width=setupentrywidth,justify="right")
pixelentry.grid(column=2,row=setuprow+5)
explabel = ttk.Label(tab1,text="Exposure time [s]",justify="left",anchor="w")
explabel.grid(column=1,row=setuprow+6)
expentry = ttk.Entry(tab1,width=setupentrywidth,justify="right")
expentry.grid(column=2,row=setuprow+6)
ttk.Separator(tab1, orient=HORIZONTAL).grid(row=setuprow+7, columnspan=3, sticky='W,E')

ttk.Separator(tab1, orient=HORIZONTAL).grid(row=setuprow-1, column=3, columnspan=3, sticky='W,E')
labelexc = ttk.Label(tab1,text ="Excitation beam ",font="Brown-gothic 11")
labelexc.grid(column=4, row=setuprow)
ttk.Separator(tab1, orient=HORIZONTAL).grid(row=setuprow+1, column=3, columnspan=3, sticky='W,E')
lambdalabel = ttk.Label(tab1,text="Wavelength [nm]",justify="left",anchor="w")
lambdalabel.grid(column=4,row=setuprow+2)
lambdaentry = ttk.Entry(tab1,width=setupentrywidth,justify="right")
lambdaentry.grid(column=5,row=setuprow+2)
powerlabel = ttk.Label(tab1,text="Beam power [W]",justify="left",anchor="w")
powerlabel.grid(column=4,row=setuprow+3)
powerentry = ttk.Entry(tab1,width=setupentrywidth,justify="right")
powerentry.grid(column=5,row=setuprow+3)
waistlabel = ttk.Label(tab1,text="Waist diameter [\u03bcm]",justify="left",anchor="w")
waistlabel.grid(column=4,row=setuprow+4)
waistentry = ttk.Entry(tab1,width=setupentrywidth,justify="right")
waistentry.grid(column=5,row=setuprow+4)
ttk.Separator(tab1, orient=HORIZONTAL).grid(row=setuprow+5, column=3, columnspan=3, sticky='W,E')
ttk.Separator(tab1, orient=VERTICAL).grid(column=6, row=0, rowspan=setuprow+16, sticky='N,S')

labelexc = ttk.Label(tab1,text ="Additional illumination:",justify="center",anchor="w")
labelexc.grid(column=4, row=setuprow+6)
illuminationvar = tkinter.IntVar(value=1)
illrad1 = ttk.Radiobutton(tab1,text='none', variable = illuminationvar, value=1, width=10)
illrad2 = ttk.Radiobutton(tab1,text='STED', variable = illuminationvar, value=2, width=10)
illrad1.grid(column=5, row=setuprow+6)
illrad2.grid(column=5, row=setuprow+7)
ttk.Separator(tab1, orient=HORIZONTAL).grid(row=setuprow+8, column=3, columnspan=3, sticky='W,E')

STEDlabel = ttk.Label(tab1,text="STED beam parameters",font="Brown-gothic 11")
STEDlabel.grid(column=4,row=setuprow+9)
ttk.Separator(tab1, orient=HORIZONTAL).grid(row=setuprow+10, column=3, columnspan=3, sticky='W,E')
STEDlambdalabel = ttk.Label(tab1,text="Wavelength [nm]",justify="left",anchor="w")
STEDlambdalabel.grid(column=4,row=setuprow+11)
STEDlambdaentry = ttk.Entry(tab1,width=setupentrywidth,justify="right")
STEDlambdaentry.grid(column=5,row=setuprow+11)
STEDpowerlabel = ttk.Label(tab1,text="Beam power [W]",justify="left",anchor="w")
STEDpowerlabel.grid(column=4,row=setuprow+12)
STEDpowerentry = ttk.Entry(tab1,width=setupentrywidth,justify="right")
STEDpowerentry.grid(column=5,row=setuprow+12)

labelexc = ttk.Label(tab1,text ="Profile type:",justify="center",anchor="w")
labelexc.grid(column=4, row=setuprow+13)
STEDprofilevar = tkinter.IntVar(value=1)
profrad1 = ttk.Radiobutton(tab1,text='sinusoidal', variable = STEDprofilevar, value=1, width=15)
profrad2 = ttk.Radiobutton(tab1,text='vortex plate', variable = STEDprofilevar, value=2, width=15)
profrad1.grid(column=5, row=setuprow+13)
profrad2.grid(column=5, row=setuprow+14)

NAefflabel = ttk.Label(tab1,text="effective NA multiplier",justify="right",anchor="w")
NAefflabel.grid(column=4,row=setuprow+15)
NAeffentry = ttk.Entry(tab1,width=setupentrywidth,justify="right")
NAeffentry.grid(column=5, row=setuprow+15)
ttk.Separator(tab1, orient=HORIZONTAL).grid(row=setuprow+16, column=3, columnspan=3, sticky='W,E')

def savepars():
    tab_parent.select(tab2)

def restorepars():
    return 0

savepars_btn = ttk.Button(tab1, text="Save values", command=savepars,width=18)
savepars_btn.grid(column=1, row=9)
restorepars_btn = ttk.Button(tab1, text="Restore defaults", command=restorepars,width=18)
restorepars_btn.grid(column=1, row=10)

label2 = ttk.Label(tab2,text ="placeholder for definition and display of created fluorophore field")
label2.grid(column=1, row=1)
label3 = ttk.Label(tab3,text ="placeholder for definition of filters and spectral map")
label3.grid(column=1, row=1)
label4 = ttk.Label(tab4,text ="results go here")
label4.grid(column=1, row=1)

window.mainloop()