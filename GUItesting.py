import tkinter
from tkinter import HORIZONTAL, VERTICAL, ttk
from PIL import Image, ImageTk
from os import path
import numpy

fluorlist = str(numpy.genfromtxt("dye_spectra/properties.dat", dtype="U",usecols=(1),comments="#",skip_header=1,max_rows=1000))
#fluorlist = ["placeholder1","placeholder2","placeholder3"]

windowwidth = 820
windowheight = 360

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

NAefflabel = ttk.Label(tab1,text="Effective NA multiplier",justify="right",anchor="w")
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



ttk.Separator(tab2, orient=VERTICAL).grid(column=0, row=0, rowspan=19, sticky='N,S')
ttk.Separator(tab2, orient=VERTICAL).grid(column=3, row=0, rowspan=19, sticky='N,S')

ttk.Separator(tab2, orient=HORIZONTAL).grid(row=0, columnspan=3, sticky='W,E')
fluorophoresetup = ttk.Label(tab2,text ="Fluorophore definition",font="Brown-gothic 11")
fluorophoresetup.grid(column=1, row=1)
ttk.Separator(tab2, orient=HORIZONTAL).grid(row=2, columnspan=3, sticky='W,E')


fluor1density = ttk.Label(tab2,text="Density [per m^2]",justify="left",anchor="w")
fluor1density.grid(column=1,row=10)
fluor2density = ttk.Label(tab2,text="Density [per m^2]",justify="left",anchor="w")
fluor2density.grid(column=1,row=14)
fluor3density = ttk.Label(tab2,text="Density [per m^2]",justify="left",anchor="w")
fluor3density.grid(column=1,row=18)
fieldtypevar = tkinter.IntVar(value=1)

def set2D3D():
    ft = fieldtypevar.get()
    if (ft == 1):
        fluor1density.configure(text="Density [per m^2]")
        fluor2density.configure(text="Density [per m^2]")
        fluor3density.configure(text="Density [per m^2]")
    elif (ft == 2):
        fluor1density.configure(text="Density [per m^3]")
        fluor2density.configure(text="Density [per m^3]")
        fluor3density.configure(text="Density [per m^3]")

typelabel = ttk.Label(tab2,text="Field type:",justify="left",anchor="w")
typelabel.grid(column=1,row=3)
typerad1 = ttk.Radiobutton(tab2,text='2D', variable = fieldtypevar, value=1, width=10, command=set2D3D)
typerad2 = ttk.Radiobutton(tab2,text='3D', variable = fieldtypevar, value=2, width=10, command=set2D3D)
typerad1.grid(column=1, row=4)
typerad2.grid(column=2, row=4)

axdimlabel = ttk.Label(tab2,text="Axial dimension [w0]",justify="left",anchor="w")
axdimlabel.grid(column=1,row=5)
axdimentry = ttk.Entry(tab2,width=setupentrywidth,justify="right")
axdimentry.grid(column=2,row=5)
latdimlabel = ttk.Label(tab2,text="Lateral dimension [w0]",justify="left",anchor="w")
latdimlabel.grid(column=1,row=6)
latdimentry = ttk.Entry(tab2,width=setupentrywidth,justify="right")
latdimentry.grid(column=2,row=6)
ttk.Separator(tab2, orient=HORIZONTAL).grid(row=7, columnspan=3, sticky='W,E')

fluor1label = ttk.Label(tab2,text ="Type 1",font="Brown-gothic 11")
fluor1label.grid(column=1, row=8)
fluor1type = ttk.Combobox(tab2, values = fluorlist,width=22)
fluor1type.set("pick fluorophore type")
fluor1type.grid(column=1,row=9)

fluor1densityentry = ttk.Entry(tab2,width=setupentrywidth,justify="right")
fluor1densityentry.grid(column=2,row=10)

ttk.Separator(tab2, orient=HORIZONTAL).grid(row=11, columnspan=3, sticky='W,E')
fluor2label = ttk.Label(tab2,text ="Type 2",font="Brown-gothic 11")
fluor2label.grid(column=1, row=12)
fluor2type = ttk.Combobox(tab2, values = fluorlist,width=22)
fluor2type.set("pick fluorophore type")
fluor2type.grid(column=1,row=13)

fluor2densityentry = ttk.Entry(tab2,width=setupentrywidth,justify="right")
fluor2densityentry.grid(column=2,row=14)

ttk.Separator(tab2, orient=HORIZONTAL).grid(row=15, columnspan=3, sticky='W,E')
fluor3label = ttk.Label(tab2,text ="Type 3",font="Brown-gothic 11")
fluor3label.grid(column=1, row=16)
fluor3type = ttk.Combobox(tab2, values = fluorlist,width=22)
fluor3type.set("pick fluorophore type")
fluor3type.grid(column=1,row=17)

fluor3densityentry = ttk.Entry(tab2,width=setupentrywidth,justify="right")
fluor3densityentry.grid(column=2,row=18)

ttk.Separator(tab2, orient=HORIZONTAL).grid(row=19, columnspan=3, sticky='W,E')

def populate():
    return 0

populate_btn = ttk.Button(tab2, text="Generate fluorophore distribution", command=populate,width=35)
populate_btn.grid(column=1, row=20)

img_placeholder_native = Image.open("GUIelements/phoremap_placeholder.png")
img_placeholder_native.thumbnail((450,300))
img_placeholder_map = ImageTk.PhotoImage(img_placeholder_native)
imgmaplabel = ttk.Label(tab2,image=img_placeholder_map)
imgmaplabel.grid(column=4, row=0,rowspan=21)

label3 = ttk.Label(tab3,text ="placeholder for definition of filters and spectral map")
label3.grid(column=1, row=1)
label4 = ttk.Label(tab4,text ="results go here")
label4.grid(column=1, row=1)

window.mainloop()