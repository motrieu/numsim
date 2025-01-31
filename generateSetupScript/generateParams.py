import tkinter as tk
from tkinter import filedialog as fd
import numpy as np
from PIL import Image
import os

global beginOfLineCoord
beginOfLineCoord = [0,0]
global beginOfLine
beginOfLine = True

def setBC(bcIndex):
    global currentBCIndex
    global currentColor
    currentBCIndex = bcIndex
    currentColor = from_rgb(colors[bcIndex])
    updateEntries(0)

def from_rgb(rgb):
    return "#%02x%02x%02x" % rgb   

def updateEntries(event):
    entryNxBegin.config(state='disabled')
    entryNyBegin.config(state='disabled')
    entryNxEnd.config(state='disabled')
    entryNyEnd.config(state='disabled')
    buttonDrawFromEntry.config(state='disabled')
    entryNxCenter.config(state='disabled')
    entryNyCenter.config(state='disabled')
    entryRadius.config(state='disabled')
    buttonDrawCircleFromEntry.config(state='disabled')
    entryUFunction.config(state='disabled')
    entryVFunction.config(state='disabled')
    entryPFunction.config(state='disabled')
    
    if (paintOptionsVar.get() in [paintOptions[1], paintOptions[2]]):
        entryNxBegin.config(state='normal')
        entryNyBegin.config(state='normal')
        entryNxEnd.config(state='normal')
        entryNyEnd.config(state='normal')
        buttonDrawFromEntry.config(state='normal')
        if (currentBCIndex == 2): #INFLOW
            entryUFunction.config(state='normal')
            entryVFunction.config(state='normal')
        elif (currentBCIndex == 4): #PRESSURE
            entryPFunction.config(state='normal')
    elif (paintOptionsVar.get() == paintOptions[3]): #paint circle
        entryNxCenter.config(state='normal')
        entryNyCenter.config(state='normal')
        entryRadius.config(state='normal')
        buttonDrawCircleFromEntry.config(state='normal')

# paint free-hand
def paint(event):
    x = (event.x // pixelSize) * pixelSize
    y = (event.y // pixelSize) * pixelSize
    if (paintOptionsVar.get() == paintOptions[0]):
        canvas.create_rectangle(x, y, x + pixelSize, y + pixelSize, fill=currentColor, outline="")

# paint line or rectangle by mouse action
def paintItems(event):
    global beginOfLineCoord
    global beginOfLine
    x = (event.x // pixelSize) * pixelSize
    y = (event.y // pixelSize) * pixelSize

    ### draw line ###
    # set starting point of line
    if beginOfLine and (paintOptionsVar.get() == paintOptions[1]):
        canvas.create_rectangle(x, y, x + pixelSize, y + pixelSize, fill=currentColor, outline="")
        beginOfLineCoord[0] = x
        beginOfLineCoord[1] = y
        beginOfLine = False
        nxBeginVar.set(event.x // pixelSize)
        nyBeginVar.set(numCellsY-1 - event.y // pixelSize)
    # set end pint of line and draw line
    elif not beginOfLine and (paintOptionsVar.get() == paintOptions[1]) and ((beginOfLineCoord[0] == x) ^ (beginOfLineCoord[1] == y)):
        canvas.create_rectangle(min(beginOfLineCoord[0],x), min(beginOfLineCoord[1],y), max(beginOfLineCoord[0],x) + pixelSize, max(beginOfLineCoord[1],y) + pixelSize, fill=currentColor, outline="")
        beginOfLine = True
        nxEndVar.set(event.x // pixelSize)
        nyEndVar.set(numCellsY-1 - event.y // pixelSize)
    
    ### draw rectangle ###
    # set first corner of rectangle
    elif beginOfLine and (paintOptionsVar.get() == paintOptions[2]):
        canvas.create_rectangle(x, y, x + pixelSize, y + pixelSize, fill=currentColor, outline="")
        beginOfLineCoord[0] = x
        beginOfLineCoord[1] = y
        beginOfLine = False
        nxBeginVar.set(event.x // pixelSize)
        nyBeginVar.set(numCellsY-1 - event.y // pixelSize)
    # set second corner of rectangle and draw rectangle
    elif not beginOfLine and (paintOptionsVar.get() == paintOptions[2]):
        canvas.create_rectangle(min(beginOfLineCoord[0],x), min(beginOfLineCoord[1],y), max(beginOfLineCoord[0],x) + pixelSize, max(beginOfLineCoord[1],y) + pixelSize, fill=currentColor, outline="")
        beginOfLine = True
        nxEndVar.set(event.x // pixelSize)
        nyEndVar.set(numCellsY-1 - event.y // pixelSize)

# draw line or rectangle based on values provided in the input entries
def drawFromEntry():
    nxBegin = int(nxBeginVar.get())
    nxEnd = int(nxEndVar.get())
    nyBegin = numCellsY-1 - int(nyBeginVar.get())
    nyEnd = numCellsY-1 - int(nyEndVar.get())
    
    # draw line
    if (paintOptionsVar.get() == paintOptions[1]) and ((nxBegin == nxEnd) ^ (nyBegin == nyEnd)):
        # draw constant line without any function applied (for all BC)
        if ((currentBCIndex not in [2,4]) or ((currentBCIndex==2) and (uFunctionVar.get()=="") and (vFunctionVar.get()=="")) or ((currentBCIndex==4) and (pFunctionVar.get()==""))):
            canvas.create_rectangle(min(nxBegin,nxEnd)*pixelSize, min(nyBegin,nyEnd)*pixelSize, (max(nxBegin,nxEnd)+1)*pixelSize, (max(nyBegin,nyEnd)+1)*pixelSize, fill=currentColor, outline="")
        # draw INFLOW condition based on function input
        elif currentBCIndex == 2:
            nxyRange = max(abs(nxEnd-nxBegin), abs(nyEnd-nyBegin))
            for i in range(min(nxBegin,nxEnd), max(nxBegin,nxEnd)+1):
                for j in range(min(nyBegin,nyEnd), max(nyBegin,nyEnd)+1):
                    x = (i-min(nxBegin,nxEnd))/nxyRange
                    y = (j-min(nyBegin,nyEnd))/nxyRange
                    uString = uFunctionVar.get()
                    u = 0
                    if uString != "":
                        u = eval(uString)
                    vString = vFunctionVar.get()
                    v = 0
                    if vString != "":
                        v = eval(vString)
                    colorCode = (int(200*u+55), 0, int(200*v+55))
                    canvas.create_rectangle(i*pixelSize, j*pixelSize, (i+1)*pixelSize, (j+1)*pixelSize, fill=from_rgb(colorCode), outline="")
        # draw PRESSURE BC based on function input
        elif currentBCIndex == 4:
            nxyRange = max(abs(nxEnd-nxBegin), abs(nyEnd-nyBegin))
            for i in range(min(nxBegin,nxEnd), max(nxBegin,nxEnd)+1):
                for j in range(min(nyBegin,nyEnd), max(nyBegin,nyEnd)+1):
                    x = (i-min(nxBegin,nxEnd))/nxyRange
                    y = (j-min(nyBegin,nyEnd))/nxyRange
                    pString = pFunctionVar.get()
                    p = 0
                    if pString != "":
                        p = eval(pString)
                    colorCode = (0, int(200*p+55), 0)
                    canvas.create_rectangle(i*pixelSize, j*pixelSize, (i+1)*pixelSize, (j+1)*pixelSize, fill=from_rgb(colorCode), outline="")
                
    # draw rectangle from input entries for all BCs
    elif (paintOptionsVar.get() == paintOptions[2]):
        canvas.create_rectangle(min(nxBegin,nxEnd)*pixelSize, min(nyBegin,nyEnd)*pixelSize, (max(nxBegin,nxEnd)+1)*pixelSize, (max(nyBegin,nyEnd)+1)*pixelSize, fill=currentColor, outline="")

# draw circle based on values provided in the input entries
def drawCircleFromEntry():
    nxCenter = int(nxCenterVar.get())
    nyCenter = numCellsY-1 - int(nyCenterVar.get())
    nradius = int(radiusVar.get())
    for i in range(nxCenter-nradius, nxCenter+nradius):
        for j in range(nyCenter-nradius, nyCenter+nradius):
            if ((nxCenter-i)**2+(nyCenter-j)**2)**0.5 < nradius:
                canvas.create_rectangle(i*pixelSize, j*pixelSize, (i+1)*pixelSize, (j+1)*pixelSize, fill=currentColor, outline="")

# generate setup image and write to selected params file
def generateImageAndParams():
    canvas.update()
    canvas.postscript(file="temp.ps", colormode='color')
    psimage = Image.open("temp.ps")
    psimage.save("setup.png")
    psimage.close()
    os.remove("temp.ps")
    generateParams("setup.png")

def generateParams(imageFileName):
    img = Image.open(imageFileName)
    npimg = np.asarray(img)
    pixelSize = np.size(npimg,1) / numCellsX
    params = open(paramsFileName, "a")
    params.write("\n\n#SETUP\n")
    for j in range(numCellsY):
        for i in range(numCellsX):
            color = npimg[int((j+0.5)*pixelSize), int((i+0.5)*pixelSize)]
            if color[0]==255 and color[1]==255 and color[2]==255:
                params.write("0 ") # Fluid cell
            elif color[0]==0 and color[1]==0 and color[2]==0:
                params.write("1 ") # NOSLIP
            elif color[0]==255 and color[1]==255 and color[2]==0:
                params.write("2 ") # SLIP
            elif color[0]>=55 and color[1]==0 and color[2]>=55:
                params.write("3 ") # INFLOW
            elif color[0]==0 and color[1]==255 and color[2]==255:
                params.write("4 ") # OUTFLOW
            elif color[0]==0 and color[1]>=55 and color[2]==0:
                params.write("5 ") # PRESSURE
            else:
                print("ERROR in setup")
        params.write("\n")
    
    params.write("\n#UIN\n")
    for j in range(numCellsY):
        for i in range(numCellsX):
            color = npimg[int((j+0.5)*pixelSize), int((i+0.5)*pixelSize)]
            if color[0]==255 and color[1]==255 and color[2]==255:
                params.write("0 ") # Fluid cell
            elif color[0]==0 and color[1]==0 and color[2]==0:
                params.write("0 ") # NOSLIP
            elif color[0]==255 and color[1]==255 and color[2]==0:
                params.write("0 ") # SLIP
            elif color[0]>=55 and color[1]==0 and color[2]>=55:
                u = (color[0]-55)/200*(uMax-uMin) + uMin
                params.write(str(u) + " ") # INFLOW
            elif color[0]==0 and color[1]==255 and color[2]==255:
                params.write("0 ") # OUTFLOW
            elif color[0]==0 and color[1]>=55 and color[2]==0:
                params.write("0 ") # PRESSURE
            else:
                print("ERROR in uin")
        params.write("\n")
    
    params.write("\n#VIN\n")
    for j in range(numCellsY):
        for i in range(numCellsX):
            color = npimg[int((j+0.5)*pixelSize), int((i+0.5)*pixelSize)]
            if color[0]==255 and color[1]==255 and color[2]==255:
                params.write("0 ") # Fluid cell
            elif color[0]==0 and color[1]==0 and color[2]==0:
                params.write("0 ") # NOSLIP
            elif color[0]==255 and color[1]==255 and color[2]==0:
                params.write("0 ") # SLIP
            elif color[0]>=55 and color[1]==0 and color[2]>=55:
                v = (color[2]-55)/200*(vMax-vMin) + vMin
                params.write(str(v) + " ") # INFLOW
            elif color[0]==0 and color[1]==255 and color[2]==255:
                params.write("0 ") # OUTFLOW
            elif color[0]==0 and color[1]>=55 and color[2]==0:
                params.write("0 ") # PRESSURE
            else:
                print("ERROR in vin")
        params.write("\n")
    
    params.write("\n#PRB\n")
    for j in range(numCellsY):
        for i in range(numCellsX):
            color = npimg[int((j+0.5)*pixelSize), int((i+0.5)*pixelSize)]
            if color[0]==255 and color[1]==255 and color[2]==255:
                params.write("0 ") # Fluid cell
            elif color[0]==0 and color[1]==0 and color[2]==0:
                params.write("0 ") # NOSLIP
            elif color[0]==255 and color[1]==255 and color[2]==0:
                params.write("0 ") # SLIP
            elif color[0]>=55 and color[1]==0 and color[2]>=55:
                params.write("0 ") # INFLOW
            elif color[0]==0 and color[1]==255 and color[2]==255:
                params.write("0 ") # OUTFLOW
            elif color[0]==0 and color[1]>=55 and color[2]==0:
                p = (color[1]-55)/200*(pMax-pMin) + pMin
                params.write(str(p) + " ") # PRESSURE
            else:
                print("ERROR in prb")
        params.write("\n")

    params.close()

def selectParamsFile():
    global paramsFileName
    paramsFileName = fd.askopenfilename()

def generateFromImage():
    global numCellsX
    numCellsX = int(numCellsXVar.get()) + 2
    global numCellsY
    numCellsY = int(numCellsYVar.get()) + 2

    global uMin
    uMin = float(uMinVar.get())
    global uMax
    uMax = float(uMaxVar.get())
    global vMin
    vMin = float(vMinVar.get())
    global vMax
    vMax = float(vMaxVar.get())
    global pMin
    pMin = float(pMinVar.get())
    global pMax
    pMax = float(pMaxVar.get())

    imageFileName = fd.askopenfilename()

    img = Image.open(imageFileName)
    npimg = np.asarray(img)
    imagePixelSize = np.size(npimg,1) / (numCellsX+1)

    buildDrawEnvironment()
    for j in range(numCellsY):
        for i in range(numCellsX):
            color = npimg[int((j+0.5)*imagePixelSize), int((i+0.5)*imagePixelSize)]
            canvas.create_rectangle(i*pixelSize, j*pixelSize, (i+1)*pixelSize, (j+1)*pixelSize, fill=from_rgb((color[0],color[1],color[2])), outline="")


    #generateParams(imageFileName)

def buildDrawEnvironment():
    global numCellsX
    numCellsX = int(numCellsXVar.get()) + 2
    global numCellsY
    numCellsY = int(numCellsYVar.get()) + 2

    global uMin
    uMin = float(uMinVar.get())
    global uMax
    uMax = float(uMaxVar.get())
    global vMin
    vMin = float(vMinVar.get())
    global vMax
    vMax = float(vMaxVar.get())
    global pMin
    pMin = float(pMinVar.get())
    global pMax
    pMax = float(pMaxVar.get())
    maxWindowWidth = 1300
    maxWindowHeight = 400
    global pixelSize
    pixelSize = 1
    while (numCellsX*pixelSize < maxWindowWidth) and (numCellsY*pixelSize < maxWindowHeight):
        pixelSize += 1

    canvas.config(width=numCellsX*pixelSize-4, height=numCellsY*pixelSize-4)
    canvas.pack()
    subbox.pack()

    column1.pack(side=tk.LEFT)
    column2.pack(side=tk.LEFT)
    column3.pack(side=tk.LEFT)
    column4.pack(side=tk.LEFT)

    for i in range(len(rbButtons)):
        rbFrames[i].pack(anchor="w")
        rbButtons[i].pack(side=tk.LEFT)
        rbLabels[i].pack(side=tk.LEFT)
    rbFrames[-1].pack()
    generateButton.pack()

    dropdown.pack(side=tk.LEFT)

    nxBeginFrame.pack()
    labelNxBegin.pack(side=tk.LEFT)
    entryNxBegin.pack(side=tk.LEFT)
    nyBeginFrame.pack()
    labelNyBegin.pack(side=tk.LEFT)
    entryNyBegin.pack(side=tk.LEFT)
    nxEndFrame.pack()
    labelNxEnd.pack(side=tk.LEFT)
    entryNxEnd.pack(side=tk.LEFT)
    nyEndFrame.pack()
    labelNyEnd.pack(side=tk.LEFT)
    entryNyEnd.pack(side=tk.LEFT)

    uFunctionFrame.pack()
    labelUFunction.pack(side=tk.LEFT)
    entryUFunction.pack(side=tk.LEFT)
    vFunctionFrame.pack()
    labelVFunction.pack(side=tk.LEFT)
    entryVFunction.pack(side=tk.LEFT)
    pFunctionFrame.pack()
    labelPFunction.pack(side=tk.LEFT)
    entryPFunction.pack(side=tk.LEFT)
    buttonDrawFromEntry.pack()

    nxCenterFrame.pack()
    labelNxCenter.pack(side=tk.LEFT)
    entryNxCenter.pack(side=tk.LEFT)
    nyCenterFrame.pack()
    labelNyCenter.pack(side=tk.LEFT)
    entryNyCenter.pack(side=tk.LEFT)
    radiusFrame.pack()
    labelRadius.pack(side=tk.LEFT)
    entryRadius.pack(side=tk.LEFT)
    buttonDrawCircleFromEntry.pack()

    buildDrawButton.config(text="Update")


currentColor = "black"
currentBCIndex = 0

window = tk.Tk()
window.title("NSE BCs")
#window.state('zoomed')

domainSizeFrame = tk.Frame(window)
domainSizeFrame.pack()
numCellsXLabel = tk.Label(domainSizeFrame, text="Number of cells in x: ")
numCellsXLabel.pack(side=tk.LEFT)
numCellsXVar = tk.StringVar()
numCellsXEntry = tk.Entry(domainSizeFrame, textvariable=numCellsXVar)
numCellsXEntry.pack(side=tk.LEFT)
numCellsYLabel = tk.Label(domainSizeFrame, text=", y: ")
numCellsYLabel.pack(side=tk.LEFT)
numCellsYVar = tk.StringVar()
numCellsYEntry = tk.Entry(domainSizeFrame, textvariable=numCellsYVar)
numCellsYEntry.pack(side=tk.LEFT)

uMinMaxFrame = tk.Frame(window)
uMinMaxFrame.pack()
uMinLabel = tk.Label(uMinMaxFrame, text="u_min: ")
uMinLabel.pack(side=tk.LEFT)
uMinVar = tk.StringVar()
uMinVar.set("0")
uMinEntry = tk.Entry(uMinMaxFrame, textvariable=uMinVar)
uMinEntry.pack(side=tk.LEFT)
uMaxLabel = tk.Label(uMinMaxFrame, text=", u_max: ")
uMaxLabel.pack(side=tk.LEFT)
uMaxVar = tk.StringVar()
uMaxVar.set("0")
uMaxEntry = tk.Entry(uMinMaxFrame, textvariable=uMaxVar)
uMaxEntry.pack(side=tk.LEFT)

vMinMaxFrame = tk.Frame(window)
vMinMaxFrame.pack()
vMinLabel = tk.Label(vMinMaxFrame, text="v_min: ")
vMinLabel.pack(side=tk.LEFT)
vMinVar = tk.StringVar()
vMinVar.set("0")
vMinEntry = tk.Entry(vMinMaxFrame, textvariable=vMinVar)
vMinEntry.pack(side=tk.LEFT)
vMaxLabel = tk.Label(vMinMaxFrame, text=", v_max: ")
vMaxLabel.pack(side=tk.LEFT)
vMaxVar = tk.StringVar()
vMaxVar.set("0")
vMaxEntry = tk.Entry(vMinMaxFrame, textvariable=vMaxVar)
vMaxEntry.pack(side=tk.LEFT)

pMinMaxFrame = tk.Frame(window)
pMinMaxFrame.pack()
pMinLabel = tk.Label(pMinMaxFrame, text="p_min: ")
pMinLabel.pack(side=tk.LEFT)
pMinVar = tk.StringVar()
pMinVar.set("0")
pMinEntry = tk.Entry(pMinMaxFrame, textvariable=pMinVar)
pMinEntry.pack(side=tk.LEFT)
pMaxLabel = tk.Label(pMinMaxFrame, text=", p_max: ")
pMaxLabel.pack(side=tk.LEFT)
pMaxVar = tk.StringVar()
pMaxVar.set("0")
pMaxEntry = tk.Entry(pMinMaxFrame, textvariable=pMaxVar)
pMaxEntry.pack(side=tk.LEFT)

selectParamsButton = tk.Button(window, text="Select Parameter File", command=selectParamsFile)
selectParamsButton.pack()

generateFrame = tk.Frame(window)
generateFrame.pack()
generateFromImageButton = tk.Button(generateFrame, text="Generate from Image", command=generateFromImage)
generateFromImageButton.pack(side=tk.LEFT)
buildDrawButton = tk.Button(generateFrame, text="Draw Setup", command=buildDrawEnvironment)
buildDrawButton.pack(side=tk.LEFT)


canvas = tk.Canvas(window, bg="white")

canvas.bind("<B1-Motion>", paint)
canvas.bind("<ButtonRelease-1>", paintItems)


subbox = tk.Frame(window)

column1 = tk.Frame(subbox)
column2 = tk.Frame(subbox)
column3 = tk.Frame(subbox)
column4 = tk.Frame(subbox)

colors = [(0,0,0), (255,255,0), (255,0,255), (0,255,255), (0,255,0), (255,255,255)]
bcs = ["NOSLIP", "SLIP", "INFLOW", "OUTFLOW", "PRESSURE", "erase"]

rbFrames = []
rbButtons = []
rbLabels = []
for i in range(6):
    rbFrames.append(tk.Frame(column1))
    rbButtons.append(tk.Button(rbFrames[i], bg=from_rgb(colors[i]), width=2, height=1, command=lambda bcIndex=i: setBC(bcIndex)))
    rbLabels.append(tk.Label(rbFrames[i], text=bcs[i], anchor="w"))

rbFrames.append(tk.Frame(column1))
generateButton = tk.Button(rbFrames[-1], text="Generate", command=generateImageAndParams)

paintOptions = ["paint", "draw line", "draw rectangle", "draw circle"]
paintOptionsVar = tk.StringVar()
paintOptionsVar.set(paintOptions[0])
dropdown = tk.OptionMenu(column2, paintOptionsVar, *paintOptions, command=updateEntries)

nxBeginFrame = tk.Frame(column3)
labelNxBegin = tk.Label(nxBeginFrame, text="n_x_begin: ")
nxBeginVar = tk.StringVar()
entryNxBegin = tk.Entry(nxBeginFrame, textvariable=nxBeginVar)
entryNxBegin.config(state='disabled')

nyBeginFrame = tk.Frame(column3)
labelNyBegin = tk.Label(nyBeginFrame, text="n_y_begin: ")
nyBeginVar = tk.StringVar()
entryNyBegin = tk.Entry(nyBeginFrame, textvariable=nyBeginVar)
entryNyBegin.config(state='disabled')

nxEndFrame = tk.Frame(column3)
labelNxEnd = tk.Label(nxEndFrame, text="n_x_end: ")
nxEndVar = tk.StringVar()
entryNxEnd = tk.Entry(nxEndFrame, textvariable=nxEndVar)
entryNxEnd.config(state='disabled')

nyEndFrame = tk.Frame(column3)
labelNyEnd = tk.Label(nyEndFrame, text="n_y_end: ")
nyEndVar = tk.StringVar()
entryNyEnd = tk.Entry(nyEndFrame, textvariable=nyEndVar)
entryNyEnd.config(state='disabled')

uFunctionFrame = tk.Frame(column3)
labelUFunction = tk.Label(uFunctionFrame, text="u(x,y): ")
uFunctionVar = tk.StringVar()
entryUFunction = tk.Entry(uFunctionFrame, textvariable=uFunctionVar)
entryUFunction.config(state='disabled')

vFunctionFrame = tk.Frame(column3)
labelVFunction = tk.Label(vFunctionFrame, text="v(x,y): ")
vFunctionVar = tk.StringVar()
entryVFunction = tk.Entry(vFunctionFrame, textvariable=vFunctionVar)
entryVFunction.config(state='disabled')

pFunctionFrame = tk.Frame(column3)
labelPFunction = tk.Label(pFunctionFrame, text="p(x,y): ")
pFunctionVar = tk.StringVar()
entryPFunction = tk.Entry(pFunctionFrame, textvariable=pFunctionVar)
entryPFunction.config(state='disabled')

buttonDrawFromEntry = tk.Button(column3, text="Apply", command=drawFromEntry)
buttonDrawFromEntry.config(state='disabled')


nxCenterFrame = tk.Frame(column4)
labelNxCenter = tk.Label(nxCenterFrame, text="n_x_center: ")
nxCenterVar = tk.StringVar()
entryNxCenter = tk.Entry(nxCenterFrame, textvariable=nxCenterVar)
entryNxCenter.config(state='disabled')

nyCenterFrame = tk.Frame(column4)
labelNyCenter = tk.Label(nyCenterFrame, text="n_y_center: ")
nyCenterVar = tk.StringVar()
entryNyCenter = tk.Entry(nyCenterFrame, textvariable=nyCenterVar)
entryNyCenter.config(state='disabled')

radiusFrame = tk.Frame(column4)
labelRadius = tk.Label(radiusFrame, text="radius: ")
radiusVar = tk.StringVar()
entryRadius = tk.Entry(radiusFrame, textvariable=radiusVar)
entryRadius.config(state='disabled')

buttonDrawCircleFromEntry = tk.Button(column4, text="Apply", command=drawCircleFromEntry)
buttonDrawCircleFromEntry.config(state='disabled')


window.mainloop()