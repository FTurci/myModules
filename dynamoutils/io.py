import os
import numpy as np
import os
import xml.etree.ElementTree as ET
import string

def loadXMLFile(filename):
    #Check if the file is compressed or not, and 
    if (os.path.splitext(filename)[1][1:].strip() == "bz2"):
        import bz2
        f = bz2.BZ2File(filename)
        doc = ET.parse(f)
        f.close()
        return doc
    else:
        return ET.parse(filename)

def read_box(filename):
  doc = loadXMLFile(filename)
  simulation = doc.find(".//Simulation")
  lx= simulation.find("SimulationSize").get("x")
  ly= simulation.find("SimulationSize").get("y")
  lz= simulation.find("SimulationSize").get("z")
  return np.array([lx,ly,lz]).astype(float)
def read_pressure(filename):
  doc = loadXMLFile(filename)
  p=doc.find(".//Misc")
  return float(p.find("Pressure").get("Avg"))

def read_temperature(filename):
  doc = loadXMLFile(filename)
  p=doc.find(".//Misc")
  return float(p.find("Temperature").get("Current") )

def read_uconfig(filename):
  doc = loadXMLFile(filename)
  p=doc.find(".//Misc")
  return float(p.find("UConfigurational").get("Current") )

def read_time(filename):
  doc = loadXMLFile(filename)
  p=doc.find(".//Misc")
  return float(p.find("Duration").get("Time") )


def read_quantity(filename,name):
  doc = loadXMLFile(filename)
  p=doc.find(".//Misc")
  return float(p.find(name).get("val") )

def conf_to_xyz(filename):
  XMLDoc=loadXMLFile(filename)

  #We can create a list of all particle tags using an xpath expression
  #(xpath expressions always return lists)
  PtTags = XMLDoc.findall(".//Pt")
  #print the x, y, and z positions of each particle
  with open(filename+".xyz", 'w') as fw:
    fw.write(f"{len(PtTags)}\nAtoms\n")
    for PtElement in PtTags:
        PosTag = PtElement.find("P")
        x,y,z=PosTag.get("x"),PosTag.get("y"),PosTag.get("z")
        fw.write(f"A {x} {y} {z}\n")

def conf_to_atom(filename,timestep):
  XMLDoc=loadXMLFile(filename)
  types=list(string.ascii_uppercase)

  SizeTag = XMLDoc.find(".//SimulationSize")
  sizes = float(SizeTag.get("x")),float(SizeTag.get("y")),float(SizeTag.get("z"))

  box =np.array([[-l*0.5, l*0.5] for l in sizes])
  xlo,xhi = box[0,0],box[0,1]
  ylo,yhi = box[1,0],box[1,1]
  zlo,zhi = box[2,0],box[2,1]

  PtTags = XMLDoc.findall(".//Pt")
  N=len(PtTags)
  #print the x, y, and z positions of each particle
  with open(filename+".atom", 'w') as fw:
    fw.write(f"""ITEM: TIMESTEP
{timestep}
ITEM: NUMBER OF ATOMS
{N}
ITEM: BOX BOUNDS pp pp pp
{xlo} {xhi}
{ylo} {yhi}
{zlo} {zhi}\n""")

    # check if there are diameters
    p0=PtTags[0]

    if p0.get("D") is not None:
      print(":: Found a diameter. Saving radius information.")
      radius = True
      fw.write("ITEM: ATOMS type x y z radius\n")
    else:
      print(":: No diameters. Assuming single species.",flush=True)
      radius = False
      fw.write("ITEM: ATOMS type x y z\n")
  
    for PtElement in PtTags:
      PosTag = PtElement.find("P")
      x,y,z=PosTag.get("x"),PosTag.get("y"),PosTag.get("z")
      
      if radius:
        diam = float(PtElement.get("D"))
        fw.write(f"A {x} {y} {z} {diam/2}\n")
      else:
        fw.write(f"A {x} {y} {z}\n")