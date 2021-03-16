from __future__ import print_function
import numpy as np
import pickle
import os

class Colors:
  # use https://coolors.co/app to get more colorschemes
  def __init__(self,num, mode=None):
    self.num=num
    self.british= np.array(["#94bfac","#5b9291","#3b6879","#264d7e","#1f3057","#2a283d","#3a73a9","#173679","#1c5680","#2c3e75","#8cc5bb","#78adc2","#3f687d","#1f4b61","#5f88c1","#2458af","#135b75","#a7c6eb","#64a0aa","#4f81c5","#bbc9a5","#bcd890","#96bf65","#698b47","#757639","#4b5729","#33533b","#254432","#428b64","#4f5241","#44945e","#476a4c","#8fc693","#2e4c1e","#364a20","#87965a","#3b3629","#68ab77","#506b52","#7e8f6e","#6b6f5a","#5f5c4b","#4f5138","#feec04","#fef963","#fef96a","#9e7339","#4c4a3c","#7b6b4f","#fced96","#fdf07a","#e9bb43","#fdd906","#fcc808","#f6c870","#dbac50","#d4b97d","#ac7c42","#fde706","#cec093","#f4f0bd","#f5e7a1","#fef6bf","#dd7b00","#feeba8","#bba38a","#eedfa5","#e8c88f","#e6c18d","#cfb48a","#e4cf93","#b2a788","#f3d163","#93612b","#74542f","#5c422e","#402d21","#a86c29","#61361e","#a89177","#845b4d","#564b47","#464334","#753b1e","#c98a71","#a65341","#83422b","#774430","#f3b28b","#67403a","#693b3f","#745f46","#613339","#fbded6","#e8a1a2","#bd8f56","#793932","#8d5b41","#573320","#59493e","#bb3016","#dd3420","#c41c22","#d21e2b","#8b1a32","#471b21","#982d57","#ef841e","#dd3524","#fb9c06","#a83c19","#d04e09","#e45523","#f24816","#a0a9aa","#bec0b8","#9d9d7e","#7a838b","#a5ad98","#9aaa9f","#6b7477","#424c53","#6f7264","#525b55","#5f7682","#8e9b9c","#6c7377","#667563","#566164","#484837","#282b2f","#4e5355","#a9b7b9","#676f76","#7b93a3","#88918d","#909a92","#b6d3cc","#6e4a75","#c9a8ce"])
    self.crayolas = np.array(["#ED0A3F","#C32148","#FD0E35","#C62D42","#CC474B","#CC3336","#E12C2C","#D92121","#B94E48","#FF5349","#FE4C40","#FE6F5E","#B33B24","#CC553D","#E6735C","#FF9980","#E58E73","#FF7F49","#FF681F","#FF8833","#FFB97B","#ECB176","#E77200","#FFAE42","#F2BA49","#FBE7B2","#F2C649","#F8D568","#FCD667","#FED85D","#FBE870","#F1E788","#FFEB00","#B5B35C","#ECEBBD","#FAFA37","#FFFF99","#FFFF9F","#D9E650","#ACBF60","#AFE313","#BEE64B","#C5E17A","#5E8C31","#7BA05B","#9DE093","#63B76C","#4D8C57","#3AA655","#6CA67C","#5FA777","#93DFB8","#33CC99","#1AB385","#29AB87","#00CC99","#00755E","#8DD9CC","#01786F","#30BFBF","#00CCCC","#008080","#8FD8D8","#95E0E8","#6CDAE7","#2D383A","#76D7EA","#7ED4E6","#0095B7","#009DC4","#02A4D3","#47ABCC","#4997D0","#339ACC","#93CCEA","#2887C8","#00468C","#0066CC","#1560BD","#0066FF","#A9B2C3","#C3CDE6","#4570E6","#7A89B8","#4F69C6","#8D90A1","#8C90C8","#7070CC","#9999CC","#ACACE6","#766EC8","#6456B7","#3F26BF","#8B72BE","#652DC1","#6B3FA0","#8359A3","#8F47B3","#C9A0DC","#BF8FCC","#803790","#733380","#D6AEDD","#C154C1","#FC74FD","#732E6C","#E667CE","#E29CD2","#8E3179","#D96CBE","#EBB0D7","#C8509B","#BB3385","#D982B5","#A63A79","#A50B5E","#614051","#F653A6","#DA3287","#FF3399","#FBAED2","#FFB7D5","#FFA6C9","#F7468A","#E30B5C","#FDD7E4","#E62E6B","#DB5079","#FC80A5","#F091A9","#FF91A4","#A55353","#CA3435","#FEBAAD","#F7A38E","#E97451","#AF593E","#9E5B40","#87421F","#926F5B","#DEA681","#D27D46","#664228","#D99A6C","#EDC9AF","#FFCBA4","#805533","#FDD5B1","#EED9C4","#665233","#837050","#E6BC5C","#D9D6CF","#92926E","#E6BE8A","#C9C0BB","#DA8A67","#C88A65","#000000","#736A62","#8B8680","#C8C8CD"])

    if num <= 8 and mode ==None:
      step=int(8/num)
      print("here")
      self.palette=np.array(["#003f5c","#2f4b7c","#665191","#a05195","#d45087","#f95d6a","#ff7c43","#ffa600"])[::step]
    elif mode=="crayola":
      step=int(len(self.crayolas)/num)
      self.palette = self.crayolas[::step]
    elif mode =="british" or mode ==None:
       step=int(len(self.british)/num)
       self.palette = self.british[::step]



    # if num
def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]

def find_nearest_index(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx

def decompress(filename,mode='rb',delete=False):
  import gzip
  filename = os.path.abspath(filename)
  if filename[-3:]==".gz":
    filenameout = filename+".decompressed"
    if os.path.isfile(filenameout):
        print("!!! Skipping decompression: file already exists.")
        return filenameout
    else:
      with gzip.open(filename, mode) as fin:
        data = fin.read()
        with open(filenameout ,'wb') as fout:
          fout.write(data)
        return filenameout
  else:
    print ("!!! Decompression failed: only .gz files are supported.")

def clean(filename,extension = ".decompressed"):
  import os
  os.system("rm "+filename+extension)
def cast_to_array(a):
  b = tuple([np.array(i) for i in a])
  return b


def grep(command):
  import subprocess
  return subprocess.call(['/bin/grep', command], shell=True)