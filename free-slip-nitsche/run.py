import shutil
import glob, os
from nitsche import *

d_t = 0.00001 # change in nitsche.py
t = 0.0
T = 0.05

shutil.copy2('../pipe/pipe_coarse.xml', './ale_mesh.xml')

counter = 1
while t < T + DOLFIN_EPS:
    
    ale_stokes(counter)       
    t = t + d_t
    print("time so far:", t)   
        
    folder = './output'    
    for filename in glob.iglob(os.path.join(folder, '*000000.vtu')):
        os.rename(filename, filename[:-10] + '.vtu')
    
    #for filename in os.listdir("./output"):
    #    if filename.endswith("000000.vtu"):
    #        os.rename(filename, filename[:-10] + '.vtu')
            
    #print(os.listdir('./output'))
        
    counter += 1