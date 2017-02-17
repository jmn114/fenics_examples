import shutil
import glob, os

t = 0
T = 3

while t < T:
    
    os.system('python conv-diff.py')
    print("TIME SO FAR:", t)           
    t +=  1
    