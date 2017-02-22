import shutil
import glob, os
from convdiff import *

i = 1
I = 50

c = 1
while i <= I:
    
    if i <= 5:
        
        s = 0.000024
    else:
        
        s = 0.00024        
        
    iter(c, s)    
    i +=  1
    c += 1
    