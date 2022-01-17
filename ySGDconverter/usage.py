#Download the python script 'ySGDconverter.py' and import
#using your notebook or shell

#Saccharomyces Genome Database (SGD) ID to Gene name converter and vice versa.
### Option 1: Convert SGDID to Gene name
### Option 2: Convert Gene name to SGDID
#version 1.0

import ySGDconverter

# Option 1
genes = ['YAL069W', 'YAL068C', 'YAL067C', 'YAL062W', 'YAL061W','YAR031W']

clist = ySGDconverter.convert(genes,1)
print(clist)

# Option 2
genes = ['YAL069W', 'PAU8', 'SEO1','GDH3','FUN51','FUN11','ERP1','YAR044W']

clist = ySGDconverter.convert(genes,2)
print(clist)
