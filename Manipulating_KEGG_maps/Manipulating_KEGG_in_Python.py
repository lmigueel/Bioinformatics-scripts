#!/usr/bin/env python
# coding: utf-8

# # KEGG access
# 
# Let´s access the KEGG information through two approaches: Biopython and Bioservices
# 
# It´s usefull if you´re a python user. 
# 
# I will teach you how to fetch genes, pathways, reactions and other information from KEGG pathways in python.
#     

# ## 1. Biopython
# 
# Biopython is a set of freely available tools for biological computation written in Python by an international team of developers.
# 
# https://biopython.org/

# In[1]:


#!pip install biopython
#!pip install reportlab


# In[2]:


#Biopython modules to interact with KEGG
from Bio import SeqIO
from Bio.KEGG import REST
from Bio.KEGG.KGML import KGML_parser
from Bio.Graphics.KGML_vis import KGMLCanvas

# Import Pandas for dataframes
import pandas as pd
import io

get_ipython().run_line_magic('matplotlib', 'inline')
from IPython.display import Image


# ### Use
# 
# The available functions are:
# 
#     kegg_conv() - convert identifiers from KEGG to those for other databases
#     kegg_find() - find KEGG entries with matching query data
#     kegg_get() - retrieve data for a specific entry from KEGG
#     kegg_info() - get information about a KEGG database
#     kegg_link() - find entries in KEGG using a database cross-reference
#     kegg_list() - list entries in a a database
# 
# The generic form of using these functions to get information from KEGG and place the output in the variable myvar is:
# 
# myvar = REST.\<function>(\<query>, \<arg1>, \<arg2>, `...`).read()
# 

# In[3]:


# Perform the query
result = REST.kegg_info("kegg").read()

# Print the result
print(result)


# In[4]:


# Print information about the PATHWAY database
result = REST.kegg_info("pathway").read()
print(result)


# In[5]:


# Print information about S.cerevisiae
result = REST.kegg_info("sce").read()
print(result)


# In[6]:


# Get all entries in the PATHWAY database as a dataframe
result = REST.kegg_list("pathway").read()
pd.read_table(io.StringIO(result), header=None)


# In[7]:


# Get all entries in the PATHWAY database for S. cerevisiae as a dataframe
result = REST.kegg_list("pathway", "sce").read()
pd.read_table(io.StringIO(result), header=None)


# In[8]:


# Get all genes from S. cerevisiae as a dataframe
result = REST.kegg_list("sce").read()
pd.read_table(io.StringIO(result), header=None)


# In[9]:


# Find a specific entry with a precise search term
result = REST.kegg_find("genes", "sce:YAL067C").read()
pd.read_table(io.StringIO(result), header=None)


# In[10]:


# Find all ethanol genes
# You may search by a specific term. Will return for all organisms
result = REST.kegg_find("genes", "ethanol").read()
pd.read_table(io.StringIO(result), header=None)


# In[11]:


# Find all ethanol genes
# You may search by an especific term. Will return for all S. cerevisiae
result = REST.kegg_find("sce", "ethanol").read()
pd.read_table(io.StringIO(result), header=None)


# In[12]:


# Get the entry information for cpd:C00469 (ethanol)
result = REST.kegg_get("cpd:C00469").read()
print(result)


# In[13]:


# Get entry information for KSE_17560
result = REST.kegg_get("sce:YAL067C").read()
print(result)


# In[14]:


# Get coding sequence for YAL067C
result = REST.kegg_get("sce:YAL067C", "ntseq").read()
print(result)



# In[15]:


# Get protein sequence for YAL067C
result = REST.kegg_get("sce:YAL067C", "aaseq").read()
print(result)


# In[16]:


# Get map of glycolysis
result = REST.kegg_get("sce00010", "image").read()
Image(result)


# In[17]:


# Get data for glycolysis in S. ceevisiae
result = REST.kegg_get("sce00010").read()
print(result)


# In[18]:


# Get genes involved with glycolysis map
result = REST.kegg_link("compound", "map00010").read()
pd.read_table(io.StringIO(result), header=None)


# In[19]:


# Get reactions involved with glycolysis
# find maps here: https://www.genome.jp/kegg/pathway.html
result = REST.kegg_link("rn", "map00010").read()
pd.read_table(io.StringIO(result), header=None)


# In[20]:


# Get reactions R09085
result = REST.kegg_get("R09085").read()
print(result)


# In[21]:


# Get EC numbers involved with glycolysis
# find maps here: https://www.genome.jp/kegg/pathway.html
result = REST.kegg_link("ec", "map00010").read()
pd.read_table(io.StringIO(result), header=None)


# In[22]:


# KO version (KEGG orthologues)
result = REST.kegg_get("ko00061").read()   
pd.read_table(io.StringIO(result), header=None)


# In[23]:


res=REST.kegg_get("sce00010").read()
print(res)


# In[24]:


# for obtain the genes corresponding to the pathways through text analysis.
# we need to get the results for a specific pathway and run over it

pathways = ["sce00010"]
genes = []

for pathway in pathways:
    #pathways[pathway]['geneid'] = set(); 
    #pathways[pathway]['gene_symbol'] = set()
    pathway_file = REST.kegg_get(pathway).read()  # query and read each pathway
    # iterate through each KEGG pathway file, keeping track of which section
    # of the file we're in, only read the gene in each pathway
    current_section = None
    for line in pathway_file.rstrip().split("\n"):
        section = line[:12].strip()  # section names are within 12 columns
        if not section == "":
            current_section = section
        if current_section == "GENE":
            try:
                gene_identifiers, _ = line[12:].split("; ")[:2]
                geneid, gene_symbol = gene_identifiers.split()
                #pathways[pathway]['geneid'].add(int(geneid))
                #pathways[pathway]['gene_symbol'].add(gene_symbol)
                genes.append(geneid)
            except: pass 


# In[25]:


genes


# ## 2. Bioservices
# 
# Bioservices is a Python package that provides access to many Bioinformatices Web Services (e.g., UniProt) and a framework to easily implement Web Services wrappers (based on WSDL/SOAP or REST protocols).
# 
# https://github.com/cokelaer/bioservices
# 
# 

# In[26]:


#!pip install bioservices


# In[27]:


from bioservices.kegg import KEGG
k = KEGG()


# In[28]:


#show organisms in KEGG and the id symbol (ex. Homo sapiens = hsa)

print(k.list("organism"))


# In[29]:


#sey an organism
k.organism = "sce"


# In[30]:


#show all pathways from S. cerevisiae
k.pathwayIds


# In[31]:


# print a pathway sce00010 (glycolysis) from S. cerevisiae
print(k.get("sce00010"))


# In[32]:


# All organisms with sacc in their names
k.lookfor_organism("sacc")


# In[33]:


# find pathways by name (words B + cell)
print(k.find("pathway", "B+cell"))


# In[34]:


# get all pathways that contain the gene YKR097W in S. cerevisiae
res = k.get_pathway_by_gene("YKR097W","sce")
res


# In[ ]:


#jupyter nbconvert --to webpdf --allow-chromium-download NAME.ipynb

