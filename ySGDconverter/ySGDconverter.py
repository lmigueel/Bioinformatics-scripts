import os
import pandas as pd

def convert(input_list,op):

	## reading structured file
	converter = pd.read_csv("saccharomyces_cerevisiae_all_conversao_GO.txt", sep='\t',comment='#',engine='python')

	## creating dictionaries
	sgd = converter[['SGDID','GeneID']]
	sgd = sgd.set_index('SGDID')['GeneID'].to_dict()

	alias = converter[['SGDID','Alias']]
	alias = alias.set_index('SGDID')['Alias'].to_dict()

	#go = converter.set_index('SGDID')['GO Term'].to_dict()

	## starting the conversion
	my_genes = input_list

	# Option 1: SGDID to Gene Name

	if op == 1:
		export =[]
		for g in my_genes:
			export.append(sgd[g])
		
		return export

	# Option 2: Gene name to SGDID, considering Aliases

	if op == 2:

		key_list = list(alias.keys())
		val_list = list(alias.values())

		alias_final= {}
		
		for val in range(0,len(val_list)):
			aux = val_list[val].split(',')
			key1 = key_list[val]
    
			if val_list[val] == '-':
				alias_final[key1] = key1
			
			for g in aux:
				alias_final[g] = key1

		export =[]
		for g in my_genes:
		    export.append(alias_final[g])

		return export

