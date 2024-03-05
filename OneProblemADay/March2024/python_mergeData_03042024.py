import pandas as pd


##	I am trying to merge two dataframe by a common columns and found this example
df1 = pd.DataFrame.from_dict( {'company':["tata", "cts", "dell"],'standard':["A1", "A2", "A3"]})
print(df1)

df2 = pd.DataFrame.from_dict({'company':["tata", "dell", "cts", "hcl"], 'return':[71, 78, 27, 23]})
print(df2)

print(pd.merge(df1, df2, on='company'))


##	but, it turn out that I need to concatenate more columns to another dataframe 
##	therefore, really I need use pandas.concat function in the following 
##	function: appendMeta

def appendMeta (mafFile, metaData):
	meta = pd.read_csv(metaData)
	print(meta.shape)
	maf  = pd.read_csv(mafFile, delimiter = "\t")
	numberOfVariants = maf.shape[0]
	TumorSampleBarcode = maf.loc[0, "Tumor_Sample_Barcode"]
	selectRow = meta.loc[meta['Tumor_Sample_Barcode'] == TumorSampleBarcode]

	## old  version of pandas,
	## pandas 2.0 has concat instead
	# meta_dm = selectRow.append([selectRow]*(numberOfVariants -1 ), ignore_index = True)

	meta_dm = pd.concat([selectRow]*numberOfVariants, ignore_index = True)
	print("make sure variant file and metadata has same number of rows")
	print(maf.shape)
	print(meta_dm.shape)
	meta_dm = meta_dm.iloc[:, 1:]
	print(meta_dm.shape)
	appended_maf = pd.concat ([maf, meta_dm], axis= 1)
	print (appended_maf.shape)
	return(appended_maf)


if __name__=="__main__":
	modifiedFile = appendMeta ("sample.maf.txt", "mouseMetaTable.csv")
	modifiedFile.to_csv("sample.maf.w.meta.txt", sep="\t", encoding='utf-8', index=False) 
	print(modifiedFile.head(2))





