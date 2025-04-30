import pandas as pd
import numpy as np


# Load the MTmutDB.csv dataframe with all possible mtDNA mutations
db = pd.read_csv('MTmutDB/MTmutDB.csv',sep='\t')
db['Strand'] = db['Strand'].astype(str)
# Create a function to get the complementary nucleotide change for Gs and As
def getCompGA(pair):
    return(pair.replace('GA', 'CT').replace('GT', 'CA').replace('GC', 'CG').replace('AG', 'TC').replace('AT', 'TA').replace('AC', 'TG'))
db['pair'] = db['reference'] + db['Allele']
db['compPair'] = db.apply(lambda x: getCompGA(x['pair']), axis = 1)

# Taking into account each substition type (SCx), calculate the number of all possible synonymous and non-synonymous variants in the mtDNA genome and infer P(syn|SCx)
mtvarscounts = db.loc[db['Consequence'] == 'synonymous_variant'].groupby(['Strand','compPair']).count()['#Uploaded_variation'].reset_index().rename( columns = {'#Uploaded_variation' : 'synonymous_count'})
nonsyscount = db.loc[db['Consequence'] != 'synonymous_variant'].groupby(['Strand','compPair']).count()['#Uploaded_variation'].reset_index().rename( columns = {'#Uploaded_variation' : 'nonsynonymous_count'})
mtvarscounts['nonsynonymous_count'] = nonsyscount['nonsynonymous_count']
mtvarscounts['P_syn-SCx'] = mtvarscounts['synonymous_count'] / (mtvarscounts['synonymous_count'] + mtvarscounts['nonsynonymous_count'])
mtvarscounts['index'] = mtvarscounts['compPair'] + '_' + mtvarscounts['Strand']

# Read the dataframe with unique positions that contains position, trinucleotide information, Gene name and other attributes (symbol, strand, ETC complex...)
mtgenome = pd.read_csv('MTmutDB/mtgenome-bp-resolution.csv', sep = '\t')

# Create a mtDNA genome reference to pick up random ttrinucleotides, which only includes positions found in db dataframe
mtgenome_coding = mtgenome.loc[mtgenome['Location'].isin(list(set(list(db['Location']))))].copy()
mtgenome['Location_'] = mtgenome['Location'].str.replace(':','_') 

# This function inputs a dataframe with observed mutations and outputs dNdS ratio and counts table
def getObserved_dNdS(df):
    # Calculate SCx counts
    df['pair'] = df['reference'] + df['Allele']
    df['compPair'] = df.apply(lambda x: getCompGA(x['pair']), axis = 1)
    df['Strand'] = df['Strand'].astype(str)
    aux = df.loc[df['Consequence'] == 'synonymous_variant'].groupby(['Strand', 'compPair']).count()['#Uploaded_variation'].reset_index().rename( columns = {'#Uploaded_variation' : 'synonymous_count'})
    aux['index'] = aux['compPair'] + '_' + aux['Strand']
    obscount = mtvarscounts[['Strand','compPair']].copy()
    obscount['index'] = obscount['compPair'] + '_' + obscount['Strand']
    obscount['synonymous_count'] = obscount['index'].map(dict(zip(list(aux['index']), list(aux['synonymous_count']))))
    obscount.fillna(0, inplace = True)
    obscount['P_SCx'] = obscount['synonymous_count'] / obscount['synonymous_count'].sum()
    obscount['P_syn-SCx'] = mtvarscounts['P_syn-SCx']
    obscount['P_syn-SCx#P_SCx'] = obscount['P_SCx'] * obscount['P_syn-SCx']
    P = obscount['P_syn-SCx#P_SCx'].sum()
    dNdS = (1 - P) / P
    return (obscount, dNdS)

# Load here your observed mtDNA variants dataset
#df = pd.read_csv('')
df = db.sample(n=1000) # This is a test dataframe with 10 mutations. Modify this line by the one above



##### Load the PCAWG data #####
pca = pd.read_csv('PCAWGS/data_mutations_pcawg.txt', sep= '\t')
pcaaux = pca.copy()
# Filter out "Noncoding" positions according to our mtgenome_coding reference to work only with coding mutations
pca = pca.loc[pca['Start_Position'].isin(list(mtgenome_coding['position']))]
# Filter out indels to calculate dNdS
pca = pca.loc[pca['Variant_Type'] == 'SNP']
# Create a dictionary to modify some cols names
pcacolsdict = {'Start_Position': 'position', 'Variant_Classification': 'Consequence', 'ShortVariantID' : '#Uploaded_variation', 'Reference_Allele' : 'reference', 'Tumor_Seq_Allele2' : 'Allele' }
# Rename some columns to match column names in the functions written below
pca.rename(columns=pcacolsdict, inplace = True)
# Modify some columns to match column names in functions below
pca['Chromosome'] = pca['Chromosome'].str.replace('chrM', 'MT')
pca['#Uploaded_variation'] = pca['Chromosome'] + '_' + pca['position'].astype(str)  + '_' + pca['reference'] + '/' +  pca['Allele']
pca['Location'] = pca['Chromosome'] + ':' + pca['position'].astype(str)
pca['Strand'] = pca['Location'].map(dict(zip(list(mtgenome['Location']), list(mtgenome['Strand']))))
pca['Consequence'] = pca['Consequence'].str.replace('Silent', 'synonymous_variant')

# Get dNdS for a specific dataset
datasetname = 'PCAWG'
obscount, dNdS = getObserved_dNdS(pca)
print('dN/dS ratio (' + datasetname + '): ' + str(dNdS))
#dN/dS ratio (PCAWG): 1.7288826759533134
#################################



print('dN/dS ratio is: ' + str(dNdS))



def getSynonymCountTable(dfaux, mtvarscounts, db, mtgenome_coding):
    # Given any mutation in the input dataset, this steps randomly assigns a similar mutaiton (same trinucleotide context) in the mtDNA genome and outputs the predicted consequence, as well as its Strand and REF/ALT (complementary if needed to get SCx)
    dfaux['RandomVariant'] = dfaux.apply(lambda x: getRandomSynonymSCx(db, x['position'], x['Allele'], x['trinucleotide'], mtgenome_coding), axis = 1)
    # Split the output column into the three desired values: Consequence, Strand and complementary pair
    dfaux['RandomConsequence'] = dfaux['RandomVariant'].str.split('#').str[0]
    dfaux['RandomStrand'] = dfaux['RandomVariant'].str.split('#').str[1]
    dfaux['RandomcompPair'] = dfaux['RandomVariant'].str.split('#').str[2]
    # Compute the SCx from random data obtained, count only synonymous mutations
    aux = dfaux.loc[dfaux['RandomConsequence'] == 'synonymous_variant'].groupby(['RandomStrand', 'RandomcompPair']).count()['#Uploaded_variation'].reset_index().rename( columns = {'#Uploaded_variation' : 'synonymous_count'})
    aux['index'] = aux['RandomcompPair'] + '_' + aux['RandomStrand']
    # Merge expected and observed data to compute the final dNdS ratio
    obscount = mtvarscounts[['Strand','compPair', 'index']].copy()
    obscount['synonymous_count'] = obscount['index'].map(dict(zip(list(aux['index']), list(aux['synonymous_count']))))
    obscount.fillna(0, inplace = True)
    obscount['P_SCx'] = obscount['synonymous_count'] / obscount['synonymous_count'].sum()
    obscount['P_syn-SCx'] = obscount['index'].map(dict(zip(list(mtvarscounts['index']), list(mtvarscounts['P_syn-SCx']))))
    obscount['P_syn-SCx#P_SCx'] = obscount['P_SCx'] * obscount['P_syn-SCx']
    P = obscount['P_syn-SCx#P_SCx'].sum()
    return((1 - P) / P)


def getRandomSynonymSCx(db, position, allele, intrin, mtgenome_coding):
    return('#'.join(list(db[['Consequence', 'Strand', 'compPair']].loc[(db['position'] == getRandomTrinuc_getConsequence(mtgenome_coding, position, intrin))].values[0]))) # & (db['Consequence'] == 'synonymous_variant')

# This function takes a mtgenome dataframe, the position of any input variant and the trinucleotide context and randomly outputs another mtDNA position with the same trinucleotide context 
def getRandomTrinuc_getConsequence(mtgenome_coding, position, intrin):
    return (mtgenome_coding['position'].loc[(mtgenome_coding['trinucleotide'] == intrin) & (mtgenome_coding['position'] != position)].sample(n=1).values[0])


# Calculate dN/dS ratio from N simulated datasets
# Define the number of iterations
N = 1000
simdata = pd.DataFrame({'n' : np.arange(1,N)})
from timeit import default_timer as timer
start = timer()
simdata['random_dNdS'] = simdata.apply(lambda x: getSynonymCountTable(df, mtvarscounts, db, mtgenome_coding), axis = 1)
end = timer()
print(end - start) # With N = 1000, it takes 30'
simdata.to_csv('MTmutDB/simDataset.csv', sep = '\t', index = None)