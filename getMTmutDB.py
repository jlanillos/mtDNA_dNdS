## THIS SCRIPT AIMS AT CREATING A mtDNA db WITH ALL POSSIBLE MUTATIONS IN THE MITOCHONDRIAL GENOME.
## IT HAS TWO PARTS (A and B).

import pandas as pd
from Bio import SeqIO
import numpy as np

# PART A

# You need the mtDNA reference genome. Given the whole reference genome, I got it like:
refgenfilepath = 'refgen/MT_Homo_sapiens.GRCh37.GATK.illumina.fasta'
fasta_sequences = SeqIO.parse(open(refgenfilepath),'fasta')
d = dict()
for fasta in fasta_sequences:
    name, sequence = fasta.id, str(fasta.seq)
    d[name] = sequence

# Create a dataframe with all possible reference>alternative variants
variants = list('ACGT')*len(d['MT'])
reference = list(d['MT'])*4
df = pd.DataFrame({'chrom': ['MT']*len(variants), 'pos': list(np.arange(1,len(d['MT'])+1))*4})
df['ref'] = reference
df.sort_values(by='pos', ascending =True, inplace = True)
df['alt'] = variants
# Remove synonymos changes (i.e., A>A)
df = df[df['ref'] != df['alt']]

#df['vep_format'] = df['ref'] + '/' + df['alt']
#df['id'] = df['chrom] + '_' + df['pos'].astype(str) + '_' + df['vep_format']
#df[['chrom', 'pos','pos','vep_format', 'id']].to_csv('allMTmutations.vepinput.txt',sep='\t', index = None, header = None)


# ooor in VCF format

df['id'] = '.'
df['qual'] = '.'
df['filter'] = 'PASS'
df['info'] = 'SOMATIC;VT=SNP'
df[['chrom', 'pos','id','ref','alt','qual','filter','info']].to_csv('MTmutDB/allMTmutations.vepinput.vcf',sep='\t', index = None, header = None)


###################################################################################################################################################
'''# Run in Bash;
cat header.txt MTmutDB/allMTmutations.vepinput.vcf > tmp.vcf # I got the header from here: https://github.com/mskcc/vcf2maf/issues/78
vep -i tmp.vcf -o MTmutDB/vep_allMTmutations.vepinput.vcf --fork 16 -format vcf  -v --force_overwrite --assembly GRCh38 --cache --dir_cache /Users/lanillj/HMF_mtdna/refgen/vep/cache/

grep -v ^## MTmutDB/vep_allMTmutations.vepinput.vcf > MTmutDB/2_vep_allMTmutations.vepinput.vcf'''
###################################################################################################################################################

# PART B
# This part takes a VEP annotated VCF with all possible mtDNA mutations, 
# annotates trinucleotide context and prioritizes annotated coding mutations over noncoding duplicates (other genes).

# Read the VEP-annotated VCF file with all possible mtDNA mutations
df = pd.read_csv('MTmutDB/2_vep_allMTmutations.vepinput.vcf', sep = '\t')
# Get REFERENCE nucleotide and mtDNA position from annotation
df['reference'] = df['#Uploaded_variation'].str.split('_').str[2].str.split('/').str[0]
df['position'] = df['#Uploaded_variation'].str.split('_').str[1].astype(int)

# Create dfuniq dataframe, which will keep just one row for each genomic position. It will be helpful to get trinucleotide context
dfuniq = df.drop_duplicates('Location', keep = 'first')
# Sort by location, if they were not sorted already
dfuniq.sort_values(by='position', ascending = True, inplace = True)
dfuniq.reset_index(inplace=True)
dfuniq.drop(columns='index', inplace = True)
# getTrinucleotide function extracts the nucleotides at +/-1 positions to create the trinucleotide
def getTriNucContext(position, dfuniq):
    if (position == 1):
        return (''.join(list(dfuniq['reference'].loc[[dfuniq.tail(1).index.values[0], position -1, position]].values)))
    elif (position == 16569):
        return (''.join(list(dfuniq['reference'].loc[[position - 2, position - 1, 0]].values)))
    else:
        return (''.join(list(dfuniq['reference'].loc[dfuniq['position'].isin([position - 1, position, position + 1])].values)))
dfuniq['trinucleotide'] = dfuniq.apply(lambda x: getTriNucContext(x['position'], dfuniq), axis = 1)

# Map trinucleotides to initial dataframe
df['trinucleotide'] = df['Location'].map(dict(zip(list(dfuniq['Location']), list(dfuniq['trinucleotide']))))

# Here, let's create a bp-resolution mtDNA reference (mtgenome df) with useful info (mt df) such as gene name, strand, gene type and complex. We downloaded such information from biomart (mtDNA genes) and manually added control regions and extra annotations.
mtgenome = dfuniq[['Location', 'reference', 'position', 'trinucleotide']]
mt = pd.read_csv('MTmutDB/mtDNA_regions.tsv', sep = '\t')
def annotateMTfeatures(mt, position, col):
    aux = mt[col].loc[(mt['Gene_start'] <= position) & (mt['Gene_end'] >= position)].values
    if len(aux) > 1:
        return(','.join(list([str(x) for x in aux])))
    elif len(aux) == 1:
        return(str(aux[0]))
    else:
        return('Check_position')
cols = ['Gene_name', 'Gene_table_ID_version','Strand','Gene_type','Complex']
# Annotate mtgenome df and use mtgenome dataframe to annotate our original df
for col in cols:
    mtgenome[col] = mtgenome.apply(lambda x: annotateMTfeatures(mt, x['position'], col), axis = 1)
    df[col] = df['position'].map(dict(zip(list(mtgenome['position']), list(mtgenome[col]))))


# Save mtgenome for future use
mtgenome.to_csv('MTmutDB/mtgenome-bp-resolution.csv', sep = '\t', index = None)


# Let's filter out regions that are "not coding"
df = df.loc[df['Gene_type'] == 'coding']



# Binarize Consequence values into "coding" or "noncoding" value. Helpful to filter duplicates later. "synonymous_variant" is considered as coding, since it happend in a coding region
cons = list(df.groupby([ 'Consequence']).count().index)  # "Cons" stands for Consequence
manual_binary_label = ['noncoding','coding','coding','coding','coding','coding','coding','coding','coding','noncoding']
d = dict(zip(cons, manual_binary_label))
'''d = {'downstream_gene_variant': 'noncoding',
 'incomplete_terminal_codon_variant,coding_sequence_variant': 'coding',
 'missense_variant': 'coding',
 'start_lost': 'coding',
 'stop_gained': 'coding',
 'stop_gained,start_lost': 'coding',
 'stop_lost': 'coding',
 'stop_retained_variant': 'coding',
 'synonymous_variant': 'coding',
 'upstream_gene_variant': 'noncoding'}
  # 'non_coding_transcript_exon_variant': 'noncoding','''

# "d" dictionary is useful to label any coding/noncoding variant, which can be further used to calculate dNdS ratio
df['simplified_Consequence'] = df['Consequence'].map(d)

### INCLUDE HERE FILTERING BY CODING/NONCODING REGION: MITOMAP

# These steps first sorts rows by this new "binary" consequence (it puts "coding" values first) and then they filter duplicates based on '#Uploaded_variation' and "simplified_Consequence". This way we keep just one value per variant and prioritize coding variants. In other words, if a specific variant (MT_XXX_A/G were labeled as coding (let's say, a missense) for one gene but algo noncoding (upstream variant) for another, it will filter out the noncoding.
df.sort_values(by='simplified_Consequence', inplace=True)
filt = df.drop_duplicates(subset='#Uploaded_variation', keep = 'first')
filt.sort_values(by='position', ascending = True, inplace=True)

# 'filt' is a dataframe with all possible unique mtDNA variants which prioritized coding versus noncoding when removind duplicated variants.
filt.to_csv('MTmutDB/MTmutDB.csv',sep='\t', index = None)