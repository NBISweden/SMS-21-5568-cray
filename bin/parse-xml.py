import xml.etree.ElementTree as ET
import os

# Input and output files
input_file = 'data/Trinity_SoF3I50bpF5_paired_mod.blast.xml'
output_file = 'results/tx2gene.tsv'

# Place in project root directory
os.chdir('/Users/erik.fasterius/projects/5568-cray')

# Open output file connection for writing
output = open(output_file, 'w')

# Read XML data from file
tree = ET.parse(input_file)
root = tree.getroot()

# Loop over each transcript ID
for transcript in root.iter('Iteration'):

    # Get transcript ID
    iter_def = transcript.find('Iteration_query-def').text

    # Loop over each Blast hit
    for hit in transcript.iter('Hit'):

        # Grab the gene name for the hit
        hit_def = hit.find('Hit_def').text
        gene = [s for s in hit_def.split(' ') if 'GN=' in s]

        # Skip to next hit if no gene is present
        if len(gene) == 0:
            continue

        # Strip gene name of unneeded string
        gene = ''.join(gene).replace('GN=', '')

        # Print current [transcript, gene] pair to output file
        # (The assignment is only for muting the output of number of characters
        # written)
        _ = output.write(iter_def + '\t' + gene + '\n')

        # Only use the first Blast hit that includes a gene
        break

# Close the file connection
output.close()
