import sys; import os
from Bio import SeqIO
from optparse import OptionParser

usage = "python " + sys.argv[0] + " -f <genbank file>";
parser = OptionParser(usage = usage);
parser.add_option("-f", "--file", dest="input_file", help="The path to the genbank file you wish to parse");
(options, args) = parser.parse_args()


if(len(sys.argv) < 2):
    print "\nUsage: " + usage;
    print "\nUse the option -h for HELP\n";
    sys.exit();
else:
    if(not os.path.exists(options.input_file)) :
        print("Error: '%s' file doesn't exist, please check and try again" %(options.input_file))
        sys.exit();


gbk_file = options.input_file


for genome in SeqIO.parse(gbk_file,'genbank'):
    for gene in genome.features:
        if gene.type != "rRNA":
            continue
        if '16S' not in gene.qualifiers['product'][0]:
			if 'l-rRNA' not in gene.qualifiers['product'][0]:
				continue
        gene_seq = gene.extract(genome.seq)

        genome_dict = vars(genome);                 # store the object's attributes into a dictionary
        gi = genome_dict['annotations']['gi'];      # python equilivent to a hash of a hash in perl

        print ">gi|%s|16S rRNA|%s| %s\n%s" % (gi, genome.id, genome.description, gene_seq)

        #exit();        # for debugging purposes