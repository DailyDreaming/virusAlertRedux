import os
import sys
import subprocess
import time
from subprocess import check_call
from os.path import dirname, realpath

# BioPython
from Bio import Entrez
from Bio.Blast import NCBIWWW # NCBI API wrapper
from Bio.Blast import NCBIXML # parse output BLAST XML
from Bio import SeqIO # ingests fastas


# SRR5383888
# Human w/ HBV (PacBio): https://www.ncbi.nlm.nih.gov/sra/SRX2678953[accn]

# SRR5383891
# Human w/ HBV (PacBio; Control): https://www.ncbi.nlm.nih.gov/sra/SRX2678956[accn]

# SRR6172653
# PacBio Soil Samples: https://www.ncbi.nlm.nih.gov/sra/SRX3283433[accn]

# SRR6172655
# PacBio Soil Samples: https://www.ncbi.nlm.nih.gov/sra/SRX3283431[accn]

def datapath(path):
    return os.path.abspath(os.path.join(os.path.dirname(realpath(__file__)), '..', 'data', path))

def toolpath(path):
    return os.path.abspath(os.path.join(os.path.dirname(realpath(__file__)), '..', 'tools', path))

def homepath(path):
    return os.path.abspath(os.path.join(dirname(os.path.realpath(__file__)), path))

class VirLib(object):
    def __init__(self):
        self.magicblast = toolpath('ncbi-magicblast-1.3.0/bin/magicblast')
        self.sratools_fastqdump = toolpath('sratoolkit.2.9.0-centos_linux64/bin/fastq-dump')

    def process_input(self, filetype, inputfile):
        """
        Returns the filepath to a local fasta file.

        Fetches fasta by SRR code if necessary.

        :param filetype:
        :param inputfile:
        :return:
        """
        print('Processing Inputs.')
        if filetype == 'fastq':
            assert (inputfile.endswith('.fastq') or
                    inputfile.endswith('.fq')), 'Your fastq file does not end in ".fq" or ".fastq".'
            return os.path.abspath(inputfile)
        if filetype == 'fasta':
            assert (inputfile.endswith('.fasta') or
                    inputfile.endswith('.fa')), 'Your fasta file does not end in ".fa" or ".fasta".'
            return os.path.abspath(inputfile)
        if filetype == 'srr':
            assert inputfile.startswith('SRR'), 'SRA codes begin with "SRR", please check your filetype.'
            return self.getSRR(SRR_code=inputfile)

    def magicblast(self, input_file, input_type, db_file, output_dir):
        output_file = os.path.join(output_dir, "magicblast_hits.gz")

        cmd = ('{magicblast} -query {in_file} -db {db_file} -out {out_file} -infmt {in_type} -gzo -outfmt tabular'
               ''.format(magicblast=self.magicblast,
                         in_file=input_file,
                         in_type=input_type,
                         db_file=db_file,
                         out_file=output_file))

        print('Now running: {cmd}'.format(cmd=cmd))
        process = subprocess.Popen(cmd, shell=True)
        process.communicate()

        assert os.path.exists(output_file)
        return output_file

    def virfinder(self):
        virfinder_script = homepath('runVirFinder.sh')
        cmd = [virfinder_script, '--vanilla', datapath('SRR6172653.fastq'), '0.95', toolpath('R')]
        check_call(cmd)

    def getSRR(self, SRR_code, output_dir=datapath('')):
        """
        Takes an SRR code, downloads it as a fastq, and returns the path to that fastq file.

        API Documentation: https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=fastq-dump

        Example: fastq-dump -X 5 -Z SRR5150787

        :param SRR_code:
        :return: filepath of downloaded fasta
        """
        print('Fetching the fastq for: ' + str(SRR_code))
        output_file = os.path.join(output_dir, SRR_code + '.fastq')

        cmd = ('{fastq_dump_path} -Z {SRR_code} > {output_file}'
               ''.format(fastq_dump_path=self.sratools_fastqdump,
                         SRR_code=SRR_code,
                         output_file=output_file))

        print('Now running: {cmd}'.format(cmd=cmd))
        process = subprocess.Popen(cmd, shell=True)
        process.communicate()

        # throttle requests so that we don't contact the server more than once every 10 seconds
        # More Info: https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=DeveloperInfo
        time.sleep(10)

        assert os.path.exists(output_file)
        return output_file

    def run(self, filetype='srr', inputfile='SRR6172653'):
        inputfilepath = self.process_input(filetype, inputfile)

class blastAPILib(object):
    """
    Fetches BLAST queries using the REST API, but it is slow,
    so it has been replaced by blast+ and magic-blast.
    """
    def __init__(self, verbose):
        self.verbose = verbose

    def processInput(self, type, input):
        """
        Returns the filepath to a local fasta file.

        Fetches fasta by SRR code if necessary.

        :param type:
        :param input:
        :return:
        """
        print('Processing Inputs.')
        if type == 'fasta' or type == 'fastq':
            assert (input.endswith('.fasta') or
                    input.endswith('.fa') or
                    input.endswith('.fastq') or
                    input.endswith('.fq'))
            return os.path.abspath(input)
        if type == 'srr':
            return self.getSRR(SRR_code=input)

    def getSRR(self, SRR_code):
        """
        Takes an SRR code, downloads it as a fastq, and returns the path to that fastq file.

        API Documentation: https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=fastq-dump

        :param SRR_code:
        :return: filepath of downloaded fasta
        """
        print('Fetching the SRR.')
        # basic example cmd: "./fastq-dump -X 5 -Z SRR5150787"
        sratools_API_cmd = "{fastq_dump_path} -Z {SRR_code} > {SRR_code}.fastq" \
                           "".format(fastq_dump_path=self.fastq_dump_path,
                                     SRR_code=SRR_code)
        # run the command
        print('Now running: {cmd}'.format(cmd=sratools_API_cmd))
        process = subprocess.Popen(sratools_API_cmd, shell=True)
        process.communicate()
        # assert the output file was created and return
        SRR_filepath = os.path.join(os.getcwd(), SRR_code + '.fastq')
        assert os.path.exists(SRR_filepath)

        # throttle requests so that we don't contact the server more than once every 10 seconds
        # More Info: https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=DeveloperInfo
        time.sleep(11)

        return os.path.join(SRR_filepath)

    def blastAPISearch(self, record, output_xml_path='default.xml'):
        """
        Takes a fasta/fastq, runs BLAST, saves the XML outputs, and returns the xml path..

        :param fastq_path:
        :param output_xml_path:
        :return:
        """
        print('BLASTing the fastq entry.')
        # process fasta as a biopython sequence record object
        result_handle = NCBIWWW.qblast("blastn", "nt", record.seq)
        print('Writing the BLAST results.')

        # save the output as an xml file
        with open(output_xml_path, "w") as out_handle:
            out_handle.write(result_handle.read())
        return os.path.abspath(output_xml_path)

    def fetchXMLResults(self, output_xml_path):
        """
        Looks for the words 'virus or 'viral' in the xml headers and returns those that match.  If
        none match, (slowly) fetch the GenBank accession files and determine the organism.

        :param output_xml_path:
        :return:
        """
        print('Parsing and fetching the xml search results.')
        matches = []
        accessions = []
        with open(output_xml_path, "r") as out_handle:
            blast_records = NCBIXML.parse(out_handle)
            for blast_record in blast_records:
                for alignment in blast_record.alignments:
                    accessions.append(self.parseAccession(alignment.title))
                    # fetch the organism by the accession number instead
                    # possibly write if we care about the e value (hsp.expect)
                    if 'virus' in alignment.title or 'viral' in alignment.title:
                        matches.append(alignment.title)
        if matches:
            print(matches[0])
            return matches[0] # only return the first match
        else:
            accessionMatches = self.checkAccessions(accessions_list=accessions)
            print(accessionMatches[0])
            return accessionMatches[0] # only return the first match

    def parseAccession(self, alignment_title):
        """"
        Given a SeqRecord, return the accession number as a string.

        Modified from biopython.org.

        e.g. "gi|2765613|emb|Z78488.1|PTZ78488" -> "Z78488.1"
        """
        parts = alignment_title.split("|")
        assert parts[0] == "gi"
        # print('Using accession: {}'.format(parts[3]))
        return parts[3]

    def organismNamefromGenBankAccession(self, accessionCode):
        """
        Returns an organism's English name from a GenBank accession number.

        :param accessionCode:
        :return: String
        """
        print('Fetching Organism from GenBank Accession.')

        # throttle requests so that we don't contact the server more than once every 10 seconds
        # More Info: https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=DeveloperInfo
        time.sleep(11)

        Entrez.email = 'lblauvel@ucsc.edu'
        records_handle = Entrez.efetch(db="nuccore", rettype="gb", id=accessionCode)
        for record in records_handle:
            record = record.strip()
            if record.startswith('ORGANISM'):
                print(record)
                return record # [len('ORGANISM'):].strip()

    def checkAccessions(self, accessions_list):
        print('Checking Accessions.')
        resulting_organisms = []
        for accession in accessions_list:
            species = self.organismNamefromGenBankAccession(accession)
            resulting_organisms.append(species)
            break # only fetch the first one
        return resulting_organisms

    def printXMLresults(self, output_xml_path):
        """
        Only for debugging/verifying output.

        :param result_handle:
        :return:
        """
        # read the xml
        with open(output_xml_path, "r") as out_handle:
            blast_records = NCBIXML.parse(out_handle)
            for blast_record in blast_records:
                for alignment in blast_record.alignments:
                    for hsp in alignment.hsps:
                        print('****Alignment****')
                        print('sequence:', alignment.title)
                        print('length:', alignment.length)
                        print('e value:', hsp.expect)
                        print(hsp.query[0:75] + '...')
                        print(hsp.match[0:75] + '...')
                        print(hsp.sbjct[0:75] + '...')

    def entriesFromFastq(self, fastq):
        """
        fastq files may have multiple entries.  This returns a list of BioPython record objects
        associated with each fastq entry.

        :param fastq:
        :return:
        """
        print('Fetching Entries from Fastq.')
        entries = []
        with open(fastq, "rU") as handle:
            for record in SeqIO.parse(handle, "fastq"):
                entries.append(record)
        return entries

    def entriesFromFasta(self, fasta):
        """
        fasta files may have multiple entries.  This returns a list of BioPython record objects
        associated with each fasta entry.

        :param fasta:
        :return:
        """
        print('Fetching Entries from Fasta.')
        entries = []
        with open(fasta, "rU") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                entries.append(record)
        return entries

    def run(self, type='srr', input='SRR6172653'):
        """

        :return:
        """
        # only fetch the first few entries in the fastq
        max_entries = 1

        # fetch a fastq local path
        fastx_path = self.processInput(type=type, input=input)

        if type == 'fastq' or type == 'srr':
            # return a list of BioPython objects representing each fastq entry in the file
            entries = self.entriesFromFastq(fastq=fastx_path)
        elif type == 'fasta':
            # return a list of BioPython objects representing each fastq entry in the file
            entries = self.entriesFromFasta(fasta=fastx_path)

        # BLAST, generate an XML, and parse to obtain the first species hit
        list_of_species = []
        for record in entries[:max_entries]:
            output_xml_path = self.blastAPISearch(record, output_xml_path='default.xml')
            species = self.fetchXMLResults(output_xml_path)
            list_of_species.append(species)

        # generate the output results in a file
        with open('species_results.txt', 'w') as f:
            print('\n\nTop {} hits by p-value in this sample are: '.format(str(max_entries)))
            for s in list_of_species:
                f.write(s)
                f.write('\n')
                print(s)
