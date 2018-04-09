import argparse

def main():
    '''
    Just a stub atm.

    :return:
    '''
    parser = argparse.ArgumentParser(description='Description of your program')
    parser.add_argument('-i', '--input', dest="input",
                        help='Either a local filepath to a fasta/fastq or an SRR code.  '
                             'The filetype option ("-t") must be specified appropriately.', required=True)
    parser.add_argument('-t', '--type', dest="type",
                        help='Specifies the type of input file.  '
                             'Options are "fasta", "fastq", and "srr".', required=True)

if __name__ == '__main__':
    main()
