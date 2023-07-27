#!/usr/bin/env python
''' RepeatMasker search against custom database
input:
- archive with sequencing data
- custom repeat database
'''
import zipfile
import shutil
import os
import subprocess
from parallel import parallel2
import glob

print(os.environ)
def extract_sequences(f):
    # check archive:
    try:
        z=zipfile.ZipFile(f)
        # extract only dirCLXXXX/reads.fas
        seq_list = []
        for filein in z.namelist():
            c1 = filein.lower().startswith("seqclust/clustering/clusters/dir_cl")
            c2 = filein.endswith("reads.fas")
            c3 = filein.endswith("reads.fasta")  # in newer RE2 versions
            if c1 and (c2 or c3):
                outdir = filein.split("/")[3]
                outfile = outdir +"/reads.fas"
                source = z.open(filein)
                os.mkdir(outdir)
                target = open(outfile, "wb")
                shutil.copyfileobj(source, target)
                seq_list.append(outfile)
        if  len(seq_list) == 0:
            raise ValueError()
                
    except zipfile.BadZipfile as e:
        print("unable to extract sequences from archive!")
        raise e
    
    except IOError as e:
        print("unable to extract sequences from archive!")
        raise e

    except ValueError as e:
        print("No sequences found in archive!")
        raise e
    
    seq_list.sort()
    return seq_list

    
def RepeatMasker(reads,database):
    args = ["RepeatMasker", reads, "-q", "-lib", database, "-pa", "1" , "-nolow", "-dir", os.path.dirname(reads)]
    print(args)
    status=subprocess.call(args , stderr = open(os.devnull, 'wb'))
    return status

def summarizeRepeatMaskerOutput(htmlout = "summary.html"):
    cmd = os.path.dirname(os.path.abspath(__file__))+"/rmsk_summary_table_multiple.r"
    args = [ cmd, "dir_CL*/reads.fas", "dir_CL*/reads.fas.out", "RM-custom_output_table"  ]
    status=subprocess.call(args)
    cmd = cmd = os.path.dirname(os.path.abspath(__file__))+"/RM_html_report.R"
    args = [cmd, htmlout]
    status=subprocess.call(args)
    return status
    

def main():
    from optparse import OptionParser
    
    parser = OptionParser()
    parser.add_option("-i", "--input_file", dest="input_file", help="seqclust zip archive")
    parser.add_option("-d", "--database", dest="database", help="custom repeatmasker database")
    parser.add_option("-g", "--galaxy_dir", dest="galaxy_dir", help="Galaxy home directory")
    parser.add_option("-r", "--report", dest="report", help="output html file with report summary",default='report.html')
    
    options, args = parser.parse_args()
    config_file = os.path.dirname(os.path.abspath(__file__))+"/seqclust.config"
    
    
    seq_files = extract_sequences(options.input_file)  ### REMOVE - TESTING
    parallel2(RepeatMasker, seq_files, [options.database])
    
    status = summarizeRepeatMaskerOutput(options.report)
    
     
        
if __name__== "__main__":
    main()
  
    
