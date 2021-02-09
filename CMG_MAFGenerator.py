import os
from os.path import expanduser, join, sep, isfile
import subprocess
from configparser import ConfigParser
import time, datetime
import logging
import socket
import pandas as pd
import concurrent.futures
logger = None

# Create log directory
if not os.path.isdir('log'):
    os.mkdir('log')

# Before we do anything, set up logging.
epoch_time = time.mktime(datetime.datetime.now().timetuple())
log_filename = os.path.join(f"log/germlineMafConverter_{str(int(epoch_time))}.log")

print(f'Writing to {log_filename}')
logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s %(message)s", filename=log_filename)
logger = logging.getLogger()

print(f"Starting download on {socket.gethostname()} with process ID: {os.getpid()}")

# Read vcf metadata file into memory
vcf_metadata: pd.DataFrame = pd.read_csv('/N/project/phi_ri/cmg_dev_ri/germline_cmg_saliva_join_vcfs/mm_cmg_vcf_germline_files.csv', delimiter='|',
                           names=['vendor', 'sample_id', 'filename', 'patient_id', 'disease', 'uuid']).astype(str)

vcfroot = '/N/project/phi_ri/cmg_dev_ri/germline_cmg_saliva_join_vcfs/cmg/ingested-WGS-vcf/'

vcf_metadata['vcftype'] = 'germline'
vcf_metadata['fileexists'] = False

# Confirm that the files that exist on disk are present in the metadata
# Gety normal_id from within the vcf file
for file in os.listdir(vcfroot):
    if file in vcf_metadata.filename.values:
        vcf_metadata.loc[vcf_metadata.filename == file, 'fileexists'] = True
        with open(join(vcfroot, file), 'r') as content:
            for line in content:
                line = line.strip()
                if line.startswith("#CHROM\t"):
                    file_sample = line.split("\t")[-1]
                    vcf_metadata.loc[vcf_metadata.filename == file, 'normal_id'] = file_sample
                    break
    else:
        logger.info(f"File {file} present on disk")

# Only keep the rows for which
vcf_metadata = vcf_metadata.loc[vcf_metadata.fileexists == True]

# logger.info(f'Normalized VCF files found:{len(vcfs)} ')
logger.info(f'VCF files found:{len(vcf_metadata)}')

# Load samtools, bcftools
config = ConfigParser()
config.read('download.cfg')
libraries = dict(config.items('libraries'))

os.environ['LOADEDMODULES'] = libraries['loadedmodules']
os.environ['PATH'] = libraries['path']
os.environ['_LMFILES_'] = libraries['lmfiles']
os.environ['LD_LIBRARY_PATH'] = libraries['ld_library_path']

toolspath = join(expanduser("~"), "cbio_tools", "vcf2maf")
vcfpath = "/N/slate/abhmalat/cmg_mm/vcf2mafConversion/vcf"
mafpath = "/N/slate/abhmalat/cmg_mm/vcf2mafConversion/maf"
#enhancedpath = "/N/slate/abhmalat/cmg_mm/vcf2mafConversion/enhanced"
tmppath = "/N/slate/abhmalat/cmg_mm/vcf2mafConversion/tmp"

[os.makedirs(path, exist_ok=True) for path in [vcfpath, mafpath, tmppath]]

logger.info(f'MAF files before processing: {len(os.listdir(mafpath))}')

vcf_metadata.to_csv(join(mafpath, "vcf_metadata_germline.csv"), sep=",")

logger.info(f"Starting the process. Total VCF files found: {vcf_metadata.shape[0]}")

success = 0
tmpdir = {}



def vcf2maf_call(sample_id, filename, vcftype, normal_id):
    run = str(int(time.mktime(datetime.datetime.now().timetuple())))
    tmpdir = join(tmppath, ("MG_" + sample_id), 'germline')
    logger.info(f'Epoch: {run}, Sample: {sample_id}, Type: {vcftype}, tmpdir: {tmpdir}')

    # Create tmp directories
    os.makedirs(tmpdir, exist_ok=True)

    vcffile = join(vcfroot, filename)
    maffile = join(mafpath, (sample_id + '.' + vcftype + '.maf'))

    if os.path.isfile(maffile):
        logger.info(f"{maffile} already exists. Skipping.")
        return True

    # if vcf.vcftype == 'somatic':
    #     vcf_id, vcf_id_val, new_id, new_id_val = '--vcf-tumor-id', vcf.tumor, '--new-normal-id', vcf.normal
    # else:
    #vcf_id, vcf_id_val, new_id, new_id_val = '--vcf-normal-id', vcf.normal_id, '--new-normal-id', vcf.sample_id

    vcf2maf = [
        'perl', join(toolspath, 'vcf2maf.pl'),
        '--input-vcf', vcffile,
        '--output-maf', maffile,
        '--vep-path', join(toolspath, 'vep'),
        '--vep-data', join(expanduser("~"), '.vep'),
        '--vep-forks', '5',
        '--tmp-dir', tmpdir,
        '--normal-id', sample_id,
        '--tumor-id', sample_id,
        '--vcf-normal-id', normal_id,
        #'--vcf-tumor-id', sample_id,
        '--custom-enst', join(toolspath, 'data', 'isoform_overrides_uniprot'),
        '--ref-fasta',
        join(expanduser("~"),
             '.vep/homo_sapiens/broad_hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta')
    ]

    #print("VCF2MAF: " + ' '.join(vcf2maf))
    logger.info("VCF2MAF: " + ' '.join(vcf2maf))

    subprocess.call(vcf2maf)
    retval = False

    if os.path.isfile(maffile):
        logger.info(f"VCF2MAF succeeded. {maffile} exists!")
        retval = True
    else:
        logger.info(
            f"VCF2MAF FAILED FOR {vcffile}")

    # for f in os.listdir(tmpdir):
    #     logger.info(f"Removing {f}")
    #     os.unlink(f)

    return(retval)


with concurrent.futures.ThreadPoolExecutor(max_workers=5) as f:
    vcf_metadata['mafresult'] = list(f.map(vcf2maf_call,
                         vcf_metadata.sample_id, vcf_metadata.filename,
                         vcf_metadata.vcftype,vcf_metadata.normal_id))
    for index, vcf in vcf_metadata.iterrows():
        logger.info(f"{vcf.filename}: {vcf.mafresult}")




logger.info(f'MAF files found: {len(os.listdir(mafpath))}')  # {success}')
for root, dirs, files in os.walk(mafpath):
    for file in files:
        logger.info(file)

print("Done with CMG Germline files")
