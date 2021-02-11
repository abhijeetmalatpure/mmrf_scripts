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
log_filename = os.path.join(f"log/MMSomaticMafConverter_{str(int(epoch_time))}.log")

print(f'Writing to {log_filename}')
logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s %(message)s", filename=log_filename)
logger = logging.getLogger()

print(f"Starting download on {socket.gethostname()} with process ID: {os.getpid()}")

# Read vcf metadata file into memory
# vcf_metadata: pd.DataFrame = pd.read_csv('/N/project/phi_ri/cmg_dev_ri/germline_cmg_saliva_join_vcfs/mm_cmg_vcf_germline_files.csv', delimiter='|',
#                            names=['vendor', 'sample_id', 'filename', 'patient_id', 'disease', 'uuid']).astype(str)

vcf_metadata = pd.DataFrame(columns=['filename', 'vcf_normal', 'vcf_tumor', 'normal_id', 'tumor_id', 'vcf_type'])
vcfroot = '/N/project/phi_ingest_mm/Targeted_Panel/Celgene_TP_IU_samples_vcf_files_hg38'

# Confirm that the files that exist on disk are present in the metadata
# Gety normal_id from within the vcf file
for file in os.listdir(vcfroot):
    if file.endswith(".vcf"):
        with open(join(vcfroot, file), 'r') as content:
            vcf_type = file.split('_')[3]
            vcf_tumor = "TUMOR"
            vcf_normal = "NORMAL"
            tumor_id = file.split('_')[1]
            normal_id = ""
            for line in content:
                line = line.strip()
                if line.startswith("#CHROM\t"):
                    if vcf_type == 'somaticSV':
                        vcf_tumor = line.split("\t")[-1]
                        vcf_normal = line.split("\t")[-2]
                        tumor_id = file.split('_')[1]
                        normal_id = vcf_normal.split("_")[1]
                    vcf_metadata = vcf_metadata.append({'filename':file,
                                                        'vcf_normal':vcf_normal, 'vcf_tumor':vcf_tumor,
                                                        'normal_id':normal_id, 'tumor_id':tumor_id,
                                                        'vcf_type':vcf_type},
                                                       ignore_index=True)
                    break


vcf_metadata = pd.merge(vcf_metadata, vcf_metadata[['normal_id', 'tumor_id']].loc[vcf_metadata.vcf_type == 'somaticSV'],
                        how='left', left_on='tumor_id', right_on='tumor_id', suffixes=('', '_y'))
vcf_metadata.normal_id.loc[vcf_metadata.normal_id == ""] = vcf_metadata.normal_id_y
vcf_metadata.drop('normal_id_y', axis=1, inplace=True)

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
tmppath = "/N/slate/abhmalat/cmg_mm/vcf2mafConversion/tmp"
mafpath = "/N/slate/abhmalat/cmg_mm/vcf2mafConversion/maf/somatic"

[os.makedirs(path, exist_ok=True) for path in [mafpath, tmppath]]

run = str(int(time.mktime(datetime.datetime.now().timetuple())))
vcf_metadata.to_csv(join(mafpath, "vcf_metadata_mm_somatic" + run + ".csv"), sep=",", index=False)

logger.info(f'MAF files before processing: {len(os.listdir(mafpath))}')
logger.info(f"Starting the process. Total VCF files found: {vcf_metadata.shape[0]}")


def vcf2maf_call(filename, vcf_normal, vcf_tumor, normal_id, tumor_id, vcf_type):
    logger.info(f"Starting with {filename}")
    tmpdir = join(tmppath, ("SM_" + tumor_id), vcf_type)
    logger.info(f'Sample: {tumor_id}, Type: {vcf_type}, tmpdir: {tmpdir}')

    # Create tmp directories
    os.makedirs(tmpdir, exist_ok=True)

    vcffile = join(vcfroot, filename)
    maffile = join(mafpath, (tumor_id + '.' + vcf_type + '.maf'))

    fasta = '.vep/homo_sapiens/102_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz' if vcf_type == 'somaticSV' \
        else '.vep/homo_sapiens/broad_hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta'

    if os.path.isfile(maffile):
        logger.info(f"{maffile} already exists. Skipping.")
        return True

    vcf2maf = [
        'perl', join(toolspath, 'vcf2maf.pl'),
        '--input-vcf', vcffile,
        '--output-maf', maffile,
        '--vep-path', join(toolspath, 'vep'),
        '--vep-data', join(expanduser("~"), '.vep'),
        '--vep-forks', '5',
        '--tmp-dir', tmpdir,
        '--normal-id', normal_id,
        '--tumor-id', tumor_id,
        '--vcf-normal-id', vcf_normal,
        '--vcf-tumor-id', vcf_tumor,
        '--ncbi-build', 'GRCh38',
        '--custom-enst', join(toolspath, 'data', 'isoform_overrides_uniprot'),
        '--ref-fasta',
        join(expanduser("~"), fasta)
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

    for fl in os.listdir(tmpdir):
        if os.path.isfile(fl) and fl.endswith('.vcf'):
            logger.info(f"Removing {fl}")
            os.unlink(f)

    return retval


with concurrent.futures.ThreadPoolExecutor(max_workers=5) as f:
    vcf_metadata['maf_result'] = list(f.map(vcf2maf_call,
                                           vcf_metadata.filename, vcf_metadata.vcf_normal,
                                           vcf_metadata.vcf_tumor,vcf_metadata.normal_id,
                                           vcf_metadata.tumor_id, vcf_metadata.vcf_type ))
    for index, vcf in vcf_metadata.iterrows():
        logger.info(f"{vcf.filename}: {vcf.maf_result}")




logger.info(f'MAF files found: {len(os.listdir(mafpath))}')  # {success}')
for root, dirs, files in os.walk(mafpath):
    for file in files:
        logger.info(file)

logger.info("Done with CMG SomaticSV, SNV, indels files")
