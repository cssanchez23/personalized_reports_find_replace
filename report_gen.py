#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys 
import os 
import zipfile
import datetime
import subprocess


SAMPLE = sys.argv[1]
template = sys.argv[2]
#copy the template to the current dir 
#os.system("cp /Users/diversigen/Desktop/astro_z/template .")
#os.system("aws s3 sync s3://garoutte-prod-temp/MQ2422_CUSTOMERBUCKET/Processed_reads/{}/ .".format(SAMPLE))
os.system("aws s3 sync s3://garoutte-prod-temp/MQ2422_CUSTOMERBUCKET/Processed_reads/{}/bbmap bbmap/ --exclude '*bam'".format(SAMPLE))
os.system("aws s3 sync s3://garoutte-prod-temp/MQ2422_CUSTOMERBUCKET/Processed_reads/{}/rast rast/ ".format(SAMPLE))
os.system("aws s3 cp s3://garoutte-prod-temp/MQ2422_CUSTOMERBUCKET/Processed_reads/{}/QUAST/report.txt QUAST/ ".format(SAMPLE))
os.system("aws s3 sync s3://garoutte-prod-temp/MQ2422_CUSTOMERBUCKET/Processed_reads/{}/busco busco/ ".format(SAMPLE))
os.system("aws s3 sync s3://garoutte-prod-temp/MQ2422_CUSTOMERBUCKET/Processed_reads/{}/prokka prokka/ ".format(SAMPLE))
os.system("aws s3 cp s3://garoutte-prod-temp/MQ2422_CUSTOMERBUCKET/Processed_reads/{}/assembly.fasta . ".format(SAMPLE))





out_unmapped = subprocess.check_output("cat bbmap/*unmapped.fq | wc -l",shell=True)
out_mapped = subprocess.check_output("cat bbmap/*mapped.fq | wc -l",shell=True)
UNMAP_COUNT= float(out_unmapped)/4
MAP_COUNT= float(out_mapped)/4
TOT_BB = UNMAP_COUNT + MAP_COUNT

unmap_p = UNMAP_COUNT * 100 / TOT_BB
map_p = MAP_COUNT * 100 /TOT_BB
map_p = str(map_p)
unmap_p = str(unmap_p)
MAP_PER =map_p[:6]
UNMAP_PER= unmap_p[:5]

 
COV = subprocess.check_output("cat bbmap/bincov.txt | grep Mean | cut -f2 | cut -d'.' -f1",shell=True)
CONTIG = subprocess.check_output("cat QUAST/report.txt | grep 'contigs (>= 0 bp)' | sed -n 2p | cut -d' ' -f14  ",shell=True)
CON1000 = subprocess.check_output("cat QUAST/report.txt | grep '# contigs (>= 1000 bp)' | cut -d' ' -f11",shell=True)
CON5000 = subprocess.check_output("cat QUAST/report.txt | grep '# contigs (>= 5000 bp)' | cut -d' ' -f11",shell=True)
CON10000 = subprocess.check_output("cat QUAST/report.txt | grep '# contigs (>= 10000 bp)' | cut -d' ' -f10",shell=True)
CON25000 = subprocess.check_output("cat QUAST/report.txt | grep '# contigs (>= 25000 bp)' | cut -d' ' -f10",shell=True)
CON50000 = subprocess.check_output("cat QUAST/report.txt | grep '# contigs (>= 50000 bp)' | cut -d' ' -f10",shell=True)
LEN0 = subprocess.check_output("cat QUAST/report.txt | grep 'Total length (>= 0 bp)' | sed -n 2p | cut -d' ' -f11",shell=True)
LEN1000 = subprocess.check_output("cat QUAST/report.txt | grep 'Total length (>= 1000 bp)' | cut -d' ' -f8",shell=True)
LEN5000 = subprocess.check_output("cat QUAST/report.txt | grep 'Total length (>= 5000 bp)' | cut -d' ' -f8",shell=True)
LEN10000 = subprocess.check_output("cat QUAST/report.txt | grep 'Total length (>= 10000 bp)' | cut -d' ' -f7",shell=True)
LEN25000 = subprocess.check_output("cat QUAST/report.txt | grep 'Total length (>= 25000 bp)' | cut -d' ' -f7",shell=True)
LEN50000 = subprocess.check_output("cat QUAST/report.txt | grep 'Total length (>= 50000 bp)' | cut -d' ' -f7",shell=True)
NUM_CON = subprocess.check_output("cat QUAST/report.txt | grep '# contigs' | sed -n 8p | cut -d' ' -f21",shell=True)
LARGEST = subprocess.check_output("cat QUAST/report.txt | grep 'Largest contig' | cut -d' ' -f16",shell=True)
TOTAL_LEN = subprocess.check_output("cat QUAST/report.txt | grep 'Total length' | sed -n 8p | cut -d' ' -f18",shell=True)
GC = subprocess.check_output("cat QUAST/report.txt | grep 'GC (%) ' | cut -d' ' -f24 # need this two times",shell=True)
N50 = subprocess.check_output("cat QUAST/report.txt | grep 'N50' | cut -d' ' -f26",shell=True)
N75 = subprocess.check_output("cat QUAST/report.txt | grep 'N75' | cut -d' ' -f26",shell=True)
L50 = subprocess.check_output("cat QUAST/report.txt | grep 'L50' | cut -d' ' -f26",shell=True)
L75 = subprocess.check_output("cat QUAST/report.txt | grep 'L75' | cut -d' ' -f26",shell=True)
NPER = subprocess.check_output("cat QUAST/report.txt | grep 'per 100 kbp' | cut -d' ' -f16",shell=True)


GEN_COM = subprocess.check_output("cat busco/short_summary_busco_*.txt | grep 'C:' | cut -d':' -f2 | cut -d'%' -f1",shell=True)
BUS_TOT = subprocess.check_output("cat busco/short_summary_busco_*.txt | grep 'Complete BUSCOs' | cut -f2",shell=True)
BUS_SING = subprocess.check_output("cat busco/short_summary_busco_*.txt | grep 'Complete and single-copy BUSCOs ' | cut -f2",shell=True)
BUS_DUP = subprocess.check_output("cat busco/short_summary_busco_*.txt | grep 'Complete and duplicated BUSCOs' | cut -f2",shell=True)
BUS_FRAG = subprocess.check_output("cat busco/short_summary_busco_*.txt | grep 'Fragmented BUSCOs (F)' | cut -f2",shell=True)
BUS_MISS = subprocess.check_output("cat busco/short_summary_busco_*.txt | grep 'Missing BUSCOs (M)' | cut -f2",shell=True)

#TOT_BB = subprocess.check_output("cat counts.txt | grep total | cut -d' ' -f5",shell=True)
#MAP_COUNT = subprocess.check_output("cat counts.txt | grep 'map counts' | cut -d' ' -f5",shell=True)
#UNMAP_COUNT = subprocess.check_output("cat counts.txt | grep 'unmapp counts' | cut -d' ' -f5",shell=True)
#UNMAP_PER = subprocess.check_output("cat counts.txt | grep 'unmap percentage' | cut -d' ' -f5 | cut -c1-5",shell=True)
#MAP_PER = subprocess.check_output("cat counts.txt | grep 'map percentage' | sed -n 2p | cut -d' ' -f6 | cut -c1-6",shell=True)

TRNA = subprocess.check_output("cat prokka/PROKKA_*.txt | grep 'tRNA' | cut -d' ' -f2",shell=True)
CDS = subprocess.check_output("cat prokka/PROKKA_*.txt | grep 'CDS' | cut -d' ' -f2",shell=True)
TMRNA = subprocess.check_output("cat prokka/PROKKA_*.txt | grep 'tmRNA' | cut -d' ' -f2",shell=True)

print (type(TRNA))
print ((CON25000))

def create_doc(template_file, out_file, replaceText):
    templateDocx = zipfile.ZipFile(template_file)
    #outdir = '/'.join(out_file.split('/')[:-1])
    #if not os.path.exists(outdir):
    #    os.makedirs(outdir)
    newDocx = zipfile.ZipFile(out_file, "w")
    for file in templateDocx.filelist:
        content = templateDocx.read(file)
        for key in replaceText.keys():
            content = content.replace((key), str(replaceText[key]))
        newDocx.writestr(file.filename, content)
    templateDocx.close()
    newDocx.close()
    
#basepath_in = ''
basepath_out = 'outputDir/'
    

replaceText = {
                "{{SAMPLE}}" : SAMPLE,
                "COV" : COV,
                "{{CONTIG}}" : CONTIG,
                "{{CON1000}}" : CON1000,
                "con_5" : CON5000,
                "con_10" : CON10000,
                "con_25" : CON25000,
                "con_50" : CON50000,
                "{{LEN0}}" : LEN0,
                "{{LEN1000}}" : LEN1000,
                "{{LEN5000}}" : LEN5000,
                "{{LEN10000}}" : LEN10000,
                "{{LEN25000}}" : LEN25000,
                "LEN50000" : LEN50000,
                "{{NUM_CON}}" : NUM_CON,
                "{{LARGEST}}" : LARGEST,
                "{{TOTAL_LEN}}" : TOTAL_LEN,
                "{{GC}}" : GC,
                "n50" : N50,
                "{{N75}}" : N75,
                "{{L50}}" : L50,
                "{{L75}}" : L75,
                "NPER" : NPER,
                "{{GEN_COM}}" : GEN_COM,
                "{{BUS_TOT}}" : BUS_TOT,
                "{{BUS_SING}}" : BUS_SING,
                "{{BUS_DUP}}" : BUS_DUP,
                "{{BUS_FRAG}}" : BUS_FRAG,
                "{{BUS_MISS}}" : BUS_MISS,
                "{{TOT_BB}}" : TOT_BB,
                "{{MAP_COUNT}}" : MAP_COUNT,
                "{{UNMAP_COUNT}}" : UNMAP_COUNT,
                "un_per" : UNMAP_PER,
                "MAP_PER" : MAP_PER,
                "{{TRNA}}" : TRNA,
                "{{CDS}}" : CDS,
                "{{TMRNA}}" : TMRNA
                }
#out_file = basepath_out+template_file+SAMPLE
out_file = "{}_Genome_report.docx".format(SAMPLE)
for key in replaceText.keys():
    out_file = out_file.replace(str(key), str(replaceText[key]))
create_doc(template, out_file,replaceText)


os.system("rm template_Genome_Report_CUSTOMER.docx")

os.system("mkdir deliverable ")
os.system("cd deliverable/")
os.system("cp *docx deliverable/")
os.system("cp assembly.fasta deliverable/")
os.system("cp prokka/PROKKA_*.faa deliverable/")
os.system("cp prokka/PROKKA_*.ffn deliverable/")
os.system("cp prokka/PROKKA_*.fsa deliverable/")
os.system("cp prokka/PROKKA_*.gbf deliverable/")
os.system("cp prokka/PROKKA_*.gff deliverable/")
os.system("cp rast/*txt deliverable/")


















