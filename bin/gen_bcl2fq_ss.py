import os
import sys


n_lanes = int(sys.argv[1]) # num lanes in run
xl_uid = sys.argv[2]  # the UID of the pool of seq libraries
xl = ['liba','libb']  # load a list of libraries reflecting the pool
center_name=sys.argv[3]  # YOUR LAB NAME

ret_txt = """
[Header]
Assay,TruSeq LT
Description,
Workflow,GenerateFASTQ
Project Name,na
Investigator Name,"""+center_name+""":::SomeLM#
Experiment Name,RU_uid
IEMFileVersion,4
Date,2001 01/01/01 01:01:01
Application,FASTQ Only
Chemistry,Amplicon

[Reads]
151
151

[Settings]

[Data]
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample
_Project,Description"""


for i in xl.libraries:
    sq = i.seq_index.unique_id
    seqa = i.seq_index.seq_a
    seqb = i.seq_index.seq_b
    ctr = 1
    while ctr <= n_lanes:
        print("{0},{1},{1},,,,{2},,{3},,".format(ctr, sq, seqa, seqb))
        ctr += 1
