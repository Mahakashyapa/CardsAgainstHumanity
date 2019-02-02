from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SearchIO
import matplotlib.pyplot as plt
import numpy as np

def start_blast(path, seq="8332115"):
    result_handle = NCBIWWW.qblast("blastn", "nt", seq)
    with open(path, "w") as out_handle:
        out_handle.write(result_handle.read())
    result_handle.close()
def show_blast(path):
    handle = open(path)
    blast_record=NCBIXML.read(handle)
    E_VALUE_THRESH = 0.01
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < E_VALUE_THRESH:
                print("****Alignment****")
                print("sequence:", alignment.title)
                print("length:", alignment.length)
                print("e value:", hsp.expect)
                print(hsp.query[0:75] + "...")
                print(hsp.match[0:75] + "...")
                print(hsp.sbjct[0:75] + "..."+"\n")


def draw_blast(path="C:/Users/arvid/Documents/arbeit/Blast.xml", max_entry=10, yMax=1000, xMax=1000):
    handle = open(path)
    blast_record = NCBIXML.read(handle)
    E_VALUE_THRESH = 0.01
    plot_werte = []
    dy=yMax/max_entry
    y=yMax-dy

    index=0
    for alignment in blast_record.alignments:
        y = y - dy/2
        index=index+1
        for hsp in alignment.hsps:
            if hsp.expect < E_VALUE_THRESH and index<max_entry:
                plt.text(hsp.query_end, y,"   "+str(index)+" "+alignment.title[alignment.title.find("PREDICTED: ") + 11:alignment.title.find("PREDICTED: ")+31])
                plt.plot([hsp.query_start, hsp.query_end], [y, y], 'g-')
    plt.plot([0,xMax],[yMax-dy,yMax-dy],'r-')


    plt.ylim((0,yMax))
    cur_axes = plt.gca()
    cur_axes.axes.get_yaxis().set_visible(False)
    plt.show()

path="c:/users/arvid/documents/arbeit/blast.xml"
show_blast(path)
draw_blast(path)
