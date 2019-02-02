import gzip
import urllib.request
import os
from Bio import Entrez

def get_cds(organismus = "homo sapiens", dateipfad = "c:/users/arvid/documents/arbeit/organismen/", mail = "arvid.hutzfeldt@student.uni-luebeck.de"):
    Entrez.email = mail  # Always tell NCBI who you are
    # Suche nach Assembly IDs in der NCBI Datenbank für den Organismus. Dabei hat jede Version eines Organismus eine eigene ID
    # (zum Beispiel jeder Stamm von E.Coli). Für höhere Organismen existiert meistens nur ein Eintrag.

    assembly_handle = Entrez.esearch(db="assembly", term=organismus+"[Organismus] AND \"complete genome\"[Assembly Level] AND \"representative genome\"[RefSeq Category]")
    record=Entrez.read(assembly_handle)
    assembly_handle.close()

    if record["Count"] == "0":
        assembly_handle = Entrez.esearch(db="assembly",term=organismus + "[Organismus] AND \"complete genome\"[Assembly Level] AND \"reference genome\"[RefSeq Category]")
        record = Entrez.read(assembly_handle)
        if record["Count"] == "1":
            assembly_level = "referenz genom"
        else:
            assembly_level = "referenz genome"
        assembly_handle.close()

        if record["Count"] == "0":
            assembly_handle = Entrez.esearch(db="assembly",term=organismus + "[Organismus] AND \"chromosome\"[Assembly Level] AND \"reference genome\"[RefSeq Category]")
            record = Entrez.read(assembly_handle)
            if record["Count"] == "1":
                assembly_level = "referenz chromosom"
            else:
                assembly_level = "referenz chromosome"
            assembly_handle.close()
            if record["Count"] == "0":
                gefunden = "vollständige Genome"
                assembly_handle = Entrez.esearch(db="assembly",term=organismus + "[Organismus] AND \"full genome representation\"[Filter] ")
                record = Entrez.read(assembly_handle)
                if record["Count"] == "1":
                    assembly_level = "genom"
                elif record["Count"]=="0":
                    print("Keine Genome gefunden für "+organismus)
                    return
                else:
                    assembly_level = "genome"
                assembly_handle.close()
    print(record["Count"] + " " + assembly_level + " für " + organismus + " gefunden. Suche nach CDS Datei...")


    # Suche für alle IDs eine Zusammenfassung, aus der die AccessionID ausgelesen werden kann. Dann starte Download der Sequenz.
    for ids in record["IdList"]:

        efetch_handle = Entrez.esummary(db="assembly", id=ids,retmode="fasta")
        record_efetch = Entrez.read(efetch_handle)
        gcf=record_efetch['DocumentSummarySet']['DocumentSummary'][0]['AssemblyAccession']
        assembly_name = record_efetch['DocumentSummarySet']['DocumentSummary'][0]['AssemblyName'] # Der Link zu der Sequenz auf dem FTP-Server ist immer gleich zusammengesetzt. Die AccessionID enthält die Ordnerstruktur:
        link="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/" + gcf[4:7] +"/"+ gcf[7:10] + "/" + gcf[10:13]+"/"+gcf+"_"+assembly_name+"/"+gcf+"_"+assembly_name+"_cds_from_genomic.fna.gz"
        if gcf[0:3]=="GCF":
            try:
                urllib.request.urlretrieve(link, dateipfad+organismus+"_"+ids+"_temp.fna.gz")
                with gzip.open(dateipfad + organismus + "_" + ids + "_temp.fna.gz", 'rb') as f:
                    formatted_file = open(dateipfad + organismus + "_" + ids + "_" + assembly_level + ".fasta", "w+")
                    formatted_file.write(str(f.read()).replace('\\n', '\n').replace('b"',''))  # die heruntergeladene Datei enthält Umbrüche, die von Python nicht erkannt werden und ersetzt werden müssen.
                    formatted_file.close()
                    f.close()
                    os.remove(dateipfad + organismus + "_" + ids + "_temp.fna.gz")
                    return
            except urllib.error.URLError:
                print("Nichts gefunden")
            # Dekomprimieren des fna.gz files und speichern:

    print("Es wurden keine Referenz Genome hinterlegt.")

def get_cds_from_file(document="c:/users/arvid/documents/arbeit/organism_list.txt"):
    try:
        organisms = open(document, "r")
        for lines in organisms:
            get_cds(lines)
    except FileNotFoundError:
        print("Kein Dokument gefunden")

get_cds("latimeria chalumnae")