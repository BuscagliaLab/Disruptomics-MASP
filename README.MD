
# Annotation and classification of Mucin Associated Surface Proteins (MASP) in Trypanozoma cruzi 

# What is Disruptomics? 

We herein developed a bioinformatic protocol for the accurate identification, characterization, quantitation and annotation of MASP molecules.

# NOTE
From GENOME
To start you must run your genome sequence in getorf EMBOSS: (we recommend these parameters)

getorf -minsize 120 -maxsize 20100 -find 1

The output obtained in EMBOSS will be the input of the algorithm. 

From PROTEOME
In case you want to run a protein fasta file,  you must ensure that the header of the sequences have the same format as EMBOOS. 
Example: >XX_2517_1_1_1 [224 - 382]

# How to run Disruptomics?
Make sure you download the ZIP folder containing dependent files and the algorithm inside. Unzip the file, and NOT DELETE anything that is inside of the folder.

Then you run:
python Disruptomics-MASP.py 
You must clarify the PATH where you saved the folder containing the algorithm, or directly run it inside the same folder. 

# Dependencies 
You should know that this tool, in order to operate correctly, needs some libraries installed previously.

--> Pandas version 2.8.0 and upwards

--> PyHMMER from  https://pyhmmer.readthedocs.io/en/stable/

--> tkinter (is the default Python interface to the Tk GUI toolkit)
