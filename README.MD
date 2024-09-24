# Annotation and classification of Mucin Associated Surface Proteins (MASP) in Trypanozoma cruzi 

# What is Disruptomics? 

We herein developed a bioinformatic protocol for the accurate identification, characterization, quantitation and annotation of MASP molecules.

# NOTE
To start you must run your genome sequence in getorf EMBOSS: (we recommend these parameters)

getorf -minsize 120 -maxsize 20100 -find 1

The output obtained in EMBOSS will be the input of the algorithm. 

# How to run Disruptomics?
Inside the Disruptomics folder:

python xxxx.py Configuration_Files/basename_configuration_file.txt
You must change basename in basename_configuration_file.txt, for whatever name you have put in the Basename field in the configuration file
