# -*- coding: utf-8 -*-
"""
Created on Sun Jan  8 13:42:18 2023

@author: Aldana A. Cepeda Dean
"""
import os
import sys
import re
import pandas as pd 
import re 
import pyhmmer
import collections
import pyhmmer.easel
from IPython.display import display
from IPython.display import display, HTML
import tkinter as tk
from tkinter import filedialog


curr_directory = os.getcwd()
print(curr_directory)
# Mostrar el diálogo de selección de directorio
directory = filedialog.askdirectory(title="Select a directory to save the outputs")

#Este script usa como in un multifasta output de getorf
#Busca las metioninas dentro de cada ORF solo si es maslargo que 40 aa
#Escribe todos los ORFs interiores que sean mayores a 40 aa.
allORFsRA_todos={}
def ORFsinORFS(x):
    if len(x) > 40:
        n=1
        for match in re.finditer(r'M+',x):
            s=match.start()  #guarda todos los indexes que corresponden con una M
            pieza=x[s:] #Recorto desde la M en adelante
            if len(pieza)>40:
                allORFsRA_todos[RA_pred][n]={}
                allORFsRA_todos[RA_pred][n]['descripcion_aprox']=allORFsRA_largos[RA_pred]['0']['descripcion']
                allORFsRA_todos[RA_pred][n]['seq']=pieza
                n=n+1
            else:
                continue

a= filedialog.askopenfilename(title= "Select .fa input file (obtained by GETORF/EMBOSS):")
RAORFpred=open(a)
print(a)

allORFsRA_largos= {}
for line2 in RAORFpred:
    if line2.startswith('>'):
        line2= line2.strip()
        line3= line2.split()
        RA_pred= line3[0][1:]
        allORFsRA_largos[RA_pred]= {}
        allORFsRA_largos[RA_pred]['0']={}
        allORFsRA_largos[RA_pred]['0']['seq']=''
        allORFsRA_largos[RA_pred]['0']['descripcion']= line2
    else:
        line2=line2.strip()
        allORFsRA_largos[RA_pred]['0']['seq']= allORFsRA_largos[RA_pred]['0']['seq']+line2

for RA_pred,v in allORFsRA_largos.items():
    #print(RA_pred)
    allORFsRA_todos[RA_pred]={}
    for kk,vv in v.items():
        #print(kk,vv)
        for kkk,seq in vv.items():
            if kkk=='descripcion':
                if 'REVERSE' in seq:
                    #print(seq)
                    range= re.findall(r'\[[^\]]*\]' ,seq)[0][1:-1]
                    bases= range.split('-')
                    pb0= int(bases[0])
                    pb1= int(bases [1])
                    pbi= pb1
                    pbf= pb0
                    strand = '-'
                    final=pbf
                else:
                    range= re.findall(r'\[[^\]]*\]' ,seq)[0][1:-1]
                    bases= range.split('-')
                    pb0= int(bases[0])
                    pb1= int(bases [1])
                    pbi= pb0
                    pbf= pb1
                    strand = '+'
                    final=pbf
            elif kkk=='seq':
                ORFsinORFS(seq)

out1= 'temp_1'
out=out1
out= open(out,'w+')

for RA_pred,v in allORFsRA_todos.items():
    for kk,vv in v.items():
        kk=str(kk)
        for clave,detalle in vv.items():
            nombre_CDS='>'+RA_pred+'_'+kk+'  '
            if clave=='descripcion_aprox':
                nombre_final = nombre_CDS+detalle[1:]
                out.write(nombre_final)
                out.write("\n")
            if clave=='seq':
                out.write(detalle)
                out.write("\n")
out.close()
RAORFpred2= out1
try:
    RAORFpred2=open(RAORFpred2)
except:
    print('No existe el archivo: %s' % RAORFpred2)

out3 = input('Name of strain to analyze: ')
out2= open(f'{directory}/{out3}_holo-ORFs.fa','w+')

allORFsRA= {}
# crear diccionarios (key es el RA_pred) value es> lista :RA_contig_##, bpi:##, bpf:##, source: ., strand: +o-, bla: .
for line2 in RAORFpred2:
    if line2.startswith('>'):
        line2=line2.strip()
        line3= line2.split()
        RA_pred= line3[0][1:]
        num_contig= RA_pred.split('_')[1]
        contig_RA= 'RA_'+ num_contig
        allORFsRA[RA_pred]= []
        range= re.findall(r'\[[^\]]*\]' ,line2)[0][1:-1]
        bases= range.split('-')
        pb0= int(bases[0])
        pb1= int(bases [1])
        sequence=""
        if pb0 < pb1:
            pbi= pb0
            pbf= pb1
            strand = '+'
        else:
            pbi= pb1
            pbf= pb0
            strand = '-'
        largo= pbf-pbi+1
        slargo = str(largo)
        allORFsRA[RA_pred]= [contig_RA,'getorf','CDS',pbi, pbf,'.', strand, '.', 'id='+ RA_pred, 'len_ORF_Original='+ slargo,sequence]
    else:
        line2=line2.strip()
        sequence=line2
        allORFsRA[RA_pred]= [contig_RA,'getorf','CDS',pbi, pbf,'.', strand, '.', 'id='+ RA_pred, 'len_ORF_Original='+ slargo,sequence]
        largo_suborf=len(line2)
        largo_suborf_3=largo_suborf*3
        largo_a_restar=largo-largo_suborf_3
        orientacion= allORFsRA[RA_pred][6]
        if orientacion=='+':
            nueva_pbi= pbi+largo_a_restar
            allORFsRA[RA_pred][3]=nueva_pbi
            out2.write('>')
            out2.write(RA_pred)
            out2.write("\t")
            out2.write("[")
            out2.write(str(allORFsRA[RA_pred][3]))
            out2.write(" - ")
            out2.write(str(allORFsRA[RA_pred][4]))
            out2.write("]")
            out2.write("\n")
            out2.write(line2)
            out2.write("\n")
        if orientacion=='-':
            nueva_pbf= pbf-largo_a_restar
            allORFsRA[RA_pred][4]=nueva_pbf
            out2.write('>')
            out2.write(RA_pred)
            out2.write("\t")
            out2.write("[")
            out2.write(str(allORFsRA[RA_pred][4]))
            out2.write(" - ")
            out2.write(str(allORFsRA[RA_pred][3]))
            out2.write("]")
            out2.write("\n")
            out2.write(line2)
            out2.write("\n")

MASPra= pd.DataFrame.from_dict(allORFsRA, orient='index')
pd.set_option('display.max_rows', None) 
pd.set_option('display.max_columns', None) 
pd.set_option('display.width', 1000) 
pd.set_option('display.colheader_justify', 'center') 
pd.set_option('display.precision', 2) 
MASPra.columns=['contig', 'getorf', 'CDS', 'Pbi', 'Pbf', '.' ,'Strand', '.', 'id','Largo', 'Sequence']

print("All internal Methionines were calculated")

# Crear el diccionario final a partir del DataFrame.
input_file = {}
final_clasification={}
for index, row in MASPra.iterrows():
    input_file[index] = {
        'Pbi': row['Pbi'],
        'Pbf': row['Pbf'],
        'Sequence': row['Sequence']
    }
    final_clasification[index] ={"gene_ID":[index] ,"N-terminal_type":"","C-terminal_type":"","N-terminal_Evalue":"", "N-terminal_bitscore":"","C-terminal_Evalue":"", "C-terminal_bitscore":"", "Chimera_RegEx_Match":"", "Annotation":[]}    

#ACA TENGO QUE CAMBIAR EL ARCHIVO PORQUE TIENE QUE SER OS.PATHLIKE Y NO TEXTIOWRAPPER
#ARMADO DE HMMS PARA SP Y GPI 

os.remove(out1)


alphabet = pyhmmer.easel.Alphabet.amino()

with pyhmmer.plan7.HMMFile(f"{curr_directory}/SP_pyhmmer.hmm") as hmm_file:
    hmm_SP = hmm_file.read()

SP_file= {}
GPI_file={}
for k,v in input_file.items():
        peptidesignal=v["Sequence"][:30]
        gpianchor=v["Sequence"][-40:]
        SP_file[k]=[v["Pbi"],v["Pbf"],peptidesignal]
        GPI_file[k]=[v["Pbi"],v["Pbf"],gpianchor]

# Escribir los datos en sp_file.txt
with open("temp2.fasta", "w") as sp_filee:
    for k, v in SP_file.items():
        sp_filee.write(f">{k}\n{v[2]}\n")
    sp_filee.close()
# Escribir los datos en gpi_file.txt
with open("temp3.fasta", "w") as gpi_filee:
    for k, v in GPI_file.items():
        gpi_filee.write(f">{k}\n{v[2]}\n")
    gpi_filee.close()

#PARA PEPTIDO SENAL

with pyhmmer.easel.SequenceFile(sp_filee.name, digital=True) as seqs_file:
    SPproteins = seqs_file.read_block()

with pyhmmer.easel.SequenceFile(gpi_filee.name, digital=True) as seqs_file:
    GPIproteins = seqs_file.read_block()

pipeline = pyhmmer.plan7.Pipeline(hmm_SP.alphabet)
hits_SP = pipeline.search_hmm(hmm_SP, SPproteins)
    
Result = collections.namedtuple("Result", ["query", "cog", "bitscore", "evalue"])

results_SP = []
result_sp_dicc={}
alineado_SP=[]
for hits_SP in pyhmmer.hmmsearch(hmm_SP, SPproteins):
    cog = hits_SP.query
    for hit in hits_SP:
        if hit.reported:
                 idd=hit.name.decode()
                 result_sp_dicc[idd]=[idd,cog,hit.score,hit.evalue,alineado_SP]
                 results_SP.append(Result(hit.name.decode(), cog, hit.score, hit.evalue))
        for domain in hit.domains:
                  ali = domain.alignment
                  ali2= str(ali)
                  ali2=ali2.replace("\n","")
                  ali2=ali2.replace(" ","")
                  find=re.search(idd, ali2)
                  if find:
                   result_sp_dicc[idd][4]=ali2
     
best_results_SP = {}
keep_query = set()
for result in results_SP:
    if result.query in best_results_SP:
        previous_bitscore = best_results_SP[result.query].bitscore
        if result.bitscore > previous_bitscore:
            best_results_SP[result.query] = result
            keep_query.add(result.query)
        elif result.bitscore == previous_bitscore:
            if best_results_SP[result.query].cog != hit.cog:
                keep_query.remove(result.query)
    else:
        best_results_SP[result.query] = result
        keep_query.add(result.query)
        
filtered_results = [best_results_SP[k] for k in sorted(best_results_SP) if k in keep_query]
            
for k,v  in result_sp_dicc.items():
     final_clasification[k]["N-terminal_Evalue"]=v[3]
     final_clasification[k]["N-terminal_bitscore"]=v[2]

#for i in filtered_results:
  # print(i,file=open('/home/aldanacepeda/Documentos/PhD/Algoritmo/resultados_MASP_allseq_HMMER_valores_SP1.txt','a'))

tabla_SP = pd.DataFrame.from_dict(result_sp_dicc, orient='index')
pd.set_option('display.max_rows', None) 
pd.set_option('display.max_columns', None) 
pd.set_option('display.max_colwidth', None)
pd.set_option('display.colheader_justify', 'center') 
pd.set_option('display.precision', 2)    
tabla_SP.style.set_properties(**{
    'text-align': 'center',
    'white-space': 'pre-wrap',
})
  
###################################### ARMADO DEL HMM PROFILE DEL GPI 1 #############################################

print("All N-terminal of the entered sequences were analyzed. ")

with pyhmmer.plan7.HMMFile(f"{curr_directory}/GPI_pyhmmer.hmm") as hmm_file_gpi:
    hmm_GPI1 = hmm_file_gpi.read()

hmm_GPI1.consensus

pipeline = pyhmmer.plan7.Pipeline(hmm_GPI1.alphabet)
hits_GPI1 = pipeline.search_hmm(hmm_GPI1, GPIproteins)
Result = collections.namedtuple("Result", ["query", "cog", "bitscore", "evalue"])
result_gpi_dicc= {}

hmm_gpi_profile=hmm_GPI1.to_profile(background=None, L=20, multihit=False, local=True)

results_GPI = []
alineado=[]
for hits_GPI1 in pyhmmer.hmmsearch(hmm_GPI1, GPIproteins):
    cog = hits_GPI1.query
    for hitgpi in hits_GPI1:
             if hitgpi.reported:
                      idd=hitgpi.name.decode()
                      result_gpi_dicc[idd]=[idd,cog,hitgpi.score,hitgpi.evalue,alineado]
                      results_GPI.append(Result(hitgpi.name.decode(), cog, hitgpi.score, hitgpi.evalue))
             for domain in hitgpi.domains:
                       ali = domain.alignment
                       ali2= str(ali)
                       ali2=ali2.replace("\n","")
                       ali2=ali2.replace(" ","")
                       find=re.search(idd, ali2)
                       if find:
                        result_gpi_dicc[idd][4]=ali2

tabla_GPI1 = pd.DataFrame.from_dict(result_gpi_dicc, orient='index')
pd.set_option('display.max_rows', None) 
pd.set_option('display.max_columns', None) 
pd.set_option('display.max_colwidth', None)
pd.set_option('display.colheader_justify', 'center') 
pd.set_option('display.precision', 2)    
tabla_GPI1.style.set_properties(**{
    'text-align': 'center',
    'white-space': 'pre-wrap',
})
       
best_results_GPI = {}
keep_query_GPI = set()
for result in results_GPI:
    if result.query in best_results_GPI:
        previous_bitscore = best_results_GPI[result.query].bitscore
        if result.bitscore > previous_bitscore:
            best_results_GPI[result.query] = result
            keep_query_GPI.add(result.query)
        elif result.bitscore == previous_bitscore:
            if best_results_GPI[result.query].cog != hit.cog:
                keep_query_GPI.remove(result.query)
    else:
        best_results_GPI[result.query] = result
        keep_query_GPI.add(result.query)
        
filtered_results_GPI = [best_results_GPI[k] for k in sorted(best_results_GPI) if k in keep_query_GPI]

#for i in filtered_results_GPI:
 #  print(i,file=open('/home/aldanacepeda/Documentos/PhD/Algoritmo/resultados_MASP_allseq_HMMER_valores_GPI1.txt','a'))

for k,v  in result_gpi_dicc.items():
     final_clasification[k]["C-terminal_Evalue"]=v[3]
     final_clasification[k]["C-terminal_bitscore"]=v[2]

os.remove('temp2.fasta')
os.remove('temp3.fasta')

print("All C-terminal of the entered sequences were analyzed. ")


################################################################################
#EVALUACION DE SECUENCIAS DESDE BASE DE DATOS 

print("Analyzing and annotating MASP sequences") 

for k1, v1 in final_clasification.items():
 for k,v in result_sp_dicc.items():
    if k1 in result_sp_dicc.keys() and v[3]<=1E-5 and v[2]>=25:
         final_clasification[k]['N-terminal_type']="MASP"
    else:
         final_clasification[k1]['N-terminal_type']="Non-MASP"
 for k2,v2 in result_gpi_dicc.items():
    if k1 in result_gpi_dicc.keys() and v2[3]<=1E-7 and v2[2]>=30:
         final_clasification[k2]['C-terminal_type']="MASP"
    else:
         final_clasification[k1]['C-terminal_type']="Non-MASP"

for k,v in final_clasification.items():
    if v['N-terminal_type']==v['C-terminal_type']:
        v['Annotation']=v['N-terminal_type']
    elif v['N-terminal_type']=='MASP' or v['C-terminal_type']=='MASP':
        v['Annotation']='MASP-related'
    else:
         v['Annotation']="Non-MASP"

#CLASIFICACION DE SECUENCIAS QUIMERAS MASP SEGUN REGEX 
print("Searching for chimeric sequences")
quimeras={}

def SP_type(valor,cadena,tipo,regex):
    quimeras[valor]=[cadena]
    final_clasification[valor]['N-terminal_type']=tipo
    final_clasification[valor]['Chimera_RegEx_Match']=[]
    final_clasification[valor]['Chimera_RegEx_Match'].append(str(regex))
       
def GPI_type(valor,cadena,tipo,regex):
    quimeras[valor]=[cadena]
    final_clasification[valor]['C-terminal_type']=tipo
    final_clasification[valor]['Chimera_RegEx_Match']=[]
    final_clasification[valor]['Chimera_RegEx_Match'].append(str(regex))
    
def match_regex(regex, string, position):
    return re.search(regex, string[:position]) if position > 0 else re.search(regex, string[position:])

SP_patterns = [
    ('TcMuc', '(M|T|R|V|I)?(M|T|V|I|K|L)(M|T|K|S|N|A)(R|C|Y|W|S)(R|C|H)(V|L|P)(L|Q|M)(R|C|Y|S)(T|V|S|A|G|I)(L|P|M|V)(L|F|M|W)(M|V|I|L|A)(V|I|L|F)(A|G|P|S|V|T|I)(M|L|F|P|V|A)(C|Y|G|W|F|S|V)(C|W|F|Y|S|L)?(C|Y|G|F|R|S|T)'),
    ('TcMuc', 'M?(M|I)(T|S)(A|G|W|L|C)R(L|V)LCVLL(L|V)(P|S|L)(S|A|L)L(M|T)C(C|F)(S|L|F)(L|F|I)C(V|M)'),
    ('TS', 'M(P|S|F|L|Y)(L|R|Q|W)(S|P|L|Y|R|Q|H|D|N|C)(M|F|P|L|Y|H|V|A|E)(F|S|P|L|Y|C|V|I|T)(I|T|A|H|D|N|C|Y|R|F|S|P|G)(F|S|P|A|Y|T)(G|S|V|A|E|T)(M|L|F|V|A|I|G)(L|P|S|F|M|H|V|A|I)(L|P|F|Y|H|V|I)?(F|P|L|H|C|V|A|I)'),
    ('TS', '(F|T)(V|S)T?A?L?(P|S)I?LL(L|F)(L|A)L(C|R)(P|F)(S|G|N)(S|E)P?A?H?'),
    ('TS', 'MP(R|Q|L|S)(L|Y)LL(Y|F)(V|G)(L|F|V)V?L?(F|L)(F|L)(L|P)(L|F)(S|N)FC(P|A)(P|N)'),
]

GPI_patterns = [
    ('TcMuc', '(D|H|E|T)(G|S|V|D|A)(S|R|G|T|C|N|A)(L|F|P|S|V|R|H)(S|G|A|V|R|I|C|N)(S|G|I|N|A|T)(T|A|S|P|F)(A|V|S|P|T)(W|C|L|S|R)(A|R|V|L|M)(C|F|G|R|V|Y)(A|T|V|G|S|P)(P|L|S|R|T)(L|Q|P)(L|Q|V|A)?(L|P)(T|V|A|S|G)(T|A|V|G|L)(S|T|A|F|C|Y)(A|V|T|S|G)(L|M|V)(A|V|S|T|M)(Y|C|H)?'),
    ('TS', '(A|V|L|P|R|M|I|W)(L|V|P|Q|M|W)?(S|L|P|Q)?(V|L|P|R|M|W|S)?(S|L|M|P|I|F)?(F|H|L|P|S)?(L|S|M|I|F)?(I|F|V|L|P|S)?(F|M|L|V|G)?(V|L|P|R|Q|M)(S|L|M|I|F|W)(G|E|V|R|W)(L|P|M|I|F)(W|C|L|R|S)(G|E|V|R|M|W|C)(I|F|Y|V|L|S)(S|A|V|L|R|T|W)'),
    ('TS', '(D|E)(S|G|T|I|N)(T|S|A)(V|M)(H|Q|R)(A|G|R|E|L|V|W)(C|Y)(V|L|A)(F|S)(R|P|Q)(L|V|M)(W|L|S|F)LL?L?(L|M)?(L|P)LGLW(S|G|V)(T|I|A|S)(A|V|T)(A|T)(L|I)(C|R|S)?'),
    ('TS', '(D|E|V)(S|D|G|T)(S|F|T)(T|M|V)(R|P|H)(G|E)(G|S|D|A|N|H)(A|I|V|M|L)(S|F|P|C)(W|R|Q|G)(A|V|L|P)(L|F)?(M|V|L|P|S)?L?(L|V)?L?L?(L|M|F)?(G|L)LG(M|L|V|F)(W|C)(G|A|V)(I|S|F|V|L)(A|S|V|L|T)(A|E|V)(I|F|L)'),
    ('TS', 'KLILH(D|E)L(T|S|R)(V|G)DSTA(R|C|F)(V|L|A)(C|R|H)(V|A)S(R|C)(V|L|A)LL(L|R)(L|P)(L|V)GL(C|F)A(L|F)VA(V|F|I)'),
    ('TS', '(L|F)(V|P)DG(G|A)A(L|F|I)S(A|T)FSGG(G|R)(V|F|L)(F|L)LCA(C|W)ALLLH(P|V|L)(L|F)(L|F)(T|M)(I|A)(V|A)(F|S)(F|V)'),
    ('TS', 'SSEASG(S|P)DSSAHECVY(L|F)V(P|L)LVSLGLWA(L|F)A(V|L|F)FYGV'),
    ('TS', '(D|S)(S|G)(D|S)(F|D)(G|N)(D|G)(D|G)ASS(L|V)S(A|T)GGRNL(V|M)YV(C|S)(A|V)LVLHLLVLVYF'),
    ('TS', '(V|A)IQ(K|R)MMLDASERGS(R|I)FS(R|L)I(D|N)TLFLLTLLLL'),
    ('TS', 'PPKGQD(T|I)GDVSLTG(S|P)LVGVLWL(A|P)L(L|A)(G|W)(L|A)(L|W)GLSVLL'),
    ('TS', 'KG(K|S)GDGKKG(S|K)GDGSMRGGVSRLLLLGFCAFVVLY'),
]

for k, v in final_clasification.items():
    if v['Annotation'] == 'MASP-related':
        string = input_file[k]['Sequence']
        sp_match = next((SP_type(k, string, tipo, match.group()) for tipo, pattern in SP_patterns if (match := match_regex(pattern, string, 30))), None)
        gpi_match = next((GPI_type(k, string, tipo, match.group()) for tipo, pattern in GPI_patterns if (match := match_regex(pattern, string, -40))), None)
   
for k,v in final_clasification.items():
    if v['N-terminal_type']==v['C-terminal_type']:
        v['Annotation']=v['N-terminal_type']
    elif (v['N-terminal_type']=='MASP' and (v['C-terminal_type']=='TS' or  v['C-terminal_type']=='TcMuc')) or ((v['N-terminal_type']=='TS' or v['N-terminal_type']=='TcMuc') and v['C-terminal_type']=='MASP'):
        v['Annotation']='MASP-chimera'
    elif v['N-terminal_type']=='MASP' or v['C-terminal_type']=='MASP':
        v['Annotation']='MASP-related'
    else:
         v['Annotation']="Non-MASP"
         
Tabla_final_clasification= pd.DataFrame.from_dict(final_clasification, orient='index')
pd.set_option('display.max_rows', None) 
pd.set_option('display.max_columns', None) 
pd.set_option('display.width', 1000) 
pd.set_option('display.colheader_justify', 'center') 
pd.set_option('display.precision', 2) 
Tabla_final_clasification.columns=["ID","N-terminal_type","C-terminal_type","N-terminal_Evalue","N-terminal_bitscore","C-terminal_Evalue","C-terminal_bitscore","Chimera_RegEx_Match","Annotation"]
Tabla_final_clasification["ID"]=Tabla_final_clasification["ID"].apply(lambda x: ', '.join(map(str,x)))

print("All sequences were classified")
print("Selecting sequences according to hierarchical ranking")
noindatabase = {}
mydatabase = {}

# Iterar sobre cada línea de TCCcds una sola vez
for index, row in Tabla_final_clasification.iterrows():
    id_mydatabase = row['ID'].strip()  # Eliminar espacios y tabulaciones
    contig=id_mydatabase.split("_")[1]
    sp = row['N-terminal_type']
    gpi = row['C-terminal_type']
    anotation = row['Annotation'].strip()  # Limpiar anotación
    n_evalue= row['N-terminal_Evalue']
    n_bitscore=row['N-terminal_bitscore']
    c_evalue= row['C-terminal_Evalue']
    c_bitscore=row['C-terminal_bitscore']
    regex_match=row['Chimera_RegEx_Match']
    description=""

    # Obtener idfinal usando rsplit
    id_parts = id_mydatabase.rsplit('_', 1)
    idfinal = id_parts[0]
    sequence = input_file[id_mydatabase]['Sequence']
    
    # Limpiar id_mydatabase
    id_mydatabase = id_mydatabase.replace('\t', '').replace(' ', '')

    # Obtener position, positfinal, y sequence usando input_file
    position = int(input_file[id_mydatabase]['Pbi'])
    positfinal = int(input_file[id_mydatabase]['Pbf'])
    # Determinar strand
    if position > positfinal:
        strand = '-'
    elif position < positfinal:
        strand = '+'
    else:
        strand = ''  # Manejar otros casos según sea necesario

    # Filtrar según condiciones de anotación
    if "MASP" in anotation and not "Non" in anotation and not "related" in anotation and (idfinal not in mydatabase.keys()):
            mydatabase[idfinal] = [id_mydatabase,contig ,anotation, sp, gpi,n_evalue,n_bitscore,c_evalue,c_bitscore,regex_match, position, positfinal, strand, description,sequence]

# Iterar sobre cada línea de TCCcds una sola vez
for index, row in Tabla_final_clasification.iterrows():
    id_mydatabase = row['ID'].strip()  # Eliminar espacios y tabulaciones
    contig=id_mydatabase.split("_")[1]
    sp = row['N-terminal_type']
    gpi = row['C-terminal_type']
    anotation = row['Annotation'].strip()  # Limpiar anotación
    n_evalue= row['N-terminal_Evalue']
    n_bitscore=row['N-terminal_bitscore']
    c_evalue= row['C-terminal_Evalue']
    c_bitscore=row['C-terminal_bitscore']
    regex_match=row['Chimera_RegEx_Match']
    description=""

    # Obtener idfinal usando rsplit
    id_parts = id_mydatabase.rsplit('_', 1)
    idfinal = id_parts[0]
    sequence = input_file[id_mydatabase]['Sequence']

    
    # Limpiar id_mydatabase
    id_mydatabase = id_mydatabase.replace('\t', '').replace(' ', '')

    # Obtener position, positfinal, y sequence usando input_file
    position = int(input_file[id_mydatabase]['Pbi'])
    positfinal = int(input_file[id_mydatabase]['Pbf'])

    # Determinar strand
    if position > positfinal:
        strand = '-'
    elif position < positfinal:
        strand = '+'
    else:
        strand = ''  # Manejar otros casos según sea necesario

    # Filtrar según condiciones de anotación
    if "chimera" in anotation and (idfinal not in mydatabase.keys()):
            mydatabase[idfinal] = [id_mydatabase,contig,anotation, sp, gpi,n_evalue,n_bitscore,c_evalue,c_bitscore,regex_match, position, positfinal, strand,description,sequence]

# Iterar sobre cada línea de TCCcds una sola vez
for index, row in Tabla_final_clasification.iterrows():
    id_mydatabase = row['ID'].strip()  # Eliminar espacios y tabulaciones
    contig=id_mydatabase.split("_")[1]
    sp = row['N-terminal_type']
    gpi = row['C-terminal_type']
    anotation = row['Annotation'].strip()  # Limpiar anotación
    n_evalue= row['N-terminal_Evalue']
    n_bitscore=row['N-terminal_bitscore']
    c_evalue= row['C-terminal_Evalue']
    c_bitscore=row['C-terminal_bitscore']
    regex_match=row['Chimera_RegEx_Match']
    description=""

    # Obtener idfinal usando rsplit
    id_parts = id_mydatabase.rsplit('_', 1)
    idfinal = id_parts[0]
    
    # Limpiar id_mydatabase
    id_mydatabase = id_mydatabase.replace('\t', '').replace(' ', '')

    # Obtener position, positfinal, y sequence usando input_file
    position = int(input_file[id_mydatabase]['Pbi'])
    positfinal = int(input_file[id_mydatabase]['Pbf'])
    sequence = input_file[id_mydatabase]['Sequence']

    # Determinar strand
    if position > positfinal:
        strand = '-'
    elif position < positfinal:
        strand = '+'
    else:
        strand = ''  # Manejar otros casos según sea necesario

    # Filtrar según condiciones de anotación
    if "related" in anotation and (idfinal not in mydatabase.keys()):
             mydatabase[idfinal] = [id_mydatabase,contig,anotation, sp, gpi,n_evalue,n_bitscore,c_evalue,c_bitscore,regex_match, position, positfinal, strand,description,sequence]

# Iterar sobre cada línea de TCCcds una sola vez
for index, row in Tabla_final_clasification.iterrows():
    id_mydatabase = row['ID'].strip()  # Eliminar espacios y tabulaciones
    contig=id_mydatabase.split("_")[1]
    sp = row['N-terminal_type']
    gpi = row['C-terminal_type']
    anotation = row['Annotation'].strip()  # Limpiar anotación
    n_evalue= row['N-terminal_Evalue']
    n_bitscore=row['N-terminal_bitscore']
    c_evalue= row['C-terminal_Evalue']
    c_bitscore=row['C-terminal_bitscore']
    regex_match=row['Chimera_RegEx_Match']
    description=""

    # Obtener idfinal usando rsplit
    id_parts = id_mydatabase.rsplit('_', 1)
    idfinal = id_parts[0]
    sequence = input_file[id_mydatabase]['Sequence']
    
    # Limpiar id_mydatabase
    id_mydatabase = id_mydatabase.replace('\t', '').replace(' ', '')

    # Obtener position, positfinal, y sequence usando input_file
    position = int(input_file[id_mydatabase]['Pbi'])
    positfinal = int(input_file[id_mydatabase]['Pbf'])

    # Determinar strand
    if position > positfinal:
        strand = '-'
    elif position < positfinal:
        strand = '+'
    else:
        strand = ''  # Manejar otros casos según sea necesario

    # Filtrar según condiciones de anotación
    if "Non" in anotation and (idfinal not in mydatabase.keys()):
            mydatabase[idfinal] = [id_mydatabase,contig,anotation, sp, gpi,n_evalue,n_bitscore,c_evalue,c_bitscore,regex_match, position, positfinal, strand,description,sequence]



# Reiniciar el puntero de TCCcds al inicio si es un archivo
key_to_delete = []
# Función para calcular la longitud y hacer la corrección MASP_related_corrected
def corregir_masp(orf, next_orfs):
    for next_orf in next_orfs:
        try:
         if mydatabase[orf][3] == "MASP" and mydatabase[orf][4] == "Non-MASP":
            if mydatabase[next_orf][3] == "Non-MASP" and mydatabase[next_orf][4] == "MASP":
                long = abs(mydatabase[next_orf][10] - mydatabase[orf][11])
                if long < 600:
                    mydatabase[orf][11] = mydatabase[next_orf][11]
                    mydatabase[orf][2] = "MASP-related"
                    mydatabase[orf][13]= "Concatenated"
                    key_to_delete.append(next_orf)
        except:
            pass

# Iterar sobre los ítems de mydatabase
for k, v in mydatabase.items():
    # Extraer la parte inicial del nombre y el número
    initial = (v[0].rsplit('_', 2))[0]
    name = (v[0].rsplit('_', 2))[1]
    kreal = int(name)  # Convertir el nombre a un entero
    orf = f"{initial}_{kreal}"  # Crear el ORF base
    next_orfs = [f"{initial}_{kreal + i}" for i in [1,2,3,4,5,6]]  # Crear una lista de ORFs subsecuentes
    # Corregir MASP_related_corrected para cada ORF y sus siguientes
    corregir_masp(orf, next_orfs)  # Llamada a la función para corregir los ORFs


# Eliminar keys de mydatabase según key_to_delete
mydatabase = {k: v for k, v in mydatabase.items() if k not in key_to_delete}

print("Finished selecting sequences according to hierarchical ranking")

tablamydatabase= pd.DataFrame.from_dict(mydatabase, orient='index')
pd.set_option('display.max_rows', None) 
pd.set_option('display.max_columns', None) 
pd.set_option('display.width', 1000) 
pd.set_option('display.colheader_justify', 'center') 
pd.set_option('display.precision', 2)
tablamydatabase.columns=['ID','Contig','Prediction_algorithm', 'N-terminal_type', 'C-terminal_type','N-terminal_Evalue','N-terminal_bitscore','C-terminal_Evalue','C-terminal_bitscore','Chimera_RegEx_Match' ,'Pbi', 'Pbf', 'Strand','Description','Sequence']

# Crear una copia del DataFrame y eliminar la última columna
tablamydatabase_copy = tablamydatabase.iloc[:, :-1]
tablamydatabase_copy.to_csv(f'{directory}/Table_{out3}_Classification.csv', index=False)

# Contar las distintas ocurrencias de cada 'Prediction_algorithm'
algorithm_counts = tablamydatabase['Prediction_algorithm'].value_counts()

# Agrupar por el algoritmo de predicción
grouped = tablamydatabase.groupby('Prediction_algorithm')

# Crear archivos FASTA separados por cada grupo de algoritmo de predicción
for algorithm, group in grouped:
    # Crear un archivo FASTA para cada algoritmo
    filename = f"{directory}/{algorithm}_{out3}_sequences.fasta"
    file_path = os.path.join(directory, filename)
    with open(filename, 'w') as fasta_file:
        for index, row in group.iterrows():
            k=row['ID']
            sequence= input_file[k]['Sequence']
            # Crear el encabezado de la secuencia
            header = f">{row['ID']}|Prediction_algorithm={row['Prediction_algorithm']}|N-terminal_type={row['N-terminal_type']}|C-terminal_type={row['C-terminal_type']}|Pbi={row['Pbi']}|Pbf={row['Pbf']}|Strand={row['Strand']}"
            # Escribir el encabezado y la secuencia al archivo FASTA
            fasta_file.write(f"{header}\n{row['Sequence']}\n")

print("FASTA files generated by prediction of the algorithm.")

# Filtrar las filas con Prediction_algorithm igual a 'NO_MASP'
filtered_table = tablamydatabase[tablamydatabase['Prediction_algorithm'] != 'Non-MASP']
# Preparar los datos en formato GFF
gff_lines = []
for index, row in filtered_table.iterrows():
    # Crear una línea en formato GFF para cada fila
    attributes = [f'ID={row["ID"]}', f'Prediction_algorithm={row["Prediction_algorithm"]}']
    
    # Agregar Description solo si tiene un valor
    if pd.notna(row['Description']) and row['Description'] != '':
        attributes.append(f'Description={row["Description"]}')
    
    # Crear la línea GFF
    gff_line = '\t'.join([
        row['ID'],                        # Nombre del cromosoma o secuencia
        row['Prediction_algorithm'],      # Origen (fuente) del feature
        'gene',                           # Tipo de feature (en este caso, gene)
        str(row['Pbi']),                  # Posición de inicio del feature
        str(row['Pbf']),                  # Posición de fin del feature
        '.',                              # Puntuación del feature (puede ser '.')
        row['Strand'],                    # Cadena del feature (puede ser '+' o '-')
        '.',                              # Fase del feature (puede ser '.' si no aplica)
        ';'.join(attributes)              # Atributos adicionales
    ])
    gff_lines.append(gff_line)
# Escribir las líneas GFF en un archivo
output_file = f"{directory}/MASP_{out3}_sequences.gff"
with open(output_file, 'w') as f:
    for line in gff_lines:
        f.write(line + '\n')

print(f"GFF file generated: {output_file}")

# Definir el nombre del archivo de salida
output_file = f'{directory}/README_{out3}.txt'

# Abrir el archivo de salida en modo escritura (esto crea el archivo si no existe o lo sobrescribe si ya existe)
with open(output_file, 'w', encoding='UTF-8') as f:
    # Escribir el nombre del archivo proporcionado por el usuario en el archivo de salida
    f.write(f'============Disruptomics v.1.0============\n\nDirectory selected:{directory}\nName of the file provided: {a}\nStrain anlyzed:{out3}\nAlgorithm Prediction:\n{algorithm_counts}\n')

print(f'The information has been stored in {output_file}')
