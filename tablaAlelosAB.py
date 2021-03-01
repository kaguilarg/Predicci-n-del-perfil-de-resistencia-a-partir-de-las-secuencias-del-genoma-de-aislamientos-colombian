import pymysql, os
import numpy as np
np.set_printoptions(formatter='int')
#-------------------------creado de tabla general real y en tipo lista
filenames = os.listdir('/vault2/homehpc/kaguilarg/AB/RGI _AB/')
open_table = open('tablageneral.txt', 'w')#tabla resumen
tablafinal= open('tablafinal.txt','w')
#limites filtrado de variables loose
limbs,limid,limcov=50,50,80
#variable inicial con nombresepa, adn, prot,%id, cov
general=[]
filtrados,c,sp=0,1,0
#----------------------------toma de datos de interes [nombre, gen, bitscore, clasificacion,ADN,PROT,%id,cov]
for filename in filenames:
    if filename[-3:] == "txt":
        #print("Reading file " + filename)#revisar los archivos de lectura
        file = os.path.abspath('/vault2/homehpc/kaguilarg/AB/RGI _AB') + "/" + filename
        #print(file)
        arch = open(file,'r')
        line = arch.readlines()
        #print(line)
        arch.close()
        c=1    
        for l in line:
            if c != 1 :
                columns = l.split("\t")
                if str(columns[5])=="Loose":
                    if float(columns[7]) > limbs and float(columns[9]) > limid and float(columns[20]) > limcov:
                        linea =filename[:-4] + "\t"+columns[8] + "\t"+columns[7]+ "\t"+columns[5]+"\t"+ columns[17] +"\t"+columns[18] +"\t"+columns[9]+"\t"+columns[20]+"\n"
                        #print(len(columns[17]))
                        open_table.write(str(linea))#txt con tabla general
                        linea=[filename[:-4],columns[8],columns[7],columns[5],columns[17],columns[18],columns[9],columns[20],columns[6]]
                        general.append(linea)
                    else:
                        filtrados+=1
                else:
                    linea =filename[:-4] + "\t"+columns[8] + "\t"+columns[7]+ "\t"+columns[5]+"\t"+ columns[17] +"\t"+columns[18] +"\t"+columns[9]+"\t"+columns[20]+"\n"
                    #print(len(columns[17]))
                    open_table.write(str(linea))#txt con tabla general
                    linea=[filename[:-4],columns[8],columns[7],columns[5],columns[17],columns[18],columns[9],columns[20],columns[6]]
                    general.append(linea)
                    sp+=1
            c +=1
        #print(c)#al descomentarlo muestra la cantidad de lineas extraidas de cada archivo menos 1 del encabezado.
print(str(len(general))+' de reportes filtrando '+str(filtrados)+' strict y perfect:  '+str(sp))#tamano de variable general
open_table.close()#txt con la tabla general.
#-----------------------------rango de los loose, %ID, coverage, pasbs
#3=clasificacion, 2=bit score,6=%id, 7=coverage
loose=[]
for gene in general:
    if gene[3]== "Loose":
        loose.append([gene[3],gene[2],gene[6],gene[7],gene[8]])
print(str(len(loose))+' cantidad de looses usando filtros de puntuacion')#tamano de los loos
Mid,Mbs,Mpbs,Mcov=0,0,0,0
mid,mbs,mcov,mpbs=100,1000,1000,1000
for l in loose:
    #print(l[1])
    if float(l[1])>float(Mbs):
        Mbs=l[1]
    if float(l[1])<float(mbs):
        mbs=l[1]
    if float(l[2])>float(Mid):
        Mid=l[2]
    if float(l[2])<float(mid):
        mid=l[2]
    if float(l[3])>float(Mcov):
        Mcov=l[3]
    if float(l[3])<float(mcov):
        mcov=l[3]
    if float(l[4])>float(Mpbs):
        Mpbs=l[4]
    if float(l[4])<float(mpbs):
        mpbs=l[4]
print("maximo bs "+str(Mbs)+'\t'+"minimo BS "+str(mbs)+'\n'+"maximo id "+str(Mid)+'\t'+'\t'+"minimo id "+str(mid)+'\n'+"maximo cov "+str(Mcov)+'\t'+"minimo cov "+str(mcov)+'\n'+"maximo BpbS "+str(Mpbs)+'\t'+"minimo BpbS "+str(mpbs))#rangos maximo y minimo
#-------------------------------------------------------crear la lista de aislamientos
anterior,actual='',''
aislamientos=[]
for gene in general:
    actual=gene[0]
    if actual != anterior:
        aislamientos.append(actual)
        #print(actual)
    anterior=actual
print('aislamientos = '+str(len(aislamientos)))#imprimir la cantidad de aislamientos
#--------------------------------------------------------------# lista de genes sin repeticion
genes=[]
cuenta,c=0,1
for gene in general:
    c=1
    for ge in genes:
        if gene[1]==ge:
            c=0
    if c == 1:
        cuenta+=1
        genes.append(gene[1])
print(str(cuenta)+' genes entre todas las sepas.')#cantidad de genes sin repeticion
#print(len(genes))
#---------------------------------------------creacion de lista[aislamiento,gen,repeticion]
aisgenesrep=[]
cuenta,cuentagen=0,0
for ais in aislamientos:
    for g in genes:
        cuenta+=1
        cuentagen=0
        for gene in general:
            if ais == gene[0] and g == gene[1]:
                cuentagen+=1
        aisgenesrep.append([ais,g,cuentagen])
#----------------------------impresion variable genes rep
#for rep in aisgenesrep:
#   print(rep)
#------------------------------unificacion en variable [aislamiento,genes,repeticion]
listagenes=[]#variable que incluye las repeticiones maximas de los genes 
cuentagen,prueba=0,0
for g in genes:
    cuentagen=0
    for agr in aisgenesrep:
        if g == agr[1] and agr[2] > cuentagen:
            cuentagen=agr[2]
    #print(cuentagen)
    for cg in range(0,cuentagen):
        prueba+=1
        if cg == 0:
            listagenes.append(g)
        else:
            listagenes.append(g+'_cp'+str(cg))
#----------------------------impresion variable lista de genes
#for lg in listagenes:
#   print(lg)
#print(prueba)
print(str(len(listagenes))+' genes incluyendo las copias.')
#print(cuenta)
# -------------------------------------------------------------# generacion de archivos txt con las secuencias de los genes y de encabezado el nombre de la sepa por cada gen.
gen=[]
comilla=0
string=''
for ge in genes:
    gen=[]
    string=str(ge)
    for gene in general:
        if ge == gene[1]:
            gen.append('>'+str(gene[0])+'\n')
            gen.append(str(gene[5])+'\n')
    comilla=string.find("'")
    string=string.replace("'","%",comilla)
    string=string.replace("(","#")    
    string=string.replace(")","*")
    string=string.replace(" ","_")
    #print(string)
    txt = open('gprot/'+string+'.txt','w')
    txt.writelines(gen)
    txt.close()
#----------------------------------------------ordenar por tamano las secuencias generando archivos.orde
#cuenta=0
protfiles=os.listdir('gprot')
for ge in protfiles:
   if ge[-3:]=='txt'and ge[1:]!='.':
       os.system('uclust --sort gprot/'+ge+' --output gpord/'+ge[:-4]+'.orde')
#print(cuenta) 
#---------------------------------------------generar .alel de los genes alelizados
ordfiles=os.listdir('gpord')
for ge in ordfiles:#limitar el .dir
   if ge[-4:]=='orde' and ge[1:]!='.':
       os.system('uclust --input gpord/'+ge+' --uc gpalel/'+ge[:-5]+'.alel --id 0.95')
#---------------------------------------------obtener datos de los archivos alelizados en variable alelos
alelfiles=os.listdir('gpalel')
#for al in alelfiles:
#    print(al)
alelos=[]
linea=[]
string=''
comilla=0
for ge in alelfiles:
    if ge[-4:]=='alel' and ge[1:]!='.':
        string=str(ge[:-5])
        comilla=string.find("%")
        string=string.replace("%","'",comilla)
        string=string.replace("#","(")    
        string=string.replace("*",")")
        string=string.replace("_"," ")
        linea=[string]
        alel=open('gpalel/'+ge,'r')
        lineas=alel.readlines()
        for l in lineas:
            columnas=l.split('\t')
            col=l.split()
            if col[0]!='#' and columnas[0]!='C':
                linea.append([columnas[8],(int(columnas[1])+1)])
        alelos.append(linea)#variable lista de genes en la que cada gen es una lista de los aislamientos y cada aislamiento es una lista con el nombre del alelo y el numero de este.
print(str(len(alelos))+'tamano de la lista alelos que debe ser la misma que la cantidad de genes sin repeticion')
#for i in alelos:#imprimir variable alelos
#    print(i[1][1])
#---------------------------------------------prueba la realizacion de la variable alelos
cuenta =0
#for alel in alelos:
   #for j in range(1,len(alel)):
       #cuenta+=1
#print(cuenta)#debe ser igual al numero de reportes 
#--------------------------------------------llenado tabla final-------------------------
alelcop=[]#contiene en cada linea primero el nombre del gen con presencia unica de los aislamientos,si hay copias entonces hay otra linea con el nombre al final _cp# que contiene el nombre de las sepas que tienen este gen copiado.
vecrep=[]
cuentasepa=0
rep=0
linea=[]
ncop=0
for alel in alelos:
    vecrep=[]
    for alei in alel[1:]:
        cuentasepa=0
        rep=0
        for alej in alel[1:]:
            #print('comparando '+str(alej[0])+' '+str(alei[0]))
            if str(alej[0])==str(alei[0]):
                cuentasepa+=1
        for acp in vecrep:
            if alei[0]==acp[0]:
                rep=1
        if rep ==0:
            vecrep.append([alei[0],cuentasepa])#vector temporal que por cada gen almacena la cantidad de veces que se repite cada aislamiento
#vecrept=[]
    #print('sepa '+str(alel[0]))
    #print(len(alel))
    #print(str(len(vecrep))+'tamano vector repeticiones')#informacion de cuantas sepas tienen cada gen
    maximo=1
    ncop=0
    while maximo!=0:
        if ncop == 0:
            linea=[alel[0]]
        else:
            linea=[str(alel[0])+'_cp'+str(ncop)]
        ncop+=1
        for vr in vecrep:
            if vr[1]>=1:
                linea.append([str(vr[0])])
                vr[1]-=1
        if maximo >0:            
            alelcop.append(linea)
            linea=[]
        maximo=0
        for vr in vecrep:
            if vr[1]>maximo:
                maximo=vr[1]        
        #print(maximo)#cuando un gen tiene copias debe aparecer 1 y luego 0, si no tiene copias enctonces solo 0
        #print(alelcop)
        #print(vecrep)#para ver el cambio de las cantidades de aislamientos relacionadas con el gen
#-----------------------------------retomado de la informacion alelica
prueba=0
for ac in alelcop:
    for als in alelos:
        if ac[0]==als[0] or als[0]==ac[0][:-4]:
            for a in ac[1:]:
                for al in als[1:]:
                    if str(a[0])==str(al[0]):
                        a.append(al[1])
                        #print(str(a))#para ver como es agregado el alelo a la variable alelcop                        
                        al[0]='*'
                        prueba+=1
                        break
print(str(prueba)+' comparaciones acertadas, debe ser igual al numero de reportes')
#print(str(alelos))#no debe quedar ningun nombre de sepa en esta variable
#-----------------------------------impresion genes incluyendo las copias, el 
#for ac in alelcop:
    #print(ac[0])
print(str(len(alelcop))+' largo de la variable alelos en las que se incluyen las copias')# debe ser igual genes incluyendo copias de arriba
#-----------------------------creado de la tabla en ceros
#-----filtrado de las letras para los nombres de los aislamientos
#****************************funcion para filtrar letras de los nombres de los aislamientos
def f(x):
    num=['0','1','2','3','4','5','6','7','8','9']
    if(x in num):
        return True
    else:
        return False
listagenes.insert(0,'aislamiento')
tabla = np.zeros([len(aislamientos),len(alelcop)+1])
tabla = tabla-1
tabla= list(tabla)
tabla.insert(0,listagenes)
#print(tabla)
#---------------------------------------llenado de la columna de aislamientos
for i in range(1,len(tabla)):
    tabla[i][0]=int(filter(f,aislamientos[i-1]))
#print(tabla)
#----------------------------------
#------------------------------------- llenado de la tabla inicializada en ceros
#------------------------- verificacion de los nombres de los genes
#for i in range(0,len(alelcop)):
    #print(str(alelcop[i][0])+'......'+str(listagenes[i]))
#--------------------------
prueba=0
for gen in range(1,len(alelcop)+1):
    for ais in range(1,len(aislamientos)+1):
        for alel in alelcop:#debe ser alelcop no alelos
            if str(alel[0])==str(listagenes[gen]) or str(alel[0])==(str(listagenes[gen])[:-4]):
                for alei in alel[1:]:
                    if alei[0]==aislamientos[ais-1]:
                        tabla[ais][gen]=int(alei[1])
                        alei[0]='*'
                        #print(alel)
                        prueba+=1
print(str(prueba)+'cantidad de datos escritos en la tabla: debe ser igual a la cantidad de reportes')
salida=''
l=''

for i in range(0,len(tabla)):# colocar nombre completo
    l=''
    for j in range(0,len(tabla[1])):
        if j!=(len(tabla[1])-1):
            if i == 0:
                l+=(str(tabla[i][j])+'\t')
            else:
                if j == 0:
                    l+=(str(aislamientos[i-1])+'\t')
                else:
                    l+=(str(tabla[i][j])+'\t')
        else:
            l+=(str(tabla[i][j])+'\n')
        #print(int(tabla[i][j]))
    salida+=l
tablafinal.writelines(str(salida))
tablafinal.close()
