import numpy as np
import random
import matplotlib.pyplot as plt
import math
import pandas as pd
import os
import xlsxwriter
from os import system
import PyInquirer as inquirer
from time import perf_counter



######## CREACION DE LOS GENERADORES ########
def big_bang():
    for i in range(p):
        # Genero cromosomas en binario y en decimal
        cromosomas_bin.append(genera_parque())
        cromosomas.append(evalua_parque(cromosomas_bin[i]))


######## GENERA CROMOSOMAS EN BINARIO ########
def genera_parque():
    parque = [[0 for i in range(celdas)] for i in range(celdas)]    #Cambio
    for i in range(cant_generadores_inicial):                           #Cambio
        while True:                                             #Cambio
            numX = np.random.randint(0, celdas)                 #Cambio
            numY = np.random.randint(0, celdas)                 #Cambio
            if parque[numY][numX] == 0:                         #Cambio
                parque[numY][numX] = 1                          #Cambio
                break                                           #Cambio

    return parque


######## CALCULA POTENCIA DE PARQUE ########
def evalua_parque(parque):
    parque_aux = [[0 for i in range(celdas)] for i in range(celdas)]
    for i in range(celdas):
        for j in range(celdas):
            if parque[i][j] == 1:
                if i == 0:
                    for k in range(len(velocidades)):
                        if u0[j] == velocidades[k]:
                            parque_aux[i][j] = potencias[k]
                else:
                    for k in range(i-1, -1, -1):
                        if parque[k][j] == 1:
                            x = (i - k)*dist_min      # Distancia entre molinos en metros.
                            ux[j] = round(u0[j]*(1 - 2*a/(1 + alfa*x/r1)**2))
                            for l in range(len(velocidades)):
                                if ux[j] == velocidades[l]:
                                    parque_aux[i][j] = potencias[l]
                                #else:
                                #   parque_aux[i][j] = 0
                            break
                    else:
                        for l in range(len(velocidades)):
                            if u0[j] == velocidades[l]:
                                parque_aux[i][j] = potencias[l]

    return parque_aux


######## CALCULA FUNCION OBJETIVO ########
def calcula_f_obj():
    for i in range(p):
        aux = 0
        for j in range(celdas):
            aux += sum(cromosomas[i][j])
        f_obj[i] = aux


######## CALCULA FITNESS ########
def calcula_fitness():
    for i in range(p):
        fitness[i] = f_obj[i]/sum(f_obj)


######## SELECCION Y CROSSOVER ########
def selec_cross():
    ruleta = calcula_ruleta()

    #Convierto la lista en un array
    aux = np.array(f_obj)

    #Busco los dos cromosomas maximos
    max1 = aux.argsort()[-1]
    max2 = aux.argsort()[-2]

    #Coloco los maximos en las primeras posiciones
    cromosomas_bin[0] = cromosomas_bin[max1]
    cromosomas_bin[1] = cromosomas_bin[max2]

    for i in range(2, p, 2):
        c1,c2 = tiradas(ruleta)
        c1,c2 = crossover(c1,c2)
        cromosomas_bin[i] = c1
        cromosomas_bin[i+1] = c2


######## RULETA ########
def calcula_ruleta():
    #Calculo la frecuencia acumulada de cada cromosoma
    frec_acum = []
    frec_acum.append(fitness[0])
    for i in range(1,p):
        acumulado = frec_acum[i - 1] + fitness[i]
        frec_acum.append(acumulado)

    return frec_acum


######## TIRADA DE RULETA ########
def tiradas(ruleta):
    padres = []
    for m in range(2): #Ciclo de 2 porque necesito dos cromosomas
        frec = random.uniform(0,1)
        #Para hacer Crossover utilizamos la frecuencia acumulada,
        # basada en los fitness de los cromosomas.
        for i in range(p):
            if ruleta[i] > frec:
                padres.append(cromosomas_bin[i])
                break
            
    return padres[0], padres[1]


######## CROSSOVER ########
def crossover(c1, c2):
    c = np.random.randint(0, 101)
    if c <= cr:
        c1 = evalua_parque(c1)      #Cambio
        c2 = evalua_parque(c2)      #Cambio
        aux_fila_1 = mejores_filas(c1)
        aux_fila_2 = mejores_filas(c2)
        aux_columna_1 = mejores_columnas(c1)
        aux_columna_2 = mejores_columnas(c2)
        c1 = aux_fila_1 + aux_fila_2
        cant_generadores_c1 = contar_generadores(c1)     #Cambio
        if cant_generadores_c1 > cant_generadores_max:       #Cambio
            c1 = corregir_parque(c1, cant_generadores_c1)     #Cambio
        c2 = np.transpose(aux_columna_1 + aux_columna_2)
        cant_generadores_c2 = contar_generadores(c2)     #Cambio
        if cant_generadores_c2 > cant_generadores_max:       #Cambio
            c2 = corregir_parque(c2, cant_generadores_c2)     #Cambio
        c1 = potencia_to_binario(c1)    #Cambio
        c2 = potencia_to_binario(c2)    #Cambio

    return c1, c2


######## SELECCIONAR FILAS ########
def mejores_filas(cromosoma):
    lista_aux = []
    cromosoma_aux = []
    for i in range(celdas):
        lista_aux.append([i, sum(cromosoma[i])])        #Cambio
    cant_filas = int(celdas/2)                          #Cambio
    lista_aux = sorted(lista_aux, reverse = True, key = lambda w: w[1])[:cant_filas]    #Cambio
    for i in range(len(lista_aux)):        #Cambio
        for j in range(celdas):
            if lista_aux[i][0] == j:    #Cambio
                cromosoma_aux.append(cromosoma[j])  #Cambio
        
    return cromosoma_aux


######## SELECCIONAR COLUMNAS ########
def  mejores_columnas(cromosoma):
    cromosoma_transpuesta = np.transpose(cromosoma)

    return mejores_filas(cromosoma_transpuesta)


######## MUTACION ########
def mutacion():
    for i in range(p):
        num = np.random.randint(0, 101)
        if num <= m:
            corte_y = np.random.randint(0, celdas)
            corte_x = np.random.randint(0, celdas)
            if i != 0 or i != 1:
                if cromosomas_bin[i][corte_y][corte_x] == 0:
                    if contar_generadores(cromosomas_bin[i]) < cant_generadores_max:
                        cromosomas_bin[i][corte_y][corte_x] = 1

                else:
                    cromosomas_bin[i][corte_y][corte_x] = 0


######## POTENCIA A BINARIO ######## 
def potencia_to_binario(cromosoma):         # Nueva funcion!!!!!
    for i in range(celdas):
        for j in range(celdas):
            if cromosoma[i][j] != 0:
                cromosoma[i][j] = 1

    return cromosoma


######## BINARIO A POTENCIA ######## 
def binario_to_potencia():
    for i in range(p):
        cromosomas[i] = evalua_parque(cromosomas_bin[i])


######## CONTAR GENERADORES ########
def contar_generadores(parque):             # Nueva funcion!!!!!
    nro_generadores = 0
    for i in range(celdas):
        for j in range(celdas):
            if parque[i][j] != 0:
                nro_generadores += 1

    return nro_generadores


######## CORREGIR PARQUES ########
def corregir_parque(parque, cant_generadores_parque):                # Nueva funcion!!!!!
    nro_generadores_a_borrar = cant_generadores_parque - cant_generadores_max
    for k in range(nro_generadores_a_borrar):
        aux = [f_obj[0]*2, 0, 0]
        for i in range(celdas):
            for j in range(celdas):
                if parque[i][j] > 0 and parque[i][j] < aux[0]:
                    aux = [parque[i][j], i, j]
        parque[aux[1]][aux[2]] = 0

    return parque


######## BUSCAR MEJOR PARQUE ########
def mejor_parque():
    lista_aux = []
    for i in range(len(lista_mejores_parques)):
        aux = 0
        for j in range(celdas):
            aux += sum(lista_mejores_parques[i][j])
        lista_aux.append(aux)
    
    return lista_mejores_parques[lista_aux.index(max(lista_aux))], max(lista_aux)




############# PROGRAMA PRINCIPAL #############
while True:
    #Ingreso de paraceldas
    system("cls")
    print()
    print('---------------------------------------------')
    print('              ALGORITMO GENETICO             ')
    print('---------------------------------------------')
    print()

    p =  50
    g =  100
    m =  20
    cr = 75
    #p =  int(input("Ingrese tamaño de la poblacion     : "))    # 50
    #g =  int(input("Ingrese cantidad de generaciones   : "))    # 1500
    #m =  int(input("Ingrese tasa de mutacion (%)       : "))    # 20
    #cr = int(input("Ingrese tasa de crossover (%)      : "))    # 75
    
    celdas = 10         # Celdas de largo y ancho del parque. Es simetrico.
    dist_min = 94           # Distancia entre generadores. Tamaño de celdas
    cant_generadores_max = 25           # Cantidad maxima de generadores.           #Cambio
    cant_generadores_inicial = np.random.randint(1, cant_generadores_max + 1)        #Cambio
    cromosomas = []
    cromosomas_bin = []

    viento_promedio = 20
    u0 = list(np.random.choice(np.arange(viento_promedio - 1, viento_promedio + 2), p=[0.2, 0.6, 0.2]) for i in range(10))    # Velocidad del viento sin turbulencia.
    ux = list(np.zeros(celdas))         # Velocidad del viento con turbulencia.
    a = 1/3                             # Coeficiente de inducción axial.
    alfa = 0.05                         # Coeficiente de arrastre.
    gamma = 2                           # Constante de proporcionalidad.
    rr = 23.5                           # Radio de la turbina.
    r1 = rr*gamma                       # Radio de la estela.

    velocidades = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]                        # Posibles velocidades del viento.
    potencias = [0, 0, 0, 0, 0, 53, 106, 166, 252, 350, 464, 560, 630, 660, 660, 660, 660, 660, 660, 660, 660, 660, 660, 660, 660, 660] # Posibles potencias que se pueden generar.

    f_obj = list(np.zeros(p))
    fitness = list(np.zeros(p))

    #Declaracion de listas
    lista_min = []
    lista_max = []
    lista_prom = []
    lista_mejores_parques = []

    #Se crea la poblacion inicial
    big_bang()

    #Se calcula la funcion objetivo y el fitness de dicha poblacion
    calcula_f_obj()
    calcula_fitness()
    
    print("Poblacion:", p)
    print("Cantidad de generaciones:", g)
    print()
    print("Viento promedio:", viento_promedio, "m/s")
    print("Cantidad de generadores iniciales:", cant_generadores_inicial)
    print()
    print("Coeficiente de induccion axial 'a':", float("{0:.4f}".format(a)) )
    print("Coeficiente de arrastre 'alfa':", alfa)
    print("Radio de las turbinas 'rr':", rr, "metros")
    print("Radio de las estelas 'r1':", r1, "metros")
    print("Constante de proporcionalidad 'gamma':", gamma)
    print()


    # Por cada generación se ejecuta...
    start_time = perf_counter()
    for i in range(g):
        print(i)
        #Se llenan las listas con los diferentes resultados de la ejecucion
        lista_min.append(min(f_obj)) #Todas las funciones objetivo minimas
        lista_max.append(max(f_obj)) #Todas las funciones objetivo maximas
        lista_prom.append(np.mean(f_obj)) #Todos los promedios de las funciones objetivo
        lista_mejores_parques.append(cromosomas[f_obj.index(max(f_obj))] )     #Cambio

        #Se hace la evaluacion, seleccion y crossover
        selec_cross()
        mutacion()
        binario_to_potencia()

        #Se vuelven a inicializar las listas para la nueva generacion
        f_obj = list(np.zeros(p))
        fitness = list(np.zeros(p))

        #Calculo de la funcion objetivo y el fitness de los cromosomas hijos
        calcula_f_obj()
        calcula_fitness()


    parque_optimo, f_obj_parque_optimo = mejor_parque()

    print()
    print(parque_optimo)
    print()
    for i in parque_optimo:
        print(i)
    print()
    print(f_obj_parque_optimo)
    print()

    print(perf_counter() - start_time)





    ################ Salida del sistema ################

    #Lista con la cantidad de generaciones
    generacion = np.arange(1, g + 1)

    ## GRÁFICAS
    plt.subplots()
    plt.title("Evolucion de cromosomas")
    plt.axhline(y = max(lista_max), color = 'r', label = "FObj del Cromosoma Optimo")
    plt.plot(generacion, lista_min, color = 'k', label = "Min")
    plt.plot(generacion, lista_max, color = 'b', label = "Max")
    plt.plot(generacion, lista_prom, color = 'g', label = "Prom")
    plt.grid(True)
    plt.xlabel("Cantidad de ciclos")
    plt.ylabel("Funcion objetivo (FO)")
    plt.legend(loc = "lower right")
    plt.tight_layout()
    plt.show()

    
    ## TABLA DE EXCEL
    Datos = pd.DataFrame({"Generacion": generacion, "Minimo FO": lista_min, "Maximo FO": lista_max, "Promedio FO": lista_prom}) 
    Datos_parque = pd.DataFrame(parque_optimo)

    Tabla = pd.ExcelWriter('D:/Descargas/Facultad/Python/TP Investigacion/tabla.xlsx', engine='xlsxwriter')
    Parque = pd.ExcelWriter('D:/Descargas/Facultad/Python/TP Investigacion/parque.xlsx', engine='xlsxwriter')

    Datos.to_excel(Tabla, sheet_name='Valores', index = False)
    Datos_parque.to_excel(Parque, sheet_name='Valores', index = False)

    ## DISEÑO TABLA
    workbook = Tabla.book
    workbook2 = Parque.book
    worksheet = Tabla.sheets["Valores"]
    worksheet2 = Parque.sheets["Valores"]

    formato = workbook.add_format({"align": "center"})
    formato2 = workbook2.add_format({"align": "center"})

    worksheet.set_column("A:D", 15, formato)                       # Creo varios worksheet para dar formato a distintas columnas por separado.
    worksheet2.set_column("A:J", 15, formato2)
    worksheet.conditional_format("C1:C"+str(len(lista_prom)+1), {"type": "3_color_scale"})
    worksheet2.conditional_format(1, 0, 11, 11, {"type": "3_color_scale"})
    Tabla.save()
    Parque.save()

    #Se borra la tabla de EXCEL
    input()
    os.remove('D:/Descargas/Facultad/Python/TP Investigacion/tabla.xlsx')
    os.remove('D:/Descargas/Facultad/Python/TPI nvestigacion/parque.xlsx')
    print("Tabla borrada. Fin de programa")
    input()