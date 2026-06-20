import numpy as np
import matplotlib.pyplot as plt
# Cargamos solo la tercera columna (índice 2)
# Si los datos están separados por espacios, no necesitas especificar delimiter.
# Si estuvieran separados por comas, añadirías: delimiter=','
iteraciones = np.loadtxt('resultados.txt', usecols=2)

# Ahora la variable 'columna_3' contiene solo esos datos
#print(iteraciones)

valores_unicos, conteos = np.unique(iteraciones, return_counts=True)
# 3. Mostramos los resultados ordenados

#print("Número -> Cantidad de veces que aparece")
#print("-" * 40)
total = 1432
P_n = conteos/total
#print('P_n',P_n)
suma=0
#for valor, p in zip(valores_unicos, P_n):
#    print(int(valor),p)
suma = 0

for i in range(len(valores_unicos)):
    suma += conteos[i]
    print(valores_unicos[i],conteos[i],suma/total)
    plt.plot(valores_unicos[i],suma/total,'rx')    
 
#plt.plot(valores_unicos,conteos)
plt.show()
