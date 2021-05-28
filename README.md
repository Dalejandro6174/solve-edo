#objetivo: Comparacion de solucion de ode
#          con la funcion odeint de scipy
#          y las Euler y Runge-Kuta

#date: 2016-10-26
#autor: duvian alejandro

import  matplotlib.pyplot as pl
import numpy as np
from scipy.integrate import odeint

def My_euler1(fn, yo, Time):#, to, yo, tf, dt):
    """ fn: funcion f(t,y) tal que y' = f(t,y)
        yo: valor de y(to)
        Time: vector con instantes de tiempo
              igualmente espaciados

        Y : listado de valores de Y respectivos
        rev. 2016-10-26

        """

    Y = [yo]
    for q in range(1,len(Time)):
        dt = Time[q]-Time[q-1]        
        y = Y[-1] + dt*fn(Y[-1],Time[q-1])
        Y.append(y)
        
    return Y
def My_rko4(fn, yo,Time):
    """ fn: funcion f(y,t)escalar tal que y' = f(y,t)
        yo: valor de y(to)
        Time: vector con instantes de tiempo
              igualmente espaciados
              
        Se calculan
        tf: tiempo final
        to: tiempo inicial
        dt: intervalo de simulacion o delta de tiempo

        Y : listado de valores de Y respectivos
            de la dimension segun la funcion
            
        Metodo Runge-Kutta orden 4
        rev. 2016-10-26
        EVH&JAVV

        """

    Y=[yo]

    for q in range(1,len(Time)):
        dt = Time[q]-Time[q-1]
        ft = fn(Y[-1], Time[q-1])
        k1 = dt*np.array(fn(Y[-1],Time[q-1]))
        k2 = dt*np.array(fn(Y[-1] + k1/2, Time[q-1]+ dt/2))
        k3 = dt*np.array(fn(Y[-1] + k2/2, Time[q-1]+ dt/2))
        k4 = dt*np.array(fn(Y[-1] + k3, Time[q-1]+ dt))
        
        y = Y[-1] + (1./6)*(k1 +2*k2 + 2*k3 + k4)
        
        Y.append(y)
        
    return Y

def Fun_dI(y, t):
    """y: funcion dependiente de t
    t: variable independiente
    Funcion para solucion de un circuito
    vf = R i + L di/dt
    """
    v = 100
    R = 4
    L = 0.7
    di = (v - R*y)/L
    return di

##codigo ppal CASO circuito RL fuente constante

Io = 0.0
print ('Corriente inicial (cir RL): ', Io)

##definir rango de solucion

t_inicial = 0
t_final= 1
t = np.linspace(t_inicial, t_final, 8)

print ("Solucion entre" ,"t_inicial", "y" ,"t_final", "segundos")

##solucion exacta de la ED

R=4.
L=0.7
V=100.
I_exacta=(V/R)*(1 - np.exp(-(R/L)*t))

##Solucionar la ED

y_cor = odeint(Fun_dI, Io, t)

y_euler = My_euler1(Fun_dI, Io, t)

y_rk = My_rko4(Fun_dI, Io, t)


##Imprimir y graficar resultados usando Matplotlib

pl.plot(t, y_cor, 'ro', t, y_euler, 'b*',
        t, y_rk, 'g+',t, I_exacta, 'k-')
pl.title('Corriente en circuito R=4, L=0.7')
pl.xlabel('Eje de tiempo')
pl.ylabel('Eje de corriente')
pl.legend( ('odeint','My_euler1',"My_rk4",'Sol exacta') )
pl.grid(True)
pl.show()

y_rk = My_rko4(Fun_dI, Io, t)
print(y_euler,t,)
