import numpy as np
import matplotlib.pyplot as plt

nu0 = 4*np.pi*10**-7

#Esta clase generaliza un diferencial de intensidad en Python
#con una posición, intensidad y dirección.

class Intensidad(object):
    
    def __init__(self, p, I, d):
        
        self.I = I
        self.x, self.y, self.z = p
        self.xd, self.yd, self.zd = d
        
    def campoB(self, p):
        
        r = [p[0]-self.x, p[1]-self.y, p[2]-self.z]
        dlr = [r[2]*self.yd - self.zd*r[1],
               r[0]*self.zd - self.xd*r[2],
               r[1]*self.xd - self.yd*r[0]]
        nu = nu0/(4*np.pi)
        r3 = (r[0]**2 + r[1]**2 + r[2]**2)**1.5
        
        return [nu*dlr[0]/r3, nu*dlr[1]/r3, nu*dlr[2]/r3]
 
#Esta clase generaliza un momento dipolar magnético con una posición
#y un vector magnetización    
 
class DipoloM(object):
    
    def __init__(self, p, m):
        
        self.x, self.y, self.z = p
        self.mx, self.my, self.mz = m
        
    def campoB(self, p):
        
        rx, ry, rz = p
        r = np.sqrt(rx**2 + ry**2 +rz**2)
        mr = self.mx*rx + self.my*ry + self.mz*rz
        nu = nu0/(4*np.pi*r**3)
        
        return [nu*((3*rx*mr)/r**2 - self.mx), nu*((3*ry*mr)/r**2 - self.my),
                nu*((3*rz*mr)/r**2 - self.mz)]
        
#Calcula el campo debido a una distribución de intensidades o momentos
    
def biotSav(p, intens):
    
    Bx, By, Bz = 0, 0, 0
    
    for I in intens:
        
        B = I.campoB(p)
        Bx += B[0]
        By += B[1]
        Bz += B[2]
        
    return Bx, By, Bz

#Genera el campo en todo el espacio a partir de una distribución de
#intensidades o dipolos

def campB(X, Y, Z, I):
    
    n = len(X)
    Bx, By, Bz = np.zeros((n, n, n)), np.zeros((n, n, n)), np.zeros((n, n, n))
    
    for i in range(n):
        
        for j in range(n):
            
            for k in range(n):
                
                p = [X[i], Y[j], Z[k]]
                B = biotSav(p, I)
                Bx[i, j, k] = B[0]
                By[i, j, k] = B[1]
                Bz[i, j, k] = B[2]
                
    return Bx, By, Bz

def main():
    
    #Definimos el espacio
    
    L = 4
    n = 20
    
    X = np.linspace(-L, L, n)
    Y = np.linspace(-L, L, n)
    Z = np.linspace(-L, L, n)
    
    #Definimos una distribución de intensidades y calculamos su campo
    
    I = 1
    r = 1
    Is = []
    ns = 20
    thetas = np.linspace(0, 2*np.pi, ns)
    
    for i in range(ns):
        
        p = [r*np.cos(thetas[i]), r*np.sin(thetas[i]), 0]
        dl = [-r*np.sin(thetas[i]), r*np.cos(thetas[i]), 0]
        Ia = Intensidad(p,I,dl)
        Is.append(Ia)
    
    B = campB(X, Y, Z, Is)
    Bx, By, Bz = B
    
    x, y, z = np.meshgrid(X, Y, Z)
    
    fig = plt.figure()
    ax = fig.gca(projection="3d")
    ax.quiver(x, y, z, B[0], B[1], B[2], normalize=True, length=1, color="b")
    #plt.savefig("Espira - Campo 3D.png", dpi=1200)
    
    plt.figure(figsize=(9, 9))
    plt.streamplot(X, Y, By[:,:,n//2-1], Bx[:,:,n//2-1])
    plt.xlabel("X (m)")
    plt.ylabel("Y (m)")
    #plt.savefig("Mini espira - Corte Z.png", dpi=1200)
    
    plt.figure(figsize=(9, 9))
    plt.streamplot(Y, Z, Bz[n//2,:,:], By[n//2,:,:])
    plt.xlabel("Y (m)")
    plt.ylabel("Z (m)")
    #plt.savefig("Mini espira - Corte X.png", dpi=1200)
    
    B1 = np.sqrt(Bx[n//2-1,n//2-1,:]**2 + By[n//2-1,n//2-1,:]**2 + 
                 Bz[n//2-1,n//2-1,:]**2)
    
    plt.figure(figsize=(9, 6))
    plt.plot(Z, B1)
    plt.xlabel("Z (m)")
    #plt.savefig("Espira - Módulo Z.png", dpi=1200)
    
    B2 = np.sqrt(Bx[:,n//2-1,n//2-1]**2 + By[:,n//2-1,n//2-1]**2 + 
                 Bz[:,n//2-1,n//2-1]**2)
    
    plt.figure(figsize=(9, 6))
    plt.plot(X, B2)
    plt.xlabel("X (m)")
    #plt.savefig("Espira - Módulo X.png", dpi=1200)
    
    #Idem pero con dipolos
    
    M = DipoloM([0, 0, 0], [0, 0, 1])
    
    B = campB(X, Y, Z, [M])
    Bx, By, Bz = B
    
    fig = plt.figure()
    ax = fig.gca(projection="3d")
    ax.quiver(x, y, z, B[0], B[1], B[2], normalize=True, length=1, color="b")
    #plt.savefig("Momento dipolar - 3D.png", dpi=1200)
    
    plt.figure(figsize=(9, 9))
    plt.streamplot(X, Y, By[:,:,n//2-1], Bx[:,:,n//2-1])
    plt.xlabel("X (m)")
    plt.ylabel("Y (m)")
    #plt.savefig("Momento dipolar - Corte Z.png", dpi=1200)
    
    plt.figure(figsize=(9, 9))
    plt.streamplot(Y, Z, Bz[n//2,:,:], By[n//2,:,:])
    plt.xlabel("X (m)")
    plt.ylabel("Z (m)")
    #plt.savefig("Momento dipolar - Corte Y.png", dpi=1200)
    
    #Distribución cuadrada
    
    nm = 20
    t = np.linspace(-L, L, nm)
    Ms = []
    
    for i in t:
        
        M = DipoloM([i, 0, 0], [1, 0, 0])
        Ms.append(M)
        
    B = campB(X, Y, Z, Ms)
    Bx, By, Bz = B
    
    fig = plt.figure()
    ax = fig.gca(projection="3d")
    ax.quiver(x, y, z, B[0], B[1], B[2], normalize=True, length=1, color="b")
    #plt.savefig("Momento distribución - 3D", dpi=1200)
    
    plt.figure(figsize=(9, 9))
    plt.streamplot(X, Y, By[:,:,n//2-1], Bx[:,:,n//2-1])
    plt.xlabel("X (m)")
    plt.ylabel("Y (m)")
    plt.savefig("Momento distribución - Corte Z", dpi=1200)
    
    plt.figure(figsize=(9, 9))
    plt.streamplot(Y, Z, By[n//2-1,:,:], Bz[n//2-1,:,:])
    plt.xlabel("Y (m)")
    plt.ylabel("Z (m)")
    plt.savefig("Momento distribución - Corte X", dpi=1200)
    
    #Dos representaciones distancia-módulo para I
    
if __name__ == "__main__":
    
    main()