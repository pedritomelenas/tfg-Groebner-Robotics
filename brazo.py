from matplotlib.widgets import Slider, RadioButtons
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from sympy import *

def main():
    myRobot=Draw_Robot()

class Draw_Robot():
    def __init__(self):
        '''VARIABLES'''
        ARTICULACIONES = int(3) #numero de ARTICULACIONES
        #lo dejo asi con la idea de que se pueda elejir el numero en el futuro

        self.l = np.array([0, 100, 100])  # longitudes
        self.x = np.array([0]*ARTICULACIONES,dtype=float) #componente x
        self.y = np.array([0]*ARTICULACIONES,dtype=float) #componente y
        self.z = np.array([0]*ARTICULACIONES,dtype=float) #componente z

        #Variables temporales para calcular groebner para guardar los datos de los slidders
        self.tx=0
        self.ty=100
        self.tz=100

        self.estado_posible = True
        self.modo = 0

        #Calculos de Groebner previos para obtener los angulos en funcion de el punto de destino
        a,b,c1,c2,s1,s2,l2,l3 = var("a,b,c1,c2,s1,s2,l2,l3")
        eq= [a-l3*(c1*c2-s1*s2)-l2*c1, b-l3*(c1*s2+c2*s1)-l2*s1, c1**2+s1**2-1,c2**2+s2**2-1]
        R = QQ.frac_field(a,b,l2,l3)
        G = groebner(eq,c2,s2,c1,s1,order='lex',domain=R)
        self.c2_var = solve(G[0].subs([(l2,self.l[1]),(l3,self.l[2])]),c2)
        self.s1_var = solve(G[3].subs([(l2,self.l[1]),(l3,self.l[2]),(c2,self.c2_var)]),s1)
        self.s2_var = []
        self.s2_var.append(solve(G[1].subs([(l2,self.l[1]),(l3,self.l[2]),(c2,self.c2_var),(s1,self.s1_var[0])]),s2)[0])
        self.s2_var.append(solve(G[1].subs([(l2,self.l[1]),(l3,self.l[2]),(c2,self.c2_var),(s1,self.s1_var[1])]),s2)[0])
        self.c1_var = []
        self.c1_var.append(solve(G[2].subs([(l2,self.l[1]),(l3,self.l[2]),(c2,self.c2_var),(s1,self.s1_var[0]),(s2,self.s2_var[0])]),c1)[0])
        self.c1_var.append(solve(G[2].subs([(l2,self.l[1]),(l3,self.l[2]),(c2,self.c2_var),(s1,self.s1_var[1]),(s2,self.s2_var[1])]),c1)[0])


        self.fig = plt.figure("Brazo Simple")  #create the frame
        #self.fig.patch.set_facecolor('white')
        self.ax = plt.axes([0, 0.2, 1, .8], projection='3d') #3d ax panel
        self.axerror = plt.axes([0.40, 0.92, 0.001, .001])#panel for error message

        '''Dibujamos los widgets'''
        #dibujamos el slider panel de la X
        axxval = plt.axes([0.35, 0.12, 0.45, 0.03])
        self.x_val = Slider(axxval, 'X', -150, 150, valinit=self.tx)
        self.x_val.set_active(False)

        #dibujamos el slider panel de la Y
        axyval = plt.axes([0.35, 0.0775, 0.45, 0.03])
        y_val = Slider(axyval, 'Y', -150, 150, valinit=self.ty)

        #dibujamos el slider panel de la Z
        axzval = plt.axes([0.35, 0.035, 0.45, 0.03])
        z_val = Slider(axzval, 'Z', 0, 200, valinit=self.tz)

        #radio buttons para seleccion de modo
        rax = plt.axes([0.05, 0.02, 0.22, 0.12])#, axisbg=axcolor)
        rax.set_title('Modo', fontsize=12)
        set_modo = RadioButtons(rax, ('2D (X disable)', '3D (X enable)'), active=0)


        '''Manejadores de eventos para los widgets'''
        def actualiza_val_x(val): # x slider event
            self.tx = val
            self.draw_robot()
        self.x_val.on_changed(actualiza_val_x)

        def actualiza_val_y(val):#y slider event
            self.ty= val
            self.draw_robot()
        y_val.on_changed(actualiza_val_y)

        def actualiza_val_z(val):#z slider event
            self.tz = val
            self.draw_robot()
        z_val.on_changed(actualiza_val_z)

        def seleccion_modo(label):#Radio Button Modo
            if label==('2D (X disable)'):
                self.modo = 0
                self.x_val.set_active(False)
                self.x_val.reset()
            if label==('3D (X enable)'):
                self.modo = 1
                self.x_val.set_active(True)
        set_modo.on_clicked(seleccion_modo)

        self.display_error()#draw the error and hide it
        self.draw_robot()#draw function to draw robot

        plt.show()#end of constructor

    '''Funciones de clase'''
    def display_error(self):
        self.axerror.set_visible(False)
        self.axerror.set_yticks([])
        self.axerror.set_xticks([])
        self.axerror.set_navigate(False)
        self.axerror.text(0, 0, 'El brazo no alcanza esa posicion!', style='oblique',
                      bbox={'facecolor':'red', 'alpha':0.5, 'pad':10}, size=20, va = 'baseline')


    #Hace todos los calculos y cambia las variables, ademas, si no es posible alcanzar el punto lo indica
    def calcula_groebner3D(self):
        dist = sqrt(self.tx**2+self.ty**2)
        s1,s2,c1= [],[],[]

        s1.append(self.s1_var[0].subs([(a, dist), (b,self.tz)]))
        s1.append(self.s1_var[1].subs([(a, dist), (b,self.tz)]))
        s2.append(self.s2_var[0].subs([(a, dist), (b,self.tz)]))
        s2.append(self.s2_var[1].subs([(a, dist), (b,self.tz)]))
        c1.append(self.c1_var[0].subs([(a, dist), (b,self.tz)]))
        c1.append(self.c1_var[1].subs([(a, dist), (b,self.tz)]))

        if(im(c1[0])==0 and im(s1[0])==0 and im(s2[0])==0):
            self.estado_posible = True

            if(self.tx>0 and self.ty>0):
                fi = atan(self.ty/self.tx)
            elif(self.tx>0 and self.ty<0):
                fi = atan(self.ty/self.tx) + 2*pi
            elif(self.tx<0):
                fi = atan(self.ty/self.tx) + pi
            elif(self.tx==0 and self.ty>0):
                fi = pi/2
            else:
                fi = -pi/2

            if ( s1[0]*self.l[1] > 0 ):
                self.x[1] = self.l[1]*c1[0]*cos(fi)
                self.y[1] = self.l[1]*c1[0]*sin(fi)
                self.z[1] = self.l[1]*s1[0]
            else:
                self.x[1] = self.l[1]*c1[1]*cos(fi)
                self.y[1] = self.l[1]*c1[1]*sin(fi)
                self.z[1] = self.l[1]*s1[1]

            self.x[2] = self.tx
            self.y[2] = self.ty
            self.z[2] = self.tz

        else :
            self.estado_posible = False

    #Hace todos los calculos y cambia las variables, ademas, si no es posible alcanzar el punto lo indica
    def calcula_groebner2D(self):
        s1,s2,c1= [],[],[]

        s1.append(self.s1_var[0].subs([(a, self.ty), (b,self.tz)]))
        s1.append(self.s1_var[1].subs([(a, self.ty), (b,self.tz)]))
        s2.append(self.s2_var[0].subs([(a, self.ty), (b,self.tz)]))
        s2.append(self.s2_var[1].subs([(a, self.ty), (b,self.tz)]))
        c1.append(self.c1_var[0].subs([(a, self.ty), (b,self.tz)]))
        c1.append(self.c1_var[1].subs([(a, self.ty), (b,self.tz)]))

        if(im(c1[0])==0 and im(s1[0])==0 and im(s2[0])==0):
            self.estado_posible = True

            self.x[1] = 0
            if (s1[0]*self.l[1]>0):
                self.y[1] = c1[0]*self.l[1]
                self.z[1] = s1[0]*self.l[1]
            else :
                self.y[1] = c1[1]*self.l[1]
                self.z[1] = s1[1]*self.l[1]

            self.x[2] = 0
            self.y[2] = self.ty
            self.z[2] = self.tz

        else :
            self.estado_posible = False



    def set_positions(self):#gets the x,y,z values for the line.
        #convert arrays to lists for drawing the line
        xs = np.array(self.x).tolist()
        ys = np.array(self.y).tolist()
        zs = np.array(self.z).tolist()
        self.ax.cla() #clear current axis
        #draw new lines,  two lines for "fancy" looks
        self.ax.plot(xs, ys, zs, 'o-', markersize=20,
                     markerfacecolor="orange", linewidth = 8, color="blue")
        self.ax.plot(xs, ys, zs, 'o-', markersize=4,
                     markerfacecolor="blue", linewidth = 1, color="silver")

    def set_ax(self):#ax panel set up
        self.ax.set_xlim3d(-200, 200)
        self.ax.set_ylim3d(-200, 200)
        self.ax.set_zlim3d(-5, 200)
        self.ax.set_xlabel('X axis')
        self.ax.set_ylabel('Y axis')
        self.ax.set_zlabel('Z axis')
        for j in self.ax.get_xticklabels() + self.ax.get_yticklabels(): #hide ticks
            j.set_visible(False)
        self.ax.set_axisbelow(True) #send grid lines to the background

    def draw_robot(self):#draw and update the 3D panel
        if self.modo==0:
            self.calcula_groebner2D()
        else:
            self.calcula_groebner3D()

        if self.estado_posible:#check boundaries
            self.axerror.set_visible(False)#turn off error message panel
            self.set_positions()
            self.set_ax()
        else:
            self.axerror.set_visible(True)#display error message panel
        plt.draw()


if __name__ == '__main__':
    main()
