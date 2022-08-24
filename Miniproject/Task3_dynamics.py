from render import Renderer, scaleAndShow
import cv2
import numpy as np
from numpy import *
import matplotlib.pyplot as plt
import time


class TwoR(Renderer):
    def __init__(self, height=600, width=600, recordLocation=None):
        super().__init__(recordLocation=recordLocation)
        self.x=300
        self.y=300
        self.q1=np.pi/6
        self.q2=np.pi/3
        self.l1=200
        self.l2=200
        self.i = 0
        self.corr=[]
        self.angle=[]
        self.theta=0
        self.m1=5
        self.m2=5
        self.g= 9.8
        self.tau_1=0
        self.tau_2=0
        self.Fx=0
        self.Fy=0
        self.tau_ext_1 =0
        self.tau_ext_2 =0
        self.x_origin=200
        self.y_origin=200
        self.x_initial=[self.x_origin-50]
        self.y_initial=[self.y_origin+50]

        # spring Constant
        self.k=2


    def getInfo(self):
        info = {
            'Fx(N)'    : round(self.Fx),
            'Fy(N)'    : round(self.Fy),
            'Tau_1 (N-m)'    : round(self.tau_1),
            'Tau_2(N-m)'    : round(self.tau_2)
        }
        return info


    def dynamics(self):
        # give trajectory of x and y 
        m1 , m2, g, l1,l2= self.m1, self.m2, self.g ,self.l1/200, self.l2/200

        x_final =  (self.x_origin-self.x_initial[-1])
        x_final = x_final+200
        y_final =  (self.y_origin-self.y_initial[-1])
        y_final= y_final+200

        x =  (1-self.i)*self.x_initial[-1] + self.i*x_final
        # print(self.x_initial)
        y =  (1-self.i)*self.y_initial[-1] + self.i*y_final
    
        if self.i>=1:
            self.x_initial.append(x_final)
            self.y_initial.append(y_final)
            self.i=0
            
        



        # x = 150*np.cos(self.i/100) +300
        # y = 200*np.sin(self.i/100) +300
        self.x= x-300
        self.y= y-300

        if (self.x**2+self.y**2) >= (self.l1+self.l2)**2:
            raise Exception("Just Be in The LIMIT :)")

        self.theta= np.arccos((self.x**2+self.y**2-self.l1**2-self.l2**2)/(2*self.l1*self.l2))
        self.q1 = np.arctan(self.y/self.x) - np.arctan((self.l2*np.sin(self.theta))/(self.l1+self.l2*np.cos(self.theta)))
        self.q2 = self.q1+self.theta
        if x<300:
            self.q1= -np.pi+self.q1
            self.q2= -np.pi+self.q2
        
        self.Fx= -self.k*(x-self.x_origin)
        self.Fy= -self.k*(y-self.y_origin)
        self.tau_ext_1 = -l1*sin(self.q1)*self.Fx + l1*cos(self.q1)*self.Fx
        self.tau_ext_2 = -l2*sin(self.q2)*self.Fx + l2*cos(self.q2)*self.Fy

        self.i +=self.k*0.01

        # if self.i % 100 ==0:
        self.corr.append((int(x), int(y), time.time()))
        self.angle.append((self.q1,self.theta, time.time()))

        # Langrangian Equations 

        for i in range(len(self.angle)):
            if i>2:
                q1_1=self.angle[i-2][0]
                q1_2=self.angle[i-1][0]
                q1_3=self.angle[i][0]
                q2_1=self.angle[i-2][1]
                q2_2=self.angle[i-1][1]
                q2_3=self.angle[i][1]
                t1=self.angle[i-2][2]
                t2=self.angle[i-1][2]
                t=self.angle[i][2]


                dt21=(t2-t1)/2
                dt32=(t-t2)/2

                dq1_prev= (q1_2-q1_1)/dt21
                dq2_prev=(q2_2-q2_1)/dt21

                dq1=(q1_3-q1_2)/dt32
                dq2=(q2_3-q2_2)/dt32

                ddq1= (2*(dq1-dq1_prev))/(t-t1)
                ddq2=(2*(dq2-dq2_prev))/(t-t1) 


                self.tau_1=  (m1*l1**2 + m2*(l1**2+2*l1*l2*cos(q2_3)+l2**2))*ddq1 + m2*(l1*l2*cos(q2_3)+l2**2)*ddq2 -m2*l1*l2*sin(q2_3)*(2*dq1*dq2+ddq2**2) + (m1+m2)*l1*g*cos(q1_3) + m2*g*l2*cos(q1_3+q2_3)
 
                self.tau_2 = m2*(l1*l2*cos(q2_3)+l2**2)*ddq1 + m2*l2**2*q2_3**2 +  m2*l1*l2*dq1**2*sin(q2_3) + m2*g*l2*cos(q1_3+q2_3)
                
                # Constant force

                self.tau_1 = self.tau_1+ self.tau_ext_1

                self.tau_2 = self.tau_2+ self.tau_ext_2
                    
            
        

        


    def draw(self, image):

        cv2.line(image, (300, 300), (300 + int(self.l1*np.cos(self.q1)), 300 + int(self.l1*np.sin(self.q1))), (0, 255, 0), 1)

        cv2.circle(image, (300 + int(self.l1*np.cos(self.q1)), 300 +  int(self.l1*np.sin(self.q1))), 5, (0, 0, 255), -1)
        
        ln1f=(300 + int(self.l1*np.cos(self.q1)), 300 +  int(self.l1*np.sin(self.q1)))

        cv2.line(image, (300 + int(self.l1*np.cos(self.q1)), 300 +  int(self.l1*np.sin(self.q1))) , (ln1f[0]+int(self.l2*np.cos(self.q2)),ln1f[1]+int(self.l2*np.sin(self.q2))), (0, 255, 0), 1)
    
        cv2.circle(image, (ln1f[0]+int(self.l2*np.cos(self.q2)),ln1f[1]+int(self.l2*np.sin(self.q2))), 7, (0, 0, 255), -1)
    
        cv2.circle(image, (self.x_origin,self.y_origin), 7, (0, 0, 0), -1)


        for i in range(len(self.corr)):
            x= self.corr[i][0]
            y= self.corr[i][1]
            cv2.circle(image, (x, y), 1, (0, 0, 120), -1)

        return image




obj = TwoR()    

for i in range(10000):
    obj.dynamics()
    if i % 1 == 0:
        obj.render(height= 600, pause = 10)

    

    
    

