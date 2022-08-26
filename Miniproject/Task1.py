from render import Renderer, scaleAndShow
import cv2
import numpy as np
import matplotlib.pyplot as plt


class TwoR(Renderer):
    def __init__(self, height=600, width=600, recordLocation=None):
        super().__init__(recordLocation=recordLocation)
        self.x=300
        self.y=300
        self.q1=np.pi/6
        self.q2=np.pi/3
        self.l1=200
        self.l2=200
        self.i =1
        self.corr=[]


    def getInfo(self):
        info = {
            'x2' : round(self.x, 4),
            'y2' : round(self.y, 4),
            'q1 (deg)' : round(self.q1*57.2958, 4),
            'q2 (deg)' : round(self.q2*57.2958, 4),

        }
        return info


    def dynamics(self):
        # Given trajectory of x and y  with self.i as the iterator
        x = 200*np.cos(self.i/100) +300
        y = 200*np.sin(self.i/100) +300
        self.x= x-300
        self.y= y-300


        # Range or workspace with singularities filtered out
        if (self.x**2+self.y**2) > (self.l1+self.l2)**2 or self.x==0 or self.y==0:
            raise Exception("Just Be in The LIMIT :)")

        theta= np.arccos((self.x**2+self.y**2-self.l1**2-self.l2**2)/(2*self.l1*self.l2))
        self.q1 = np.arctan(self.y/self.x) - np.arctan((self.l2*np.sin(theta))/(self.l1+self.l2*np.cos(theta)))
        self.q2 = self.q1+theta
        if x<300:
            self.q1= -np.pi+self.q1
            self.q2= -np.pi+self.q2

        self.i +=1
        if self.i %10 ==0:
            self.corr.append((int(x), int(y)))
        


    def draw(self, image):

        cv2.line(image, (300, 300), (300 + int(self.l1*np.cos(self.q1)), 300 + int(self.l1*np.sin(self.q1))), (0, 255, 0), 1)

        cv2.circle(image, (300 + int(self.l1*np.cos(self.q1)), 300 +  int(self.l1*np.sin(self.q1))), 5, (0, 0, 255), -1)
        
        ln1f=(300 + int(self.l1*np.cos(self.q1)), 300 +  int(self.l1*np.sin(self.q1)))

        cv2.line(image, (300 + int(self.l1*np.cos(self.q1)), 300 +  int(self.l1*np.sin(self.q1))) , (ln1f[0]+int(self.l2*np.cos(self.q2)),ln1f[1]+int(self.l2*np.sin(self.q2))), (0, 255, 0), 1)
    
        cv2.circle(image, (ln1f[0]+int(self.l2*np.cos(self.q2)),ln1f[1]+int(self.l2*np.sin(self.q2))), 7, (0, 0, 255), -1)
        
        for x, y  in self.corr:
            cv2.circle(image, (x, y), 1, (0, 0, 120), -1)

        return image


obj = TwoR()    

for i in range(8000):
    obj.dynamics() #give trajectory in this function with i as the iterator
    print("Press 'q' to exit")
    if i % 1 == 0:
        obj.render(height= 600, pause = 10)

    

    
    

