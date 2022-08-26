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
        self.l1=100
        self.l2=100
        self.i =1
        self.corr=[]


    def getInfo(self):
        info = {
            'x2' : round(self.x, 4),
            'y2' : round(self.y, 4),
            'q1 (deg)' : round(-self.q1*57.2958, 4),
            'q2 (deg)' : round(-self.q2*57.2958, 4),

        }
        return info


    def dynamics(self,min,max,speed):
        for q1 in range(min,max+1, speed):
            for q2 in range(min,max+1, speed):
                self.q1=-q1*np.pi/180
                self.q2=-q2*np.pi/180
                self.render()

        
        


    def draw(self, image):

        cv2.line(image, (300, 300), (300 + int(self.l1*np.cos(self.q1)), 300 + int(self.l1*np.sin(self.q1))), (0, 255, 0), 1)

        cv2.circle(image, (300 + int(self.l1*np.cos(self.q1)), 300 +  int(self.l1*np.sin(self.q1))), 3, (0, 0, 255), -1)
        
        ln1f=(300 + int(self.l1*np.cos(self.q1)), 300 +  int(self.l1*np.sin(self.q1)))

        cv2.line(image, (300 + int(self.l1*np.cos(self.q1)), 300 +  int(self.l1*np.sin(self.q1))) , (ln1f[0]+int(self.l2*np.cos(self.q2)),ln1f[1]+int(self.l2*np.sin(self.q2))), (0, 255, 0), 1)
    
        cv2.circle(image, (ln1f[0]+int(self.l2*np.cos(self.q2)),ln1f[1]+int(self.l2*np.sin(self.q2))), 4, (0, 0, 255), -1)
        
        for x, y  in self.corr:
            cv2.circle(image, (x, y), 1, (0, 0, 120), -1)
        
        self.i +=1
        # if self.i %10 ==0:
        self.corr.append((ln1f[0]+int(self.l2*np.cos(self.q2)),ln1f[1]+int(self.l2*np.sin(self.q2))))

        return image


obj = TwoR()    

for i in range(8000):
    obj.dynamics(min=35,max=145,speed=3)
    print("Press 'q' to exit")
    if i % 1 == 0:
        obj.render(height= 600, pause = 10)

    

    
    

