import numpy as np

def Scara_type(l1,l2,d,theta_1,theta_2):

    theta_1=theta_1*np.pi/180
    theta_2=theta_2*np.pi/180
           
    H_01=np.array([[np.cos(theta_1),-np.sin(theta_1),0,0],
                    [np.sin(theta_1),np.cos(theta_1),0,0],
                    [0,0,1,0],
                    [0,0,0,1]])


    H_12=np.array([[np.cos(theta_2),-np.sin(theta_2),0,0],
                    [np.sin(theta_2),np.cos(theta_2),0,l1],
                    [0,0,1,0],
                    [0,0,0,1]])

    H_23=np.array([[1,0,0,0],
                    [0,1,0,l2],
                    [0,0,1,0],
                    [0,0,0,1]])

    p_3=np.array([[0],
             [0],
             [-d],
             [1]])

    P_0=(H_01@H_12@H_23@p_3)[:-1]   # End effector coordinates
  


    return print(P_0)

#Enter the arguments below 
theta_1=80     # angle in degree
theta_2=50     # angle isn degree
l1=4             # length of arm 1
l2=5            # length of arm 2
d=2         # instantaneous distance from prismatic joint

Scara_type(l1,l2,d,theta_1,theta_2)