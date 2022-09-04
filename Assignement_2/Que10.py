import numpy as np

def RRR_elbow(l1,l2,l3,d,theta_1,theta_2,theta_3):

    theta_1=theta_1*np.pi/180
    theta_2=theta_2*np.pi/180
    theta_3=theta_3*np.pi/180
           
    H_01=np.array([[np.cos(theta_1),-np.sin(theta_1),0,0],
                    [np.sin(theta_1),np.cos(theta_1),0,0],
                    [0,0,1,0],
                    [0,0,0,1]])


    H_12=np.array([[np.cos(theta_2),-np.sin(theta_2),0,l1],
                    [np.sin(theta_2),np.cos(theta_2),0,0],
                    [0,0,1,0],
                    [0,0,0,1]])
    
    H_23=np.array([[np.cos(theta_3),-np.sin(theta_3),0,l2],
                    [np.sin(theta_3),np.cos(theta_3),0,0],
                    [0,0,1,0],
                    [0,0,0,1]])

    p_1=np.array([[l1],[0],[0],[1]])
    p_2=np.array([[l2],[0],[0],[1]])
    p_3=np.array([[l3],[0],[0],[1]])
    
    O0=np.array([[0],[0],[0]])
    O1=(H_01@p_1)[:-1]
    O2=(H_01@H_12@p_2)[:-1]
    O3=(H_01@H_12@H_23@p_3)[:-1]

    z0=np.array([[0],[0],[1]])
    z1=np.array([[0],[0],[1]])
    z2=np.array([[0],[0],[1]])


    jacobian=np.concatenate((np.concatenate((np.cross(z0.T,(O3-O0).T).T,np.cross(z1.T,(O3-O1).T).T,np.cross(z2.T,(O3-O2).T).T),axis=1),
            np.concatenate((z0,z1,z2),axis=1)),axis=0)

  

    return print(jacobian.real)

#Enter the arguments below 
theta_1=80     # angle in degree
theta_2=50     # angle isn degree
theta_3=40

l1=4             # length of arm 1
l2=5            # length of arm 2
l3=5
d=2         # instantaneous distance from prismatic joint

RRR_elbow(l1,l2,l3,d,theta_1,theta_2,theta_3)