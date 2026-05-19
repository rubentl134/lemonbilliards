'''
This is the script for lemon billiards simulations.
We are doing it here step-by-step.
'''
import numpy as np
import matplotlib.pyplot as plt
from lemon1 import *

print('WELCOME TO LEMON BILLIARDS')
print('This a version that shows the trajectory of the particle after n iterations')
plt.hlines(y=0,xmin=-2,xmax=2,color='red')
plt.vlines(x=0,ymin=-2,ymax=2,color='red')

#circunferences
center1 = np.array([0,0]) # this is set to (0,0)
center2 = np.array([0.4,0]) # need to be changed

r1 = 1
r2 = 0.8
C1 = circunference(center1,r1)
C2 = circunference(center2,r2)

print('center1:',center1,'radius1:',r1)
print('center2:',center2,'radius2:',r2)

print('parameter B=',np.abs(center1[0])+np.abs(center2[0]))

plt.plot(C1[0],C1[1],color='blue')
plt.plot(C2[0],C2[1],color='blue')

plt.plot(center1[0],center1[1],'ro')
plt.plot(center2[0],center2[1],'ro')

# Now we delimit the lemon region
x_intersection = asymmetric_x_intersection(center1,r1,center2,r2)
y_intersection = (r1**2 - (x_intersection - center1[0])**2)**0.5  #this has two values and we are going to use de + and the -

plt.vlines(x_intersection,ymin=-1,ymax=1,color='gray',linestyle='--')

UpLemon = [x_intersection, y_intersection]
DownLemon = [x_intersection, -y_intersection]
RightLemon = right(center1,r1)
LeftLemon = left(center2,r2)

plt.plot(UpLemon[0],UpLemon[1],'go')
plt.text(UpLemon[0],UpLemon[1]+0.1,'UpLemon')
plt.plot(DownLemon[0],DownLemon[1],'go')
plt.text(DownLemon[0],DownLemon[1]-0.1,'DownLemon')
plt.plot(RightLemon[0],RightLemon[1],'go')
plt.text(RightLemon[0]+0.01,RightLemon[1]-0.1,'RightLemon')
plt.plot(LeftLemon[0],LeftLemon[1],'go')
plt.text(LeftLemon[0]-0.2,LeftLemon[1]-0.1,'LeftLemon')

s_total = r1*angle_bw_vectors(DownLemon-center1,UpLemon-center1) + r2*angle_bw_vectors(UpLemon-center2,DownLemon-center2)
s_angle_right =angle_bw_vectors(DownLemon-center1,UpLemon-center1)
s_angle_left =angle_bw_vectors(UpLemon-center2,DownLemon-center2)

s_right = r1*s_angle_right
s_left = r2*s_angle_left

#####################################################################################
### initial conditions
#s and alpha
s=0.1 #from 0 to 1
alpha=np.deg2rad(10) #in degrees
P = s2xy_asym(s,DownLemon,RightLemon,UpLemon,LeftLemon,s_total,s_right,s_left,r2,center2)
plt.plot(P[0],P[1],'kX')
#velocity
v = velocity_lemon(P,center1,center2,x_intersection,alpha)
s_initial = xy2s_asymmetrical(P,center1,center2,DownLemon - center1,UpLemon - center2,x_intersection,r1,r2,DownLemon,UpLemon,LeftLemon)

###########################################################################################
#exit
h = 0.01
s_exit = 0.5
P_exit0 = s2xy_asym(s_exit,DownLemon,RightLemon,UpLemon,LeftLemon,s_total,s_right,s_left,r2,center2)
P_exit1 = s2xy_asym(s_exit+h,DownLemon,RightLemon,UpLemon,LeftLemon,s_total,s_right,s_left,r2,center2)
P_exit2 = s2xy_asym(s_exit-h,DownLemon,RightLemon,UpLemon,LeftLemon,s_total,s_right,s_left,r2,center2)
plt.plot(P_exit0[0],P_exit0[1],'kd')
plt.plot(P_exit1[0],P_exit1[1],'kd')
plt.plot(P_exit2[0],P_exit2[1],'kd')

############################################################################################



plt.arrow(P[0],P[1],v[0]*0.2,v[1]*0.2,color='green',width=0.01)


for i in range(30):
    #intersecion
    c1_intersection = line_and_circle(P,v,center1,r1)
    c2_intersection = line_and_circle(P,v,center2,r2)
    #which is the correct intersection point?
    #print(20*'---')
    candidates = [c1_intersection[0],c1_intersection[1],c2_intersection[0],c2_intersection[1]]
    val_count = 0
    for k in range(len(candidates)):
        val = belongs_asymmetrical_lemon(candidates[k],P,UpLemon,DownLemon,LeftLemon,RightLemon,x_intersection,center1,r1,center2,r2)
        if val == True:
            val_count =+ 1
            hit = candidates[k]
            #print('hit at:',hit)
    if val_count != 1:
        #checking there is one solution
        print('DIFFERENTS SOLUCIOTN POINTS! CHECK IT')

    #velocity
    vReflex = reflexVelocity_asymmetric(hit,v,center1,center2,x_intersection)
    #print('reflexion',vReflex)
    #plt.arrow(hit[0],hit[1],vReflex[0],vReflex[1])
    plt.plot([P[0],hit[0]],[P[1],hit[1]],color='skyblue',alpha = 0.6)
    #update
    P = hit
    v = vReflex
    
    #print('hit',hit)
    s = xy2s_asymmetric_2(P,r1,r2,center1,center2,DownLemon,UpLemon,LeftLemon,RightLemon,x_intersection)
    alpha_final = alpha_angle_asymmetrical(hit,vReflex,center1,center2,x_intersection)
    print(s,np.rad2deg(alpha_final),hit)
    if s > s_exit-h and s < s_exit + h:
        print(i,'exit at____________________________',s)
        #print(float(sys.argv[1]),float(sys.argv[2]),i,'OUT')
        #break

plt.show()
