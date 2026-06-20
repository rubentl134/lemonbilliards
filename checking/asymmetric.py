'''
THIS IS ASYMMETRIC LEMON BILLIARDS
CHECKING FOR THE R1=R2=1 case
'''
import sys
import numpy as np
import matplotlib.pyplot as plt
from lemon1 import *

#print(40*'=')
#print('WELCOME TO ASYMMETRIC LEMON BILLIARDS')
#print(40*'=')

#scenario
#plt.hlines(y=0,xmin=-2,xmax=2,color='red')
#plt.vlines(x=0,ymin=-2,ymax=2,color='red')

#circunferences
center1 = np.array([0,0])
center2 = np.array([0.1,0])
#print('B=',np.abs(center1[0])+np.abs(center2[0]))
r1 = 1
r2 = 1
C1 = circunference(center1,r1)
C2 = circunference(center2,r2)

#print('center1:',center1)
#print('radio1:',r1)
#print('center2:',center2)
#print('radio2:',r2)
#plt.plot(C1[0],C1[1],color='blue')
#plt.plot(C2[0],C2[1],color='blue')

x_intersection = asymmetric_x_intersection(center1,r1,center2,r2)
y_intersection = (r1**2 - (x_intersection - center1[0])**2)**0.5  #this has two values and we are going to use de + and the -
#print(x_intersection)

UpLemon = [x_intersection, y_intersection]
DownLemon = [x_intersection, -y_intersection]

RightLemon = right(center1,r1)
#print('Right lemon',RightLemon)

LeftLemon = left(center2,r2)
#print('Left lemon',LeftLemon)

#plt.plot(UpLemon[0],UpLemon[1],'rX')
#plt.plot(DownLemon[0],DownLemon[1],'rX')

fix1 = DownLemon - center1
fix2 = UpLemon - center2
#plt.arrow(center1[0],center1[1],fix1[0],fix1[1],color='red')
#plt.arrow(center2[0],center2[1],fix2[0],fix2[1],color='red')
s_angle_right =angle_bw_vectors(DownLemon-center1,UpLemon-center1)
s_angle_left =angle_bw_vectors(UpLemon-center2,DownLemon-center2)
#print('****************************************************')
#print('angle_right',np.rad2deg(s_angle_right))
#print('angle_left',np.rad2deg(s_angle_left))
s_right = r1*s_angle_right
s_left = r2*s_angle_left
#print('s_right',s_right)
#print('s_left',s_left)
s_total = s_right + s_left
#print('s_total',s_total)
#plt.plot(RightLemon[0],RightLemon[1],'rX')
#plt.plot(LeftLemon[0],LeftLemon[1],'rX')

#position
#x = 0.85
#P = initial_position2(x,center1,center2,r1,r2,LeftLemon[0],x_intersection,RightLemon[0],'-')
s=float(sys.argv[1])
P = s2xy_asym(s,DownLemon,RightLemon,UpLemon,LeftLemon,s_total,s_right,s_left,r2,center2)
#print('P:',P)
#plt.plot(P[0],P[1],'kX')

#tangent
#T = tangent_lemon(P,center1,center2,x_intersection)
#print('tangent',T)
#plt.arrow(P[0],P[1],T[0],T[1])
#velocity
alpha=np.deg2rad(float(sys.argv[2]))
v = velocity_lemon(P,center1,center2,x_intersection,alpha)
#s_initial = xy2s_asymmetrical(P,center1,center2,fix1,fix2,x_intersection,r1,r2,DownLemon,UpLemon,LeftLemon)
#print('s_initial',s_initial)
#print('v',v)
#plt.arrow(P[0],P[1],v[0],v[1],color='green')

#survival
h = 0.050
s_exit = 0.50

MAX_ITER = 10000

#plot hole??????
P_exit1 = s2xy_asym(s_exit+h,DownLemon,RightLemon,UpLemon,LeftLemon,s_total,s_right,s_left,r2,center2)
P_exit2 = s2xy_asym(s_exit-h,DownLemon,RightLemon,UpLemon,LeftLemon,s_total,s_right,s_left,r2,center2)
#plt.plot(P_exit1[0],P_exit1[1],'ko')
#plt.plot(P_exit2[0],P_exit2[1],'ko')


for i in range(MAX_ITER):
    #print(10*'----')
    #print(i)
    
    #intersecion
    c1_intersection = line_and_circle(P,v,center1,r1)
    #print('c1 intersection',c1_intersection)
    #plt.plot(c1_intersection[0][0],c1_intersection[0][1],'kd')
    #plt.plot(c1_intersection[1][0],c1_intersection[1][1],'kd')
    #plt.plot([c1_intersection[0][0],c1_intersection[1][0]],[c1_intersection[0][1],c1_intersection[1][1]],linestyle='dashed')

    c2_intersection = line_and_circle(P,v,center2,r2)
    #print('c2 intersection',c2_intersection)
    #plt.plot(c2_intersection[0][0],c2_intersection[0][1],'kd')
    #plt.plot(c2_intersection[1][0],c2_intersection[1][1],'kd')
    #plt.plot([c2_intersection[0][0],c2_intersection[1][0]],[c2_intersection[0][1],c2_intersection[1][1]],linestyle='dashed')

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
    #plt.plot([P[0],hit[0]],[P[1],hit[1]],color='skyblue',alpha = 0.6)
    #update
    P = hit
    v = vReflex
    
    #print('hit',hit)
    s = xy2s_asymmetric_2(P,r1,r2,center1,center2,DownLemon,UpLemon,LeftLemon,RightLemon,x_intersection)
    alpha_final = alpha_angle_asymmetrical(hit,vReflex,center1,center2,x_intersection)
    #print(s,alpha_final)

    #SURVIVE STUFF
    if s > s_exit-h and s < s_exit + h:
        #print('exit at',s)
        print(float(sys.argv[1]),float(sys.argv[2]),i+1)
        break


if i == MAX_ITER-1 :
    print(float(sys.argv[1]),float(sys.argv[2]),i+1)
    

#plt.show()


