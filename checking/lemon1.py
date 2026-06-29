import numpy as np
import matplotlib.pyplot as plt

def module(v):
    '''
    Returns the module of a 2d-vector v
    INPUT: 2d-vector
    l v\=/sqrt{v_x^2 + v_y^2}
    OUT: module of the vector
    '''
    
    return np.sqrt(v[0]**2 + v[1]**2)

def unitary(v):
    '''
    Returns the unitary vector of a 2d-vector
    IN: 2d-vector
    OUT: unitary 2d-vector
    '''
    return v/(np.sqrt(v[0]**2 + v[1]**2))

def rotation(v,theta):
    '''
    Rotation of a 2d-vector v an angle theta (clockwise)
    This is the rotation matrix
    R = \ cos(\theta)   -sin(\theta) \
        \ sin(\theta)    cos(theta)  \
    IN: 2d-vector and angle theta
    OUT: 2d rotated vector
    '''
    R = np.matrix([[np.cos(theta),-np.sin(theta) ],
                   [np.sin(theta), np.cos(theta)]])
    vv = np.matmul(R,v)
    w = np.array([vv[0,0] ,vv[0,1] ])
    return w


def solve2equation(A,B,C):
    '''
    Solves 2nd degree equation
    A*x^2 + B*x + C = 0
    IN: three coeficients of the 2nd degree equation
    OUT: two solutions
    '''
    discriminant = B**2 - 4*A*C
    if discriminant == 0:
        #just one solution
        x1 =  -B/(2*A)
        x2 = x1
    elif discriminant < 0:
        #no real solutions
        x1 = None
        x2 = None
    elif discriminant > 0:
        x1 = (-B+np.sqrt(discriminant))/(2*A)
        x2 = (-B-np.sqrt(discriminant))/(2*A)
    else:
        print('PROBLEM ON solve2equation() FUNCTION')
        print('No real solutions for the equation')
        print('Please contact complex numbers or your closest physicist')

    return np.array([x1,x2])

def line_and_circle(p,v,c,R):
    '''
    Line and circunference intersection. This is an improved function.
    IN:
    p: (x,y) coordinates of the lemon
    v: (vx,vy) velocity vector
    c: (h,k) center coordinates of the circunference
    R: radii of the circunference.
    OUT: Two solutions
    returns sol1,sol2
    two solutions, you need to choose the correct one.
    '''
    x0,y0 = p[0],p[1] #initial points
    m = v[1]/v[0] #slope
    h,k = c[0],c[1]

    A = 1 + m**2
    B = 2*(-h- x0*m**2 +m*y0 -m*k)
    C = h**2+m**2*x0**2-2*m*x0*y0+y0**2+2*m*x0*k-2*y0*k+k**2-R**2

    xsol = solve2equation(A,B,C)
    ysol = m*(xsol-x0)+y0
    sol1 = xsol[0],ysol[0]
    sol2 = xsol[1],ysol[1]
    return sol1,sol2


def asymmetric_x_intersection(center1,r1,center2,r2):
    '''
    Given the two circunference, this functions returns the intersection of the
    circuferences
    IN:
    center1: center of circunference 1
    r1: radii of circunference 1
    center2: center of circunference 2
    r2: radii of circunference 2
    OUT:
    x: x.coordinate of the vertical line
    '''
    h1,h2 = center1[0],center2[0]
    x = (1/2) * ((r1**2 - r2**2)/(h2-h1) + h1 + h2)
    return x

def right(center1,r1):
    '''
    Coordinates of the right side of the lemon
    IN:
    center1: circunference1 center
    r1: circunference1 radii
    OUT:
    [r,0] : coordinates of the lemon right side
    '''
    h1 = center1[0]
    x = solve2equation(1,-2*h1,(h1**2-r1**2))
    if x[0] > x[1]:
        r = x[0]
    else:
        r = x[1]
    return [r,0]

def left(center2,r2):
    '''
    Coordinates of the left side of the lemon
    IN:
    center2: circunference2 center
    r2: circunference2 radii
    OUT:
    [l,0] : coordinates of the lemon left side
    '''
    h2 = center2[0]
    x = solve2equation(1,-2*h2,h2**2 - r2**2)
    if x[0] < x[1]:
        l = x[0]
    else:
        l = x[1]
    return [l,0]

		
def circunference(c,r,num=5000):
    '''
    Returns the points for plot a circunference
    IN:
    c: centre (h,k)
    r: radius
    num: number of points
    OUT:
    [x,y]: points of the circunferences
    '''
    
    theta = np.linspace(0,2*np.pi,num)
    x = c[0] + r*np.cos(theta)
    y = c[1] + r*np.sin(theta)
    return np.array([x,y])

    


def belongsC1(point,B,R=1,tol=10e-12):
    '''
    checking if a point belongs to circunference1  of the form
    (x-h)^2+y^2=1
    with a tolerance tol set to 10e-12
    IN:
    point: checking point
    B: distante from center to center
    R: radii (set to 1)
    tol: tolerance (set to 10e-12)
    OUT:
    val: Boolean that indicates if the point belongs to C1
    '''
    x,y=point[0],point[1]
    d = (x+B)**2 + y**2 -1
    if d<=tol:
        val = True
    else:
        val = False
    return val

def belongsC2(point,B,R=1,tol=10e-12):
    '''
    checking if a point belongs to a circunferences of the form
    (x-h)^2+y^2=1
    with a tolerance tol
    IN:
    point: checking point
    B: distante from center to center
    R: radii (set to 1)
    tol: tolerance (set to 10e-12)
    OUT:
    val: Boolean that indicates if the point belongs to C2
    '''
    x,y=point[0],point[1]
    d = (x-B)**2 + y**2 -1
    if d<=tol:
        val = True
    else:
        val = False
    return val



def tangent_vector_circunference(P,center):
    '''
    Given a point of a circunference, it returns the anticlockwise unitary tangent vector
    IN:
    P: point of circunference (x,y)
    center: x-coordinate of the center of the circunferences
    OUT:
    tanget_vector: tangent vector
    '''
    if P[0]>0:
        radio_vector = P - np.array([-center,0])
        tangent_vector = rotation(radio_vector,np.pi/2)
    else:
        radio_vector = P - np.array([center,0])
        tangent_vector = rotation(radio_vector,np.pi/2)
        
    return unitary(tangent_vector)

def tangent_lemon(P,center1,center2,x_intersection):
    '''
    Returns the unitary tangent vector on the lemon
    IN:
    P: point
    center1: center of circunference1
    center2: center of circunference2
    x_intersection: vertical line on lemon
    OUT:
    unitary(T): unitary tangent vector
    '''
    x = P[0]
    if x > x_intersection:
        R1 = P - center1
        T = rotation(R1,np.pi/2)
    else: #x<x_intersection
        R2 = P - center2
        T = rotation(R2,np.pi/2)
    return unitary(T)

def velocity_lemon(P,center1,center2,x_intersection,alpha):
    '''
    Gives the velocity vector of a particle afeer the colission with the lemon wall
    IN:
    P:
    center1:
    center2:
    x_intersection: vertical line of the lemon
    OUT:
    v: reflexion velocity vector
    '''
    T = tangent_lemon(P,center1,center2,x_intersection)
    v = rotation(T,alpha)
    return v

def alpha_angle(point,v,B):
    '''
    Alpha angle is the angle between the velocity and the tangent on that point
    IN:
    point: point of collision
    v: velocity vector
    B: distante between centers
    OUT:
    alpha: angle between velocity and tangent
    '''
    if point[0] < 0: #left circunference
        T = tangent_vector_circunference(point,B)
    else: #right circunference
        T = tangent_vector_circunference(point,B)
    alpha = angle_bw_vectors(v,T)
    
    return alpha

def alpha_angle_asymmetrical(P,v,center1,center2,x_intersection):
    '''
    Measuring the angle between the velocity vector and the tangent on the lemon
    IN:
    P: point of collision
    v: velocity vector
    center1: center of circunference1
    center2: center of circunference1
    x_intersection: vertical lemon axis
    OUT:
    alpha: angle alpha
    
    '''
    T = tangent_lemon(P,center1,center2,x_intersection)
    return angle_bw_vectors(v,T)
    



def angle_bw_vectors(v1,v2):
    '''
    Measure the angle between two vectors
    IN:
    v1: vector 1
    v2: vector 2
    OUT:
    angle: angle between vectors
    '''
    angle = np.arccos((v1[0]*v2[0]+v1[1]*v2[1])/(module(v1)*module(v2)))
    return angle



def reflexVelocity_asymmetric(hit,velocity,center1,center2,x_intersection):
    '''
    reflexion velocity on the asymmetrical lemon
    v-(2v.n)n
    IN:
    hit: collision point
    velocity: velocity vector
    center1: center of circunference1
    center2: center of circunference2
    x_intersection: vertical axis on lemon
    OUT:
    reflex: reflexion velocity vector
    '''
    if hit[0]>x_intersection: #lemon right side
        normal =unitary(center1 - hit)
    elif hit[0]<x_intersection: #lemon left side
        normal =unitary(center2 - hit)
    else:
        print('problems on reflexVelocity function')
        
    reflex = np.float128(velocity - (np.dot(2*velocity,normal))*normal) #vector velocidad reflexion
        
    return reflex	



    
def xy2s_asymmetric_2(P,r1,r2,center1,center2,DownLemon,UpLemon,LeftLemon,RightLemon,x_intersection):
    '''
    Function to convert (x,y) coordinates to s coordinates. Version 2 of this function
    IN:
    P: point on lemon
    r1: radii 1
    r2: radii2
    center1: center 1
    center2: center 2
    DownLemon: down lemon point
    UpLemon: up lemon point
    LeftLemon: left lemon point
    RightLemon: right lemon point
    x_intersection: vertical axis of the lemon
    OUT:
    s: s coordinates
    '''
    #print(20*'---')
    #print('We are on xy2s_asymetric_2() function')
    center_lemon = np.array([x_intersection,0])
    beta1 = angle_bw_vectors( RightLemon - center1, UpLemon - center1)
    beta2 = angle_bw_vectors(UpLemon-center2,LeftLemon-center2)
    stotal = r1*2*beta1 + r2*2*beta2
    #print('stotal',stotal)
    #print('beta1',np.rad2deg(beta1))
    #print('beta2',np.rad2deg(beta2))
    if P[0] > x_intersection:
        if P[1] <0: #region I
          #print('This is region I ')
          a = angle_bw_vectors(DownLemon-center1,P-center1)
          OD = angle_bw_vectors(DownLemon-center1,P-center1)
          s = r1 * OD
          #print('a:',np.rad2deg(a))
          #print('angle is',np.rad2deg(OD),'s is',s/stotal)
        else: #region II
          #print('This is region II')
          a=angle_bw_vectors(RightLemon-center1,P-center1)
          #print('a:',np.rad2deg(a))
          OD = angle_bw_vectors(RightLemon-center1,P-center1) + beta1
          s = r1 * OD
          #print('angle is',np.rad2deg(OD),'s is',s/stotal)
            
    else:
        if P[1] > 0: #region III
          #print('This is region III')
          a =  angle_bw_vectors(UpLemon - center2,P-center2)
          #print('a:',np.rad2deg(a))
          OD = angle_bw_vectors(UpLemon - center2,P-center2) 
          s = r1*2*beta1 + r2*OD
          #print('angle is',np.rad2deg(OD),'s is',s/stotal)

        else: #region IV
          #print('This is region IV')
          a = angle_bw_vectors(LeftLemon - center2,P-center2)
          #print('a:',np.rad2deg(a))
          OD = angle_bw_vectors(LeftLemon - center2,P-center2)
          s = r1*2*beta1 + r2*beta2 + r2*OD
          #print('angle is',np.rad2deg(OD),'s is',s/stotal)
    return s/stotal



def belongs_asymmetrical_lemon(point,position,UpLemon,DownLemon,LeftLemon,RightLemon,x_intersection,c1,r1,c2,r2):
    '''
    checking if points belongs on asymmetrical lemon
    IN:
    point: checking point
    position: position of particle
    UpLemon: up lemon point
    DownLemon: down lemon point
    LeftLemon: left lemon point
    RightLemon: right lemon point
    x_intersection: vertical lemon axis
    c1: center1
    r1: radii 1
    c2: center 2
    r2: radii 2
    OUT:
    is_hit: boolean T particle collides F particle did not collide
    '''
    #is this point the position?
    if (np.abs(point[0] - position[0]) < 10e-10) and (np.abs(point[1] -position[1]) < 10e-10):
        #print('i am the position')
        is_hit = False
    else:
        #print('i am possible a collision point')
        # am i on the x- range domain?
        if point[0] > LeftLemon[0] and point[0] < RightLemon[0]:
            #print('i am on the x domain ')
            if point[0] < x_intersection:  #but i am the correct one?
                #print('I need to be on C2 circunference')
                if np.abs((point[0]-c2[0])**2 + point[1]**2 - r2**2) < 10e-10:
                    #print('i am on C2 ===> :) collision point DETECTED ')
                    is_hit = True
                else:
                    is_hit = False
            elif point[0] > x_intersection:
                #print('I need to be on C1 circunference ')
                if np.abs((point[0]-c1[0])**2 + point[1]**2 - r1**2) < 10e-10:
                    #print('i am on C1 ===> :) collision point DETECTED ')
                    is_hit = True
                else:
                    is_hit = False
            else:
                print('PROBLEMS!!!!')
                #hit = [None,None]
        else:
            #print('i am not a collision points')
            is_hit = False
    return is_hit


def s2xy_asym(s,DownLemon,RightLemon,UpLemon,LeftLemon,s_total,s_right,s_left,R2,center2,R1=1):
    '''
    LAST TRY TO s2xy
    :(
    
    '''
    right_rate = s_right/s_total
    left_rate = s_left/s_total
    if s < right_rate: #irght side
        #print('RIGHT SIDE')
        angle = (s*s_total)/R1
        pos = rotation(DownLemon,angle)
    else: #left side
        #print('LEFT SIDE')
        angle = ((s-right_rate)*s_total)/R2
        pos = rotation(UpLemon-center2,angle) + np.array([center2[0],0])
    return pos
    
    
def step(s,alpha,s_exit,h,MAX_ITER):
    '''
    THIS IS ASYMMETRIC LEMON BILLIARDS
    CHECKING FOR THE R1=R2=1 case
    '''

    #scenario
    #plt.hlines(y=0,xmin=-2,xmax=2,color='red')
    #plt.vlines(x=0,ymin=-2,ymax=2,color='red')

    #circunferences
    center1 = np.array([0,0])
    center2 = np.array([0.2,0])
    #print('B=',np.abs(center1[0])+np.abs(center2[0]))
    r1 = 1
    r2 = 1
    C1 = circunference(center1,r1)
    C2 = circunference(center2,r2)

    #plt.plot(C1[0],C1[1],color='blue')
    #plt.plot(C2[0],C2[1],color='blue')

    x_intersection = asymmetric_x_intersection(center1,r1,center2,r2)
    y_intersection = (r1**2 - (x_intersection - center1[0])**2)**0.5  #this has two values and we are going to use de + and the -

    UpLemon = [x_intersection, y_intersection]
    DownLemon = [x_intersection, -y_intersection]

    RightLemon = right(center1,r1)

    LeftLemon = left(center2,r2)

    fix1 = DownLemon - center1
    fix2 = UpLemon - center2
    s_angle_right =angle_bw_vectors(DownLemon-center1,UpLemon-center1)
    s_angle_left =angle_bw_vectors(UpLemon-center2,DownLemon-center2)
    s_right = r1*s_angle_right
    s_left = r2*s_angle_left
    s_total = s_right + s_left

    #position

    s_initial = s
    P = s2xy_asym(s,DownLemon,RightLemon,UpLemon,LeftLemon,s_total,s_right,s_left,r2,center2)

    #plt.plot(P[0],P[1],'kX')

    #velocity
    alpha_initial = alpha
    v = velocity_lemon(P,center1,center2,x_intersection,alpha)


    #plot hole
    P_exit1 = s2xy_asym(s_exit+h,DownLemon,RightLemon,UpLemon,LeftLemon,s_total,s_right,s_left,r2,center2)
    P_exit2 = s2xy_asym(s_exit-h,DownLemon,RightLemon,UpLemon,LeftLemon,s_total,s_right,s_left,r2,center2)
    #plt.plot(P_exit1[0],P_exit1[1],'ko')
    #plt.plot(P_exit2[0],P_exit2[1],'ko')


    for i in range(MAX_ITER):

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

        plt.plot([P[0],hit[0]],[P[1],hit[1]],color='skyblue',alpha = 0.6)
        #update
        P = hit
        v = vReflex

        #print('hit',hit)
        s = xy2s_asymmetric_2(P,r1,r2,center1,center2,DownLemon,UpLemon,LeftLemon,RightLemon,x_intersection)
        alpha_final = alpha_angle_asymmetrical(hit,vReflex,center1,center2,x_intersection)
        #print(s,alpha_final)

        #SURVIVE STUFF
        if s > s_exit-h and s < s_exit + h:
            print(s_initial,alpha_initial,i+1,'SCAPED')
            break


    if i == MAX_ITER-1 :
        print(s_initial,alpha_initial,i+1,'SURVIVOR')

