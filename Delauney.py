import matplotlib
import numpy as np
import numpy.ma as ma
import random
import math


def main():
    TestArray=np.random.uniform(0,30,(100,2))
    isX=True
    Tris=Delauney2D(np.array([]),TestArray,np.array([]),isX)


    plot(TestArray,Tris)

def Delauney2D(AFL,P,TriList,isX):
    simplex_list=[]

    P1,P2,spliting_plane,x_range,y_range=PointsetPartion(P,isX)

    if len(AFL)==0:
        tri=MakeInitialSimplex(P1,P2,isX,spliting_plane,x_range,y_range)
        return tri
        P=p_updated
        if tri!=None:
            return None
        else:
            AFL.append(GetFaces(tri))
            TriList.append(tri)

    AFL_alpha,AFL_1,AFL_2=SplitFL(AFL)
    while(len(AFL_alpha>0)):
        wallFace=AFL_alpha[0]
        AFL_alpha=AFL_alpha[1:]
        tri_2,p_updated=MakeSimplex(wallFace,P)
        P=p_updated

        if tri_2!=None:
            TriList.append(tri_2)
            temp_faces=GetFaces(tri_2)
            temp_fl_alpha,temp_fl_1,temp_fl_2=SplitFL(temp_faces)

            if len(temp_fl_alpha)>0:
                AFL_alpha.append(temp_fl_alpha)
            if len(temp_fl_1>0):
                AFL_1.append(temp_fl_1)
            if len(temp_fl_2>0):
                AFL_2.append(temp_fl_2)

    if len(AFL_1)>0:
        TriList.append(Delauney2D(AFL_1,P1,TriList,not isX))
    if len(AFL_2>0):
        TriList.append(Delauney2D(AFL_2,P2,TriList,not isX))

    return TriList

def PointsetPartion(P,isX):
    P1=None
    P2=None
    min_x=10000
    min_y=10000
    max_x=-10000
    max_y=-10000

    if isX:
        line_mid=np.median(P[:len(P),0])
    else:
        line_mid=np.median(P[:len(P),1])

    for point in P:
        x=point[0]
        y=point[1]

        if x<min_x:
            min_x=x
        elif y<min_y:
            min_y=y
        elif x>max_x:
            max_x=x
        elif y>max_y:
            max_y=y

        if isX:
            if x<line_mid:
                if type(P1) is not np.ndarray:
                    P1=np.array([point])
                else:
                    P1=np.append(P1,[point],axis=0)
            if x>line_mid:
                if type(P2) is not np.ndarray:
                    P2=np.array([point])
                else:
                    P2=np.append(P2,[point],axis=0)
        else:
            if y<line_mid:
                if type(P1) is not np.ndarray:
                    P1=np.array([point])
                else:
                    P1=np.append(P1,[point],axis=0)
            if y>line_mid:
                if type(P2) is not np.ndarray:
                    P2=np.array([point])
                else:
                    P2=np.append(P2,[point],axis=0)

    x_range=(min_x,max_x)
    y_range=(min_y,max_y)

    return P1,P2,line_mid,x_range,y_range

def MakeInitialSimplex(P1,P2,isX,line_mid,x_range,y_range):
    SearchPoints_1=ma.array(P1.copy())
    SearchPoints_2=ma.array(P2.copy())


    if len(P1)+len(P2)<=2:

        return None

    min_distance_1=1000
    min_point_1=None
    min_distance_2=1000
    min_point_2=None
    min_point_index_1=0
    min_point_index_2=0

    reject_threshold=(abs(x_range[1]-x_range[0])/4,abs(y_range[1]-y_range[0])/4)


    for pointIndex in range(0,len(SearchPoints_1)):

        if isX:
            distance=SearchPoints_1[pointIndex][0]-line_mid
            if distance>reject_threshold[0]:
                SearchPoints_1[pointIndex]=ma.masked
        else:
            distance=SearchPoints_1[pointIndex][1]-line_mid
            if distance>reject_threshold[1]:
                SearchPoints_1[pointIndex]=ma.masked

        if distance<min_distance_1:
            min_distance_1=distance
            min_point_1=SearchPoints_1[pointIndex].copy()
            min_point_index_1=pointIndex

    SearchPoints_1[min_point_index_1]=ma.masked

    for pointIndex_2 in range(0,len(SearchPoints_2)):

        if isX:
            distance=SearchPoints_2[pointIndex_2][0]-line_mid
            if distance>reject_threshold[0]:
                SearchPoints_2[pointIndex_2]=ma.masked
        else:
            distance=SearchPoints_2[pointIndex_2][1]-line_mid
            if distance>reject_threshold[1]:
                SearchPoints_2[pointIndex_2]=ma.masked

        if distance<min_distance_2:
            min_distance_2=distance
            min_point_2=SearchPoints_2[pointIndex_2].copy()
            min_point_index_2=pointIndex_2

    SearchPoints_2[min_point_index_2]=ma.masked
    P_SS=ma.array([])

    if abs(min_distance_1)<=abs(min_distance_2):
        P_SS=SearchPoints_2.compressed().reshape(-1,2)
        SS_point=min_point_1
    else:
        P_SS=SearchPoints_1.compressed().reshape(-1,2)
        SS_point=min_point_2

    SS_min_distance=10000
    SS_min_point=np.array([])
    SS_min_index=0
    
    for pointIndex in range(0,len(P_SS)):
        distance=PointDistance(SS_point,P_SS[pointIndex])
        if distance<SS_min_distance:
            SS_min_point=P_SS[pointIndex]
            SS_min_distance=distance
            SS_min_index=0

    if min_distance_1<=min_distance_2:
        SearchPoints_1[SS_min_index]=ma.masked
    else:
        SearchPoints_2[SS_min_index]=ma.masked
    TS_min_radius=10000
    TS_min_point=np.array([])

    for point in np.append(SearchPoints_1,SearchPoints_2,axis=0).compressed().reshape(-1,2):

        center_point,radius=CircumCenter(SS_min_point,SS_point,point)

        if type(center_point)==np.ndarray and radius!=None:

            if radius<TS_min_radius:
                TS_min_radius=radius
                TS_min_point=point



    return np.append(np.append(SS_min_point,SS_point,axis=0),TS_min_point,axis=0)

def SplitFL():
    return None
def CircumCenter(P1,P2,P3):

    midpoint_1=np.true_divide((P1+P2),2)
    midpoint_2=np.true_divide(P1+P3,2)


    x=np.subtract(P2,P1)


    if x[0]==0:
        return None,None


    inv_slope_1=-1/(x[1]/x[0])


    x2=np.subtract(P3,P1)

    if x2[0]==0:
        return None,None
    inv_slope_2=-1/(x2[1]/x2[0])

    center_point=LineIntersection(x,x2,inv_slope_1,inv_slope_2)

    if type(center_point)!=np.ndarray:
        return None,None
    radius=PointDistance(x,center_point)

    return center_point,radius
def GetFaces():

    return None

def LineIntersection(P_1,P_2,slope_1,slope_2):
    yintercept_1=((slope_1*-1)*P_1[1])
    yintercept_2=((slope_2*-1)*P_2[1])

    if slope_1==slope_2:
        return None
    x=(yintercept_1+yintercept_2)/(slope_1-slope_2)
    y=(slope_2*x)-(slope_2*P_2[0])+P_2[1]

    return(np.array([x,y]))

def MakeSimplex():
    return None
def PointDistance(p1,p2):

    return math.sqrt( ((p1[0]-p2[0])**2)+((p1[1]-p2[1])**2) )

def plot(points,tris):

    fig, ax = plt.subplots()
    ax.set_xlim((0, 30))
    ax.set_ylim((0, 30))
    ax.plot(points,'bo')

    print(tris)
    ax.plot(tris)

    plt.show()

def circleTest(x,y,center_x,center_y):
    dx = abs(x-center_x)
    dy = abs(y-center_y)
    R = radius

    if dx+dy <= R:
        return True
    else:
        return False


main()
