import math
import random
import copy
from perm import *
import matplotlib.pyplot as plt
import matplotlib.patches as ptc
import itertools
from time import time
import numpy as np

def generate_points(N):
    #generate a random (uniformly distributed) set of points of size N
    point_list=[(0,0)]
    for i in range(N-1):
        point_list.append((random.random(),random.random()))
    return point_list

def sort_by_x(point_list): #simple bubblesort on the x-coordinate
    sorted_x = point_list.copy()
    for i in range(len(sorted_x)-1):
        for j in range(len(sorted_x)-1):
            if sorted_x[j][0]>sorted_x[j+1][0]:
                sorted_x[j],sorted_x[j+1]=sorted_x[j+1],sorted_x[j]
    return sorted_x

def factorial_algorithm(setz):
    #optimum algorithm, based on just listing all possible permutations, very slow
    #returns maximum area, optimum ordering
    sets=list(itertools.permutations(setz))
    sorted_x=sort_by_x(setz)
    bestval=0
    bestorder=[]
    for s in sets:
        permz=list(s)
        val=greedy_algorithm(permz,sorted_x)[0]
        if val>bestval:
            bestval=val
            bestorder=permz
    return bestval, bestorder

def optimum_algorithm(S):
    #dynamic programming optimum algorithm
    #returns maximum area, optimum ordering
    setz=S.copy()
    setz.pop(0)
    righttop=[]
    for p in setz:
        prighttop=True
        for q in setz:
            if q[0]>p[0] and q[1]>p[1]:
                prighttop=False
        if prighttop:
            righttop.append(p)
    for p in righttop:
        setz.remove(p)
    N=len(setz)
    sets=[[{p},[p],greedy_algorithm([p]+righttop,sort_by_x([p]+righttop))[0]] for p in setz]
    for i in range(2,N+1):
        newsets=list(itertools.combinations(setz,i)) #all subsets of size i
        sets2=[]
        for s in newsets: #range over all subsets
            s=set(s)
            bestval=0
            bestorder=[]
            for t in s: #remove an element from the subset and retrieve its order and value
                tapproved=True
                setminust=s.copy()
                setminust.remove(t)
                for u in setminust: #check if t dominates a point in s
                    if u[0]<t[0] and u[1]<t[1]:
                        tapproved=False
                        break
                if tapproved:
                    for v in sets: #find the result of the set minus t
                        if v[0]==setminust:
                            orderz=v[1]+[t]
                            val=greedy_algorithm(righttop+orderz,sort_by_x(righttop+orderz))[0]
                            if val>bestval:
                                bestval=val
                                bestorder=orderz
                            break
            sets2.append([set(s),bestorder,bestval])
        sets=sets2
    if len(sets)>0:
        return greedy_algorithm(righttop+sets[0][1]+[(0,0)],sort_by_x(righttop+sets[0][1]+[(0,0)]))[0],righttop+sets[0][1]+[(0,0)]
    else:
        return greedy_algorithm(righttop+[(0, 0)], sort_by_x(righttop+[(0, 0)]))[0],righttop+[(0, 0)]




def greedy_algorithm(permutation,sorted_x):
    #greedy algorithm, input: list of points in some order, and sorted list of points by x
    #output: maximum area, list of rectangles
    rectangles=[] #rectangle is of form [(xlb,ylb),(xrt,yrt)]
    A=0
    for p in permutation:
        yprime=1.0 #maximum y-coordinate of right top candidate rectangle
        xprime=1.0 #"" x-coordinate ""
        for r in rectangles:
            if r[0][0]<=p[0] and r[1][0]>p[0] and r[0][1]<yprime and r[0][1]>p[1]:
                yprime=r[0][1]
            if r[0][1]<=p[1] and r[1][1]>p[1] and r[0][0]<xprime and r[0][0]>p[0]:
                xprime=r[0][0]
        AR=0
        righttop=(p[0],yprime)
        for q in sorted_x:
            if q[0]>p[0] and q[0]<xprime and q[1]<=yprime and q[1]>p[1]:
                candidate_area=(yprime-p[1])*(q[0]-p[0])
                if candidate_area>AR:
                    AR=candidate_area
                    righttop=(q[0],yprime)
                yprime=q[1]
            elif q[0]==xprime:
                candidate_area = (yprime - p[1]) * (q[0] - p[0])
                if candidate_area > AR:
                    AR = candidate_area
                    righttop = (q[0], yprime)
                    break
            #print(p,q,righttop,xprime,yprime)
        if xprime==1.0:
            candidate_area=(yprime - p[1]) * (1.0 - p[0])
            if candidate_area > AR:
                AR = candidate_area
                righttop = (1.0, yprime)
        rectangles.append([p,righttop])
        A+=AR
    return A,rectangles

def fancyplot(S,rectangles1,rectangles2):
    #plot the points and rectangles
    plt.figure()
    plt.subplot(121)
    currentAxis = plt.gca()
    plt.plot([s[0] for s in S],[s[1] for s in S],'o',color='black')
    plt.plot([0,1,1,0,0],[0,0,1,1,0],color='grey')
    plt.title('Random ordering')
    for r in rectangles1:
        currentAxis.add_patch(ptc.Rectangle(r[0], r[1][0] - r[0][0], r[1][1] - r[0][1], facecolor='red'))
        currentAxis.add_patch(ptc.Rectangle(r[0],r[1][0]-r[0][0],r[1][1]-r[0][1],fill=None))
    plt.subplot(122)
    currentAxis = plt.gca()
    plt.plot([s[0] for s in S], [s[1] for s in S], 'o', color='black')
    plt.plot([0, 1, 1, 0, 0], [0, 0, 1, 1, 0], color='grey')
    plt.title('Optimum ordering')
    for r in rectangles2:
        currentAxis.add_patch(ptc.Rectangle(r[0], r[1][0] - r[0][0], r[1][1] - r[0][1], facecolor='red'))
        currentAxis.add_patch(ptc.Rectangle(r[0], r[1][0] - r[0][0], r[1][1] - r[0][1], fill=None))
    plt.show()

def fancyplot2(S,orderzz,rectangles2):
    #plot the rectangles with the optimum ordering
    plt.figure()
    currentAxis = plt.gca()
    plt.plot([s[0] for s in S], [s[1] for s in S], 'o', color='black')
    plt.plot([0, 1, 1, 0, 0], [0, 0, 1, 1, 0], color='grey')
    plt.title('Optimum ordering')
    for i in range(len(S)):
        plt.text(orderzz[i][0],orderzz[i][1],i+1)
    for r in rectangles2:
        currentAxis.add_patch(ptc.Rectangle(r[0], r[1][0] - r[0][0], r[1][1] - r[0][1], facecolor='red'))
        currentAxis.add_patch(ptc.Rectangle(r[0], r[1][0] - r[0][0], r[1][1] - r[0][1], fill=None))
    plt.show()

def pointnorm(p,a):
    #norm of point
    if a==-math.inf:
        return min(p[0],p[1])
    elif a==0:
        return math.sqrt(p[0]*p[1])
    elif a==math.inf:
        return max(p[0],p[1])
    else:
        return pow((pow(p[0],a)+pow(p[1],a))/2,1/a)

def compare_norms(N,n):
    for i in range(n):
        S=generate_points(N)
        B,orderzz=optimum_algorithm(S)

def investigate_order(N,n):
    boxcounter=[[0.0 for i in range(100)] for j in range(100)]
    boxtotal=[[0.0 for i in range(100)] for j in range(100)]
    start=time()
    for ii in range(n):
        S=generate_points(N)
        B, orderzz = optimum_algorithm(S)
        for jj in range(len(orderzz)):
            xbox=math.floor(100*orderzz[jj][0])
            ybox=math.floor(100*orderzz[jj][1])
            boxcounter[xbox][ybox]+=1
            boxtotal[xbox][ybox]+=jj+1
        if ii%1000==0:
            finish=time()
            print(ii,'sets have been processed. It took', finish-start,'seconds.')
            start=time()
    boxaverage=[[0.0 for i in range(100)] for j in range(100)]
    for xx in range(100):
        for yy in range(100):
            if boxcounter[xx][yy]>0:
                boxaverage[xx][yy]=boxtotal[xx][yy]/boxcounter[xx][yy]
    xv, yv = np.meshgrid(np.linspace(0.005, 0.995, 100), np.linspace(0.005, 0.995, 100))
    fig, ax = plt.subplots()
    CS = ax.contourf(xv, yv, boxaverage,np.arange(0,N,0.5))
    ax.set_title('Contour plot of average number in optimum ordering per position.')
    cbar = fig.colorbar(CS)
    #cbar.add_lines(CS)
    plt.show()
    return boxaverage

if __name__=="__main__":

    if False:
        S=generate_points(16)
        sorted_x = sort_by_x(S)
        A,rectangles=greedy_algorithm(S,sorted_x)
        start=time()
        B,orderzz=optimum_algorithm(S)
        finish=time()
        print('The maximum area is',B)
        print('The normal order yields an area of',A)
        print('It took',finish-start,'seconds for',len(orderzz),'points.')
        fancyplot(S,rectangles,greedy_algorithm(orderzz,sorted_x)[1]) #compare random and optimum
        fancyplot2(S,orderzz,greedy_algorithm(orderzz,sorted_x)[1]) #look at optimum order
    print(investigate_order(10,150000))
