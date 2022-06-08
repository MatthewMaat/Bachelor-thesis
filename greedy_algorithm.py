import math
import random
import copy
from perm import *
import matplotlib.pyplot as plt
import matplotlib.patches as ptc
import itertools
from time import time
import numpy as np

### This file contains a lot of functions concerning greedy and optimal algorithms for anchored rectangle packings
### Things that can be done by running this file:
MAKEPLOT=False #run the greedy and optimum algorithm on a random set, make a plot, for 100 times
INVESTIGATEORDER=True #run the optimum algorithm a lot of times, make a contour plot of f_X,N
TESTOPTIMUM=False #run the DP against an algorithm that tries all cominations to test if it works



###########Point generating functions
def generate_points(N):
    #generate a random (uniformly distributed) set of points of size N
    point_list=[(0,0)]
    for i in range(N-1):
        point_list.append((random.random(),random.random()))
    return point_list

def generate_points_triangle(N):
    #generate a random (triangularly distributed) set of points of size N
    point_list=[(0,0)]
    for i in range(N-1):
        point_list.append((random.triangular(0,1,0),random.triangular(0,1,0)))
    return point_list

def generate_points_exp(N):
    #generate a random (triangularly distributed) set of points of size N
    point_list=[(0,0)]
    for i in range(N-1):
        xpp,ypp=2,2
        while xpp>1:
            xpp=random.expovariate(5)
        while ypp>1:
            ypp=random.expovariate(5)
        point_list.append((xpp,ypp))
    return point_list

###########Sorting functions

def sort_by_x(point_list): #simple bubblesort on the x-coordinate
    sorted_x = point_list.copy()
    for i in range(len(sorted_x)-1):
        for j in range(len(sorted_x)-1):
            if sorted_x[j][0]>sorted_x[j+1][0]:
                sorted_x[j],sorted_x[j+1]=sorted_x[j+1],sorted_x[j]
    return sorted_x

def sort_by_y(point_list): #simple bubblesort on the x-coordinate
    sorted_y = point_list.copy()
    for i in range(len(sorted_y)-1):
        for j in range(len(sorted_y)-1):
            if sorted_y[j][1]>sorted_y[j+1][1]:
                sorted_y[j],sorted_y[j+1]=sorted_y[j+1],sorted_y[j]
    return sorted_y

def pointnorm(p,a):
    #norm of point
    if p==(0,0):
        return 0
    elif a==-math.inf:
        return min(p[0],p[1])
    elif a==0:
        return math.sqrt(p[0]*p[1])
    elif a==math.inf:
        return max(p[0],p[1])
    else:
        return pow((pow(p[0],a)+pow(p[1],a))/2,1/a)

def diagonal_value(p):
    if p==(0,0):
        return 0
    return (p[0]+p[1]-math.fabs(p[0]-p[1]))/(1-math.fabs(p[0]-p[1]))

def inverse_pointnorm(p,a):
    return -pointnorm((1-p[0],1-p[1]),a)

def sort_by_norm(point_list,a):
    sorted_x = point_list.copy()
    for i in range(len(sorted_x) - 1):
        for j in range(len(sorted_x) - 1):
            if pointnorm(sorted_x[j],a) > pointnorm(sorted_x[j + 1],a):
                sorted_x[j], sorted_x[j + 1] = sorted_x[j + 1], sorted_x[j]
    return sorted_x

def sort_by_diagonal(point_list):
    sorted_x = point_list.copy()
    for i in range(len(sorted_x) - 1):
        for j in range(len(sorted_x) - 1):
            if diagonal_value(sorted_x[j]) > diagonal_value(sorted_x[j + 1]):
                sorted_x[j], sorted_x[j + 1] = sorted_x[j + 1], sorted_x[j]
    return sorted_x

def sort_by_inverse_norm(point_list,a):
    sorted_x = point_list.copy()
    for i in range(len(sorted_x) - 1):
        for j in range(len(sorted_x) - 1):
            if inverse_pointnorm(sorted_x[j],a) > inverse_pointnorm(sorted_x[j + 1],a):
                sorted_x[j], sorted_x[j + 1] = sorted_x[j + 1], sorted_x[j]
    return sorted_x

def weird_inverse_pointnorm(p,a):
    return -pointnorm((1.1-p[0],1.1-p[1]),a)


def weird_sort_by_inverse_norm(point_list,a):
    sorted_x = point_list.copy()
    for i in range(len(sorted_x) - 1):
        for j in range(len(sorted_x) - 1):
            if weird_inverse_pointnorm(sorted_x[j],a) > weird_inverse_pointnorm(sorted_x[j + 1],a):
                sorted_x[j], sorted_x[j + 1] = sorted_x[j + 1], sorted_x[j]
    return sorted_x

#########Algorithms

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

def make_string_from_set(listzz):
    listzz=list(listzz)
    listzz=sort_by_x(listzz)
    return str(listzz)


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
    sets={make_string_from_set([p]):[p] for p in setz}
    for i in range(2,N+1):
        newsets=list(itertools.combinations(setz,i)) #all subsets of size i
        sets2={}
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
                    orderz=sets[make_string_from_set(setminust)]+[t]
                    val=greedy_algorithm(righttop+orderz,sort_by_x(righttop+orderz))[0]
                    if val>bestval:
                        bestval=val
                        bestorder=orderz
            sets2[make_string_from_set(s)]=bestorder
        sets=sets2
    if len(sets)==1:
        bestorder=setz
    if len(sets)>0:
        return greedy_algorithm(righttop+sets[make_string_from_set(bestorder)]+[(0,0)],sort_by_x(righttop+sets[make_string_from_set(bestorder)]+[(0,0)]))[0],righttop+sets[make_string_from_set(bestorder)]+[(0,0)]
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

###########Plotting functions

def fancyplot(S,orderzz,rectangles1,rectangles2):
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
    for i in range(len(S)):
        plt.text(S[::-1][i][0],S[::-1][i][1],i+1)
    plt.subplot(122)
    currentAxis = plt.gca()
    plt.plot([s[0] for s in S], [s[1] for s in S], 'o', color='black')
    plt.plot([0, 1, 1, 0, 0], [0, 0, 1, 1, 0], color='grey')
    plt.title('Optimum ordering')
    for r in rectangles2:
        currentAxis.add_patch(ptc.Rectangle(r[0], r[1][0] - r[0][0], r[1][1] - r[0][1], facecolor='red'))
        currentAxis.add_patch(ptc.Rectangle(r[0], r[1][0] - r[0][0], r[1][1] - r[0][1], fill=None))
    for i in range(len(S)):
        plt.text(orderzz[i][0],orderzz[i][1],i+1)
    plt.show()

def fancyplot2(S,orderzz,rectangles2):
    #plot the rectangles with the optimum ordering
    plt.figure()
    currentAxis = plt.gca()
    plt.plot([s[0] for s in S], [s[1] for s in S], 'o', color='black')
    plt.plot([orderzz[i][0] for i in range(len(S))], [orderzz[i][1] for i in range(len(S))], '-', color='blue',linewidth=1)
    plt.plot([0, 1, 1, 0, 0], [0, 0, 1, 1, 0], color='grey')
    plt.title('Optimum ordering')
    for i in range(len(S)):
        plt.text(orderzz[i][0],orderzz[i][1],i+1)
    for r in rectangles2:
        currentAxis.add_patch(ptc.Rectangle(r[0], r[1][0] - r[0][0], r[1][1] - r[0][1], facecolor='red'))
        currentAxis.add_patch(ptc.Rectangle(r[0], r[1][0] - r[0][0], r[1][1] - r[0][1], fill=None))
    plt.show()


def investigate_order(N,n,sizet):
    boxcounter=[[0.0 for i in range(100)] for j in range(100)]
    boxtotal=[[0.0 for i in range(100)] for j in range(100)]
    start=time()
    for ii in range(n):
        """topmatch=False
        while not topmatch:
            S=generate_points(N)
            righttop=[]
            for p in S:
                prighttop = True
                for q in S:
                    if q[0] > p[0] and q[1] > p[1]:
                        prighttop = False
                if prighttop:
                    righttop.append(p)
            if len(righttop)==sizet:
                topmatch=True
        #print(len(righttop))"""
        S=generate_points_exp(N)
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
    np.savez('optimum', boxaverage=boxaverage)
    plt.show()
    return boxaverage

if __name__=="__main__":

    if MAKEPLOT:
        for ghbjn in range(100):
            S=generate_points(20)
            sorted_x = sort_by_x(S)
            A,rectangles=greedy_algorithm(S,sorted_x)
            start=time()
            B,orderzz=optimum_algorithm(S)
            finish=time()
            print('The maximum area is',B)
            print('The normal order yields an area of',A)
            print('It took',finish-start,'seconds for',len(orderzz),'points.')
            fancyplot(S,orderzz,rectangles,greedy_algorithm(orderzz,sorted_x)[1]) #compare random and optimum
            fancyplot2(S,orderzz,greedy_algorithm(orderzz,sorted_x)[1]) #look at optimum order
    if INVESTIGATEORDER:
        print(investigate_order(10,150000,2))
    if TESTOPTIMUM:
        for dsf in range(1000):
            S=generate_points(7)
            a=factorial_algorithm(S)[0]
            b=optimum_algorithm(S)[0]
            if math.fabs(a-b)>1e-14:
                print('ERRORRR')
                break
            if dsf%50==0:
                print(dsf,'factorial:',a,'dp:',b)