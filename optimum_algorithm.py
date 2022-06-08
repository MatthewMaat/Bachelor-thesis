from greedy_algorithm import *

#### File to investigate some things regarding optimal orderings
RUNSOMETHINGS=False #run the optimum algorithm, compare with some orderings
MAKECONTOURPLOT=False #load data from the experiment calculating f_X,N, and make a plot
COUNTOPTIMUM=True #run the optimum algorithm a lot of times, count how many times different orderings yield the optimum


def load_and_plot_with_hyperbola(N,showhyperbola=False):
    b=np.load('optimum-7-150000.npz')
    boxaverage=b['boxaverage']
    xv, yv = np.meshgrid(np.linspace(0.005, 0.995, 100), np.linspace(0.005, 0.995, 100))
    fig, ax = plt.subplots()
    CS = ax.contourf(xv, yv, boxaverage,np.arange(0,N,0.5))
    #ax.set_title('Contour plot of average number in optimum ordering per position.')
    plt.ylim(0.0,1.0)
    cbar = fig.colorbar(CS)
    X=np.linspace(0.0,0.99,100)
    if showhyperbola:
        for c in np.linspace(0.05,1.95,39):
            Y=[1.1-c/(1.1-x) for x in X]
            #Y=[1-1/(1/c-1/(1-x)) for x in X]
            #Y=[(c-math.sqrt(1-x)>0)*(1-(c-math.sqrt(1-x))**2) for x in X]
            plt.plot(X,Y,color='black',linewidth=1)
    #cbar.add_lines(CS)
    plt.show()

def count_how_many_optimum(N,n):
    random_order=0
    norm_order=[0,0,0,0,0,0,0] #-inf, -2,-1,0,1,2,inf
    inverse_norm_order=[0,0,0,0,0,0,0] #-inf,-2,-1,0,1,2,inf
    x_order=0
    y_order=0
    diagonal_order=0
    start = time()
    for ii in range(n):
        S=generate_points_exp(N)
        sorted_x = sort_by_x(S)
        Opt=optimum_algorithm(S)[0]
        if math.fabs(greedy_algorithm(S,sorted_x)[0]-Opt)<1e-13:
            random_order+=1
        if math.fabs(greedy_algorithm(sort_by_norm(S,-math.inf)[::-1],sorted_x)[0]-Opt)<1e-13:
            norm_order[0]+=1
        if math.fabs(greedy_algorithm(sort_by_norm(S,-2)[::-1],sorted_x)[0]-Opt)<1e-13:
            norm_order[1]+=1
        if math.fabs(greedy_algorithm(sort_by_norm(S,-1)[::-1],sorted_x)[0]-Opt)<1e-13:
            norm_order[2]+=1
        if math.fabs(greedy_algorithm(sort_by_norm(S,0)[::-1],sorted_x)[0]-Opt)<1e-13:
            norm_order[3]+=1
        if math.fabs(greedy_algorithm(sort_by_norm(S,1)[::-1],sorted_x)[0]-Opt)<1e-13:
            norm_order[4]+=1
        if math.fabs(greedy_algorithm(sort_by_norm(S,2)[::-1],sorted_x)[0]-Opt)<1e-13:
            norm_order[5]+=1
        if math.fabs(greedy_algorithm(sort_by_norm(S,math.inf)[::-1],sorted_x)[0]-Opt)<1e-13:
            norm_order[6]+=1
        if math.fabs(greedy_algorithm(sort_by_inverse_norm(S,-math.inf)[::-1],sorted_x)[0]-Opt)<1e-13:
            inverse_norm_order[0]+=1
        if math.fabs(greedy_algorithm(sort_by_inverse_norm(S,-2)[::-1],sorted_x)[0]-Opt)<1e-13:
            inverse_norm_order[1]+=1
        if math.fabs(greedy_algorithm(sort_by_inverse_norm(S,-1)[::-1],sorted_x)[0]-Opt)<1e-13:
            inverse_norm_order[2]+=1
        if math.fabs(greedy_algorithm(sort_by_inverse_norm(S,0)[::-1],sorted_x)[0]-Opt)<1e-13:
            inverse_norm_order[3]+=1
        if math.fabs(greedy_algorithm(sort_by_inverse_norm(S,1)[::-1],sorted_x)[0]-Opt)<1e-13:
            inverse_norm_order[4]+=1
        if math.fabs(greedy_algorithm(sort_by_inverse_norm(S,2)[::-1],sorted_x)[0]-Opt)<1e-13:
            inverse_norm_order[5]+=1
        if math.fabs(greedy_algorithm(sort_by_inverse_norm(S,math.inf)[::-1],sorted_x)[0]-Opt)<1e-13:
            inverse_norm_order[6]+=1
        if math.fabs(greedy_algorithm(sort_by_x(S)[::-1],sorted_x)[0]-Opt)<1e-13:
            x_order+=1
        if math.fabs(greedy_algorithm(sort_by_y(S)[::-1],sorted_x)[0]-Opt)<1e-13:
            y_order+=1
        if math.fabs(greedy_algorithm(sort_by_diagonal(S)[::-1],sorted_x)[0]-Opt)<1e-13:
            diagonal_order+=1
        if ii%1000==0:
            finish=time()
            print(ii,'sets have been processed. It took', finish-start,'seconds.')
            start=time()
    np.savez('ordercount',x_order=x_order,y_order=y_order,norm_order=norm_order,inverse_norm_order=inverse_norm_order,diagonal_order=diagonal_order,random_order=random_order)
    print('Random:',random_order)
    print('X:',x_order)
    print('Y:',y_order)
    print('Norm:',norm_order)
    print('Inverse norm:',inverse_norm_order)
    print('Diagonal:',diagonal_order)
    return


if __name__=="__main__":
    if RUNSOMETHINGS:
        S=generate_points(20)
        sorted_x=sort_by_x(S)
        print('The optimum algorithm yields',optimum_algorithm(S)[0])
        print('Random order yields',greedy_algorithm(S,sorted_x)[0])
        print('Normal hyperbolasort yields', greedy_algorithm(sort_by_inverse_norm(S, 0)[::-1], sorted_x)[0])
        print('Weird hyperbolasort yields',greedy_algorithm(weird_sort_by_inverse_norm(S,0)[::-1],sorted_x)[0])
        print('Harmonic sort yields', greedy_algorithm(sort_by_inverse_norm(S, -1)[::-1], sorted_x)[0])
    if MAKECONTOURPLOT:
        load_and_plot_with_hyperbola(7,True)
    if COUNTOPTIMUM:
        count_how_many_optimum(10,200000)
