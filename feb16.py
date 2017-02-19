
import myfunc as mf
import numpy as np
import time
import matplotlib.pyplot as plt


import myfunc
import importlib
importlib.reload(myfunc)


# test data set from lecture ###########################################################################
points = np.array([0, 1, 2, 3, 4, 5, 8, 11, 14, 17, 20, 21])
mf.plot_delaunay1d(points)
mf.plot_voronoi1d(points)
[N,cell_edges]=mf.voronoi1d(points)
A=np.diff(cell_edges)

blocksOP=mf.optinterval(N,A,.2)
blocksPELT=mf.optintervalPELT(N,A,.2)
blocksBS=mf.BinarySeg(N,A,0.2)


# second test data set from lecture ########################################################################
points = np.array([0, 1,1.6,2.4,4.4,5.2,7.4,10.8,12.4,13,14.6,18.4,19.4,20])
mf.plot_delaunay1d(points)
mf.plot_voronoi1d(points)
[N,cell_edges]=mf.voronoi1d(points)
A=np.diff(cell_edges)

blocksOP=mf.optinterval(N,A,.2)
blocksPELT=mf.optintervalPELT(N,A,.2)
blocksBS=mf.BinarySeg(N,A,0.2)




# test data set 1 ###########################################################################
points=  np.hstack((0,np.cumsum(np.arange(1,21)),250))
mf.plot_voronoi1d(points)
[N,cell_edges]=mf.voronoi1d(points)
A=np.diff(cell_edges)

blocksOP=mf.optinterval(N,A,.5)
blocksPELT=mf.optintervalPELT(N,A,.5)
blocksBS=mf.BinarySeg(N,A,0.5)










# test data set 2 ###########################################################################
points=np.array(np.cumsum(np.hstack((range(20),range(1,20)[::-1]))))
mf.plot_delaunay1d(points)
mf.plot_voronoi1d(points)
[N,cell_edges]=mf.voronoi1d(points)
A=np.diff(cell_edges)

blocksOP=mf.optinterval(N,A,.5)
blocksPELT=mf.optintervalPELT(N,A,.5)
blocksBS=mf.BinarySeg(N,A,0.5)








# test set 2 (random points) #####################################################################
points=np.sort(np.random.random(80))
mf.plot_delaunay1d(points)
mf.plot_voronoi1d(points)
[N,cell_edges]=mf.voronoi1d(points)
A=np.diff(cell_edges)

blocksOP=mf.optinterval(N,A,.5)
blocksPELT=mf.optintervalPELT(N,A,1.5)
blocksBS=mf.BinarySeg(N,A,0.5)






# test set 3 (simple simulated data) #######################################
points=mf.simplePPSim([2,10,4],[4,4,4])
mf.plot_delaunay1d(points)
mf.plot_voronoi1d(points)
[N,cell_edges]=mf.voronoi1d(points)
A=np.diff(cell_edges)

blocksOP=mf.optinterval(N,A,2)
blocksPELT=mf.optintervalPELT(N,A,2)
blocksBS=mf.BinarySeg(N,A,2)




################# Plotting ############################################################################################




#Basic Block Segmentation Plot ############

block_object = blocksPELT #PELT is used for this example
plt.close()
for i in range(len(A)):
    point,=plt.plot([cell_edges[i],cell_edges[i+1]],[1/A[i],1/A[i]],c="black",lw=2)
    plt.ylim([0,np.percentile(1/A,90)*1.2])
plt.xlim([points[0],points[-1]])
numblocks=len(block_object["changepoints"])
cp=np.hstack((block_object["changepoints"],len(A)))
ints=block_object["intensities"]
for j in range(numblocks):
    block,=plt.plot([cell_edges[cp[j]],cell_edges[cp[j+1]]],[ints[j],ints[j]],lw=5,c="green")
plt.legend([point,block], ["Voronoi cell Intensity","Block Intensity"])
plt.xlabel("Time")
plt.ylabel("Intensity")
plt.title("Block Segmentation")



## Pruned Log-Likelihood Plots #################

# Gray shaded areas represent pruned points. no likelihood computation is done for these points
rows=3
cols=4
plot_points=np.linspace(1,len(N)-1,rows*cols).astype(int)
if len(N) == 12:
    plot_points = range(12)
f, axarr = plt.subplots(rows, cols)
k = 0
for i in range(rows):
    for j in range(cols):
        LLs=blocksPELT["LLs"][plot_points[k]]
        axarr[i, j].plot(LLs)
        axarr[i, j].set_title(str(plot_points[k] + 1) + " points")
        axarr[i, j].set_xlim([0, len(A)])
        axarr[i,j].set_ylim()
        for a in range(len(LLs)):
            if LLs[a]==-np.inf:
                axarr[i, j].axvspan(a, a+1, color='grey', alpha=0.5, lw=0)
        axarr[i,j].axvline(x=plot_points[k],c="black",lw=2,ls="dashed")
        axarr[i, j].yaxis.set_major_formatter(plt.NullFormatter())
        axarr[i, j].xaxis.set_major_formatter(plt.NullFormatter())
        k += 1



##################### Testing the Computation Time for each Algorithm ##############################



#### run the tests##########################################
#set the number of cells to stop at
max_points=500

#set the number of tests you want to run
num_tests=20

tests=np.linspace(3,max_points,num_tests).astype(int)
BS_results=np.array([0,0])
OP_results=np.array([0,0])
PELT_results=np.array([0,0])
for i in range(len(tests)):
    start=time.clock()
    results=mf.one_run(tests[i],0.3)
    BS_results=np.vstack((BS_results,results["BS"]))
    OP_results = np.vstack((OP_results, results["OP"]))
    PELT_results = np.vstack((PELT_results, results["PELT"]))
    print("Iteration", i, ": ","{0:.3f}".format(time.clock()-start)," seconds")
print("done")

#plot the computation time results #################################
BSline, = plt.plot(tests,BS_results[1:,0], label='BS',lw=2)
PELTline, = plt.plot(tests,PELT_results[1:,0], label='PELT',lw=2)
OPline, = plt.plot(tests,OP_results[1:,0], label='OP',lw=2)
plt.legend([OPline,PELTline,BSline], ['Optimal Partitioning', 'OP + PELT','Binary Segmentation'],loc=2)
plt.xlabel("number of random cells")
plt.ylabel("Seconds")
plt.title("Algorithm Computation Time Comparison")


#plot the ratio of blocks to cells on each test#######################
# (to verify that we are estimated approximately the same number of blocks each time)
BSline, = plt.plot(tests,BS_results[1:,1], label='BS',lw=2)
PELTline, = plt.plot(tests,PELT_results[1:,1], label='PELT',lw=2)
OPline, = plt.plot(tests,OP_results[1:,1], label='OP',lw=2)
plt.legend([BSline,PELTline,OPline], ['Binary Segmentation', 'OP + PELT','Optimal Partitioning'],loc=2)
plt.ylim((0,0.5))
plt.xlabel("number of random cells")
plt.ylabel("proportion of blocks to cells")
plt.title("block to cell ratio")



