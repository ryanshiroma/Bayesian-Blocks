
import numpy as np
import matplotlib.pyplot as plt
import numpy.ma as ma
import time

def newobjective(N,M,c): #objective function
    with np.errstate(all='ignore'): #when N or M are zero, ignore the divide by zero warning and return 0
        return np.nan_to_num(N * np.log(np.array(N)/M)-c)





#####
# Optimal Interval function returns an object with the following elements
# blocks["changepoints"] gives you a list of changepoints using Dr Scargle' changepoint convention
# blocks["step_CPs"] gives a list of what the changepoints were after each consecutively added cell
# blocks["maxLL"] gives you the optimal log likelihood after each consecutively added cell
# blocks["LLs"] gives you more detailed information on the log likelihoods at each step(just like in the matlab function output)
# blocks["intensities"] gives you the intensities of each block
# blocks["step_Intens"] gives you the intensities of each block after each consecutively added cell
def optinterval(N,A,c): #optinterval function(should have same results as the matlab function)
    cells = len(N)
    blocks = {"changepoints": [], "step_CPs": [], "maxLL": [],"LLs":[],"step_Intens":[], "intensities": []}
    for i in np.arange(cells): #loop through cells 2 through to the end

        # create larger and larger blocks by backtracking from the last cell
        blocksum = np.cumsum(A[i::-1])[::-1]
        pointsum = np.cumsum(N[i::-1])[::-1]

        #here we create the log likelihoods from the ojbective function for each backtracked block combination
        #LLs = newobjective(pointsum, blocksum,c) + blocks["maxLL"]
        LLs = np.array([newobjective(pointsum[k], blocksum[k], c) for k in np.arange(i + 1)]) + np.hstack(
            (0, blocks["maxLL"]))

        #save the totals for maybe later viewing
        blocks["LLs"].append(LLs)

        #find where the best partition occurred
        best_partition = np.argmax(LLs)

        #add the resulting log likelihood sum to the array of best LL's
        blocks["maxLL"].append(np.amax(LLs))

        #calculate the intensity of this new block
        best_block_intensity=pointsum[best_partition]/blocksum[best_partition]

        # if the best partition includes all previous cells, assign the last changepoint to 0(first cell)
        if best_partition == 0:
            blocks["step_CPs"].append(np.array([0]))
            blocks["step_Intens"].append(np.array([best_block_intensity]))
        #if the best partition is somwhere else, find the optimal partitioning from previous calculations and add the new changepoint
        else:
            blocks["step_CPs"].append(np.hstack((blocks["step_CPs"][best_partition - 1], best_partition)))
            blocks["step_Intens"].append(np.hstack((blocks["step_Intens"][best_partition - 1],best_block_intensity)))
    blocks["changepoints"] = blocks["step_CPs"][-1]
    blocks["intensities"] = blocks["step_Intens"][-1]
    return blocks




#The PELT implementation is the same as above with PELT pruning added
# blocks["pruned"] gives you an array containing the pruned potential changepoints
def optintervalPELT(N,A,c):

    #initialize the block object
    cells = len(N)
    blocks = {"changepoints": [], "step_CPs": [], "maxLL": [], "LLs": [], "intensities": [], "step_Intens": [],"pruned": []}
    for i in np.arange(cells):

        # create larger and larger blocks by backtracking from the last cell
        blocksum = np.cumsum(A[i::-1])[::-1]
        pointsum = np.cumsum(N[i::-1])[::-1]
        # here we create the log likelihoods from the ojbective function for each backtracked block combination
        LLs=np.ones(i+1)*-np.inf
        for k in np.setdiff1d(range(i+1),blocks["pruned"]):
            LLs[k]=newobjective(pointsum[k], blocksum[k],c)
        LLs+= np.hstack((0, blocks["maxLL"]))

        # save the totals for maybe later viewing
        blocks["LLs"].append(LLs)

        # find where the best partition occurred
        best_partition = np.argmax(LLs)

        # update prune list
        blocks["pruned"]=np.where(LLs[best_partition] > c + LLs)[0]

        # add the resulting log likelihood sum to the array of best LL's
        blocks["maxLL"].append(np.amax(LLs))

        # calculate the intensity of this new block
        best_block_intensity = pointsum[best_partition] / blocksum[best_partition]

        # if the best partition includes all previous cells, assign the last changepoint to 0(first cell)
        if best_partition == 0:
            blocks["step_CPs"].append(np.array([0]))
            blocks["step_Intens"].append(np.array([best_block_intensity]))
        # if the best partition is somwhere else, find the optimal partitioning from previous calculations and add the new changepoint
        else:
            blocks["step_CPs"].append(np.hstack((blocks["step_CPs"][best_partition - 1], best_partition)))
            blocks["step_Intens"].append(np.hstack((blocks["step_Intens"][best_partition - 1], best_block_intensity)))

    blocks["changepoints"] = blocks["step_CPs"][-1]
    blocks["intensities"] = blocks["step_Intens"][-1]
    return blocks


def BinarySeg(N,A,c):
    cp=__BS_recurse(N,A,c)
    edge=np.hstack((cp,len(N)))
    return {"changepoints": cp,  "intensities": [(np.sum(N[edge[i]:edge[i+1]])/np.sum(A[edge[i]:edge[i+1]])) for i in range(len(cp))]}




def one_run(size,c):
    #create a random dataset of size "size"
    points = np.sort(np.random.random(size))
    [N, cell_edges] = voronoi1d(points)
    A = np.diff(cell_edges)

    #run the PELT algorithm
    start_time = time.clock()
    blocksPELT=optintervalPELT(N, A, c)
    time_PELT = time.clock() - start_time
    prop_PELT=len(blocksPELT["changepoints"])/size
    #print(time_PELT, "seconds PELT")

    #run the standard Optimal Partitioning algorithm
    start_time = time.clock()
    blocksOP=optinterval(N, A,c)
    time_OP = time.clock() - start_time
    prop_OP=len(blocksOP["changepoints"])/size
    #print(time_OP, "seconds OP")

    #run the binary segmentation algorithm
    start_time = time.clock()
    blocksBS=BinarySeg(N, A, c)
    time_BS = time.clock() - start_time
    prop_BS=len(blocksBS["changepoints"])/size
    #print(time_BS, "seconds BS")

    return {"OP":[time_OP,prop_OP],"PELT":[time_PELT,prop_PELT],"BS":[time_BS,prop_BS]}


#Returns an array containing the changepoints
def __BS_recurse(N,A,c):
    #create an array for the left side segmentation
    leftA = np.hstack((0,np.cumsum(A[:-1])))
    leftN = np.hstack((0,np.cumsum(N[:-1])))

    #create the array for the right side
    totalA=np.sum(A)
    rightA=totalA-leftA
    totalN = np.sum(N)
    rightN=totalN-leftN

    #select the best split point
    changepoint =np.argmax([newobjective(leftN,leftA,c)+newobjective(rightN,rightA,c)])


    if changepoint==0 :# no more splitting needed - stop recursion
        return [0]
    else: #continue down both branches
        leftCPs=__BS_recurse(N[:changepoint],A[:changepoint],c)
        rightCPs=__BS_recurse(N[changepoint:],A[changepoint:],c)+changepoint
        return np.unique(np.hstack((0,leftCPs,changepoint,rightCPs)))



#basic 1 dimensional Voronoi plot
def plot_voronoi1d(points):
    plt.close()
    A = np.hstack((points[0], (points[2:-1] + points[1:-2]) / 2, points[-1]))
    plt.vlines([points[0],points[-1]], ymin=-2, ymax=2,lw=3)
    plt.scatter(points[1:-1],np.zeros(len(points)-2))
    plt.vlines(A,ymin=-1,ymax=1,linestyle="dashed")
    plt.hlines(0,points[0],points[-1])
    plt.yticks([])
    plt.title("Voronoi Diagram")
    plt.xlim([points[0],points[-1]])

#basic 1 dimensional Delaunay Triangulation plot
def plot_delaunay1d(points):
    plt.close()
    plt.vlines([points[0],points[-1]], ymin=-2, ymax=2,lw=3)
    plt.scatter(points[1:-1],np.zeros(len(points)-2))
    plt.vlines(points[1:-1],ymin=-1,ymax=1,linestyle="dashed")
    plt.hlines(0,points[0],points[-1])
    plt.yticks([])
    plt.title("Delaunay Diagram")
    plt.xlim([points[0],points[-1]])


#plots the intensity view of the points(with or without changepoints)
def plot_intensity(N,cell_edges,changepoints=[],intensities=[]):
    plt.close()
    for i in range(len(N)):
        plt.plot([cell_edges[i], cell_edges[i + 1]], [1 / A[i], 1 / A[i]], c="black", lw=2)
        plt.ylim([0, max(1 / A) * 1.2])
    numblocks = len(changepoints)
    cp = np.hstack((changepoints, len(A)))
    ints = blocks["intensities"]
    for j in range(numblocks):
        plt.plot([cell_edges[cp[j]], cell_edges[cp[j + 1]]], [ints[j], ints[j]], lw=3, c="green")


#this function outputs the cell edges and sizes for Voronoi segmentation
def voronoi1d(points):
    return [np.ones(len(points)-2),np.hstack((points[0], (points[2:-1] + points[1:-2]) / 2, points[-1]))]






#this function outputs the cell edges and sizes for Delaunay segmentation
def delaunay1d(points):
    return [np.hstack((0.5,np.ones(len(points)),0.5)),np.diff(points)]








#this function creates some simulated data. the inputs are:
# rates: a vector of event rates wanted
# events: a vector with the number of events desired per event rate

# example: simplePPSim([3,5],[10,20]) will output a vector where the first 10 units of time with have an event rate of 3
#   and the next 20 units of time will have an event rate of 5
def simplePPSim(rates,times):
    data=[0]
    lasttime=0
    j=0
    for i in range(len(rates)):
        while data[-1]-lasttime<times[i]:
            data.append(data[j]+np.random.exponential(1/rates[i], 1))
            j+=1
        lasttime=data[-1]
    return np.array(data)


