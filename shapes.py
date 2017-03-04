import numpy as np
import matplotlib.pyplot as plt
import time


#dr scargles standard fitness function(poisson process data)
def newobjective(N, M, c):
    with np.errstate(all='ignore'):  # when N or M are zero, ignore the divide by zero warning and return 0
        return np.nan_to_num(N * np.log(np.array(N) / M) - c)

#bayesian block model object
class BBModel(object):

    #initializing the bayesian block model object with defaults
    def __init__(self, points=[], cell_type="voronoi", data_type="poisson",online=False):  # initialize the bayesian block object

        if data_type == "poisson":  # first check that the data type is valid
            self.data_type = "poisson"
        elif data_type == "binomial":
            self.data_type = "binomial"
        else:
            raise ValueError("data type must be \"poisson\" or \"binomial\"")  # throw error

        if len(points) == 0 and online==False:  # make sure points are provided
            raise ValueError("Points array is empty.")  # throw error
        elif online==False:

            if cell_type == "voronoi":  # make sure a cell type is provided
                self.N = np.ones(len(points) - 2)
                self.cell_edges = np.hstack((points[0], (points[1:-2] + points[2:-1]) / 2, points[-1]))
                self.M = np.diff(self.cell_edges)
                self.cell_type = "Voronoi"
            elif cell_type == "delaunay":
                self.N = np.hstack((0.5, np.ones(len(points) - 3), 0.5))
                self.cell_edges = points
                self.M = np.diff(self.cell_edges)
                self.cell_type = "Delaunay"
            else:
                raise ValueError("Cell type must be \"voronoi\" or \"delaunay\"")  # throw error if not

            # initialize variables
            self.points = points
            self.true_changepoints = []

        else: # online algorithm is used
            self.online=online
            self.tte = []

    # online algorithm (call this function with a new single point)
    def newpoint(self,newpoint,c):
        self.tte.append(newpoint)
        if len(self.tte) > 1:  # start algorithm when at least 2 points have been collected
            events = len(self.tte) - 1

            # create left and right blocks
            blockright = np.cumsum(np.diff(self.tte)[::-1])[::-1]
            pointright = np.arange(events, 0, -1)
            blockleft = np.hstack(([0], np.cumsum(np.diff(self.tte))[:-1]))
            pointleft = np.arange(events)

            LLs = newobjective(pointright, blockright, c) + newobjective(pointleft, blockleft, c)

            # print(LLs)
            # find where the best partition occurred
            best_partition = np.argmax(LLs)

            self.points = self.tte#update points vector for plotting purposes
            if best_partition == 0: #if untriggered, return 0
                return 0
            else:

                self.tte = [] #if trigger is found, reset data
                return (best_partition) #return trigger point
        return 0

    #run the optinterval algorithm and return the estimated changepoints -- PELT is implemented here
    def optinterval(self, c=0.5):  # optinterval function(should have same results as the matlab function)
        self.changepoints = []
        self.stepCP = []
        self.maxLL = []
        self.LLs = []
        self.stepIntensities = []
        self.intensities = []
        self.pruned = []

        # initialize the block object
        cells = len(self.N)
        blocks = {"changepoints": [], "step_CPs": [], "maxLL": [], "LLs": [], "intensities": [], "step_Intens": [],
                  "pruned": []}
        progress = False
        last = 0
        for i in np.arange(cells):
            start = time.clock()

            # create larger and larger blocks by backtracking from the last cell
            blocksum = np.cumsum(self.M[i::-1])[::-1]
            pointsum = np.cumsum(self.N[i::-1])[::-1]
            # here we create the log likelihoods from the ojbective function for each backtracked block combination
            LLs = np.ones(i + 1) * -np.inf
            for k in np.setdiff1d(range(i + 1), self.pruned):
                LLs[k] = newobjective(pointsum[k], blocksum[k], c)
            LLs += np.hstack((0, self.maxLL))

            # save the totals for maybe later viewing
            self.LLs.append(LLs)

            # find where the best partition occurred
            best_partition = np.argmax(LLs)

            # update prune list
            self.pruned = np.where(LLs[best_partition] > c + LLs)[0]

            # add the resulting log likelihood sum to the array of best LL's
            self.maxLL.append(np.amax(LLs))

            # calculate the intensity of this new block
            best_block_intensity = pointsum[best_partition] / blocksum[best_partition]

            # if the best partition includes all previous cells, assign the last changepoint to 0(first cell)
            if best_partition == 0:
                self.stepCP.append(np.array([0]))
                self.stepIntensities.append(np.array([best_block_intensity]))
            # if the best partition is somwhere else, find the optimal partitioning from previous calculations and add the new changepoint
            else:
                self.stepCP.append(np.hstack((self.stepCP[best_partition - 1], best_partition)))
                self.stepIntensities.append(
                    np.hstack((self.stepIntensities[best_partition - 1], best_block_intensity)))
            stop = time.clock()
            if progress or (stop - start) * cells > 120:
                progress = True
                if (i * 100 / cells) - last > 5:
                    print('{:2f}'.format((i * 100 / cells)) + " Percent Complete " + '{:6f}'.format(
                        (cells - i) * (stop - start)) + " Seconds Remaining")
                    last = i * 100 / cells

        self.changepoints = self.stepCP[-1]
        self.intensities = self.stepIntensities[-1]
        return self.changepoints

    #plot just the points in a Voronoi/Delaunay diagram
    def plot(self):
        plt.close()
        plt.vlines([self.points[0], self.points[-1]], ymin=-2, ymax=2, lw=3)
        plt.scatter(self.points[1:-1], np.zeros(len(self.points) - 2))
        plt.vlines(self.cell_edges, ymin=-1, ymax=1, linestyle="dashed")
        plt.hlines(0, self.points[0], self.points[-1])
        plt.yticks([])
        plt.title(self.cell_type + " Diagram")
        plt.xlim([self.points[0], self.points[-1]])

    #plot the estimated changepoints on top of a Voronoi/Delaunay diagram
    def plot_changepoints(self, include_true=False):
        self.plot()
        cp = plt.vlines(self.cell_edges[self.changepoints[1:]], ymin=-2, ymax=2, lw=5, color="red")
        if include_true and len(self.true_changepoints) != 0:
            true = plt.vlines(self.true_changepoints[1:], ymin=-2, ymax=2, lw=4, color="green")
            plt.legend([cp, true],
                       ["Estimated Changepoints", "True Changepoints"])

    #plot the data by intensity
    def plot_intensities(self, include_true=False):
        if self.online:
            self.N = np.hstack((0.5, np.ones(len(self.points) - 3), 0.5))
            self.cell_edges = self.points
            self.M = np.diff(self.cell_edges)

        plt.close()
        for i in range(len(self.M)):
            point, = plt.plot([self.cell_edges[i], self.cell_edges[i + 1]], [1 / self.M[i], 1 / self.M[i]],
                              color="black", lw=2)
            plt.ylim([0, np.percentile(1 / self.M, 90) * 1.2])
        plt.xlim([self.points[0], self.points[-1]])
        numblocks = len(self.changepoints)
        cp = np.hstack((self.changepoints, len(self.M)))
        ints = self.intensities
        for j in range(numblocks):
            block, = plt.plot([self.cell_edges[cp[j]], self.cell_edges[cp[j + 1]]], [ints[j], ints[j]], lw=4,
                              color="red")
        if include_true and len(self.true_changepoints) != 0:
            true = plt.vlines(self.true_changepoints[1:], ymin=0, ymax=np.max(1 / self.M), lw=4, color="green",linestyles="dashed")
            plt.legend([point, block, true],
                       [self.cell_type + " cell Intensity", "Block Intensity", "True Changepoints"])
        else:
            plt.legend([point, block], [self.cell_type + " cell Intensity", "Block Intensity"])
        plt.xlabel("Time")
        plt.ylabel("Intensity")
        plt.title("Block Segmentation")

    #if true changepoints are known, add them with this function to include them in the plots
    def add_true_changepoints(self, truecp=[]):
        self.true_changepoints = truecp


 # standard objective function



#####################simulated data generators##################################

#POISSON Process
#rates: vector of event rates for each respective block
#durations: time length of each block
def poissonSim(rates,durations):
    data=[0]
    lasttime=0
    j=0
    for i in range(len(rates)):
        while data[-1]-lasttime<durations[i]:
            data.append(data[j]+np.random.exponential(1/rates[i], 1))
            j+=1
        lasttime=data[-1]
    return np.array(data)


#BINOMIAL blocks
#rates: vector of proportion rates for each respective block
#durations: time length of each block
def binSim(rates,durations):
    data=np.array([])
    j=0
    for i in range(len(rates)):
        print(i,rates[i],durations[i])
        data=np.append(data,np.random.binomial(1,rates[i],durations[i]))
    return data


#Fast Rise, Early Decay Gamma Ray Burst
#rate: peak event rate
#rise: time to rise
#decay: time to decay
def FREDsim(rate,rise,decay):
    current = 0
    data = np.array([0])
    la = -np.log(0.003) / decay
    print(la)
    while current < rise+decay:
        e = np.random.exponential((1/rate), 1)
        current += e
        u = np.random.uniform()
        if u < ((la * np.power(np.e, -la * current)) / ( 1/rate)):
            data = np.append(data, current)
    return data


points=FREDsim(1000,0,10)
############ script ########################


############# REGULAR ALGORITHM ###################
#first create the Bayesian Block Model object using the data points as a parameter

points=poissonSim([100,200,100],[1,1,1])
blocks=BBModel(points,cell_type="voronoi")
blocks.plot()
blocks.optinterval(c=5)
blocks.plot_changepoints()
blocks.plot_intensities()
blocks.add_true_changepoints([0,1,2])
blocks.plot_changepoints(include_true=True)
blocks.plot_intensities(include_true=True)


########### ONLINE ALGORITHM ######################

points = poissonSim([50,100], [2,2])
mod=BBModel(online=True)
for i in range(len(points)):
    result = mod.newpoint(points[i], 5)
    if result:
        break

mod.add_true_changepoints([2])
mod.plot_changepoints()
print(result)