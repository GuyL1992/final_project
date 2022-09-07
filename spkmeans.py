import pandas as pd
import numpy as np
import sys
import spkmeans
from enum import Enum

class Goals(Enum):
    WAM = "wam"
    DDG = "ddg"
    LNORM = "lnorm"
    JACOBI = "jacobi"
    SPK = "spk"

def initial_points(file_name):

    """ call this method to read the input file and translate it to nested list represents the observations 
    Args:
        file_name (string): input file name (path)

    Returns:
        list[list[float]], int, int: 3 returned values, one represents the observed vectors (nested list) and two int values
        to represent the dimentions (n & d)
    """

    with open(file_name,'r') as file:

        def parse_row(row):
            return list(map(lambda point: float(point),row.split(',')))
        observations = list(map(lambda row: parse_row(row), file))

    return observations, len(observations), len(observations[0])

### K_means_pp process and the required help functions ###

def calc_distance(vec1, vec2): # calculate distance between two vectors 
    return sum((vec1 - vec2)**2)
          
def kmeans_pp(k,data_arr,n,d):
    """_summary_

    Args:
        k (int): number of clusters required 
        data_arr (list[list[double]]): input observations 
        n (int): first dimention
        d (int): second dimention

    Returns:
        list[int], dataFrame: list of initiate k index, and the actual init k vectors as df
    """
    np.random.seed(0)
    first_centroid = np.random.choice(n)
    centroids_index = []
    p_list =  np.zeros((n), dtype=np.float64) 
    d_list =  np.full(n, np.inf) 
    centroids_index.append(first_centroid)
    d_sum = 0
    j = 1
    while j < k:
        new_centroid = data_arr.loc[centroids_index[j-1]]
        for i in range(n):
            curr_distance = calc_distance(data_arr.loc[i], new_centroid)
            if curr_distance < d_list[i]:
                d_list[i] = curr_distance
        d_sum = d_list.sum()

        for i in range(n):
            p_list[i] = d_list[i] / d_sum
        next_centroid = np.random.choice(n, size=None, p=p_list)
        centroids_index.append(next_centroid)
        j += 1
    
    return centroids_index, data_arr.loc[centroids_index]


def recieve_input(): 
    """ recieve and validate the user input. assigned it to the appropriate variables

    Returns:
        int, str, str: k desired by the user - might not be provided, the desired flow, the input file path
    """
    arguments_size = len(sys.argv)
    if arguments_size == 4:
        k_float = float(sys.argv[1])
        k = int(k_float)   
        if k_float != k:
            print("Invalid Input")
            exit()   
        flow = sys.argv[2]
        input = sys.argv[3]      

    elif arguments_size == 3:
        k=0
        flow = sys.argv[1]
        input = sys.argv[2]

    else :
        print("Invalid Input!")
        exit()
    
    # TO DO - add a call to validate input function with all the validation required 

    return k, flow, input

def validate_input(k): 
    if(k < 1 or k >= n):
        print("Invalid Input")
        exit()
    return

if __name__ == "__main__":
    MAX_ITER = 300
    k,goal,input_path = recieve_input()
    observations, n, d = initial_points(input_path)

    match(goal):
        case (Goals.SPK.value):
            new_data = spkmeans.get_new_data(n,d,k,observations)
            k = len(new_data[0])
            df = pd.DataFrame(new_data)
            init_index, init_centroids = kmeans_pp(k,df,n,k)
            any(map(lambda i: print(init_index[i], end=",") if i < len(init_index) - 1 else print(init_index[i]),range(len(init_index))))
            init_centroids = init_centroids.values.tolist()
            spkmeans.kmeans_pp(n,k,k,MAX_ITER,new_data,init_centroids)

        case(Goals.WAM.value):
            spkmeans.wam(n,d,observations)

        case(Goals.DDG.value):
            spkmeans.ddg(n,d,observations)

        case(Goals.LNORM.value):
            spkmeans.lnorm(n,d,observations)

        case(Goals.JACOBI.value):
            spkmeans.jacobi(n,d,observations)






















    # try:
    #     vectors = initial_Points(input)
    # except:
    #     print("An Error Has Occured")
    #     exit()

    # n = len(vectors) # n = num of rows
    # d = len(vectors[0])# d=  num of columns
    # k,flow, input = recieve_input() 
    # validate_input(k,n)
    # max_iter=300

    # try:
    #     if flow=="spk":
             #spkmeans.get_new_data(n,d,k,vectors)
            # data_Points_df = pd.DataFrame(newDataPoints)
            # k=len(newDataPoints[0])
            # centroids_indexes,centroids_points  = kmeans_pp(k, data_Points_df,n,d)
            # centroids_list =centroids_points.values.tolist()
            # for i in range(k - 1):
            #     print(centroids_indexes[i], end=',')
            # print(centroids_indexes[k - 1])
            # k_centroids = spkmeans.kmeans_pp(n, k, k, max_iter, newDataPoints, centroids_list)
        # elif flow=="wam":
        #     spkmeans.wam(n,d,vectors)
        # elif flow=="ddg":
        #     spkmeans.ddg(n,d,vectors)
        # elif flow == "lnorm":
        #     spkmeans.lnorm(n,d,vectors)
        # elif flow=="jacobi":
        #     spkmeans.jacobi(n,vectors)
    # except:
    #     print("An Error Has Occurred")
    #     exit()   
