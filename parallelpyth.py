from mpi4py import MPI
from os import write
import DLS as cdls
import numpy as np
import sys

def calculate_mandelBrot(start, end):
    
    for i in range(start, end):
        
        # t1 = time.time()
        x = i // h_pixel
        y = i % h_pixel
        a = (x - (w_pixel/2)) / (w_pixel/4)
        b = (y - (h_pixel/2)) / (h_pixel/4)
   
        # c = Complex(a, b)
        c_r = a
        c_i = b
        # z = Complex(0, 0)
        z_r = 0
        z_i = 0
        it = 0
        for it in range(0,iterations):

            t = (z_r * z_r) - (z_i * z_i)
            z_i = 2.0 * z_r * z_i
            z_r = t

            t = (z_r * z_r) - (z_i * z_i)
            z_i = 2.0 * z_r * z_i
            z_r = t

            z_r = c_r + z_r
            z_i = z_i + c_i
            # z = add(z, c)
            
            mag =  (z_r*z_r) + (z_i*z_i)

            if mag > 4.0:
                pixel_arr[i] = it
                break

        if pixel_arr[i] == -1:
            pixel_arr[i] = iterations
        # t2 = time.time()

        # print(t2-t1," - ", i)

def writeFile(h, w, pixel_arr):
    f = open("py_parallel_mandel.csv", "w")
    for i in range(0, h*w):
        x = i // h
        y = i % h
        strval = str(y)+","+str(x)+","+str(pixel_arr[i])+",\n"
        f.write(strval)
    f.close()



if __name__ == "__main__":

    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
    # iterations = 10000
    h_pixel = 1024
    w_pixel = 1024
    total_pixel = h_pixel * w_pixel
    master = 0
    requestWhen = -1
    breakAfter = -1
    minChunk = 20
    h_overhead = 0.013716 
    mu = 0.1
    sigma = 0.0005
    alpha = 0.0605
    nKNL = 0 # number of KNL cores
    KNL_speed = 0.414 
    Xeon_speed = 1.85
    B = 2
    X = 4
    SWR = .7
    worktime = 0
    itersDone = 0
    # print ("Number of args: ", len(sys.argv))
    arr_Args = sys.argv

    method = int(arr_Args[1])
    iterations = int(arr_Args[2])
    pixel_arr = np.full(total_pixel, -1)
    results_arr = np.full(total_pixel, -1)
    # pixel_arr = [-1] * total_pixel
    # results_arr = [-1] * total_pixel
    t1 = MPI.Wtime()
    info = cdls.infoDLS(comm, size, rank, 0, total_pixel, master, method, requestWhen, breakAfter, minChunk, h_overhead, sigma, mu, alpha, nKNL, Xeon_speed, KNL_speed, X, B, SWR)

    cdls.StartLoop(info)
    start = 0
    chunks = 0
    end = 0
    while(not cdls.DLS_Terminated(info)):

        start, chunks = cdls.DLS_StartChunk(info)

        if start < (start+chunks):
            if(start+chunks - total_pixel == 1):
                calculate_mandelBrot(start, start+chunks-1)
            else:
                calculate_mandelBrot(start, start+chunks)
            
        cdls.DLS_EndChunk(info)    
    
    itersDone, worktime = cdls.DLS_EndLoop(info)
    t2 = MPI.Wtime()
   
    t2 = t2 - t1
    updatetime = 0
    comm.Reduce(pixel_arr, results_arr, MPI.MAX, 0)
    updatetime = comm.reduce(t2, MPI.MAX, 0)
    
    if rank == 0:
        writeFile(h_pixel,w_pixel, results_arr)
        print("Total time: ",updatetime)
   
