from src.PPM import ppm
from matplotlib import pyplot as plt
import numpy as np
import cProfile
import pstats
import time
from multiprocessing import Process
import datetime

# Dispatched to simulate part of the design space
def sim_worker(data_path,Lc,n_seeds,m_samples,sigma,X):
    print(data_path)
    # Define base PPM
    PPM = ppm.mechanism(Lc=1,H=.2,D=.4,G=60)
    with open(data_path, "w") as f:
        f.write("G,H,HD,E_mean\n")

    # Loop over all specified designs
    for i in range(np.shape(X)[0]):
        h,hd,g = X[i,0],X[i,1],X[i,2]
        PPM.update_links(Lc=Lc,H=h,D=h/hd,G=g)
        E_tot = 0
        # Test PPM (could possibly be ill-defined for given theta-max)
        try:
            # Run error sim
            for i in range(0,n_seeds):
                PPM.add_noise(sigma,i)
                e_out,w = PPM.find_random_surface(m_samples)
                E_tot += e_out/sigma
            # Find mean error
            E_mean = E_tot/n_seeds
        except:
            E_mean = np.NAN
        # Write to path
        with open(data_path, "a") as f:
            np.savetxt(f, np.array([[g,h,hd,E_mean]]),delimiter=',')

# Splits design space across multiple threads and simulates error
def sim_main(s_range,n_sweep,num_workers,Lc,n_seeds,m_samples,sigma):
    # Define our simulation bounds
    H = np.linspace(s_range[0][0],s_range[0][1],n_sweep[0])
    HD = np.linspace(s_range[1][0],s_range[1][1],n_sweep[1])
    G = np.linspace(s_range[2][0],s_range[2][1],n_sweep[2])
    print("H_size: {}, HD_size: {}, G_size: {}".format(np.size(H),np.size(HD),np.size(G)))

    # Combine into a single vector
    X = np.zeros((np.size(H)*np.size(HD)*np.size(G),3))
    idx = 0
    for h in H:
        for hd in HD:
            for g in G:
                X[idx,0],X[idx,1],X[idx,2] = h,hd,g
                idx += 1
    
    print("Sim size:",np.shape(X)[0])
    csv_list = ["./data/sim_{}".format(i) for i in range(num_workers)]
    params = [(csv_list[i],Lc,n_seeds,m_samples,sigma,np.array_split(X,num_workers)[i]) for i in range(num_workers)]
    p_list = []

    # Divide sim over workers
    for worker in range(num_workers):
        p = Process(target = sim_worker,args=params[worker])
        p_list.append(p)
    
    # Start each process
    for p in p_list: p.start()
    # Wait to finish
    for p in p_list: p.join()

    # Combine results into one big csv
    results = "./data/error_sim_"+datetime.datetime.now().strftime("%Y_%m_%d-%H_%M_%S")
    with open(results, "w") as f:
        f.write("G,H,HD,E_mean\n")
    for f_csv in csv_list:
        worker_file = np.loadtxt(f_csv,delimiter=',',skiprows=1)
        with open(results, "a") as f:
            np.savetxt(f,worker_file,delimiter=',')
    
if __name__ == '__main__':

    t_start = time.time()
    print("Starting sim: "+str(datetime.datetime.now()))
    # Define design space sweep
    sweep_ranges = ((.1,.5),(.1,3.5),(10,350)) # (H, HD, G)
    delta_sweep = (0.01,0.025,2) # For np.arange
    n_sweep = (100,100,100) # For np.linspace
    num_workers = 6
    # Simulation constants
    n_seeds = 100
    m_samples = 50
    Lc = 1
    sigma = Lc*(.05/100)
    # Run simulation
    sim_main(sweep_ranges,n_sweep,num_workers,Lc,n_seeds,m_samples,sigma)
    print("Run time: {:.3f}".format(time.time()-t_start))
