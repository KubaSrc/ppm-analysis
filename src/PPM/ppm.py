from os import link
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
from pyparsing import alphanums
from scipy.linalg import norm

class mechanism:
    """
    Docstring For PPM
    """
    def __init__(self,Rr=1,**kwargs):
        # Initialize drive link to zero state
        self.theta = 0
        self.phi = 0
        self.theta_max = np.deg2rad(25)
        # Link and node geometry
        self.link_r = 0.75
        self.node_r = 1

        # Create node pair assigments
        self.links = ["A0","A1","B0","B1","B2","B3","B4","B5","C0","C1","C2","D0","D1"]
        self.node_pairs = [(0,1),(1,2),(2,3),(2,4),(2,5),(3,6),(4,6),(5,6),(0,3),(0,4),(0,5),(3,4),(4,5)]

        # Used for initializing arrays
        self.XYZ_rand = False

        # Define transofrmation matricies
        self.T_xyz = np.vstack(([1.0,0.0,0.0,0.0],[0.0,1.0,0.0,0.0],[0.0,0.0,1.0,0.0],[0.0,0.0,0.0,1.0]))
        self.R_x = np.vstack(([1.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,1.0]))
        self.R_y = np.vstack(([0.0,0.0,0.0,0.0],[0.0,1.0,0.0,0.0],[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,1.0]))
        self.R_z = np.vstack(([0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0],[0.0,0.0,1.0,0.0],[0.0,0.0,0.0,1.0]))

        # Update geometry
        self.update_links(Rr,**kwargs)
        self.update_state()
        self.find_surface()

    def update_links(self,Rr=1,**kwargs):
        # Link length parametrization
        if ("A" in kwargs and "B" in kwargs and "C" in kwargs and "D" in kwargs):
            self.A0,self.B0,self.C0,self.D0 = kwargs["A"],kwargs["B"],kwargs["C"],kwargs["D"]
            theta_CA = self.__law_of_cos(self.B0,2*self.A0,self.C0)
            D = self.C0*np.sin(theta_CA)
            H = self.B0*np.sin(np.arccos(D/self.B0))
            Lc = 2*self.A0+2*H
            G = 360-2*self.__law_of_cos(self.D0,D,D,mode="deg")
            self.L_char = Lc
            self.Dt = D
            self.Ht = H
            self.Gt = G
        # Triangular Bipyramid and Lc
        elif ("Lc" in kwargs and "H" in kwargs and "D" in kwargs and "G" in kwargs):
            Lc,H,D,G = kwargs["Lc"],kwargs["H"],kwargs["D"],np.deg2rad(kwargs["G"])
            self.A0 = (Lc-2*H)/2
            self.B0 = np.sqrt(H**2+D**2)
            self.C0 = np.sqrt((2*self.A0+H)**2+D**2)
            self.D0 = np.sqrt((D+D*np.cos(G/2))**2+(D*np.sin(G/2))**2)
            # Save TriPyr parametrization
            self.L_char = Lc
            self.Ht = H
            self.Dt = D
            self.Gt = np.rad2deg(G)
        else: raise Exception("Incorrect PPM parameters given!")
        # Create link arrays
        self.A = np.array([self.A0]*2)
        self.B = np.array([self.B0]*6)
        self.C = np.array([self.C0]*3)
        self.D = np.array([self.D0]*2)
        self.link_lengths = np.concatenate((self.A,self.B,self.C,self.D))

    def print_params(self):
        print("A: {:.4f} B: {:.4f} C: {:.4f} D: {:.4f} Lc: {:.4f} Ht: {:.4f} Dt: {:.4f} Gt: {:.4f}"\
            .format(self.A0,self.B0,self.C0,self.D0,self.L_char,self.Ht,self.Dt,self.Gt))

    # Creates a new permutation of a noisy PPM
    def add_noise(self,sigma,ppm_seed=None):
        if ppm_seed != None: np.random.seed(ppm_seed)
        self.A = self.A0 + np.random.normal(scale=sigma,size=2)
        self.B = self.B0 + np.random.normal(scale=sigma,size=6)
        self.C = self.C0 + np.random.normal(scale=sigma,size=3)
        self.D = self.D0 + np.random.normal(scale=sigma,size=2)
        self.E_in = np.hstack((self.A-self.A0,self.B-self.B0,self.C-self.C0,self.D-self.D0))
        self.RMSE_in = np.sqrt(np.sum(self.E_in**2)/np.size(self.E_in))
        self.link_lengths = np.concatenate((self.A,self.B,self.C,self.D))

    # Updates angles for PPM
    def update_angles(self,theta=0,phi=0,debug=False):
        self.theta = theta
        self.phi = phi
        self.update_state(debug)

    # Returns angle between L2 and L3
    def __law_of_cos(self,L1,L2,L3,mode="rad"):
        ratio  = (L2**2+L3**2-L1**2)/(2*L2*L3)
        if ratio > 1 or ratio < -1:
            raise  Exception("Triangle cannot be solved in arccos domain!",L1,L2,L3,ratio)
        T = np.arccos(ratio)
        if mode == "deg": T = np.rad2deg(T)
        return T

    # Returns length across from theta
    def __law_of_cos_theta(self,L1,L2,theta,mode="rad"):
        try:
            return np.sqrt(L1**2+L2**2-2*L1*L2*np.cos(theta))
        except:
            print("Arc Cos Error:",L1,L2,theta)

    # Returns angle across from L2
    def __law_of_sin(self,L1,theta1,L2):
        return np.arcsin(L2*(np.sin(theta1)/L1))
                           
    # Find the intersection of three spheres                 
    # P1,P2,P3 are the centers, r1,r2,r3 are the radii       
    # Implementaton based on Wikipedia Trilateration article.
    # https://stackoverflow.com/questions/1406375/finding-intersection-points-between-3-spheres                              
    def __trilaterate(self,P1,P2,P3,r1,r2,r3):                      
        temp1 = P2-P1                                        
        e_x = temp1/norm(temp1)                              
        temp2 = P3-P1                                        
        i = np.dot(e_x,temp2)                                   
        temp3 = temp2 - i*e_x                                
        e_y = temp3/norm(temp3)                              
        e_z = np.cross(e_x,e_y)                                 
        d = norm(P2-P1)                                      
        j = np.dot(e_y,temp2)                                   
        x = (r1*r1 - r2*r2 + d*d) / (2*d)                    
        y = (r1*r1 - r3*r3 -2*i*x + i*i + j*j) / (2*j)       
        temp4 = r1*r1 - x*x - y*y                            
        if temp4<0:                                          
            raise Exception("The three spheres do not intersect!")
        z = np.sqrt(temp4)                                      
        p_12_a = P1 + x*e_x + y*e_y + z*e_z                  
        p_12_b = P1 + x*e_x + y*e_y - z*e_z # We can ignore this solution for the PPM                  
        return p_12_b     


    # Finds distance between nodes in the structure
    def __dist(self,V1,V2):
        return np.linalg.norm(V1-V2)

    ####################
    # Forward Kinematics
    ####################

    def __validate_geo(self):
        for i in range(len(self.node_pairs)):
            print(self.links[i],round(self.link_lengths[i]-self.__dist(self.N[self.node_pairs[i][0]],self.N[self.node_pairs[i][1]]),10))
        pass
    
    def __t_xyz(self,dX,dY,dZ):
        self.T_xyz[0,3] = dX
        self.T_xyz[1,3] = dY
        self.T_xyz[2,3] = dZ
        return np.copy(self.T_xyz)
    
    def __r_x(self,tx):
        self.R_x[1,1] = np.cos(tx)
        self.R_x[1,2] = -np.sin(tx)
        self.R_x[2,1] = np.sin(tx)
        self.R_x[2,2] = np.cos(tx)
        return np.copy(self.R_x)

    def __r_y(self,ty):
        self.R_y[0,0] = np.cos(ty)
        self.R_y[0,2] = np.sin(ty)
        self.R_y[2,0] = -np.sin(ty)
        self.R_y[2,2] = np.cos(ty)
        return np.copy(self.R_y)

    def __r_z(self,tz):
        self.R_z[0,0] = np.cos(tz)
        self.R_z[0,1] = -np.sin(tz)
        self.R_z[1,0] = np.sin(tz)
        self.R_z[1,1] = np.cos(tz)
        return np.copy(self.R_z)

    # Finds the position of each
    def forward_kinematics(self,theta,phi,debug=False):
        # N0 (origin)
        N0 = np.array([0,0,0])
        # N1
        T_N1 = self.__t_xyz(self.A[0],0,0)
        N1 = T_N1[0:3,3]
        # N2 - Find using virtual link L from N0 -> N2 to simplify kinematics
        L = self.__law_of_cos_theta(self.A[0],self.A[1],np.pi-theta)
        theta_L = self.__law_of_sin(L,np.pi-theta,self.A[1])
        T_N2 = self.__r_x(phi)@self.__r_z(theta_L)@self.__t_xyz(L,0,0)@self.__r_x(-phi)
        N2 = T_N2[0:3,3]
        # Find thetas for B links
        theta_B0 = np.pi-self.__law_of_cos(self.C[0],L,self.B[0])
        theta_B1 = np.pi-self.__law_of_cos(self.C[1],L,self.B[1])
        theta_B2 = np.pi-self.__law_of_cos(self.C[2],L,self.B[2])
        # Find projection angle for D link
        C0_L = self.__law_of_cos(self.B[0],L,self.C[0])
        C1_L = self.__law_of_cos(self.B[1],L,self.C[1])
        C2_L = self.__law_of_cos(self.B[2],L,self.C[2])
        skew_D0 = self.C[0]*np.cos(C0_L)-self.C[1]*np.cos(C1_L)
        delta_D0 = np.tan(skew_D0/self.D[0])
        skew_D1 = self.C[2]*np.cos(C2_L)-self.C[1]*np.cos(C1_L)
        delta_D1 = np.tan(skew_D1/self.D[1])
        # Find phis for B links
        phi_B1 = 0
        phi_B0 = self.__law_of_cos(self.D[0]*np.cos(delta_D0),self.B[0]*np.sin(theta_B0),self.B[1]*np.sin(theta_B1))
        phi_B2 = -self.__law_of_cos(self.D[1]*np.cos(delta_D1),self.B[1]*np.sin(theta_B1),self.B[2]*np.sin(theta_B2))
        # N3
        T_N3 = T_N2@self.__r_x(phi_B0)@self.__r_z(theta_B0)@self.__t_xyz(self.B[0],0,0)
        N3 = T_N3[0:3,3]
        # N4
        T_N4 = T_N2@self.__r_x(phi_B1)@self.__r_z(theta_B1)@self.__t_xyz(self.B[1],0,0)
        N4 = T_N4[0:3,3]
        # N5
        T_N5 = T_N2@self.__r_x(phi_B2)@self.__r_z(theta_B2)@self.__t_xyz(self.B[2],0,0)
        N5 = T_N5[0:3,3]
        # N6
        N6 = self.__trilaterate(N3,N4,N5,self.B[3],self.B[4],self.B[5])
        # Group nodes into a full list and return
        N = [N0,N1,N2,N3,N4,N5,N6]
        # Verify distances between node links
        if debug:
            self.__validate_geo()
        return N

    # Updates all nodes within mechanism
    def update_state(self,debug=False):
        self.N = self.forward_kinematics(self.theta,self.phi)
        if debug: self.__validate_geo()


    # Finds plane of best fit for a set of 3D points
    def plane_fit(self,P):
        m = P.sum(axis=0) / P.shape[0]
        # run SVD
        u, s, vh = np.linalg.svd(P - m)
        # unitary normal vector
        n = vh[2, :]
        return(n,m)

    # Returns distances from plane of best fit
    def find_error(self,X,Y,Z):
        n,m = self.plane_fit(np.hstack((X,Y,Z)))
        # Coefficients of plane
        a,b,c = n[0],n[1],n[2]
        d = a*m[0]+b*m[1]+c*m[2]
        # Calculate distance of points
        D = (a*X+b*Y+c*Z-d)/np.sqrt(a**2+b**2+c**2) # Residuals
        RMSE_out = np.sqrt(np.sum(D**2)/(np.size(D))) # RMSE
        return D,RMSE_out

    # Sweeps over input vectors to find output surface
    def find_surface(self,n=20,m=20,flip_yz=False,debug=False):
        self.X = np.zeros([n,m])
        self.Y = np.zeros([n,m])
        self.Z = np.zeros([n,m])
        for i,theta in enumerate(np.linspace(0,self.theta_max,n)):
            for j,phi in enumerate(np.linspace(0,2*np.pi,m)):
                N = self.forward_kinematics(theta,phi)
                self.X[i][j],self.Y[i][j],self.Z[i][j] = N[6][0],N[6][1],N[6][2]
        D,RMSE = self.find_error(np.reshape(self.Z,(n*m,1)),np.reshape(self.Y,(n*m,1)),np.reshape(self.X,(n*m,1)))
        return (RMSE,100*np.max(self.Y)/self.L_char)

    # Sweep of randomized configurations
    def find_random_surface(self,n,flip_yz=False,debug=False):
        if not self.XYZ_rand:
            self.X_rand = np.zeros((n,1))
            self.Y_rand = np.zeros((n,1))
            self.Z_rand = np.zeros((n,1))
            self.XYZ_rand = True # For efficient memory usage during big sims
        for i in range(n):
            theta = np.random.uniform(0,self.theta_max)
            phi = np.random.uniform(0,2*np.pi)
            N = self.forward_kinematics(theta,phi,debug)
            if flip_yz:
                self.X_rand[i,0],self.Y_rand[i,0],self.Z_rand[i,0] = N[6][0],N[6][2],N[6][1]
            else:
                self.X_rand[i,0],self.Y_rand[i,0],self.Z_rand[i,0] = N[6][0],N[6][1],N[6][2]
        D,RMSE = self.find_error(self.Z_rand,self.Y_rand,self.X_rand)
        return (RMSE,100*np.max(self.Y)/self.L_char)

    ####################
    # Visualization
    ####################

    def plot_node(self,N,i):
        N = np.flip(N)
        x,y,z = N[0],N[1],N[2]
        c = [0.5, 0.5, 0.5]
        if i == 6:
            c = [0.6350, 0.0780, 0.1840]
        self.ax.scatter(x,y,z,s=[40],color=c,zorder=3)

    def plot_link(self,N1,N2,a):
        N1 = np.flip(N1)
        N2 = np.flip(N2)
        plt.plot(*zip(N1,N2),linewidth=3,color="black",alpha=a,zorder=1)

    def draw(self,animation=False):
        # Set up for plot
        if not animation:
            self.fig = plt.figure(1)
            self.ax = self.fig.add_subplot(projection='3d',computed_zorder=False)
        # Set up axis for animation
        else:
            self.ax.cla()

        # Plot output surface
        alpha = 0.1
        self.ax.plot_surface(self.Z, self.Y, self.X, color=[0.3010, 0.7450, 0.9330, alpha],edgecolor=[0, 0.4470, 0.7410, 0.5],zorder=2)    
        
        # Draw nodes
        for i,N in enumerate(self.N):
            self.plot_node(N,i)
        # Draw links
        for pair in self.node_pairs:
            alpha = 1
            if pair == (0,1): alpha = 0.4
            self.plot_link(self.N[pair[0]],self.N[pair[1]],alpha)           
        
        # Configure plot
        self.ax.set_xlabel("Z")
        self.ax.set_ylabel("Y")
        self.ax.set_zlabel("X")
        self.ax.set_zlim(0,85)
        self.ax.set_ylim(-50,50)
        self.ax.set_xlim(-50,50)
        self.ax.view_init(0,0)

        # Show plot is not in animation
        if not animation:
            plt.axis('equal') 
            plt.show()

    def plot_surface(self):
        self.find_surface()
        self.fig = plt.figure(3)
        self.ax = self.fig.add_subplot(projection='3d',computed_zorder=False)
        self.ax.plot_surface(self.Z, self.Y, self.X,zorder=2)   
        plt.show()

    # Used to update frames ion animation
    def update_frame(self,sweep):
        self.theta = sweep[0]
        self.phi = sweep[1]
        self.update_state()
        return self.draw(animation=True)

    # Animates mechanism moving across theta range
    def animate(self,sweep,save=False,path = "./videos/ppm_animation"):
        # Plotting Variables
        self.fig = plt.figure(2)
        self.ax = self.fig.add_subplot(projection='3d',computed_zorder=False)
        self.ppm_animation = animation.FuncAnimation(self.fig, self.update_frame, frames=sweep, interval=100, repeat=True)
        if save:
            print("SAVING VIDEO")
            self.ppm_animation.save(path+".gif", writer='imagemagick', fps=60)
        
        plt.show()