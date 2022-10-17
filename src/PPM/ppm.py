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
        self.theta_max = np.deg2rad(10)
        # Link and node geometry
        self.link_r = 0.75
        self.node_r = 1

        # Create node pair assigments
        self.links = ["A0","A1","B0","B1","B2","B3","B4","B5","C0","C1","C2","D0","D1"]
        self.node_pairs = [(0,1),(1,2),(2,3),(2,4),(2,5),(3,6),(4,6),(5,6),(0,3),(0,4),(0,5),(3,5),(4,5)]

        # Update geometry
        self.update_links(Rr,**kwargs)
        self.update_state()
        self.find_surface()
        self.find_Lchar()

    def update_links(self,Rr,**kwargs):
        # Link length parametrization
        if ("A" in kwargs and "B" in kwargs and "C" in kwargs and "D" in kwargs):
            self.A = np.array([kwargs["A"]]*2)
            self.B = np.array([kwargs["B"]]*6)
            self.C = np.array([kwargs["C"]]*3)
            self.D = np.array([kwargs["D"]]*2)
        if ("A_links" in kwargs and "B_links" in kwargs and "C_links" in kwargs and "D_links" in kwargs):
            self.A = kwargs["A_links"]
            self.B = kwargs["B_links"]
            self.C = kwargs["C_links"]
            self.D = kwargs["D_links"]
        self.link_lengths = np.concatenate((self.A,self.B,self.C,self.D))
        # Rhombic pyramid parametrization

    def find_Lchar(self):
        self.Lchar = self.__dist(self.N[0],self.N[6])

    # Updates angles for PPM
    def update_angles(self,theta=0,phi=0,debug=False):
        self.theta = theta
        self.phi = phi
        self.update_state(debug)

    # Returns angle between L2 and L3
    def __law_of_cos(self,L1,L2,L3,mode="rad"):
        ratio  = round((L2**2+L3**2-L1**2)/(2*L2*L3),10)
        try:
            T = np.arccos(ratio)
            if mode == "deg": T = np.rad2deg(T)
            return T
        except:
            print("Arc Cos Error:",L1,L2,L3,ratio)

    # Returns length across from theta
    def __law_of_cos_theta(self,L1,L2,theta,mode="rad"):
        try:
            return np.sqrt(L1**2+L2**2-2*L1*L2*np.cos(theta))
        except:
            print("Arc Cos Error:",L1,L2,theta)

    # Returns angle across from L2
    def __law_of_sin(self,L1,theta1,L2):
        return np.arcsin(L2*(np.sin(theta1)/L1))

    def __is_valid(self):
        pass
                           
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
        return p_12_a     


    # Finds distance between nodes in the structure
    def __dist(self,V1,V2):
        return np.linalg.norm(V1-V2)

    ####################
    # Forward Kinematics
    ####################

    def __validate_geo(self):
        for i in range(len(self.node_pairs)):
            print(self.links[i],round(self.link_lengths[i]-self.__dist(self.N[self.node_pairs[i][0]],self.N[self.node_pairs[i][1]]),5))
        pass
    
    # Finds the position of each
    def forward_kinematics(self,theta,phi):
        # Define translation and rotation matricies
        t_xyz = lambda dX,dY,dZ: np.vstack(([1,0,0,dX],[0,1,0,dY],[0,0,1,dZ],[0,0,0,1]))
        r_x = lambda tx: np.vstack(([1,0,0,0],[0,np.cos(tx),-np.sin(tx),0],[0,np.sin(tx),np.cos(tx),0],[0,0,0,1]))
        r_y = lambda ty: np.vstack(([np.cos(ty),0,np.sin(ty),0],[0,1,0,0],[-np.sin(ty),0,np.cos(ty),0],[0,0,0,1]))
        r_z = lambda tz: np.vstack(([np.cos(tz),-np.sin(tz),0,0],[np.sin(tz),np.cos(tz),0,0],[0,0,1,0],[0,0,0,1]))
        # N0 (origin)
        N0 = np.array([0,0,0])
        # N1
        T_N1 = t_xyz(self.A[0],0,0)
        N1 = T_N1[0:3,3]
        # N2 - Find using virtual link L from N0 -> N2 to simplify kinematics
        L = self.__law_of_cos_theta(self.A[0],self.A[1],np.pi-theta)
        theta_L = self.__law_of_sin(L,np.pi-theta,self.A[1])
        T_N2 = r_x(phi)@r_z(theta_L)@t_xyz(L,0,0)@r_x(-phi)
        N2 = T_N2[0:3,3]
        # Find thetas for B links
        theta_B0 = np.pi-self.__law_of_cos(self.C[0],L,self.B[0])
        theta_B1 = np.pi-self.__law_of_cos(self.C[1],L,self.B[1])
        theta_B2 = np.pi-self.__law_of_cos(self.C[2],L,self.B[2])
        # Find phis for B links
        phi_B1 = 0
        phi_B0 = self.__law_of_cos(self.D[0],self.B[0]*np.sin(theta_B0),self.B[1]*np.sin(theta_B1))
        phi_B2 = -self.__law_of_cos(self.D[1],self.B[1]*np.sin(theta_B1),self.B[2]*np.sin(theta_B2))
        # N3
        T_N3 = T_N2@r_x(phi_B0)@r_z(theta_B0)@t_xyz(self.B[0],0,0)
        N3 = T_N3[0:3,3]
        # N4
        T_N4 = T_N2@r_x(phi_B2)@r_z(theta_B2)@t_xyz(self.B[2],0,0)
        N4 = T_N4[0:3,3]
        # N5
        T_N5 = T_N2@r_x(phi_B1)@r_z(theta_B1)@t_xyz(self.B[1],0,0)
        N5 = T_N5[0:3,3]
        # N6
        N6 = self.__trilaterate(N3,N4,N5,self.B[3],self.B[4],self.B[5])
        # Group nodes into a full list and return
        N = [N0,N1,N2,N3,N4,N5,N6]
        # Verify distances between node links
        return N

    # Updates all nodes within mechanism
    def update_state(self,debug=False):
        self.N = self.forward_kinematics(self.theta,self.phi)
        if debug: self.__validate_geo()

    # Sweeps over input vectors to find output surface
    def find_surface(self,n=10,m=24):
        self.X = np.zeros([n,m])
        self.Y = np.zeros([n,m])
        self.Z = np.zeros([n,m])
        for i,theta in enumerate(np.linspace(0,self.theta_max,n)):
            for j,phi in enumerate(np.linspace(0,2*np.pi,m)):
                N = self.forward_kinematics(theta,phi)
                self.X[i][j],self.Y[i][j],self.Z[i][j] = np.round(N[6][0],7),np.round(N[6][1],7),np.round(N[6][2],7)
        return np.std(self.X)
    ####################
    # Visualization
    ####################

    def plot_node(self,N):
        N = np.flip(N)
        x,y,z = N[0],N[1],N[2]
        self.ax.scatter(x,y,z,s=[40],color=[0.6350, 0.0780, 0.1840],zorder=3)

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
        for N in self.N:
            self.plot_node(N)
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
        self.ax.view_init(30,45)

        # Show plot is not in animation
        if not animation:
            plt.axis('equal') 
            plt.show()

    def plot_surface(self):
        self.find_surface()
        self.fig = plt.figure(3)
        self.ax = self.fig.add_subplot(projection='3d',computed_zorder=False)
        self.ax.plot_surface(self.Z, self.Y, self.X, color=[0.3010, 0.7450, 0.9330],edgecolor=[0, 0.4470, 0.7410, 0.5],zorder=2)   
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