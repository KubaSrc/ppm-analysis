from src.PPM import ppm
import numpy as np

# Different visualizations for PPM
def vis_ppm(mode = "draw"):
    # Generate mechanism
    PPM = ppm.mechanism(A=25.15,B=38,C=76,D=58.25)
    # Draw in a given position
    if mode == "draw":
        PPM.update_angles(theta=np.deg2rad(27.5),phi=np.deg2rad(30),debug=False) # 30 - 35
        PPM.draw()
        np.savetxt("nodes.csv",PPM.N,delimiter=',')
        # np.savetxt("x.csv",PPM.X,delimiter=',')
        # np.savetxt("y.csv",PPM.Y,delimiter=',')
        # np.savetxt("z.csv",PPM.Z,delimiter=',')
    # Animate across range of angles
    if mode == "surf":
        PPM.plot_surface()
    if mode == "animate":
        sweep = []
        for theta in np.linspace(0,PPM.theta_max,10):
            for phi in np.linspace(0,2*np.pi,10):
                sweep.append((theta,phi))
        PPM.animate(sweep,save=True)
    # Show surface plot of PPM

if __name__ == "__main__":
    vis_ppm("draw")
