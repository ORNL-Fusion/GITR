import io,libconf
import os
import numpy as np
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import math

def makeGeom(inFile='d3drzZ.csv',filename="gitrD3DGeometry2DWringsPy.cfg",plot=1):
    float_format = "%.6f"
    rzZ=np.loadtxt(inFile, dtype='float',delimiter=',')    
    Z = rzZ[:,2]
    Z = np.append(Z,0)
    x1 = rzZ[:,0]
    x2 = rzZ[1:,0]
    x2 = np.append(x2,rzZ[0,0])
    z1 = rzZ[:,1]
    z2 = rzZ[1:,1]
    z2 = np.append(z2,rzZ[0,1])
    slope = (z2-z1)/(x2-x1)
    tol = 1.0e-9
    vertical_lines = np.where(abs(x2-x1) < tol)
    horizontal_lines = np.where(abs(z2-z1) < tol)
    intercept = -slope*x1 + z1
    print('x1 z1 ',x1[0:3],z1[0:3])
    print('x2 z2 ',x2[0:3],z2[0:3])
    print('slope intercept', slope[0:3], intercept[0:3])
    slope[vertical_lines] = np.sign(z2[vertical_lines]-z1[vertical_lines])*1.0e12;
    intercept[vertical_lines] = 1.0e12;
    print('slope intercept', slope[0:3], intercept[0:3])
    slope[horizontal_lines]=0;

    inDir = np.ones(x1.size+1)
    reverseDir =[0,1,2,4,5,6,7,8,9,10,11,17,31,49,63]
    inDir[reverseDir] = -1;
    intercept[horizontal_lines] = z1[horizontal_lines];
    print('slope intercept', slope[0:3], intercept[0:3])
    
    line_length = np.sqrt(np.multiply(x2-x1,x2-x1) + np.multiply(z2-z1,z2-z1))
    surface = np.zeros(x1.size+1)
    w_indices = np.where(Z>0)
    w_indices = w_indices[0]
    w_indices2 = np.append(w_indices[1:len(w_indices)],w_indices[0])
    geom_indices = np.where(np.subtract(w_indices2,w_indices) == 1)
    surface[w_indices[geom_indices]]=1
    
    if plot ==1:
        plt.close() 
        plt.figure(1,figsize=(6, 10), dpi=60)
        for i in range(len(x1)):
            if (slope[i] == 0.0):
                perpSlope = 1.0e12;
            else:
                perpSlope = math.copysign(1.0,slope[i]) / abs(slope[i]);
            Bx = -1.0 /math.sqrt(perpSlope * perpSlope + 1.0);
            By = 0.0;
            Bz = math.copysign(1.0,perpSlope) * math.sqrt(1.0 - Bx * Bx);
            plt.quiver(0.5*(x1[i]+x2[i]),0.5*(z1[i]+z2[i]),Bx*inDir[i],Bz*inDir[i])
            plt.text(0.5*(x1[i]+x2[i]),0.5*(z1[i]+z2[i]),str(i))
            if(surface[i]==1):
                plot_color='green'
                surface_plot, =plt.plot([x1[i],x2[i]],[z1[i],z2[i]],linewidth=2.0,color=plot_color,label='Surface')
            else:
                plot_color='blue'
                boundary_plot, =plt.plot([x1[i],x2[i]],[z1[i],z2[i]],linewidth=2.0,color=plot_color,label='Boundary')

       
        plt.scatter(x1,z1)
        plt.autoscale(enable=True, axis='x', tight=True)
        plt.title('DIII-D W Geometry Creation',fontsize=20)
        plt.xlabel('r [m]',fontsize=16)
        plt.ylabel('z [m]',fontsize=16)
        plt.legend(handles=[surface_plot,boundary_plot])
        print('saving geometry plot')
        plt.savefig('geometry.png')
        plt.show()

    fileExists = os.path.exists(filename)
    if not fileExists:
        f=open(filename,"w")
        f.close()
    
    with io.open(filename) as f:
        config = libconf.load(f)
    
    config['geom'] = {}
    config['geom']['x1'] = x1.tolist()   
    config['geom']['z1'] = z1.tolist()
    config['geom']['x2'] = x2.tolist()   
    config['geom']['z2'] = z2.tolist()
    config['geom']['slope'] = ['%.6f' % elem for elem in slope.tolist()]
    config['geom']['intercept'] = ['%.6f' % elem for elem in intercept.tolist()]
    config['geom']['length'] = line_length.tolist()
    config['geom']['Z'] = Z.tolist()
    config['geom']['surface'] = ['%i' % elem for elem in surface.tolist()]
    config['geom']['inDir'] = ['%i' % elem for elem in inDir.tolist()]
    config['geom']['y1'] = 0.0
    config['geom']['y2'] = 0.0
    config['geom']['periodic'] = 0

    with io.open(filename,'w') as f:
        libconf.dump(config,f)
def removeQuotes(infile='this.cfg',outfile='that.cfg'):
    with open(infile, 'r') as f, open(outfile, 'w') as fo:
        for line in f:
            fo.write(line.replace('"', ''))

if __name__ == "__main__":
    makeGeom(inFile='python/d3drzZ.csv',filename="gitrD3DGeometry2DWringsPy1.cfg",plot=1)
    removeQuotes(infile="gitrD3DGeometry2DWringsPy1.cfg",outfile="gitrD3DGeometry2DWringsPy.cfg")
