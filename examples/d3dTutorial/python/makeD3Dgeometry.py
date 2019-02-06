import io,libconf
import os
import numpy as np

def makeGeom(inFile='d3drzZ.csv',filename="gitrD3DGeometry2DWringsPy.cfg"):
    float_format = "%.6f"
    int_format=1
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
    vertical_lines = np.where((x2-x1)==0)
    horizontal_lines = np.where((z2-z1)==0)
    intercept = -slope*x1 - z1
    slope[vertical_lines] = np.sign(z2[vertical_lines]-z1[vertical_lines])*1.0e12;
    intercept[vertical_lines] = 1.0e12;
    slope[horizontal_lines]=0;
    intercept[horizontal_lines] = z1[horizontal_lines];
    
    line_length = np.sqrt(np.multiply(x2-x1,x2-x1) + np.multiply(z2-z1,z2-z1))
    surface = np.zeros(x1.size+1)
    surface[np.where(Z>0)]=1

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
    makeGeom(inFile='d3drzZ.csv',filename="gitrD3DGeometry2DWringsPy1.cfg")
    removeQuotes(infile="gitrD3DGeometry2DWringsPy1.cfg",outfile="gitrD3DGeometry2DWringsPy.cfg")
