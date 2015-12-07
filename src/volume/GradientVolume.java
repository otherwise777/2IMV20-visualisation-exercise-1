/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volume;

import volvis.TransferFunction;

/**
 *
 * @author michel
 */
public class GradientVolume {
    
    TransferFunction tFunc;
    private Volume gradVolume = null;
    public boolean DoneMakingGrad = false; 
    
    private void getGradientArray() {
        
        gradVolume = new Volume(volume.getDimX(), volume.getDimY(), volume.getDimZ());
    
        //Gradient function
        for (int x = 0; x < dimX; x++) {
            for (int y = 0; y < dimY; y++) {
                for (int z = 0; z < dimZ; z++) {
                    double gv_1 = 0.5 * (volume.getVoxel(x + 1, y, z) - volume.getVoxel(x - 1, y, z));
                    double gv_2 = 0.5 * (volume.getVoxel(x, y + 1, z) - volume.getVoxel(x, y - 1, z));
                    double gv_3 = 0.5 * (volume.getVoxel(x, y, z + 1) - volume.getVoxel(x, y, z - 1));
                    
                    double gradMagVal = Math.sqrt(Math.pow(gv_1, 2) + Math.pow(gv_2, 2) + Math.pow(gv_3, 2));
                    
                    double a;
                    double aTotal = 1;
                    for (int p = 1; p < tFunc.getControlPoints().size() - 1; p++) {
                        double funcVal = ((double) tFunc.getControlPoints().get(p).value) / 255;
                        double alphaVal = tFunc.getControlPoints().get(p).color.a;
                        
                        // Equation to get a(xi)
                        double Fvp1 = ((double) tFunc.getControlPoints().get(p + 1).value) / 255;
                        short voxelVal = volume.getVoxel(x, y, z);
                        double fxi = ((double) voxelVal) / 255;
                        if (funcVal <= fxi && fxi <= Fvp1) {
                            double alphaVp1 = tFunc.getControlPoints().get(p + 1).color.a;
                            a = (gradMagVal / 255) * (alphaVp1 * ((fxi - funcVal) / (Fvp1 - funcVal)) + alphaVal * ((Fvp1 - fxi) / (Fvp1 - funcVal)));
                        } else {
                            a = 0;
                        }

                        aTotal = (1 - a) * aTotal;
                    }
                    aTotal = 1 - aTotal;

                    gradVolume.setVoxel(x, y, z, (short) (aTotal * 255));
                }
            }
        }
        DoneMakingGrad = true;
    }

    public GradientVolume(Volume vol) {
        volume = vol;
        dimX = vol.getDimX();
        dimY = vol.getDimY();
        dimZ = vol.getDimZ();
        data = new VoxelGradient[dimX * dimY * dimZ];
        compute();
        maxmag = -1.0;
    }

    public VoxelGradient getGradient(int x, int y, int z) {
        return data[x + dimX * (y + dimY * z)];
    }

    
    public void setGradient(int x, int y, int z, VoxelGradient value) {
        data[x + dimX * (y + dimY * z)] = value;
    }

    public void setVoxel(int i, VoxelGradient value) {
        data[i] = value;
    }

    public VoxelGradient getVoxel(int i) {
        return data[i];
    }

    public int getDimX() {
        return dimX;
    }

    public int getDimY() {
        return dimY;
    }

    public int getDimZ() {
        return dimZ;
    }

    private void compute() {

        // this just initializes all gradients to the vector (0,0,0)
        for (int i=0; i<data.length; i++) {
            data[i] = zero;
        }
                
    }
    
    public double getMaxGradientMagnitude() {
        if (maxmag >= 0) {
            return maxmag;
        } else {
            double magnitude = data[0].mag;
            for (int i=0; i<data.length; i++) {
                magnitude = data[i].mag > magnitude ? data[i].mag : magnitude;
            }   
            maxmag = magnitude;
            return magnitude;
        }
    }
    
    private int dimX, dimY, dimZ;
    private VoxelGradient zero = new VoxelGradient();
    VoxelGradient[] data;
    Volume volume;
    double maxmag;
}