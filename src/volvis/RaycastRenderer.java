/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volvis;

import com.jogamp.opengl.util.texture.Texture;
import com.jogamp.opengl.util.texture.awt.AWTTextureIO;
import gui.RaycastRendererPanel;
import gui.TransferFunction2DEditor;
import gui.TransferFunctionEditor;
import java.awt.image.BufferedImage;
import javax.media.opengl.GL;
import javax.media.opengl.GL2;
import util.TFChangeListener;
import util.VectorMath;
import volume.GradientVolume;
import volume.VoxelGradient;
import volume.Volume;

/***
 *
 * @author michel
 */
public class RaycastRenderer extends Renderer implements TFChangeListener {

    private Volume volume = null;
    private GradientVolume gradients = null;
    RaycastRendererPanel panel;
    TransferFunction tFunc;
    TransferFunctionEditor tfEditor;
    TransferFunction2DEditor tfEditor2D;
    public static String type = "slicer";
    public static String filename = "";
    
    public RaycastRenderer() {
        panel = new RaycastRendererPanel(this);
        panel.setSpeedLabel("0");
    }
    
    public void setType(String setType) {
        type = setType;
        changed();
    }
    
    public void currentFile(String f) {
        filename = f;
    }
     
    public void setVolume(Volume vol) {
        System.out.println("Assigning volume");
        volume = vol;

        System.out.println("Computing gradients");
        gradients = new GradientVolume(vol);

        // set up image for storing the resulting rendering
        // the image width and height are equal to the length of the volume diagonal
        int imageSize = (int) Math.floor(Math.sqrt(vol.getDimX() * vol.getDimX() + vol.getDimY() * vol.getDimY()
                + vol.getDimZ() * vol.getDimZ()));
        if (imageSize % 2 != 0) {
            imageSize = imageSize + 1;
        }
        image = new BufferedImage(imageSize, imageSize, BufferedImage.TYPE_INT_ARGB);
        // create a standard TF where lowest intensity maps to black, the highest to white, and opacity increases
        // linearly from 0.0 to 1.0 over the intensity range
        tFunc = new TransferFunction(volume.getMinimum(), volume.getMaximum());
        
        // uncomment this to initialize the TF with good starting values for the orange dataset 
        //tFunc.setTestFunc();
        
        tFunc.addTFChangeListener(this);
        tfEditor = new TransferFunctionEditor(tFunc, volume.getHistogram());
        
        tfEditor2D = new TransferFunction2DEditor(volume, gradients);
        tfEditor2D.addTFChangeListener(this);

        System.out.println("Finished initialization of RaycastRenderer");
    }

    public RaycastRendererPanel getPanel() {
        return panel;
    }

    public TransferFunction2DEditor getTF2DPanel() {
        return tfEditor2D;
    }
    
    public TransferFunctionEditor getTFPanel() {
        return tfEditor;
    }

    short getVoxel(double[] coord) {

        if (coord[0] < 0 || coord[0] > volume.getDimX() - 1 || coord[1] < 0 || coord[1] > volume.getDimY() - 1
                || coord[2] < 0 || coord[2] > volume.getDimZ() - 1) {
            return 0;
        }

        int x0 = (int) Math.floor(coord[0]);
        int y0 = (int) Math.floor(coord[1]);
        int z0 = (int) Math.floor(coord[2]);
        
        if(interactiveMode) {
            return volume.getVoxel(x0,y0,z0);
        }
        //trilinear interpolation
        int x1 = (int) x0 + 1;
        int y1 = (int) y0 + 1;
        int z1 = (int) z0 + 1;
        
        double xs = coord[0] - Math.floor(coord[0]); //segment of x, removes numbers before the komma, 11,7 becomes 0,7
        double ys = coord[1] - Math.floor(coord[1]);
        double zs = coord[2] - Math.floor(coord[2]);
        
        //making the slice, plane
        double p1 = xs * volume.getVoxel(x0, y0, z0) + (1 - xs) * volume.getVoxel(x1, y0, z0);
        double p2 = xs * volume.getVoxel(x0, y0, z1) + (1 - xs) * volume.getVoxel(x1, y0, z1);
        double p3 = xs * volume.getVoxel(x0, y1, z0) + (1 - xs) * volume.getVoxel(x1, y1, z0);
        double p4 = xs * volume.getVoxel(x0, y1, z1) + (1 - xs) * volume.getVoxel(x1, y1, z1);

        //every pn now contains the relative value, using the segment of x
        
        //making the line from b1 to b2 in the y direction
        double b1 = zs * p1 + (1 - zs) * p2;
        double b2 = zs * p3 + (1 - zs) * p4;
        
        //we now have a line from b1 to b2 over the y axis, 
        
        double p = (ys * b1 + (1 - ys) * b2);
        
        //we now have a value p with the voxel value

        
        return (short) p;
        // end of trilinear interpolation
    }


    void slicer(double[] viewMatrix) {
        //System.out.println("selected " + type);
        //type = "mip";
        
        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
                //System.out.println("viva la Nanne");
            }
        }

        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        double limit = 1.0;
        
        // Accept user input for gradient
        int gradI = tfEditor2D.triangleWidget.baseIntensity;
        double gradR = tfEditor2D.triangleWidget.radius;
        TFColor gradC = tfEditor2D.triangleWidget.color;
        
        TFColor voxelColor = new TFColor();
        
        for (int i = 0; i < image.getWidth(); i++) {
            for (int j = 0; j < image.getHeight(); j++) {
            
                
                if(type.equals("slicer")) {
                    // set gray color as default
                    tFunc.setTFcolor(filename, type);
                    
                    pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                            + volumeCenter[0];
                    pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                            + volumeCenter[1];
                    pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                            + volumeCenter[2];

                   int val = getVoxel(pixelCoord);

                    // Map the intensity to a grey value by linear scaling
                    voxelColor.r = val/max;
                    voxelColor.g = voxelColor.r;
                    voxelColor.b = voxelColor.r;
                    voxelColor.a = val > 0 ? 1.0 : 0.0;  // this makes intensity 0 completely transparent and the rest opaque
                    // Alternatively, apply the transfer function to obtain a color
                    //voxelColor = tFunc.getColor(val);


                    // BufferedImage expects a pixel color packed as ARGB in an int
                    int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                    int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                    int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                    int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                    int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                    image.setRGB(i, j, pixelColor);
                } else if(type.equals("mip")) {
                    //Maximum Intensity Projection
                    //System.out.println("selected " + type);
                    if(interactiveMode) {
                        if(i%3 == 0 && j%3 == 0) {
                            //Increase loading speed Maximum Intensity Projection
                            //System.out.println("selected " + type);

                            int maxVal = 0;
                            double maxRange = Math.abs(viewVec[0]) > (Math.abs(viewVec[1]) > Math.abs(viewVec[2]) ? volume.getDimY() : volume.getDimZ()) ? volume.getDimX() : (Math.abs(viewVec[1]) > Math.abs(viewVec[2]) ? volume.getDimY() : volume.getDimZ());
                            double minRange = maxRange*-1;
                            
                            //Loops through the pixels
                             for (double n = minRange; n < maxRange; n++) {
                                pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                                        + viewVec[0] * (n - (maxRange / 2)) + volumeCenter[0];
                                pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                                        + viewVec[1] * (n - (maxRange / 2)) + volumeCenter[1];
                                pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                                        + viewVec[2] * (n - (maxRange / 2)) + volumeCenter[2];

                                int val = getVoxel(pixelCoord);

                                if (val > maxVal) {
                                    maxVal = val;
                                }
                            }
                             // Map the intensity to a grey value by linear scaling
                                //voxelColor.r = val/max;
                                //voxelColor.g = voxelColor.r;
                                //voxelColor.b = voxelColor.r;
                                //voxelColor.a = val > minRange ? 1.0 : 0.0;  // this makes intensity 0 completely transparent and the rest opaque
                            // Apply the transfer function to obtain a color
                            voxelColor = tFunc.getColor(maxVal);

                            // BufferedImage expects a pixel color packed as ARGB in an int
                            int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                            int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                            int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                            int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                            int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                            
                            for(int ri = 0; ri < 4;ri++) {
                                for(int rj = 0; rj < 4;rj++) {
                                    if ((i + ri < image.getHeight()) && (j + rj < image.getWidth())) {
                                        image.setRGB(ri + i, rj + j, pixelColor); 
                                    }
                                }  
                            }
                        }   
                    } else {
                        //If optimization fails
                        //System.out.println("selected " + type);

                        int maxVal = 0;
                        double maxRange = Math.abs(viewVec[0]) > (Math.abs(viewVec[1]) > Math.abs(viewVec[2]) ? volume.getDimY() : volume.getDimZ()) ? volume.getDimX() : (Math.abs(viewVec[1]) > Math.abs(viewVec[2]) ? volume.getDimY() : volume.getDimZ());

                        //Loops through the pixels
                         for (int n = 0; n < maxRange; n++) {
                            pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                                    + viewVec[0] * (n - (maxRange / 2)) + volumeCenter[0];
                            pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                                    + viewVec[1] * (n - (maxRange / 2)) + volumeCenter[1];
                            pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                                    + viewVec[2] * (n - (maxRange / 2)) + volumeCenter[2];

                            int val = getVoxel(pixelCoord);

                            if (val > maxVal) {
                                maxVal = val;
                            }
                        }
                        // Apply the transfer function to obtain a color
                        voxelColor = tFunc.getColor(maxVal);

                        // BufferedImage expects a pixel color packed as ARGB in an int
                        int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                        int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                        int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                        int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                        int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;

                        image.setRGB(i,j, pixelColor); 

                    }
                } else if(type.equals("comp")) {
                //Composite rendering (Direct Volume Rendering - DVR)
                    if(interactiveMode) {
                        if(i%3 == 0 && j%3 == 0) {

                    // set colorscheme according to filename
                    tFunc.setTFcolor(filename, type);

                    TFColor compColor = new TFColor(0, 0, 0, 0);
                    double maxRange = Math.abs(viewVec[0]) > (Math.abs(viewVec[1]) > Math.abs(viewVec[2]) ? volume.getDimY() : volume.getDimZ()) ? volume.getDimX() : (Math.abs(viewVec[1]) > Math.abs(viewVec[2]) ? volume.getDimY() : volume.getDimZ());

                    //Loops through the pixels
                     for (int n = 0; n < maxRange; n++) {
                        pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                                + viewVec[0] * (n - (maxRange / 2)) + volumeCenter[0];
                        pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                                + viewVec[1] * (n - (maxRange / 2)) + volumeCenter[1];
                        pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                                + viewVec[2] * (n - (maxRange / 2)) + volumeCenter[2];

                        int val = getVoxel(pixelCoord);

                        // Apply the transfer function to obtain a color
                        voxelColor = tFunc.getColor(val);

                        compColor.a = voxelColor.a * voxelColor.a + (1 - voxelColor.a) * compColor.a;
                        compColor.r = voxelColor.r * voxelColor.a + (1 - voxelColor.a) * compColor.r;
                        compColor.g = voxelColor.g * voxelColor.a + (1 - voxelColor.a) * compColor.g;
                        compColor.b = voxelColor.b * voxelColor.a + (1 - voxelColor.a) * compColor.b;
                    }

                    // BufferedImage expects a pixel color packed as ARGB in an int;
                    int c_alpha = compColor.a <= 1.0 ? (int) Math.floor(compColor.a * 255) : 255;
                    int c_red = compColor.r <= 1.0 ? (int) Math.floor(compColor.r * 255) : 255;
                    int c_green = compColor.g <= 1.0 ? (int) Math.floor(compColor.g * 255) : 255;
                    int c_blue = compColor.b <= 1.0 ? (int) Math.floor(compColor.b * 255) : 255;

                    int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;

                    for(int ri = 0; ri < 4;ri++) {
                                for(int rj = 0; rj < 4;rj++) {
                                    if ((i + ri < image.getHeight()) && (j + rj < image.getWidth())) {
                                        image.setRGB(ri + i, rj + j, pixelColor); 
                                    }
                                }  
                            }
                        }   
                    } else {
                        //If optimization fails
                        //System.out.println("selected " + type);

                        // set colorscheme according to filename
                        tFunc.setTFcolor(filename, type);

                        TFColor compColor = new TFColor(0, 0, 0, 0);
                        double maxRange = Math.abs(viewVec[0]) > (Math.abs(viewVec[1]) > Math.abs(viewVec[2]) ? volume.getDimY() : volume.getDimZ()) ? volume.getDimX() : (Math.abs(viewVec[1]) > Math.abs(viewVec[2]) ? volume.getDimY() : volume.getDimZ());

                        //Loops through the pixels
                         for (int n = 0; n < maxRange; n++) {
                            pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                                    + viewVec[0] * (n - (maxRange / 2)) + volumeCenter[0];
                            pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                                    + viewVec[1] * (n - (maxRange / 2)) + volumeCenter[1];
                            pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                                    + viewVec[2] * (n - (maxRange / 2)) + volumeCenter[2];

                            int val = getVoxel(pixelCoord);

                            // Apply the transfer function to obtain a color
                            voxelColor = tFunc.getColor(val);

                            compColor.a = voxelColor.a * voxelColor.a + (1 - voxelColor.a) * compColor.a;
                            compColor.r = voxelColor.r * voxelColor.a + (1 - voxelColor.a) * compColor.r;
                            compColor.g = voxelColor.g * voxelColor.a + (1 - voxelColor.a) * compColor.g;
                            compColor.b = voxelColor.b * voxelColor.a + (1 - voxelColor.a) * compColor.b;
                        }

                        // BufferedImage expects a pixel color packed as ARGB in an int;
                        int c_alpha = compColor.a <= 1.0 ? (int) Math.floor(compColor.a * 255) : 255;
                        int c_red = compColor.r <= 1.0 ? (int) Math.floor(compColor.r * 255) : 255;
                        int c_green = compColor.g <= 1.0 ? (int) Math.floor(compColor.g * 255) : 255;
                        int c_blue = compColor.b <= 1.0 ? (int) Math.floor(compColor.b * 255) : 255;

                        int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                        
                        image.setRGB(i,j, pixelColor); 
                    }
                } else if(type.equals("gradient")){
                    // Gradient-based opacity weighting
                    TFColor predecessor = new TFColor();
                    TFColor successor = new TFColor(); 
                        
                    pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                            + volumeCenter[0];
                    pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                            + volumeCenter[1];
                    pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                            + volumeCenter[2];

                    int val = getVoxel(pixelCoord);
                        
                    if ( ( (pixelCoord[0] < volume.getDimX() && pixelCoord[0] >= 0) || (pixelCoord[1] < volume.getDimY() && pixelCoord[1] >= 0) || (pixelCoord[2] < volume.getDimZ() && pixelCoord[2] >= 0) ) && val > limit) {
                        // Get user defined color
                        voxelColor = gradC;
                        VoxelGradient voxGrad = gradients.getGradient((int)Math.floor(pixelCoord[0]), (int)Math.floor(pixelCoord[1]), (int)Math.floor(pixelCoord[2])); 

                        // gradient intensity check
                        if (val == gradI && voxGrad.mag == 0) {
                            voxelColor.a = 1.0;
                        }
                        else if (voxGrad.mag > 0.0 && ((val - gradR * voxGrad.mag) <= gradI) && ((val + gradR * voxGrad.mag) >= gradI)){
                            voxelColor.a = 1.0 - (1 / gradR) * (Math.abs((gradI - val)/ voxGrad.mag));
                        }
                        else 
                            voxelColor.a = 0.0;

                        // apply opacity weights to colors of voxels and their successors                
                        successor.r = voxelColor.a * voxelColor.r + (1 - voxelColor.a) * predecessor.r;
                        successor.g = voxelColor.a * voxelColor.g + (1 - voxelColor.a) * predecessor.g;
                        successor.b = voxelColor.a * voxelColor.b + (1 - voxelColor.a) * predecessor.b;

                        successor.a = (1 - voxelColor.a) * predecessor.a;

                        predecessor = successor;

                    }

                        // BufferedImage expects a pixel color packed as ARGB in an int;
                        int c_alpha = (1 - successor.a) <= 1.0 ? (int) Math.floor((1 - successor.a) * 255) : 255;
                        int c_red = successor.r <= 1.0 ? (int) Math.floor(successor.r * 255) : 255;
                        int c_green = successor.g <= 1.0 ? (int) Math.floor(successor.g * 255) : 255;
                        int c_blue = successor.b <= 1.0 ? (int) Math.floor(successor.b * 255) : 255;

                        int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;

                        // Set multiple pixels at lower resolution
                        image.setRGB(i, j, pixelColor);  
                    }
            }
        }
    }


    private void drawBoundingBox(GL2 gl) {
        gl.glPushAttrib(GL2.GL_CURRENT_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glColor4d(1.0, 1.0, 1.0, 1.0);
        gl.glLineWidth(1.5f);
        gl.glEnable(GL.GL_LINE_SMOOTH);
        gl.glHint(GL.GL_LINE_SMOOTH_HINT, GL.GL_NICEST);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glDisable(GL.GL_LINE_SMOOTH);
        gl.glDisable(GL.GL_BLEND);
        gl.glEnable(GL2.GL_LIGHTING);
        gl.glPopAttrib();

    }

    @Override
    public void visualize(GL2 gl) {


        if (volume == null) {
            return;
        }

        drawBoundingBox(gl);

        gl.glGetDoublev(GL2.GL_MODELVIEW_MATRIX, viewMatrix, 0);

        //long startTime = System.currentTimeMillis();
        slicer(viewMatrix);    
        
        //long endTime = System.currentTimeMillis();
        //double runningTime = (endTime - startTime);
        //panel.setSpeedLabel(Double.toString(runningTime));

        Texture texture = AWTTextureIO.newTexture(gl.getGLProfile(), image, false);

        gl.glPushAttrib(GL2.GL_LIGHTING_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        // draw rendered image as a billboard texture
        texture.enable(gl);
        texture.bind(gl);
        double halfWidth = image.getWidth() / 2.0;
        gl.glPushMatrix();
        gl.glLoadIdentity();
        gl.glBegin(GL2.GL_QUADS);
        gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
        gl.glTexCoord2d(0.0, 0.0);
        gl.glVertex3d(-halfWidth, -halfWidth, 0.0);
        gl.glTexCoord2d(0.0, 1.0);
        gl.glVertex3d(-halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 1.0);
        gl.glVertex3d(halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 0.0);
        gl.glVertex3d(halfWidth, -halfWidth, 0.0);
        gl.glEnd();
        texture.disable(gl);
        texture.destroy(gl);
        gl.glPopMatrix();

        gl.glPopAttrib();


        if (gl.glGetError() > 0) {
            System.out.println("some OpenGL error: " + gl.glGetError());
        }

    }
    private BufferedImage image;
    private double[] viewMatrix = new double[4 * 4];

    @Override
    public void changed() {
        for (int i=0; i < listeners.size(); i++) {
            listeners.get(i).changed();
        }
    }

    public void setFileName(String filename) {
        this.filename = filename;
    }
}
