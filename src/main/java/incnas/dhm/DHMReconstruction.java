package incnas.dhm;

import ij.IJ;
import ij.ImagePlus;
import ij.measure.Measurements;
import ij.measure.ResultsTable;
import ij.plugin.ContrastEnhancer;
import ij.plugin.filter.ParticleAnalyzer;
import ij.process.AutoThresholder;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ijopencv.ij.ImagePlusMatConverter;
import ijopencv.opencv.MatImagePlusConverter;
import incnas.dhm.utils.CurveFunction;
import net.imagej.Dataset;
import net.imagej.ImageJ;
import net.imagej.ImgPlus;
import net.imagej.ops.OpService;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.RealType;
import org.bytedeco.javacpp.opencv_core;
import org.bytedeco.javacpp.opencv_phase_unwrapping;
import org.ddogleg.optimization.FactoryOptimization;
import org.ddogleg.optimization.UnconstrainedLeastSquares;
import org.ddogleg.optimization.UtilOptimize;
import org.ddogleg.optimization.functions.FunctionNtoM;
import org.ejml.data.DMatrixRMaj;
import org.jtransforms.fft.FloatFFT_2D;
import org.scijava.command.Command;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;
import org.scijava.plugin.PluginService;
import org.scijava.ui.UIService;
import unal.od.jdiffraction.cpu.FloatAngularSpectrum;
import unal.od.jdiffraction.cpu.utils.ArrayUtils;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

//TODO: Use more efficient array elementwise operations?

@Plugin(type = Command.class, menuPath = "Plugins>DHM Reconstruction")
public class DHMReconstruction <T extends RealType<T>> implements Command {
    private boolean debug = false;

    @Parameter
    private Dataset currentData;

    @Parameter
    private PluginService pluginService;

    @Parameter
    private UIService uiService;

    @Parameter
    private OpService opService;

    private int size_x,size_y;
    private float z, lambda, dx, dy;

    private float[][] cropOrders(ImageProcessor ip, int small_const, int large_const){
        int iy = 1; //Number of objects
        float inc = 0; //Increment
        int sigma = 2;

        int complex_size_y = size_y*2;

        float[][] f_img = ip.getFloatArray();

        //Spread to allow for complex data by doubling array size
        float[][] img = new float[size_x][complex_size_y];

        for (int i = 0 ; i < size_x; i++){
            for (int j=0; j<size_y; j++){
                img[i][j*2] = (float) f_img[i][j];
            }
        }

        //FFT
        float[][] imgFFT = img.clone();
        FloatFFT_2D fft2D = new FloatFFT_2D(size_x,size_y);

        fft2D.complexForward(imgFFT);
        ArrayUtils.complexShift(imgFFT);

        //FFT Log Magnitude
        float[][] magImgFFt = ArrayUtils.modulusSq(imgFFT);
        FloatProcessor floatProcessor = new FloatProcessor(magImgFFt);
        System.out.println(floatProcessor.getMax());
        floatProcessor.add(1);
        floatProcessor.resetMinAndMax();
        System.out.println(floatProcessor.getMax());
        floatProcessor.log();
        floatProcessor.resetMinAndMax();
        System.out.println(floatProcessor.getMax());

        //Normalise
        floatProcessor.subtract(floatProcessor.getMin());
        floatProcessor.multiply(1/floatProcessor.getMax());

        FloatProcessor blurred = (FloatProcessor) floatProcessor.duplicate();

        //Gaussian blur
        blurred.blurGaussian(sigma); //TODO: Check gaussian size
        new ContrastEnhancer().equalize(blurred);

        //otsu threshold
        blurred.multiply(255.0);
        //blurred.setAutoThreshold(AutoThresholder.Method.Otsu,true);
        //int level = blurred.getAutoThreshold();
        float level = new AutoThresholder().getThreshold(AutoThresholder.Method.Otsu, blurred.getHistogram());
        float step = level / 100.0f;

        int small_img = Math.round((((float) size_x)/small_const)*(((float) size_y)/small_const));
        int big_img = Math.round((((float) size_x)/large_const)*(((float) size_y)/large_const));

        ImageProcessor mask = null;
        ResultsTable rt = null;

        //TODO: Try with float instead of int mask
        while (iy != 3){
            System.out.println(level+inc);
            float[][] float_mask = blurred.getFloatArray();
            byte[] int_mask = new byte[size_x*size_y];

            for(int i=0; i<size_x;i++){
                for(int j=0; j<size_y;j++){
                    int_mask[(j*size_x)+i] = (byte) Math.round(float_mask[i][j]);
                }
            }
            //mask = new ShortProcessor(size_x,size_y,int_mask, new java.awt.image.ColorModel(DataBuffer.TYPE_USHORT));
            mask = new ByteProcessor(size_x,size_y,int_mask);
            //mask.setThreshold(140,255,ImageProcessor.RED_LUT);
            mask.threshold((int) (level+inc));
            mask.invert();

            Img ImageMask = ImageJFunctions.wrap(new ImagePlus("Mask",mask));

            /*
            float[][] mask_array = mask.getFloatArray();
            float[] flat_mask = new float[size_x*size_y];

            //Create mask
            for(int i=0; i<size_x;i++){
                for(int j=0; j<size_y;j++){
                    flat_mask[i*size_y+j] = (mask_array[i][j] < level+inc) ? 0 : mask_array[i][j];
                }
            }

            //Remove small objects
            */

            rt = new ResultsTable();
            Double max_size = Double.POSITIVE_INFINITY;
            ParticleAnalyzer pa = new ParticleAnalyzer(ParticleAnalyzer.SHOW_MASKS + ParticleAnalyzer.IN_SITU_SHOW,
                    Measurements.CENTROID + Measurements.SHAPE_DESCRIPTORS + Measurements.AREA + Measurements.RECT, rt, small_img, max_size);
            pa.analyze(new ImagePlus("Mask",mask));

            if(debug) {
                System.out.println(rt);
                System.out.println("----------------------");
                System.out.println(rt.size());
                System.out.println(Arrays.toString(rt.getColumn(ResultsTable.AREA)));
                System.out.println(Arrays.toString(rt.getColumn(ResultsTable.X_CENTROID)));
                System.out.println(Arrays.toString(rt.getColumn(ResultsTable.Y_CENTROID)));
                System.out.println(Arrays.toString(rt.getColumn(ResultsTable.ROI_HEIGHT)));
                System.out.println(Arrays.toString(rt.getColumn(ResultsTable.ROI_WIDTH)));
                System.out.println(Arrays.toString(rt.getColumn(ResultsTable.ROI_X)));
                System.out.println(Arrays.toString(rt.getColumn(ResultsTable.ROI_Y)));
            }

            iy = rt.size();

            if (rt.size() == 3) {

                for (float x : rt.getColumn(ResultsTable.AREA)) {
                    if(x > big_img){
                        iy = 1;
                    }
                }
            }

            inc += step;

            if (inc > 250){
                IJ.handleException(new Exception("Failed to find threshold"));
                if(debug){
                    ImagePlus fft = new ImagePlus("FFT ", floatProcessor);
                    uiService.show("FFT",ImageJFunctions.wrap(fft));
                }

                return null;
            }
        }

        //Crop to region
        int crop_width = (int) rt.getColumn(ResultsTable.ROI_WIDTH)[0];
        int crop_height = (int) rt.getColumn(ResultsTable.ROI_HEIGHT)[0]*2;

        int crop_x = (int) rt.getColumn(ResultsTable.ROI_X)[0];
        int crop_y = (int) rt.getColumn(ResultsTable.ROI_Y)[0];
        float[][] cropped_area = new float[crop_width][crop_height];

        //Retrieve complex masked
        for(int i=0; i<crop_width; i++){
            for (int j=0; j<crop_height; j++){
                cropped_area[i][j] = imgFFT[i+crop_x][j+crop_y*2];
            }
        }

        //Move object into center of image
        float[][] centered_complex_img = new float[size_x][size_y*2];
        int displacement_x = (size_x-crop_width)/2;
        int displacement_y = (complex_size_y-crop_height)/2;

        for(int i=0; i<crop_width; i++){
            for(int j=0; j<crop_height; j++){
                centered_complex_img[i+displacement_x][j+displacement_y] = cropped_area[i][j];
            }
        }

        float[][] centered_angle = ArrayUtils.modulusSq(centered_complex_img);
        FloatProcessor angleProcessor = new FloatProcessor(centered_angle);
        angleProcessor.add(1);
        angleProcessor.log();

        //debug show selection
        ImagePlus fft = new ImagePlus("FFT ", floatProcessor);
        uiService.show("FFT",ImageJFunctions.wrap(fft));

        ImagePlus mask_img = new ImagePlus("Mask ", mask);
        uiService.show("Mask",ImageJFunctions.wrap(mask_img));

        ImagePlus angle_img = new ImagePlus("Cropped Mag", angleProcessor);
        uiService.show("Cropped Mag",ImageJFunctions.wrap(angle_img));

        //return centered image
        return centered_complex_img;
    }

    private ImagePlus retrieve_phase(float[][] centered_fftImg){
        //TODO: multiply by 2/pi?

        //Propagate fft - TODO
        FloatAngularSpectrum as = new FloatAngularSpectrum(this.size_x,this.size_y,this.lambda,this.z,this.dx,this.dy);
        float[][] diffracted = java.util.Arrays.stream(centered_fftImg).map(float[]::clone).toArray(float[][]::new);
        as.diffract(centered_fftImg);

        // Loop through all rows


        //IFFT
        FloatFFT_2D fft2D = new FloatFFT_2D(size_x,size_y);
        ArrayUtils.complexShift(centered_fftImg);
        fft2D.complexInverse(centered_fftImg,false);

        //retrieve phase
        float[][] phase = ArrayUtils.phase(centered_fftImg);
        FloatProcessor phase_processor = new FloatProcessor(phase);

        ImagePlus wrapped_phase = new ImagePlus("Wrapped Phase ", phase_processor);

        //unwrap phase
        ImagePlusMatConverter ic = new ImagePlusMatConverter();
        MatImagePlusConverter mip = new MatImagePlusConverter();

        opencv_core.Mat wrapped_phase_m = ic.convert(wrapped_phase, opencv_core.Mat.class);
        opencv_core.Mat unwrapped_phase_m = new opencv_core.Mat();

        opencv_phase_unwrapping.HistogramPhaseUnwrapping.Params params = new opencv_phase_unwrapping.HistogramPhaseUnwrapping.Params();
        params.width(wrapped_phase.getWidth());
        params.height(wrapped_phase.getHeight());

        opencv_phase_unwrapping.HistogramPhaseUnwrapping phaseUnwrapper = opencv_phase_unwrapping.HistogramPhaseUnwrapping.create(params);
        phaseUnwrapper.unwrapPhaseMap(wrapped_phase_m,unwrapped_phase_m);

        ImagePlus unwrapped_phase = mip.convert(unwrapped_phase_m,ImagePlus.class);

        //FunctionNtoM test = new FunctionLineDistanceEuclidean(points);

        //remove curve
        FunctionNtoM func = new CurveFunction((FloatProcessor) unwrapped_phase.getProcessor());
        UnconstrainedLeastSquares<DMatrixRMaj> optimizer = FactoryOptimization.levenbergMarquardt(null, true);
        optimizer.setFunction(func,null);

        optimizer.initialize(new double[]{1,1,1,1,1,1},1e-12,1e-12);
        UtilOptimize.process(optimizer,500);

        double optimal[] = new double[unwrapped_phase.getWidth()*unwrapped_phase.getHeight()];
        func.process(optimizer.getParameters(),optimal);

        float[][] flat_phase = new float[unwrapped_phase.getWidth()][unwrapped_phase.getHeight()];
        for(int i=0; i<unwrapped_phase.getWidth(); i++){
            for (int j=0; j<unwrapped_phase.getHeight(); j++){
                flat_phase[i][j] = (float) optimal[i*unwrapped_phase.getHeight()+j];
            }
        }

        FloatProcessor flat_processor = new FloatProcessor(flat_phase);
        ImagePlus flat_img = new ImagePlus("Flat Phase", flat_processor);

        if (debug) {
            uiService.show("Unwrapped Phase", ImageJFunctions.wrap(unwrapped_phase));

            uiService.show("Wrapped Phase", ImageJFunctions.wrap(wrapped_phase));

            uiService.show("Flat Phase", ImageJFunctions.wrap(flat_img));
        }

        return flat_img;
    }

    private RandomAccessibleInterval<T> retrieve_magnitude(float[][] centered_fftImg){
        //Propagate fft

        //IFFT

        //retrieve magnitude

        //remove curve

        //return magnitude

        return null;
    }


    private void reconstruct(ImagePlus img){
        ImageProcessor ip = img.getProcessor();

        size_x = ip.getWidth();
        size_y = ip.getHeight();

        float[][] centered_fftImg = cropOrders(ip, 25,4);

        ImagePlus phase = retrieve_phase(centered_fftImg);
        RandomAccessibleInterval magnitude = retrieve_magnitude(centered_fftImg);

        //uiService.show(phase);
        //uiService.show(magnitude);
    }

    @Override
    public void run() {
        z = 0;
        lambda = 0.638f;
        dx = 0.434f;
        dy = 0.434f;

        //nu.pattern.OpenCV.loadShared();

        final ImgPlus image = currentData.getImgPlus();

        /*float[][] float_img = new float[(int) image.dimension(0)][(int) image.dimension(1)];

        float[] list_img = new float[(int) image.dimension(0)*(int) image.dimension(1)];

        int i=0;
        for (final T type : (Img<T>)image.getImg()){

        }

        for (int i=0; i<image.dimension(0);i++){
            for(int j=0; j<image.dimension(1);j++){
                float_img[i][j] =
            }
        }

        new ImagePlus("Image",new FloatProcessor(new float[10][10]));*/

        ImagePlus img = ImageJFunctions.wrap(image,"Input Image");

        System.out.println(img==null);

        //
        // Enter image processing code here ...
        // The following is just a Gauss filtering example
        //
        final double[] sigmas = {1.0, 3.0, 5.0};

        List<RandomAccessibleInterval<T>> results = new ArrayList<>();

        for (double sigma : sigmas) {
            results.add(opService.filter().gauss(image, sigma));
        }

        RandomAccessibleInterval<T> fft_res = opService.filter().fft(image);

        reconstruct(img);

        //IJ.handleException(new Exception("Test Exception"));

        //fft_res.

        //uiService.show(fft_res);

        // display result
        //for (RandomAccessibleInterval<T> elem : results) {
        //    uiService.show(elem);
        //}
    }

    /**
     * This main function serves for development purposes.
     * It allows you to run the plugin immediately out of
     * your integrated development environment (IDE).
     *
     * @param args whatever, it's ignored
     * @throws Exception
     */
    public static void main(final String... args) throws Exception {
        // create the ImageJ application context with all available services
        final ImageJ ij = new ImageJ();
        ij.ui().showUI();

        // ask the user for a file to open
        final File file = ij.ui().chooseFile(null, "open");

        if (file != null) {
            // load the dataset
            final Dataset dataset = ij.scifio().datasetIO().open(file.getPath());

            // show the image
            ij.ui().show(dataset);

            // invoke the plugin
            ij.command().run(DHMReconstruction.class, true);
        }
    }
}
