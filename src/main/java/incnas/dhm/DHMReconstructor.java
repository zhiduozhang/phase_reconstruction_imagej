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
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import org.bytedeco.javacpp.opencv_core;
import org.bytedeco.javacpp.opencv_phase_unwrapping;
import org.ddogleg.optimization.FactoryOptimization;
import org.ddogleg.optimization.UnconstrainedLeastSquares;
import org.ddogleg.optimization.UtilOptimize;
import org.ddogleg.optimization.functions.FunctionNtoM;
import org.ejml.data.DMatrixRMaj;
import org.jtransforms.fft.FloatFFT_2D;
import unal.od.jdiffraction.cpu.FloatAngularSpectrum;
import unal.od.jdiffraction.cpu.utils.ArrayUtils;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

public class DHMReconstructor {
    public final static String FFT = "fft_img";
    public final static String MASK = "mask";
    public final static String WRAPPED = "wrapped";
    public final static String UNWRAPPED = "unwrapped";
    public final static String MAGNITUDE = "magnitude";

    int size_x,size_y;
    Config config;

    boolean debug = false;

    private ImagePlus fft_img = null;
    private ImagePlus mask_img = null;
    private ImagePlus wrapped_img = null;
    private ImagePlus unwrapped_img = null;
    private ImagePlus magnitude_img = null;

    public static class Config{
        float z, lambda, dx, dy, dz;

        int small_contraint, big_contraint;

        boolean autofocus;

        boolean show_fft, show_mask, show_wrapped_phase, show_unwrapped_phase, show_magnitude;

        public Config(float lambda, float dx, float dy, float dz){
            /*
                z = 0;
                lambda = 0.638f;
                dx = 0.434f;
                dy = 0.434f;
                dz = 1.503f;
             */

            z = 0;
            this.lambda = lambda;
            this.dx = dx;
            this.dy = dy;
            this.dz = dz;

            small_contraint = 25;
            big_contraint = 4;

            autofocus = false;
            show_fft = false;
            show_mask = false;
            show_wrapped_phase = false;
            show_unwrapped_phase = true;
            show_magnitude = false;
        }

        public void update(float z,float lambda,float dx,float dy,float dz, int small_contraint, int big_contraint,
                           boolean autofocus, boolean show_fft, boolean show_mask, boolean show_wrapped_phase,
                           boolean show_unwrapped_phase, boolean show_magnitude){
            this.z = z;
            this.lambda = lambda;
            this.dx = dx;
            this.dy = dy;
            this.dz = dz;
            this.small_contraint = small_contraint;
            this.big_contraint = big_contraint;
            this.autofocus = autofocus;
            this.show_fft = show_fft;
            this.show_mask = show_mask;
            this.show_wrapped_phase = show_wrapped_phase;
            this.show_unwrapped_phase = show_unwrapped_phase;
            this.show_magnitude = show_magnitude;
        }

        public float getZ() {
            return z;
        }

        public void setZ(float z) {
            this.z = z;
        }

        public float getLambda() {
            return lambda;
        }

        public void setLambda(float lambda) {
            this.lambda = lambda;
        }

        public float getDx() {
            return dx;
        }

        public void setDx(float dx) {
            this.dx = dx;
        }

        public float getDy() {
            return dy;
        }

        public void setDy(float dy) {
            this.dy = dy;
        }

        public float getDz() {
            return dz;
        }

        public void setDz(float dz) {
            this.dz = dz;
        }

        public int getSmall_contraint() {
            return small_contraint;
        }

        public void setSmall_contraint(int small_contraint) {
            this.small_contraint = small_contraint;
        }

        public int getBig_contraint() {
            return big_contraint;
        }

        public void setBig_contraint(int big_contraint) {
            this.big_contraint = big_contraint;
        }

        public boolean isAutofocus() {
            return autofocus;
        }

        public void setAutofocus(boolean autofocus) {
            this.autofocus = autofocus;
        }

        public boolean isShow_fft() {
            return show_fft;
        }

        public void setShow_fft(boolean show_fft) {
            this.show_fft = show_fft;
        }

        public boolean isShow_mask() {
            return show_mask;
        }

        public void setShow_mask(boolean show_mask) {
            this.show_mask = show_mask;
        }

        public boolean isShow_wrapped_phase() {
            return show_wrapped_phase;
        }

        public void setShow_wrapped_phase(boolean show_wrapped_phase) {
            this.show_wrapped_phase = show_wrapped_phase;
        }

        public boolean isShow_unwrapped_phase() {
            return show_unwrapped_phase;
        }

        public void setShow_unwrapped_phase(boolean show_unwrapped_phase) {
            this.show_unwrapped_phase = show_unwrapped_phase;
        }

        public boolean isShow_magnitude() {
            return show_magnitude;
        }

        public void setShow_magnitude(boolean show_magnitude) {
            this.show_magnitude = show_magnitude;
        }
    }

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
        floatProcessor.add(1);
        floatProcessor.resetMinAndMax();
        floatProcessor.log();
        floatProcessor.resetMinAndMax();

        //Normalise
        floatProcessor.subtract(floatProcessor.getMin());
        floatProcessor.multiply(1/floatProcessor.getMax());

        FloatProcessor blurred = (FloatProcessor) floatProcessor.duplicate();

        //Gaussian blur
        blurred.blurGaussian(sigma); //TODO: Check gaussian size
        new ContrastEnhancer().equalize(blurred);

        //otsu threshold
        blurred.multiply(255.0);
        float level = new AutoThresholder().getThreshold(AutoThresholder.Method.Otsu, blurred.getHistogram());
        float step = level / 100.0f;

        int small_img = Math.round((((float) size_x)/small_const)*(((float) size_y)/small_const));
        int big_img = Math.round((((float) size_x)/large_const)*(((float) size_y)/large_const));

        ImageProcessor mask = null;
        ResultsTable rt = null;

        //TODO: Try with float instead of int mask
        while (iy != 3){
            if(debug)
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
                //TODO: Don't rely on IJ to handle exception
                IJ.handleException(new Exception("Failed to find threshold"));
                if(debug){
                    ImagePlus fft = new ImagePlus("FFT ", floatProcessor);
                    //uiService.show("FFT",ImageJFunctions.wrap(fft_img));
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
        //uiService.show("FFT",ImageJFunctions.wrap(fft_img));
        this.fft_img = fft;

        ImagePlus mask_img = new ImagePlus("Mask ", mask);
        //uiService.show("Mask",ImageJFunctions.wrap(mask_img));
        this.mask_img = mask_img;

        //ImagePlus angle_img = new ImagePlus("Cropped Mag", angleProcessor);
        //uiService.show("Cropped Mag",ImageJFunctions.wrap(angle_img));

        //return centered image
        return centered_complex_img;
    }

    private ImagePlus retrieve_magnitude(float[][] centered_fftImg){
        //Propagate fft_img
        FloatAngularSpectrum as = new FloatAngularSpectrum(size_x,size_y,config.getLambda(),
                config.getZ(),config.getDx(),config.getDy());
        as.diffract(centered_fftImg);

        //IFFT
        FloatFFT_2D fft2D = new FloatFFT_2D(size_x,size_y);
        ArrayUtils.complexShift(centered_fftImg);
        fft2D.complexInverse(centered_fftImg,false);

        //retrieve magnitude
        FloatProcessor magnitude_processor = new FloatProcessor(ArrayUtils.modulusSq(centered_fftImg));

        this.magnitude_img = new ImagePlus("Magnitude ", magnitude_processor);

        return magnitude_img;
    }

    private ImagePlus retrieve_phase(float[][] centered_fftImg){
        //TODO: multiply by 2/pi?

        //Propagate fft_img
        FloatAngularSpectrum as = new FloatAngularSpectrum(size_x,size_y,config.getLambda(),
                config.getZ(),config.getDx(),config.getDy());
        as.diffract(centered_fftImg);

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
                flat_phase[i][j] = (float) optimal[i*unwrapped_phase.getHeight()+j] * this.config.getDz();
            }
        }

        FloatProcessor flat_processor = new FloatProcessor(flat_phase);
        ImagePlus flat_img = new ImagePlus("Flat Phase", flat_processor);

        this.unwrapped_img = flat_img;
        this.wrapped_img = wrapped_phase;

        return flat_img;
    }

    public DHMReconstructor(Config config){
        this.config = config;
    }

    /**
     *
     * @param img - Source hologram
     *
     * @return Map of requested images by name
     */
    public Map<String,ImagePlus> reconstruct(ImagePlus img){

        ImageProcessor ip = img.getProcessor();

        size_x = ip.getWidth();
        size_y = ip.getHeight();

        float[][] cropped = cropOrders(ip,config.getSmall_contraint(),config.getBig_contraint());
        retrieve_phase(cropped);
        retrieve_magnitude(cropped);

        Map<String,ImagePlus> imageMap = new HashMap<>();

        if(config.isShow_fft()){
            imageMap.put(FFT, fft_img);
        }

        if(config.isShow_mask()){
            imageMap.put(MASK,mask_img);
        }

        if(config.isShow_wrapped_phase()){
            imageMap.put(WRAPPED,wrapped_img);
        }

        if(config.isShow_unwrapped_phase()){
            imageMap.put(UNWRAPPED,unwrapped_img);
        }

        if(config.isShow_magnitude()){
            imageMap.put(MAGNITUDE,magnitude_img);
        }

        return imageMap;
    }

}
