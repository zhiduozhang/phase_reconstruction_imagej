package incnas.dhm;

/**
 * Original matlab code by Xuefei He (Australian National University).
 *
 * @author Zhiduo Zhang
 */

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
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
import incnas.dhm.utils.CurveProblem;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresBuilder;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresOptimizer;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem;
import org.apache.commons.math3.fitting.leastsquares.LevenbergMarquardtOptimizer;
import org.apache.commons.math3.linear.RealVector;
import org.bytedeco.javacpp.opencv_core;
import org.bytedeco.javacpp.opencv_phase_unwrapping;
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

    static boolean debug = false;

    private ImagePlus fft_img = null;
    private ImagePlus mask_img = null;
    private ImagePlus wrapped_img = null;
    private ImagePlus unwrapped_img = null;
    private ImagePlus magnitude_img = null;

    private static void log(String log){
        //TODO: Display in ImageJ Log
        if(debug) {
            System.out.println(log);
        }
    }

    public static class Config{
        float z, lambda, dx, dy, dz;

        int small_contraint, big_contraint;

        boolean autofocus, black_background;

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
            black_background = true;
        }

        public void update(float z,float lambda,float dx,float dy,float dz, int small_contraint, int big_contraint,
                           boolean autofocus, boolean show_fft, boolean show_mask, boolean show_wrapped_phase,
                           boolean show_unwrapped_phase, boolean show_magnitude, boolean black_background){
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
            this.black_background = black_background;
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

        int complex_size_y = this.size_y*2;

        float[][] f_img = ip.getFloatArray();

        //Spread to allow for complex data by doubling array size
        float[][] img = new float[this.size_x][complex_size_y];

        for (int i = 0 ; i < this.size_x; i++){
            for (int j=0; j<this.size_y; j++){
                img[i][j*2] = f_img[i][j];
            }
        }

        //FFT
        float[][] imgFFT = img.clone();
        FloatFFT_2D fft2D = new FloatFFT_2D(this.size_x,this.size_y);

        fft2D.complexForward(imgFFT);
        ArrayUtils.complexShift(imgFFT);

        //FFT Log Magnitude
        float[][] magImgFFt = ArrayUtils.modulus(imgFFT);
        FloatProcessor floatProcessor = new FloatProcessor(magImgFFt);
        floatProcessor.add(1);
        floatProcessor.resetMinAndMax();
        floatProcessor.log();
        floatProcessor.resetMinAndMax();

        //FloatProcessor unorm = (FloatProcessor) floatProcessor.duplicate();

        //Normalise
        floatProcessor.subtract(floatProcessor.getMin());
        floatProcessor.resetMinAndMax();
        floatProcessor.multiply(1.0f/floatProcessor.getMax());
        floatProcessor.resetMinAndMax();

        /*for(float[] row : floatProcessor.getFloatArray()){
            System.out.println(Arrays.toString(row));
        }*/

        FloatProcessor blurred = (FloatProcessor) floatProcessor.duplicate();
        blurred.multiply(255.0);
        blurred.resetMinAndMax();

        //Gaussian blur
        blurred.blurGaussian(sigma); //TODO: Check gaussian size
        new ContrastEnhancer().equalize(blurred);

        //otsu threshold

        float level = new AutoThresholder().getThreshold(AutoThresholder.Method.Otsu, blurred.getHistogram());
        float step = level / 100.0f;

        int small_img = Math.round((((float) this.size_x)/small_const)*(((float) this.size_y)/small_const));
        int big_img = Math.round((((float) this.size_x)/large_const)*(((float) this.size_y)/large_const));

        log("Small Img: " + small_img);
        log("Big Img: " + big_img);

        ImageProcessor mask = null;
        ResultsTable rt = null;

        float[][] float_mask = blurred.getFloatArray();
        log("Max: " + blurred.getMax());

        ImageStack maskStack = new ImageStack(this.size_x,this.size_y);

        //TODO: Try with float instead of int mask
        while (iy != 3){
            log("Theshold:" + (level + inc));

            byte[] int_mask = new byte[this.size_x*this.size_y];
            short max = 0;

            for(int i=0; i<this.size_x;i++){
                for(int j=0; j<this.size_y;j++){
                    //log("Float: " + float_mask[i][j]);
                    //log("Byte: " + (byte) Math.round(float_mask[i][j]));
                    int_mask[(j*this.size_x)+i] = (byte) Math.round(float_mask[i][j]);
                    max = (short) Math.max(max, Math.round(float_mask[i][j]));
                }
            }

            //mask = new ShortProcessor(size_x,size_y,int_mask, new java.awt.image.ColorModel(DataBuffer.TYPE_USHORT));
            mask = new ByteProcessor(this.size_x,this.size_y,int_mask);
            //mask.setThreshold(140,255,ImageProcessor.RED_LUT);
            mask.threshold((int) (level+inc));
            if(!config.black_background) {
                mask.invert();
            }

            ByteProcessor cleanMask = new ByteProcessor(this.size_x,this.size_y,int_mask);
            maskStack.addSlice(mask);

            Img ImageMask = ImageJFunctions.wrap(new ImagePlus("Mask",cleanMask));

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

            log(rt.toString());
            log("----------------------");
            log("Size:" + rt.size());
            log("Area:" + Arrays.toString(rt.getColumn(ResultsTable.AREA)));
            log("X Centroid: " + Arrays.toString(rt.getColumn(ResultsTable.X_CENTROID)));
            log("Y Centroid: " + Arrays.toString(rt.getColumn(ResultsTable.Y_CENTROID)));
            log("Height: " + Arrays.toString(rt.getColumn(ResultsTable.ROI_HEIGHT)));
            log("Width: " + Arrays.toString(rt.getColumn(ResultsTable.ROI_WIDTH)));
            log("ROI X: " + Arrays.toString(rt.getColumn(ResultsTable.ROI_X)));
            log("ROI Y: " + Arrays.toString(rt.getColumn(ResultsTable.ROI_Y)));

            iy = rt.size();

            if (rt.size() == 3) {

                for (float x : rt.getColumn(ResultsTable.AREA)) {
                    if(x > big_img){
                        log("Discarding threshold: Area too big");
                        iy = 1;
                        break;
                    } else {
                        mask = pa.getOutputImage().getProcessor();
                    }
                }
            }

            inc += step;

            if (inc+level > max){
                //TODO: Don't rely on IJ to handle exception
                //this.fft_img = new ImagePlus("FFT ", cleanMask);
                this.fft_img = new ImagePlus("FFT ", maskStack);
                IJ.handleException(new Exception("Failed to find threshold. Max reached " + (inc+level)));

                return null;
            }
        }

        log("Found threshold!");


        int min_x = (int) rt.getColumn(ResultsTable.X_CENTROID)[0];
        int min_index = 0;
        for (int i=1; i<3; i++){
            if (min_x > rt.getColumn(ResultsTable.X_CENTROID)[i]){
                min_x = (int) rt.getColumn(ResultsTable.X_CENTROID)[i];
                min_index = i;
            }
        }

        log("Min Index: " + min_index);

        //Crop to region
        int crop_width = (int) rt.getColumn(ResultsTable.ROI_WIDTH)[min_index];
        int crop_height_2 = (int) rt.getColumn(ResultsTable.ROI_HEIGHT)[min_index]*2;

        int crop_x = (int) rt.getColumn(ResultsTable.ROI_X)[min_index];
        int crop_y = (int) rt.getColumn(ResultsTable.ROI_Y)[min_index];
        float[][] cropped_area = new float[crop_width][crop_height_2];

        //Retrieve complex masked
        for(int i=0; i<crop_width; i++){
            for (int j=0; j<crop_height_2; j++){
                cropped_area[i][j] = imgFFT[i+crop_x][j+crop_y*2]; // * mask_vals[i+crop_x][j+crop_y];
            }
        }

        //Move object into center of image
        float[][] centered_complex_img = new float[this.size_x][this.size_y*2];
        float[][] croppedMag = ArrayUtils.modulus(cropped_area);

        //Find peak
        int peak_x = 0;
        int peak_y = 0;
        float max = 0;
        for(int i=0; i<crop_width; i++){
            for(int j=0; j<crop_height_2/2; j++){
                if(croppedMag[i][j] > max){
                    max = croppedMag[i][j];
                    peak_x = i;
                    peak_y = j;
                }
            }
        }

        //Recreate cropped area centered around the peak
        int bound_width = Math.max(peak_x,Math.abs(crop_width-peak_x))*2;
        int bound_height_2 = Math.max(peak_y,Math.abs((crop_height_2/2)-(peak_y)))*4;

        int bound_x = (crop_x+peak_x) - bound_width/2;
        int bound_y = (crop_y+peak_y)*2 - bound_height_2/2;

        float[][] new_crop = new float[bound_width][bound_height_2];
        for(int i=0; i<bound_width; i++){
            for (int j=0; j<bound_height_2; j++){
                new_crop[i][j] = imgFFT[i+bound_x][j+bound_y]; // * mask_vals[i+crop_x][j+crop_y];
            }
        }

        int new_displacement_x = (int) Math.round((this.size_x-bound_width)/2.0);
        int new_displacement_y = (int) Math.round((complex_size_y-bound_height_2)/2.0);

        for(int i=0; i<bound_width; i++){
            for(int j=0; j<bound_height_2; j++){
                centered_complex_img[i+new_displacement_x][j+new_displacement_y] = new_crop[i][j];
            }
        }

        //floatProcessor.drawRect(crop_x,crop_y,crop_width,crop_height_2/2);
        ImagePlus fft = new ImagePlus("FFT ", floatProcessor);
        this.fft_img = fft;

        ImagePlus mask_img = new ImagePlus("Mask ", mask);
        this.mask_img = mask_img;

        float[][] abs = ArrayUtils.modulus(centered_complex_img);
        FloatProcessor abs_fp = new FloatProcessor(abs);
        abs_fp.add(1);
        abs_fp.log();
        ImagePlus abs_img = new ImagePlus("Abs FFT ", abs_fp);

        this.mask_img = abs_img;

        //return centered image
        return centered_complex_img;
    }

    private ImagePlus retrieve_magnitude(float[][] centered_fftImg){
        for(int i=0; i<this.size_x; i++){
            for(int j=0; j<this.size_y;j++){
                centered_fftImg[i][j] = (float) (centered_fftImg[i][j]*Math.PI/2);
            }
        }

        //Propagate fft_img
        FloatAngularSpectrum as = new FloatAngularSpectrum(size_x,size_y,config.getLambda(),
                config.getZ(),config.getDx(),config.getDy());
        as.diffract_fft(centered_fftImg);

        //IFFT
        FloatFFT_2D fft2D = new FloatFFT_2D(size_x,size_y);
        ArrayUtils.complexShift(centered_fftImg);
        fft2D.complexInverse(centered_fftImg,true);

        //retrieve magnitude
        FloatProcessor magnitude_processor = new FloatProcessor(ArrayUtils.modulus(centered_fftImg));

        this.magnitude_img = new ImagePlus("Magnitude ", magnitude_processor);

        return magnitude_img;
    }

    private ImagePlus retrieve_phase(float[][] centered_fftImg){
        for(int i=0; i<this.size_x; i++){
            for(int j=0; j<this.size_y;j++){
                centered_fftImg[i][j] = (float) (centered_fftImg[i][j]*Math.PI/2);
            }
        }

        //Propagate fft_img
        FloatAngularSpectrum as = new FloatAngularSpectrum(size_x,size_y,config.getLambda(),
                config.getZ(),config.getDx(),config.getDy());
        as.diffract_fft(centered_fftImg);

        //IFFT
        FloatFFT_2D fft2D = new FloatFFT_2D(size_x,size_y);
        ArrayUtils.complexShift(centered_fftImg);
        fft2D.complexInverse(centered_fftImg,true);

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
        double[] flat = new double[unwrapped_phase.getWidth()*unwrapped_phase.getHeight()];
        for(int i=0; i<this.size_x; i++) {
            for (int j = 0; j < this.size_y; j++) {
                flat[i * this.size_y + j] = unwrapped_phase.getProcessor().getf(i,j);
            }
        }

        CurveProblem curveProblem = new CurveProblem((FloatProcessor) unwrapped_phase.getProcessor());
        LeastSquaresProblem problem = new LeastSquaresBuilder()
                .start(new double[]{0.5,0.5,0.5,0.5,0.5,0.5})
                .model(curveProblem)
                .target(flat)
                .lazyEvaluation(false)
                .maxEvaluations(500)
                .maxIterations(500)
                .build();

        LeastSquaresOptimizer.Optimum optimum = new LevenbergMarquardtOptimizer().optimize(problem);
        RealVector opt_params = optimum.getPoint();

        double[] optimal = curveProblem.value(opt_params).getFirst().toArray();

        log(Arrays.toString(opt_params.toArray()));

        float[][] unwrapped = unwrapped_phase.getProcessor().getFloatArray();
        float[][] flat_phase = new float[unwrapped_phase.getWidth()][unwrapped_phase.getHeight()];
        for(int i=0; i<unwrapped_phase.getWidth(); i++){
            for (int j=0; j<unwrapped_phase.getHeight(); j++){
                flat_phase[i][j] = (float) Math.abs(optimal[i*unwrapped_phase.getHeight()+j]-unwrapped[i][j]) * this.config.getDz();
                flat_phase[i][j] = Math.max(flat_phase[i][j],0);
            }
        }

        FloatProcessor flat_processor = new FloatProcessor(flat_phase);
        flat_processor.blurGaussian(3);
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
     * @return Map of requested images by name or null if we cannot find a threshold.
     */
    public Map<String,ImagePlus> reconstruct(ImagePlus img){
        Map<String,ImagePlus> imageMap = new HashMap<>();

        ImageProcessor ip = img.getProcessor();

        size_x = ip.getWidth();
        size_y = ip.getHeight();

        float[][] cropped = cropOrders(ip,config.getSmall_contraint(),config.getBig_contraint());

        float[][] mag = new float[this.size_x][this.size_y*2];
        for(int i=0; i<this.size_x; i++){
            if (this.size_y >= 0) System.arraycopy(cropped[i], 0, mag[i], 0, this.size_y);
        }

        if (cropped != null){
            retrieve_phase(cropped);
            retrieve_magnitude(mag);

            if(config.isShow_wrapped_phase()){
                imageMap.put(WRAPPED,wrapped_img);
            }

            imageMap.put(UNWRAPPED,unwrapped_img);

            if(config.isShow_magnitude()){
                imageMap.put(MAGNITUDE,magnitude_img);
            }

            if(config.isShow_mask()){
                imageMap.put(MASK,mask_img);
            }
        }

        if(config.isShow_fft()){
            imageMap.put(FFT, fft_img);
        }


        return imageMap;
    }

}
