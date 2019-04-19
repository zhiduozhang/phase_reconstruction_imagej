package incnas.dhm;

import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.plugin.SurfacePlotter;
import ij.plugin.filter.Duplicater;
import net.imagej.Dataset;
import net.imagej.ImageJ;
import net.imagej.ImgPlus;
import net.imagej.ops.OpService;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.img.imageplus.ImagePlusImg;
import net.imglib2.type.numeric.RealType;
import org.scijava.command.Command;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;
import org.scijava.plugin.PluginService;
import org.scijava.ui.UIService;

import java.awt.*;
import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

//TODO: Use more efficient array elementwise operations?

@Plugin(type = Command.class, menuPath = "Plugins>DHM Reconstruction")
public class DHMReconstruction<T extends RealType<T>> implements Command {
    @Parameter
    private Dataset currentData;

    @Parameter
    private PluginService pluginService;

    @Parameter
    private UIService uiService;

    @Parameter
    private OpService opService;

    @Override
    public void run() {
        //Load Prefs
        float wavelength = (float) Prefs.get("dhm.reconstruction.wavelength", 0.635);
        float dx = (float) Prefs.get("dhm.reconstruction.dx", 0.434);
        float dy = (float) Prefs.get("dhm.reconstruction.dy", 0.434);
        float dz = (float) Prefs.get("dhm.reconstruction.dz", 1.503);
        int small_constraint = Prefs.getInt("dhm.reconstruction.small_constraint", 25);
        int large_constraint = Prefs.getInt("dhm.reconstruction.large_constraint", 4);

        boolean autofocus = Prefs.get("dhm.recosntruction.autofocus", false);
        float focus_distance = (float) Prefs.get("dhm.reconstruction.focus_distance", 0);

        boolean show_fft = Prefs.get("dhm.reconstruction.show_fft",false);
        boolean show_mask = Prefs.get("dhm.reconstruction.show_mask",false);
        boolean show_wrapped_phase = Prefs.get("dhm.reconstruction.show_wrapped",false);
        boolean show_unwrapped_phase = Prefs.get("dhm.reconstruction.show_unwrapped",true);
        boolean show_magnitude = Prefs.get("dhm.reconstruction.show_magnitude",false);
        boolean show_3d_phase = Prefs.get("dhm.reconstruction.show_3d_phase",false);

        GenericDialog gd = new GenericDialog("DHM Reconstruction");
        //gd.addMessage("");
        gd.addNumericField("Wavelength (um)", wavelength,3);
        gd.addNumericField("dx", dx, 3);
        gd.addNumericField("dy", dy, 3);
        gd.addNumericField("dz", dz, 3);
        gd.addNumericField("Small Constraint", small_constraint, 0);
        gd.addNumericField("Large Constraint", large_constraint, 0);

        gd.addCheckbox("Enable Autofocus (pending)",autofocus);
        gd.addNumericField("Focus Distance (um)", focus_distance, 3);

        gd.addMessage("Output Images:");
        gd.addCheckbox("FFT",show_fft);
        gd.addCheckbox("Mask + Cropped Region", show_mask);
        gd.addCheckbox("Wrapped Phase", show_wrapped_phase);
        gd.addCheckbox("Unwrapped Phase", show_unwrapped_phase);
        gd.addCheckbox("Show 3D Phase", show_3d_phase);
        gd.addCheckbox("Magnitude", show_magnitude);

        gd.showDialog();

        if (gd.wasCanceled()) return;

        wavelength= (float) gd.getNextNumber();
        dx = (float) gd.getNextNumber();
        dy = (float) gd.getNextNumber();
        dz = (float) gd.getNextNumber();

        small_constraint = (int) gd.getNextNumber();
        large_constraint = (int) gd.getNextNumber();

        autofocus = gd.getNextBoolean();
        focus_distance = (float) gd.getNextNumber();

        show_fft = gd.getNextBoolean();
        show_mask = gd.getNextBoolean();
        show_wrapped_phase = gd.getNextBoolean();
        show_unwrapped_phase = gd.getNextBoolean();
        show_3d_phase = gd.getNextBoolean();
        show_magnitude = gd.getNextBoolean();

        //Save Prefs
        Prefs.set("dhm.reconstruction.wavelength", wavelength);
        Prefs.set("dhm.reconstruction.dx", dx);
        Prefs.set("dhm.reconstruction.dy", dy);
        Prefs.set("dhm.reconstruction.dz", dz);
        Prefs.set("dhm.reconstruction.small_constraint", small_constraint);
        Prefs.set("dhm.reconstruction.large_constraint", large_constraint);

        Prefs.set("dhm.recosntruction.autofocus", autofocus);
        Prefs.set("dhm.reconstruction.focus_distance", focus_distance);

        Prefs.set("dhm.reconstruction.show_fft",show_fft);
        Prefs.set("dhm.reconstruction.show_mask",show_mask);
        Prefs.set("dhm.reconstruction.show_wrapped",show_wrapped_phase);
        Prefs.set("dhm.reconstruction.show_unwrapped",show_unwrapped_phase);
        Prefs.set("dhm.reconstruction.show_magnitude",show_magnitude);
        Prefs.set("dhm.reconstruction.show_3d_phase",show_3d_phase);

        Prefs.savePreferences();

        //Reconstruct Hologram
        boolean black_background = Prefs.blackBackground; //See https://imagej.net/Troubleshooting#The_same_plugin_gives_different_results_on_different_machines.21 and https://imagej.nih.gov/ij/source/ij/plugin/filter/Binary.java

        //IJ.newImage(title, "8-bit", width, height, 1);

        DHMReconstructor.Config config = new DHMReconstructor.Config(wavelength,dx,dy,dz);
        config.update(focus_distance,wavelength,dx,dy,dz,small_constraint,large_constraint,autofocus,show_fft,show_mask,show_wrapped_phase,
            show_unwrapped_phase,show_magnitude,black_background);

        final ImgPlus image = currentData.getImgPlus();
        ImagePlus img = ImageJFunctions.wrap(image,"Input Image");
        if(img.isStack()){
            ImageStack fftStack = new ImageStack(img.getWidth(), img.getHeight());
            ImageStack maskStack = new ImageStack(img.getWidth(), img.getHeight());
            ImageStack wrappedStack = new ImageStack(img.getWidth(), img.getHeight());
            ImageStack unwrappedStack = new ImageStack(img.getWidth(), img.getHeight());
            ImageStack magnitudeStack = new ImageStack(img.getWidth(), img.getHeight());

            for(int i=1; i<=img.getImageStack().getSize(); i++){
                ImagePlus input = new ImagePlus("",img.getImageStack().getProcessor(i));

                DHMReconstructor reconstructor = new DHMReconstructor(config);
                Map<String, ImagePlus> images = reconstructor.reconstruct(input);

                if(show_fft && images.containsKey(DHMReconstructor.FFT)){
                    fftStack.addSlice(images.get(DHMReconstructor.FFT).getProcessor());
                }

                if(show_mask && images.containsKey(DHMReconstructor.MASK)){
                    maskStack.addSlice(images.get(DHMReconstructor.MASK).getProcessor());
                }

                if(show_wrapped_phase && images.containsKey(DHMReconstructor.WRAPPED)){
                    wrappedStack.addSlice(images.get(DHMReconstructor.WRAPPED).getProcessor());
                }

                if(images.containsKey(DHMReconstructor.UNWRAPPED)){
                    unwrappedStack.addSlice(images.get(DHMReconstructor.UNWRAPPED).getProcessor());
                }

                if(show_magnitude && images.containsKey(DHMReconstructor.MAGNITUDE)){
                    magnitudeStack.addSlice(images.get(DHMReconstructor.MAGNITUDE).getProcessor());
                }
            }

            if(show_unwrapped_phase && unwrappedStack.getSize() > 0){
                //uiService.show("Unwrapped Phase",ImageJFunctions.wrap(new ImagePlus("Unwrapped Phase", unwrappedStack)));
                new ImagePlus("Unwrapped Phase", unwrappedStack).show();
            }

            if(show_3d_phase){
                if(!show_unwrapped_phase && unwrappedStack.getSize() > 0){
                    //uiService.show("Unwrapped Phase", ImageJFunctions.wrap(new ImagePlus("Unwrapped Phase", unwrappedStack)));
                    new ImagePlus("Unwrapped Phase", unwrappedStack).show();
                }

                new SurfacePlotter().run("");
            }

            if(show_fft && fftStack.getSize()>0){
                //uiService.show("FFT",ImageJFunctions.wrap(new ImagePlus("FFT", fftStack)));
                new ImagePlus("FFT", fftStack).show();
            }

            if(show_mask && maskStack.getSize()>0){
                //uiService.show("Cropping Mask",ImageJFunctions.wrap(new ImagePlus("Mask", maskStack)));
                uiService.show(new ImagePlus("Mask", maskStack));
            }

            if(show_wrapped_phase && wrappedStack.getSize() > 0){
                //uiService.show("Wrapped Phase",ImageJFunctions.wrap(new ImagePlus("Wrapped Phase", wrappedStack)));
                uiService.show(new ImagePlus("Wrapped Phase", wrappedStack));
            }

            if(show_magnitude && magnitudeStack.getSize() > 0){
                //uiService.show("Magnitude",ImageJFunctions.wrap(new ImagePlus("Magnitude",magnitudeStack)));
                uiService.show(new ImagePlus("Magnitude",magnitudeStack));
            }
        }else{
            DHMReconstructor reconstructor = new DHMReconstructor(config);
            Map<String, ImagePlus> images = reconstructor.reconstruct(img);

            if(images == null){

            }

            if(show_fft && images.containsKey(DHMReconstructor.FFT)){
                uiService.show("FFT",ImageJFunctions.wrap(images.get(DHMReconstructor.FFT)));
            }

            if(show_mask && images.containsKey(DHMReconstructor.MASK)){
                uiService.show("Cropping Mask",ImageJFunctions.wrap(images.get(DHMReconstructor.MASK)));
            }

            if(show_wrapped_phase && images.containsKey(DHMReconstructor.WRAPPED)){
                uiService.show("Wrapped Phase",ImageJFunctions.wrap(images.get(DHMReconstructor.WRAPPED)));
            }

            if(show_unwrapped_phase && images.containsKey(DHMReconstructor.UNWRAPPED)){
                uiService.show("Unwrapped Phase",ImageJFunctions.wrap(images.get(DHMReconstructor.UNWRAPPED)));
            }

            if(show_3d_phase){
                if(!show_unwrapped_phase){
                    uiService.show("Unwrapped Phase",ImageJFunctions.wrap(images.get(DHMReconstructor.UNWRAPPED)));
                }

                new SurfacePlotter().run("");
            }

            if(show_magnitude && images.containsKey(DHMReconstructor.MAGNITUDE)){
                uiService.show("Magnitude",ImageJFunctions.wrap(images.get(DHMReconstructor.MAGNITUDE)));
            }
        }
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
