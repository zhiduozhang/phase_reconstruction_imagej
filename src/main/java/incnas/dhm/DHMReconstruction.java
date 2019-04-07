package incnas.dhm;

import ij.ImagePlus;
import ij.gui.GenericDialog;
import net.imagej.Dataset;
import net.imagej.ImageJ;
import net.imagej.ImgPlus;
import net.imagej.ops.OpService;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.RealType;
import org.scijava.command.Command;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;
import org.scijava.plugin.PluginService;
import org.scijava.ui.UIService;

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
        GenericDialog gd = new GenericDialog("DHM Reconstruction");
        //gd.addMessage("");
        gd.addNumericField("Wavelength (um)", 0.635,3);
        gd.addNumericField("dx", 0.434, 3);
        gd.addNumericField("dy", 0.434, 3);
        gd.addNumericField("dz", 1.503, 3);
        gd.addNumericField("Small Constraint", 25, 0);
        gd.addNumericField("Large Constraint", 4, 0);

        gd.addCheckbox("Enable Autofocus (pending)",true);
        gd.addNumericField("Focus Distance (um)", 0, 3);

        gd.addMessage("Output Images:");
        gd.addCheckbox("FFT",false);
        gd.addCheckbox("Mask + Cropped Region", false);
        gd.addCheckbox("Wrapped Phase", false);
        gd.addCheckbox("Unwrapped Phase", true);
        gd.addCheckbox("Magnitude", false);

        gd.showDialog();

        if (gd.wasCanceled()) return;

        float lambda= (float) gd.getNextNumber();
        float dx = (float) gd.getNextNumber();
        float dy = (float) gd.getNextNumber();
        float dz = (float) gd.getNextNumber();

        int small_contraint = (int) gd.getNextNumber();
        int big_contraint = (int) gd.getNextNumber();

        boolean autofocus = gd.getNextBoolean();
        float z = (float) gd.getNextNumber();

        boolean show_fft = gd.getNextBoolean();
        boolean show_mask = gd.getNextBoolean();
        boolean show_wrapped_phase = gd.getNextBoolean();
        boolean show_unwrapped_phase = gd.getNextBoolean();
        boolean show_magnitude = gd.getNextBoolean();

        //IJ.newImage(title, "8-bit", width, height, 1);

        DHMReconstructor.Config config = new DHMReconstructor.Config(lambda,dx,dy,dz);
        config.update(z,lambda,dx,dy,dz,small_contraint,big_contraint,autofocus,show_fft,show_mask,show_wrapped_phase,
            show_unwrapped_phase,show_magnitude);

        final ImgPlus image = currentData.getImgPlus();
        ImagePlus img = ImageJFunctions.wrap(image,"Input Image");

        DHMReconstructor reconstructor = new DHMReconstructor(config);
        Map<String, ImagePlus> images = reconstructor.reconstruct(img);

        if(show_fft){
            uiService.show("FFT",ImageJFunctions.wrap(images.get(DHMReconstructor.FFT)));
        }

        if(show_mask){
            uiService.show("Cropping Mask",ImageJFunctions.wrap(images.get(DHMReconstructor.MASK)));
        }

        if(show_wrapped_phase){
            uiService.show("Wrapped Phase",ImageJFunctions.wrap(images.get(DHMReconstructor.WRAPPED)));
        }

        if(show_unwrapped_phase){
            uiService.show("Unwrapped Phase",ImageJFunctions.wrap(images.get(DHMReconstructor.UNWRAPPED)));
        }

        if(show_magnitude){
            uiService.show("Magnitude",ImageJFunctions.wrap(images.get(DHMReconstructor.MAGNITUDE)));
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
