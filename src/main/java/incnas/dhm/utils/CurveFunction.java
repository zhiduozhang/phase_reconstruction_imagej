package incnas.dhm.utils;

import ij.process.FloatProcessor;
import org.ddogleg.optimization.functions.FunctionNtoM;

public class CurveFunction implements FunctionNtoM {
    protected int imgSize,size_x,size_y;
    protected float[][] data;

    public CurveFunction(FloatProcessor processor){
        this.size_x = processor.getWidth();
        this.size_y = processor.getHeight();

        this.imgSize = size_x*size_y;

        this.data = processor.getFloatArray();
    }

    @Override
    public void process(double[] input, double[] output){
        for(int i=0; i<this.size_x; i++){
            for (int j=0; j<this.size_y; j++){
                int index = i*this.size_x+j;
                double estimate = input[0] + input[1]*(i+1) + input[2]*(j+1) + input[3]*(i+1)*(j+1) + input[4]*(i+1)*(i+1) + input[5]*(j+1)*(j+1);
                output[index] = this.data[i][j] - estimate;
            }
        }
    }

    @Override
    public int getNumOfInputsN() {
        return 6;
    }

    @Override
    public int getNumOfOutputsM() {
        return this.imgSize;
    }
}
