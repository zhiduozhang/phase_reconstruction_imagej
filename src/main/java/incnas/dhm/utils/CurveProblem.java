package incnas.dhm.utils;

import ij.process.FloatProcessor;
import org.apache.commons.math3.fitting.leastsquares.MultivariateJacobianFunction;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.Pair;

public class CurveProblem implements MultivariateJacobianFunction {
    protected int imgSize,size_x,size_y;
    protected float[][] data;

    public CurveProblem(FloatProcessor processor){
        this.size_x = processor.getWidth();
        this.size_y = processor.getHeight();

        this.imgSize = size_x*size_y;

        this.data = processor.getFloatArray();
    }

    @Override
    public Pair<RealVector, RealMatrix> value(RealVector realVector) {
        double[] input = realVector.toArray();
        RealVector value = new ArrayRealVector(this.size_x*this.size_y);
        RealMatrix jacobian = new Array2DRowRealMatrix(this.size_x*this.size_y, 6);

        for(int i=0; i<this.size_x; i++){
            for (int j=0; j<this.size_y; j++){
                int index = i*this.size_x+j;
                double estimate = input[0] + input[1]*(i+1) + input[2]*(j+1) + input[3]*(i+1)*(j+1) + input[4]*(i+1)*(i+1) + input[5]*(j+1)*(j+1);
                value.setEntry(index, estimate);

                jacobian.setEntry(index,0,1);
                jacobian.setEntry(index,1,i+1);
                jacobian.setEntry(index,2,j+1);
                jacobian.setEntry(index,3,(i+1)*(j+1));
                jacobian.setEntry(index,4,(i+1)*(i+1));
                jacobian.setEntry(index,5,(j+1)*(j+1));

            }
        }

        return new Pair<>(value, jacobian);
    }
}
