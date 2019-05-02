package incnas.dhm;


import java.io.BufferedReader;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.LinkedList;

/**
 * Port of C++ Implementation of  "Fast two-dimensional phase-unwrapping algorithm based on sorting by
 * reliability following a noncontinuous path"
 * by  Miguel Arevallilo Herra´ez, David R. Burton, Michael J. Lalor, and Munther A. Gdeisat
 * published in the Applied Optics, Vol. 41, No. 35, pp. 7437, 2002.
 */
public class MigualPhaseUnwrapper {

    private class Pixel{
        int increment;
        float value;
        float reliability;
        int group;

        Pixel(float value){
            increment = 0;
            this.value = value;
            group = -1;
            reliability = (float) (9999999+Math.random()*10000);
        }
    }

    private class Edge implements Comparable<Edge>{
        float reliability;
        int increment;
        int pixPosA;
        int pixPosB;

        public Edge() {
            pixPosA = -1;
            pixPosB = -1;
        }

        @Override
        public int compareTo(Edge o) {
            return Float.compare(this.reliability, o.reliability);
        }
    }

    private float wrap(float value){
        if(value > Math.PI) {
            return (float) (value - (Math.PI*2));
        } else if(value < -Math.PI){
            return (float) (value + (Math.PI*2));
        } else {
            return value;
        }
    }

    private int find_wrap(float a, float b){
        if(a-b > Math.PI){
            return -1;
        } else if (a-b < -Math.PI){
            return 1;
        } else {
            return 0;
        }
    }

    private void calculateReliability(Pixel[] pixArray, int width, int height){
        //float H, V, D1, D2;
        int index = width+1;
        int wip = width+1;
        for (int i=1; i<height-1; i++){
            for(int j=1; j<width-1;j++){
                float H = wrap(pixArray[wip-1].value-pixArray[wip].value) - wrap(pixArray[wip].value-pixArray[wip+1].value);
                float V = wrap(pixArray[wip-width].value-pixArray[wip].value) - wrap(pixArray[wip].value-pixArray[wip+width].value);
                float D1 = wrap(pixArray[wip-width-1].value-pixArray[wip].value) - wrap(pixArray[wip].value-pixArray[wip+width+1].value);
                float D2 = wrap(pixArray[wip-width+1].value-pixArray[wip].value) - wrap(pixArray[wip].value-pixArray[wip+width-1].value);

                pixArray[index].reliability = H*H + V*V + D1*D1 + D2*D2;
                index++;
                wip++;
            }
            index+=2;
            wip+=2;
        }
    }

    private void horizontalEDGEs(Pixel[] pixArray, Edge[] edgeArray, int width, int height){
        int pixCounter = 0;
        int edgeCounter = 0;
        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width - 1; j++) {
                edgeArray[edgeCounter] = new Edge();
                edgeArray[edgeCounter].reliability = pixArray[pixCounter].reliability + pixArray[pixCounter + 1].reliability;
                edgeArray[edgeCounter].increment = find_wrap(pixArray[pixCounter].value, pixArray[pixCounter + 1].value);
                edgeArray[edgeCounter].pixPosA = pixCounter;
                edgeArray[edgeCounter].pixPosB = pixCounter + 1;

                pixCounter++;
                edgeCounter++;
            }
            pixCounter++;
        }
    }

    private void verticalEDGEs(Pixel[] pixArray, Edge[] edgeArray, int width, int height) {
        int edgeCounter = height*(width-1);
        int pixCounter = 0;
        for(int i=0; i<(height-1); i++){
            for(int j=0; j<width;j++){
                edgeArray[edgeCounter] = new Edge();
                edgeArray[edgeCounter].reliability = pixArray[pixCounter].reliability + pixArray[pixCounter+width].reliability;
                edgeArray[edgeCounter].increment = find_wrap(pixArray[pixCounter].value, pixArray[pixCounter+width].value);
                edgeArray[edgeCounter].pixPosA = pixCounter;
                edgeArray[edgeCounter].pixPosB = pixCounter+width;

                pixCounter++;
                edgeCounter++;
            }
        }
    }

    private void gatherPIXELs(Pixel[] pixArray, Edge[] edgeArray, int width, int height){
        int num_edges = (width)*(height-1) + (width-1)*(height);
        LinkedList<Integer>[] pixGroups = new LinkedList[width*height];
        int num_groups = 0;

        for(int i=0; i< num_edges; i++){
            Edge edge = edgeArray[i];

            int pixPosA = edge.pixPosA;
            int pixPosB = edge.pixPosB;

            Pixel pixelA = pixArray[pixPosA];
            Pixel pixelB = pixArray[pixPosB];

            //PIXEL 1 and PIXEL 2 belong to different groups
            //initially each pixel is a group by it self and one pixel can construct a group
            //no else or else if to this if
            if(pixelA.group==pixelB.group && pixelA.group != -1) {
                //System.out.println(edge.reliability + " - " + "Same group");
                continue;
            }

            //System.out.print("Pix A: " + pixelA.value + " - Pix B: " + pixelB.value + " - ");

            //PIXEL 2 is alone in its group
            //merge this pixel with PIXEL 1 group and find the number of 2 pi to add
            //to or subtract to unwrap it
            if(pixelB.group == -1){
                //Create new group
                if(pixelA.group == -1){
                    //System.out.println(edge.reliability + " - " + "Both new ");
                    pixelA.group = num_groups;
                    pixGroups[num_groups] = new LinkedList<>();
                    pixGroups[num_groups].add(pixPosA);
                    num_groups++;
                }

                pixelB.increment = pixelA.increment - edge.increment; //Update increment
                pixelB.group = pixelA.group; //Assign to group of pixel A
                pixGroups[pixelA.group].add(pixPosB);
            }
            //PIXEL 1 is alone in its group
            //merge this pixel with PIXEL 2 group and find the number of 2 pi to add
            //to or subtract to unwrap it
            else if(pixelA.group == -1) {
                //System.out.println(edge.reliability + " - " + "New A");
                pixelA.increment = pixelB.increment + edge.increment; //Update increment
                pixelA.group = pixelB.group; //Assign to group of pixel B
                pixGroups[pixelB.group].add(pixPosA);
            }
            //PIXEL 1 and PIXEL 2 both have groups
            else {
                //the no. of pixels in PIXEL 1 group is large than the no. of PIXELs
                //in PIXEL 2 group.   Merge PIXEL 2 group to PIXEL 1 group
                //and find the number of wraps between PIXEL 2 group and PIXEL 1 group
                //to unwrap PIXEL 2 group with respect to PIXEL 1 group.
                //the no. of wraps will be added to PIXEL 2 grop in the future
                if(pixGroups[pixelA.group].size() > pixGroups[pixelB.group].size()){ //merge PIXEL 2 with PIXEL 1 group
                    if(pixGroups[pixelA.group].contains(pixPosB)){
                        System.out.println("??????");
                    }
                    //System.out.println(edge.reliability + " - " + "Big A");
                    int inc = pixelA.increment - pixelB.increment - edge.increment;

                    for(Integer index : pixGroups[pixelB.group]){
                        pixArray[index].group = pixelA.group;
                        pixArray[index].increment += inc;
                        pixGroups[pixelA.group].add(index);
                    }

                    //pixGroups[pixelB.group] = new LinkedList<>();
                }
                //the no. of PIXELs in PIXEL 2 group is large than the no. of PIXELs
                //in PIXEL 1 group.   Merge PIXEL 1 group to PIXEL 2 group
                //and find the number of wraps between PIXEL 2 group and PIXEL 1 group
                //to unwrap PIXEL 1 group with respect to PIXEL 2 group.
                //the no. of wraps will be added to PIXEL 1 grop in the future
                else {
                    if(pixGroups[pixelB.group].contains(pixPosA)){
                        System.out.println("??????");
                    }
                    //System.out.println(edge.reliability + " - " + "Big B");
                    int inc = pixelB.increment - pixelA.increment + edge.increment;

                    for(Integer index : pixGroups[pixelA.group]){
                        pixArray[index].group = pixelB.group;
                        pixArray[index].increment += inc;
                        pixGroups[pixelB.group].add(index);
                    }
                    //pixGroups[pixelA.group] = new LinkedList<>();
                }
            }
        }
    }

    private void unwrapImage(Pixel[] pixArray, int width, int height){
        for(int i=0; i<width*height; i++){
            pixArray[i].value += (float) (Math.PI*2 * pixArray[i].increment);
        }
    }

    private float[][]returnImage(Pixel[] pixArray, int width, int height){
        float[][] retVal = new float[width][height];

        for(int i=0; i<height; i++){
            for(int j=0; j<width; j++){
                int index = width*i + j;
                retVal[j][i] = pixArray[index].value;
            }
        }

        return retVal;
    }

    public float[][] unwrap(float[][] input_img, int width, int height){
        int num_edges = (width)*(height-1) + (width-1)*(height);

        Pixel[] pixArray = new Pixel[width*height];
        Edge[] edgeArray = new Edge[num_edges];

        for(int i=0; i<height;i++){
            for(int j=0; j<width; j++) {
                pixArray[i*width+j] = new Pixel(input_img[j][i]);
            }
        }

        calculateReliability(pixArray,width,height);

        /*for(int i=0; i<width*height; i++){
            //System.out.println(pixArray[i].reliability);
        }*/

        horizontalEDGEs(pixArray,edgeArray,width,height);
        verticalEDGEs(pixArray,edgeArray,width,height);

        Arrays.sort(edgeArray);
        /*for(int i=0; i<num_edges;i++){
            //float diff = pixArray[edgeArray[i].pixPosA].value - pixArray[edgeArray[i].pixPosB].value;
            //System.out.println(edgeArray[i].increment + " : " + edgeArray[i].reliability + "---" + diff);
        }*/

        //gather PIXELs into groups
        gatherPIXELs(pixArray, edgeArray, width, height);

        /*for(int i=0; i<width*height; i++){
            //System.out.println(pixArray[i].increment);
        }*/

        //unwrap the whole image
        unwrapImage(pixArray, width, height);

        return returnImage(pixArray, width, height);
    }

    public static void main(final String... args) throws Exception {
        float[][] wrapped_img = new float[286][286];
        try (BufferedReader br = Files.newBufferedReader(Paths.get("C:\\Users\\Duo\\Desktop\\image data\\Wrapped Phase.csv"),
                StandardCharsets.US_ASCII)) {

            // read the first line from the text file
            String line = br.readLine();
            int line_count = 0;

            // loop until all lines are read
            while (line != null) {

                // use string.split to load a string array with the values from
                // each line of
                // the file, using a comma as the delimiter
                String[] elements = line.split(",");

                for(int i=0; i<elements.length;i++){
                    wrapped_img[line_count][i] = Float.parseFloat(elements[i]);
                }

                // read next line before looping
                // if end of file reached, line would be null
                line = br.readLine();
                line_count++;
            }

        } catch (IOException ioe) {
            ioe.printStackTrace();
        }

        MigualPhaseUnwrapper unwrapper = new MigualPhaseUnwrapper();
        float[][] unwrapped_phase = unwrapper.unwrap(wrapped_img,286,286);

        for(float[] row : unwrapped_phase) {
            //System.out.println(Arrays.toString(row));
        }

    }
}

