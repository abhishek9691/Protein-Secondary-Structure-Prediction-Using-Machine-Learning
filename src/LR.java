import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;

/**
 * Created by abhi on 12/7/2014.
 */
public class LR {
    private static ArrayList<Protein> trainingData = new ArrayList<Protein>();
    private static String longestSeq = "";

    public static void main(String args[]) {
        getTrainingData();
        //execute logistic regression
        //int length = longestSeq.length();
        double weightsAlpha[] = logisticReg('A');
        /*for(int i = 0; i < length; i++){
            System.out.println(weightsAlpha[i]);
        }*/
        //System.out.println("************************************************************");
        double weightsBeta[] = logisticReg('B');
        /*for(int i = 0; i < length; i++){
            System.out.println(weightsBeta[i]);
        }*/
        //System.out.println("************************************************************");
        double weightsCoil[] = logisticReg('C');
        /*for(int i = 0; i < length; i++){
            System.out.println(weightsCoil[i]);
        }*/
        //can now use the 3 weight vectors to get probabilities for each amino acid in a protein
        //rather than do a full dot product, use the sigmoid to get probability for each residue by
        //multiply the ith component of the weight vector by the ith residue
        int pickedProtein = (int) (Math.random() * trainingData.size());
        Protein protein = trainingData.get(pickedProtein);
        int proteinSize = protein.getSequence().length();
        double alphaProb[] = new double[proteinSize];
        double betaProb[] = new double[proteinSize];
        double coilProb[] = new double[proteinSize];

        for (int i = 0; i < proteinSize; i++) {
            double alpha = 1 / (1 + Math.pow(Math.E, (-1) * (weightsAlpha[i] * protein.getSequence().charAt(i))));
            double beta = 1 / (1 + Math.pow(Math.E, (-1) * (weightsBeta[i] * protein.getSequence().charAt(i))));
            double coil = 1 / (1 + Math.pow(Math.E, (-1) * (weightsCoil[i] * protein.getSequence().charAt(i))));
            alphaProb[i] = alpha;
            betaProb[i] = beta;
            coilProb[i] = coil;
        }
        protein.setAlphaProb(alphaProb);
        protein.setBetaProb(betaProb);
        protein.setCoilProb(coilProb);
        //assign the centroids, and create them
        int numCentroids = (int) (Math.random() * 10) + 7;
        //the second dimension contains the positions of each residue in the original protein
        ArrayList<ArrayList<Integer>> clusters = new ArrayList<ArrayList<Integer>>(numCentroids);
        for (int i = 0; i < numCentroids; i++) {
            clusters.add(new ArrayList<Integer>());
        }
        int centroids[] = new int[numCentroids];
        numCentroids--;
        while (numCentroids >= 0) {
            int centroid = (int) (Math.random() * protein.getSequence().length());
            centroids[numCentroids] = centroid;
            numCentroids--;
        }
        for(int x = 0; x < 5; x++){
            //assign points to closest centroid
            for (int i = 0; i < proteinSize; i++) {
                char residue = protein.getSequence().charAt(i);
                double alpha = protein.getAlphaProb()[i];
                double beta = protein.getBetaProb()[i];
                double coil = protein.getCoilProb()[i];
                //check against each centroid

                ArrayList<Integer> sameCentroids = new ArrayList<Integer>();
                for (int k = 0; k < centroids.length; k++) {
                    //get max for the current residue
                    char max = getMax(alpha, beta, coil);
                    //get max of the centroid
                    char centroidMax = getMax(protein.getAlphaProb()[centroids[k]], protein.getBetaProb()[centroids[k]], protein.getCoilProb()[centroids[k]]);
                    if (max == centroidMax) {
                        sameCentroids.add(k);
                    }
                }
                if (sameCentroids.size() > 0) {
                    int assignment = 0;
                    int minDistance = Integer.MAX_VALUE;
                    for (int k = 0; k < sameCentroids.size(); k++) {
                        int distance = Math.abs(sameCentroids.get(k) - i);
                        if (distance < minDistance) {
                            minDistance = distance;
                            assignment = k;
                        }
                    }
                    clusters.get(assignment).add(i);
                } else {
                    //add the ith residue to the first cluster
                    clusters.get(0).add(i);
                }
            }
            //compute new centroid for each cluster based on distance metric
            for (int i = 0; i < centroids.length; i++) {
                int clusterSum = 0;
                for (int k = 0; k < clusters.get(i).size(); k++) {
                    clusterSum += clusters.get(i).get(k);
                }
                //System.out.println("Before: "+centroids[i]);
                if(clusters.get(i).size() != 0) {
                    int avg = clusterSum / clusters.get(i).size();
                    centroids[i] = avg;
                }
                //System.out.println("After: "+centroids[i]);
            }
            //if any change, then repeat until convergence
            System.out.println("Centroids: ");
            for (int i = 0; i < centroids.length; i++) {
                System.out.println(clusters.get(i));
            }

            clusters.clear();
            for (int i = 0; i < centroids.length; i++) {
                clusters.add(new ArrayList<Integer>());
            }
            System.out.println("**************************************************************");
        }
    }

    private static char getMax(double a, double b, double c){
        if(a == b && a == c){
            return 'A';
        } else if(a > b && a > c){
            return 'A';
        } else if(b > a && b > c){
            return 'B';
        } else if(c > a && c > b){
            return 'C';
        } else {
            return 'A';
        }
    }

    private static double[] logisticReg(char struct){
        int longestSeqLength = longestSeq.length();
        double learningRate = 0.05;
        double weights[] = new double[longestSeqLength];
        for(int i = 0; i < longestSeqLength; i++){
            weights[i] = 0;
        }
        int iterations = 100;
        while(iterations > 0){
            double gradients[] = new double[longestSeqLength];
            for(int i = 0; i < longestSeqLength; i++){
                gradients[i] = 0;
            }
            int trainingExamples = trainingData.size();
            for(int i = 0; i < trainingExamples; i++){
                //System.out.println(dotProduct(weights, trainingData.get(i), 'A'));
                double p = 1/(1+Math.pow(Math.E, (-1)*(dotProduct(weights, trainingData.get(i), struct))));
                //System.out.println(p);
                double error = getFrequency(trainingData.get(i), struct) - p;
                //System.out.println(error);
                for(int j = 0; j < trainingData.get(i).getSequence().length(); j++){
                    //System.out.println("Before: "+gradients[j]);
                    gradients[j] = gradients[j] +
                            (error*
                                    ((int)trainingData.get(i).getStructureForSequence().charAt(j))*
                                    ((int)trainingData.get(i).getSequence().charAt(j)));
                    //System.out.println("After: "+gradients[j]);
                }
            }
            for(int i = 0; i < longestSeqLength; i++){
                weights[i] = weights[i] + (learningRate)*(gradients[i]);
                //System.out.println(weights[i]);
            }
            iterations--;
        }
        return weights;
    }

    private static double getFrequency(Protein protein, char charToComp){
        int length = protein.getSequence().length();
        int freq = 0;
        for(int i = 0; i < length; i++){
            if(protein.getStructureForSequence().charAt(i) == charToComp){
                freq++;
            }
        }
        return (freq/length);
    }

    private static double dotProduct(double weights[], Protein sequence, char charToComp){
        int length = sequence.getSequence().length();
        double dotProduct = 0;
        for(int i = 0; i < length; i++){
            int residue = (int)sequence.getSequence().charAt(i);
            char struct = sequence.getStructureForSequence().charAt(i);
            double product;
            if(struct == charToComp){
                product = weights[i]*residue*((int) charToComp);
            } else {
                product = weights[i]*residue;
            }
            dotProduct += product;
        }
        return dotProduct;
    }

    private static void getTrainingData(){
        BufferedReader br = null;

        try {

            String sCurrentLine;

            br = new BufferedReader(new FileReader("training_data.txt"));
            HashSet<String> temp = new HashSet<String>();

            while ((sCurrentLine = br.readLine()) != null) {
                if(sCurrentLine.indexOf('>') == -1) {
                    temp.add(sCurrentLine);
                }
            }
            Iterator<String> it = temp.iterator();
            String maxSeqLength = "";
            while(it.hasNext()){
                String seq = it.next();
                if(seq.length() > maxSeqLength.length())
                    maxSeqLength = seq;
                trainingData.add(new Protein(seq));
            }
            longestSeq = maxSeqLength;

        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            try {
                if (br != null)br.close();
            } catch (IOException ex) {
                ex.printStackTrace();
            }
        }
    }
}
