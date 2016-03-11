/**
 * Created by abhi on 12/7/2014.
 */
public class Protein {
    private String sequence;
    private String structureForSequence;
    private double alphaProb[];
    private double betaProb[];
    private double coilProb[];

    public Protein(String seq, String structure){
        sequence = seq;
        structureForSequence = structure;
        alphaProb = new double[sequence.length()];
        betaProb = new double[sequence.length()];
        coilProb = new double[sequence.length()];
    }

    public Protein(String seq){
        sequence = seq;
        createStructure();
        alphaProb = new double[sequence.length()];
        betaProb = new double[sequence.length()];
        coilProb = new double[sequence.length()];
    }

    private void createStructure(){
        String secondary = "";
        int length = sequence.length();
        char[] structs = {'A', 'B', 'C'};
        while(secondary.length() != length){
            int numSeq = (int)(Math.random() * (length)/2);
            int struc = (int)(Math.random() * 3);
            char secondary_struc = structs[struc];
            for(int i = 0; i < numSeq; i++){
                secondary += secondary_struc;
                if(secondary.length() == length){
                    break;
                }
            }
        }
        structureForSequence = secondary;
    }

    public String toString(){
        return "Sequence:"+sequence+"; Structure:"+structureForSequence;
    }

    public double[] getAlphaProb() {
        return alphaProb;
    }

    public void setAlphaProb(double alphaProb[]) {
        this.alphaProb = alphaProb;
    }

    public double[] getBetaProb() {
        return betaProb;
    }

    public void setBetaProb(double betaProb[]) {
        this.betaProb = betaProb;
    }

    public double[] getCoilProb() {
        return coilProb;
    }

    public void setCoilProb(double coilProb[]) {
        this.coilProb = coilProb;
    }

    public String getStructureForSequence() {
        return structureForSequence;
    }

    public void setStructureForSequence(String structureForSequence) {
        this.structureForSequence = structureForSequence;
    }

    public String getSequence() {
        return sequence;
    }

    public void setSequence(String sequence) {
        this.sequence = sequence;
    }
}
