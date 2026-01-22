package main;


public class ProdRatio {
	private int prodId;
	private double ratio;

	public void setRatio(double ratio) {
		this.ratio = ratio;
	}

	public ProdRatio(int prodId,double ratio) {
		this.prodId=prodId;
		this.ratio=ratio;
	}
	
	public int getId() {
		return prodId;
	}
	public double getRatio() {
		return ratio;
	}
	
}
