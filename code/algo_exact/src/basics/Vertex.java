package basics;

import java.io.Serializable;

public class Vertex implements Serializable{
	/**
	 * 
	 */
	private static final long serialVersionUID = 8528549193712508659L;
	private int id;
	private double xcoord;
	private double ycoord;

	public void setCoord(double xcoord, double ycoord) {
		this.xcoord = xcoord;
		this.ycoord = ycoord;
	}

	public void setId(int id) {
		this.id = id;
	}
	
	public int getId() {
		return id;
	}

	public double getXcoord() {
		return xcoord;
	}

	public double getYcoord() {
		return ycoord;
	}
	

}
