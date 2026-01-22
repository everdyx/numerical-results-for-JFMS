package basics;

import java.io.Serializable;

public class Product implements Serializable{

	/**
	 * 
	 */
	private static final long serialVersionUID = 5951897112043677532L;
	
	private int Id;
	private int demand;
	private double cycleTime;
	private double taf;//tafEarl_no_dist
	public double taf_with_dist_Earl;
	public double taf_with_dist_Late;
	public double tafHeur_Late;
	public double tafHeur_Earl;
	private double sumLambda;
	public double getSumLambda() {
		return sumLambda;
	}
	public void setSumLambda(double sumLambda) {
		this.sumLambda = sumLambda;
	}
	public double getProcessTime_without_logi() {
		return processTime_without_logi;
	}
	public void setProcessTime_without_logi(double processTime_without_logi) {
		this.processTime_without_logi = processTime_without_logi;
	}
	private double processTime_without_logi;


	
	public double getTaf() {
		return taf;
	}
	public void setTaf(double taf) {
		this.taf = taf;
	}
	public int getId() {
		return Id;
	}
	public void setId(int id) {
		Id = id;
	}
	public int getDemand() {
		return demand;
	}
	public void setDemand(int demand) {
		this.demand = demand;
	}
	public double getCycleTime() {
		return cycleTime;
	}
	public void setCycleTime(double cycleTime) {
		this.cycleTime = cycleTime;
	}
	

}
