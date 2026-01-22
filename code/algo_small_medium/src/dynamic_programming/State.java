package dynamic_programming;

public class State {
	private double totalCompTime;
	private double totalWeightedCompTime;
	private double transTime;
	private double dualValue; //对应列生成主问题的对偶变量,OR文章中用的符号是lambda
	private double lambda;//对应目标函数中任务i的完工时间的权重
	private double startTime;
	
	public void setStartTime(double startTime) {
		this.startTime=startTime;
	}
	public double getStartTime() {
		return this.startTime;
	}
	public void setTotalCompTime(double totalCompTime) {
		this.totalCompTime=totalCompTime;
	}
	public double getTotalCompTime() {
		return this.totalCompTime;
	}
	public void setTotalWeightedCompTime(double totalWeightedCompTime) {
		this.totalWeightedCompTime=totalWeightedCompTime;
	}
	public double getWeightedTotalCompTime() {
		return this.totalWeightedCompTime;
	}
	public void setLambda(double lambda) {
		this.lambda=lambda;
	}
	public double getLambda() {
		return this.lambda;
	}
	public void setTransTime(double transTime) {
		this.transTime=transTime;
	}
	public double getTransTime() {
		return this.transTime;
	}
	public void setDualValue(double dualValue) {
		this.dualValue=dualValue;
	}
	public double getDualValue() {
		return this.dualValue;
	}
}
