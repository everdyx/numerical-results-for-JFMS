package dynamic_programming;

public class ActionValue {
	private boolean action;
	private double value;
	
	public boolean getAction() {
		return action;
	}
	public void setAction(boolean action) {
		this.action = action;
	}
	public double getValue() {
		return value;
	}
	public void setValue(double value) {
		this.value = value;
	}
	
	public String toString() {
		return action + " and " + value;
	}
	
}
