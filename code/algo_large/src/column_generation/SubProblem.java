package column_generation;

import java.io.IOException;
import java.util.ArrayList;

import basics.AgvPath;
import basics.LogiTask;
import dynamic_programming.ActionValue;
import dynamic_programming.DynamicProgramming;

public class SubProblem {
	
	public double reducedCost;
	//	public Map<LogiTask,Double> logiTaskDualValues;
//	public double agvResourceDualValue;
//	public Instance para;
	public ArrayList<LogiTask> logiTaskSet=new ArrayList<LogiTask>();
	private double lambda0;
	public AgvPath optPath=new AgvPath();
	
	public SubProblem(ArrayList<LogiTask> logiTaskSet, double dualVal_AgvResource) throws ClassNotFoundException, IOException {
		this.logiTaskSet=logiTaskSet;
    	this.lambda0=dualVal_AgvResource;
	}

	
	public AgvPath solve(){
		AgvPath optPath=new AgvPath();
		DynamicProgramming dp= new DynamicProgramming(logiTaskSet);
		ActionValue opt_action_value=new ActionValue();
		opt_action_value = dp.OptValFunction(0 , 0.0);
		reducedCost= opt_action_value.getValue() - lambda0;
		optPath=dp.getOptPath(0,0.0);
		this.optPath=optPath;
		return optPath;
	}
	
	public AgvPath getPath() {
		return this.optPath;
	}
	
	public double getReducedCost() {
		return reducedCost;
	}

	public void displaySolution() {
		System.out.println(">>>Optimal path for the subproblem in iteration:\n"+ optPath.toString()+"with reduced cost "+ reducedCost);
	}

}
