package basics;

import java.util.ArrayList;

import ilog.concert.IloNumVar;

public class AgvPath{

	private double weightedTotalCompletionTime;
	private IloNumVar theta;
	private ArrayList<LogiTask> executeSequence=new ArrayList<LogiTask>();
	private ArrayList<Integer> executeIdSequence=new ArrayList<Integer>();
	
	public double getWeightedTotalCompletionTime() {
		return weightedTotalCompletionTime;
	}
	
//	public void setWeightedTotalCompletionTime(double weightedTotalCompletionTime) {
//		this.weightedTotalCompletionTime = weightedTotalCompletionTime;
//	}

	public IloNumVar getTheta() {
		return theta;
	}

	public void setTheta(IloNumVar theta) {
		this.theta = theta;
	}

	public ArrayList<LogiTask> getExecuteSequence() {
		return executeSequence;
	}

	public void setExecuteSequence(ArrayList<LogiTask> executeSequence) {
		this.executeSequence = executeSequence;
	}

	public ArrayList<Integer> getExecuteIdSequence() {
		return executeIdSequence;
	}

	public void setExecuteIdSequence(ArrayList<Integer> executeIdSequence) {
		this.executeIdSequence = executeIdSequence;
	}

	
	public void addAndUpdate(LogiTask logiTask) {
		executeSequence.add(logiTask);
		executeIdSequence.add(logiTask.getId());
		double completionTime=0;
		for (LogiTask task :executeSequence) {
			completionTime=completionTime+task.getProcessTime();
		}
		this.weightedTotalCompletionTime=this.weightedTotalCompletionTime + logiTask.getLambda() * completionTime;
	}
	
	public boolean ContainsLogiTask(LogiTask logiTask) {
		return executeIdSequence.contains(logiTask.getId())? true : false;
	}
	
	public String toString() {
		StringBuffer print = new StringBuffer();
		
		for (int i = 0; i < executeSequence.size()-1; ++i) {
			print.append(String.format("%d", executeSequence.get(i).getId())+" -> ");
		}
		print.append(String.format("%d",executeSequence.get(executeSequence.size()-1).getId()));
		print.append(": ");
		for (int i = 0; i < executeSequence.size()-1; ++i) {
			print.append(String.format("M%d", executeSequence.get(i).getMach().getId())+" -> ");
		}
		print.append(String.format("M%d",executeSequence.get(executeSequence.size()-1).getMach().getId()));
		print.append(String.format(", totalWeightCompletionTime=%f",this.weightedTotalCompletionTime));
		print.append("\n");
		return print.toString();
	}
	

}
