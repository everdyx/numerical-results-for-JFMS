package dynamic_programming;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

import basics.AgvPath;
import basics.LogiTask;
import basics.LogiTaskSortByLambdaTransTimeRatio;

public class DynamicProgramming {
	public ArrayList<State> states = new ArrayList<State>();
	public Map<StageState, ActionValue> optInformMap = new HashMap<StageState, ActionValue>();
	public ArrayList<Double> optVals = new ArrayList<Double>();
	public int numStage;
	public ArrayList<LogiTask> logiTaskSet;

	public DynamicProgramming(ArrayList<LogiTask> logiTaskSet) {
		this.logiTaskSet = logiTaskSet;
		Collections.sort(logiTaskSet, new LogiTaskSortByLambdaTransTimeRatio());

		for (LogiTask logiTask : logiTaskSet) {
			State state = new State();
			state.setTransTime(logiTask.getProcessTime());
			state.setDualValue(logiTask.getDualValue());
			state.setLambda(logiTask.getLambda());
			states.add(state);
		}
		this.numStage = this.states.size();
	}

	public ActionValue OptValFunction(int stage, double startTime) {
		double valAdd;
		double valSkip;
		StageState stage_state=new StageState();
		stage_state.setStage(stage);
		stage_state.setState((int) Math.round(startTime));

		ActionValue action_value= new ActionValue();

		if (stage == states.size())  {//最后一个阶段
			if (optInformMap.containsKey(stage_state)) {//最优解的映射表中有最优信息
				action_value=optInformMap.get(stage_state);
			}
			else {
				action_value.setAction(true);
				action_value.setValue(0.0);
				optInformMap.put(stage_state, action_value);
			}
		}
		else {
			if(optInformMap.containsKey(stage_state)) {//最优解的映射表中有最优信息
				action_value=optInformMap.get(stage_state);
			}
			else {
				// if add the current logistics task into the path 
				State state=states.get(stage);
				double stateTransTime=state.getTransTime();
				double lambda=state.getLambda();
				double dualVal=state.getDualValue();
				double val=(startTime+stateTransTime)*lambda-dualVal;
				StageState stage_state_Next= new StageState();
				ActionValue action_value_Next=new ActionValue();
				stage_state_Next.setStage(stage+1);
				stage_state_Next.setState( (int) Math.round(startTime+stateTransTime) );
				if (optInformMap.containsKey(stage_state_Next)) {
					action_value_Next=optInformMap.get(stage_state_Next);
				}
				else {
					action_value_Next=OptValFunction(stage+1, startTime+stateTransTime);
				}
				double optVal_Next= action_value_Next.getValue();
				valAdd=val+optVal_Next;
				
				// if exclude the current logistics task from the path 
				stage_state_Next.setState( (int) Math.round(startTime) );
				if (optInformMap.containsKey(stage_state_Next)) {
					action_value_Next=optInformMap.get(stage_state_Next);
				}
				else {
					action_value_Next=OptValFunction(stage+1, startTime);
				}
				optVal_Next=action_value_Next.getValue();
				valSkip=optVal_Next;
				if (valAdd<=valSkip) {
					action_value.setAction(true);
					action_value.setValue(valAdd);
				}
				else {
					action_value.setAction(false);
					action_value.setValue(valSkip);
				}
				optInformMap.put(stage_state, action_value);
			}
		}
		return action_value;
	}

	public Map<StageState, ActionValue> getOptInformMap() {
		return optInformMap;
	}

	public AgvPath getOptPath(int i, double state) {
		AgvPath optPath = new AgvPath();
		for (int stage = i; stage < logiTaskSet.size(); stage++) {
			StageState stage_state = new StageState((int) stage, (int) Math.round(state));
			ActionValue action_value = optInformMap.get(stage_state);
			if (action_value != null && action_value.getAction()) {
				state = state + logiTaskSet.get(stage).getProcessTime();
				optPath.addAndUpdate(logiTaskSet.get(stage));
			}
		}
		return optPath;
	}


}
