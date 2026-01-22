package column_generation;

import java.io.IOException;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import basics.AgvPath;
import basics.Instance;
import basics.LogiTask;
import basics.LogiTaskSortById;
import basics.LogiTaskSortByLambdaTransTimeRatio;
import ilog.concert.IloColumn;
import ilog.concert.IloConversion;
import ilog.concert.IloException;
import ilog.concert.IloNumVar;
import ilog.concert.IloNumVarType;
import ilog.concert.IloObjective;
import ilog.concert.IloRange;
import ilog.cplex.IloCplex;

public class RstMasterProblem {

	public IloCplex rmpModel;
	private IloObjective rmpObj;
	private List<IloConversion> mipConversion=new ArrayList<IloConversion>();
	private List<AgvPath> agvPaths=new ArrayList<AgvPath>();
//	private List<AgvPath> optAgvPaths=new ArrayList<AgvPath>();
	private Map<LogiTask,IloRange> rowTasksConstrains=new HashMap<LogiTask,IloRange>();
	private IloRange rowResourceConstrains;
//	private Map<LogiTask,Double> dualVals_LogiTasks=new HashMap<LogiTask,Double>();
	private double dualVal_AgvResource;
	public Instance para;
	public ArrayList<LogiTask> logiTaskSet=new ArrayList<LogiTask>();
	public double objVal;
	public List<IloNumVar> thetaList=new ArrayList<IloNumVar>();
	
	public double getAgvResourceDualValue() {
		return this.dualVal_AgvResource;
	}
	
	public RstMasterProblem(Instance para) throws ClassNotFoundException, IOException {
		this.para=para;
		this.logiTaskSet=para.logiTaskSet;
	}
	
	public void buildModel() {
		try{		
			rmpModel = new IloCplex();
			rmpModel.setOut(null);
			rmpObj = rmpModel.addMinimize();
			Collections.sort(logiTaskSet,new LogiTaskSortById());
			for (LogiTask logiTask : logiTaskSet) {
				rowTasksConstrains.put(logiTask, rmpModel.addRange( 1,Float.POSITIVE_INFINITY, "logitask "+logiTask.getId()));
//				rowTasksConstrains.put(logiTask, rmpModel.addRange( 1,1, "logitask "+logiTask.getId()));
			}
			rowResourceConstrains=rmpModel.addRange(Float.NEGATIVE_INFINITY,para.numVeh, "agv resource");
//			rowResourceConstrains=rmpModel.addRange(para.numVeh,para.numVeh, "agv resource");
		}
		catch (IloException e) {
			System.err.println("Concert exception caught: " + e);
		}
    }
	
    public void constructDefaultPath() {

    	Collections.sort(logiTaskSet,new LogiTaskSortByLambdaTransTimeRatio());
    	
    	ArrayDeque<LogiTask> logiTaskQueue=new ArrayDeque<LogiTask>();
    	//将logiTask按顺序压入队列
    	for (LogiTask logiTask:logiTaskSet) {
    		logiTaskQueue.addLast(logiTask);
    	}
    	//为各车分配初始任务
    	for (int k=0;k<para.numVeh;k++) {
    		AgvPath agvPath=new AgvPath();
    		LogiTask logiTask=logiTaskQueue.pollFirst();
    		agvPath.addAndUpdate(logiTask);
    		agvPaths.add(agvPath);
    	}
    	while(!logiTaskQueue.isEmpty()) {
    		LogiTask logiTask = logiTaskQueue.pollFirst();
    		double currentMinWeightedTotalCompletionTime=Integer.MAX_VALUE;
    		int indexWithMinWeightedTotalCompletionTime=0;
    		for (int i=0;i<agvPaths.size();i++) {//遍历各车，看当前各车的任务承担情况
    			double weightedTotalCompletionTime=agvPaths.get(i).getWeightedTotalCompletionTime();
    			if (weightedTotalCompletionTime<currentMinWeightedTotalCompletionTime){
    				currentMinWeightedTotalCompletionTime = weightedTotalCompletionTime;
    				indexWithMinWeightedTotalCompletionTime=i;
    			}
    		}
    		agvPaths.get(indexWithMinWeightedTotalCompletionTime).addAndUpdate(logiTask);
    	}
    	for (int i=0; i< agvPaths.size();i++) {
    		addNewColumn(agvPaths.get(i),i);
    	}
    	
    }
	
	public void addNewColumn(AgvPath agvPath,int pathId) {//可能存在的问题：找出了一个新的列添加到主问题中，这个新的列和现有所有列中的某一个一样
        try {
            IloColumn new_column = rmpModel.column(rmpObj, agvPath.getWeightedTotalCompletionTime());
            for (LogiTask logiTask : logiTaskSet ) {
                new_column = new_column.and( rmpModel.column(rowTasksConstrains.get(logiTask),(agvPath.ContainsLogiTask(logiTask)? 1 : 0) ) );
            }
            new_column = new_column.and(rmpModel.column(rowResourceConstrains,1) );            	
            agvPath.setTheta(rmpModel.numVar(new_column, 0, 1, "theta." + pathId)); 
            //注意这里并没有强制说theta是整数解 若强制整数解，则括号内的参数是(new_column, 0, 1,IloNumVarType.Int, "theta." + pathId)
            if(pathId>=101) {
            	agvPaths.add(agvPath);
            }
        } catch (IloException e) {
            System.err.println("Concert exception caught: " + e);
        }
	}
    
	public void solveRelaxation() {
		try {
			if (rmpModel.solve()) {
//				saveDualValues();
				objVal = rmpModel.getObjValue();
			}
		}
		catch (IloException e) {
			System.err.println("Concert exception caught: " + e);
		}
		
	}
	
    public void saveDualValues() {
        try {
            for (LogiTask logiTask : logiTaskSet ) {
            	logiTask.setDualValue(rmpModel.getDual(rowTasksConstrains.get(logiTask)));
            }
            dualVal_AgvResource=rmpModel.getDual(rowResourceConstrains);
        } catch (IloException e) {
            System.err.println("Concert exception caught: " + e);
        }
    }
    
	public void convertToMIP() {
		try {
			for (AgvPath agvPath : agvPaths) {
				mipConversion.add(rmpModel.conversion(agvPath.getTheta(), IloNumVarType.Bool)) ;
				rmpModel.add(mipConversion.get(mipConversion.size()-1));
			}
		}
		catch (IloException e) {
			System.err.println("Concert exception caught: " + e);
		}
	}
	
    public void solveMIP() {
		try {
			convertToMIP();

			if (rmpModel.solve()) {
//				saveDualValues();
				objVal = rmpModel.getObjValue();
//				displaySolution();
			}
			else {
				System.out.println("Integer solution not found");
			}
		}
		catch (IloException e) {
			System.err.println("Concert exception caught: " + e);
		}
	}
    
    public ArrayList<AgvPath> getOptAgvPaths() {
    	ArrayList<AgvPath> optAgvPaths=new ArrayList<AgvPath>();
		try {
			for (AgvPath agvPath : agvPaths) {
				if (rmpModel.getValue(agvPath.getTheta() )>0.99) {
					optAgvPaths.add(agvPath);
				}
			}
		}
		catch (IloException e) {
			System.err.println("Concert exception caught: " + e);
		}
		return optAgvPaths;
    }

	
	public void displaySolution() throws IloException {
			double totalCompletionTime_All_Paths = 0;
			System.out.println("\n>>>Solution for the restricted master problem" );

			ArrayList<AgvPath> optAgvPaths=new ArrayList<AgvPath>();
			optAgvPaths=getOptAgvPaths();
			for (AgvPath agvPath : optAgvPaths) {
				totalCompletionTime_All_Paths += agvPath.getWeightedTotalCompletionTime();
				System.out.print(agvPath.toString());
			}
	        for (LogiTask logiTask : logiTaskSet) {
	        	totalCompletionTime_All_Paths=totalCompletionTime_All_Paths-logiTask.getLambda()*(logiTask.getProcessTime()/2);
			}
			System.out.print("num of theta: "+agvPaths.size()+" and ");
			System.out.println("totalCompletionTime= "+ totalCompletionTime_All_Paths);
			
//			LogiProblem_Transformed logiProblem_t = new LogiProblem_Transformed(para);
//			logiProblem_t.buildModel();
//			logiProblem_t.getSolution();
//			System.out.println("optimal total completion time (obtained from MIP_Transformed): " + logiProblem_t.getObjVal());
//			
//			LogiProblem logiProblem = new LogiProblem(para);
//			logiProblem.buildModel();
//			logiProblem.getSolution();
//			System.out.println("optimal total completion time (obtained from MIP): " + logiProblem.getObjVal());

	}
	


}
