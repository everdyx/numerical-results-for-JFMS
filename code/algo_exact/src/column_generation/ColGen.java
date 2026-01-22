package column_generation;

import java.io.IOException;
import java.util.ArrayList;

import basics.Instance;
import basics.LogiTask;
import ilog.concert.IloException;

public class ColGen {
	public Instance para;
	public RstMasterProblem rstMasterProblem;
	public SubProblem subProblem;
	private double zero_reduced_cost_AbortColGen=-0.5;
	private int iter_no_improve=0;
	private int max_iter_no_improve=10;
	private double obj;
	private double objBest=Integer.MAX_VALUE;
	ArrayList<LogiTask> logiTaskSet=new ArrayList<LogiTask>();
	
	public int getMax_iter_no_improve() {
		return max_iter_no_improve;
	}

	public void setMax_iter_no_improve(int max_iter_no_improve) {
		this.max_iter_no_improve = max_iter_no_improve;
	}

	public double getObj() {
		return obj;
	}

	public void setObj(double obj) {
		this.obj = obj;
	}

//	public double getObjBest() {
//		return objBest;
//	}
//
//	public void setObjBest(double objBest) {
//		this.objBest = objBest;
//	}

//	public static void main(String[] args) throws ClassNotFoundException, IOException, IloException {
//		Instance para=new Instance();
//		ColGen cg=new ColGen(para);
//		cg.run();
//	}
	

	public ColGen(Instance para) throws ClassNotFoundException, IOException {
		this.para = para;
		this.logiTaskSet = para.logiTaskSet; 
		rstMasterProblem = new RstMasterProblem( this.para);
		rstMasterProblem.buildModel();
		rstMasterProblem.constructDefaultPath();
		Parameters.configureCplex(rstMasterProblem);

	}
	
	public void run() throws ClassNotFoundException, IOException, IloException {
        int iter = 0;
        do {
            iter++;
            rstMasterProblem.solveRelaxation();
//            rstMasterProblem.solveMIP();
//            rstMasterProblem.displaySolution();
            
            subProblem = new SubProblem(logiTaskSet,rstMasterProblem.getAgvResourceDualValue());
            subProblem.solve();
//            subProblem.displaySolution();

            rstMasterProblem.addNewColumn(subProblem.getPath(),iter+100);
//            displayIteration(iter);
            if (rstMasterProblem.objVal>=objBest) {
            	iter_no_improve++;
            }
            else {
            	objBest=rstMasterProblem.objVal;
            }
//        } while(iter_no_improve<=max_iter_no_improve);
        } while (subProblem.getReducedCost() < zero_reduced_cost_AbortColGen&&iter<500);
        
        rstMasterProblem.solveRelaxation();
        obj=rstMasterProblem.objVal;
        rstMasterProblem.solveMIP();
//        rstMasterProblem.displaySolution();
		
	}
	
//    private void displayIteration(int iter) {
//        if ((iter) % 20 == 0 || iter == 1) {
//            System.out.println();
//            System.out.print("Iteration");
//            System.out.print("       MP lb");
//            System.out.print("      SB int");
//            System.out.println();
//        }
//        System.out.format("%9.0f", (double) iter);
//        System.out.format("%15.2f", rstMasterProblem.objVal);//master lower bound
//        System.out.format("%12.4f", subProblem.reducedCost);
//        System.out.println();
//    }

}
