package basics;


import java.util.ArrayList;

import ilog.concert.IloException;
import ilog.concert.IloNumExpr;
import ilog.concert.IloNumVar;
import ilog.concert.IloNumVarType;
import ilog.cplex.IloCplex;

public class LogiProblem {
	public Instance para;

	IloCplex model; // 定义cplex内部类的对象
	IloNumVar[][] y;
	IloNumVar[] b;
	IloNumVar[] r;
	IloNumVar[] tau;
	double[][] yVal;
	LogiSolution sol;
	double runTime;
	ArrayList<ArrayList<Integer>> routes = new ArrayList<>();
	private double objVal; // 目标值object
//	ArrayList<ArrayList<Integer>> sequences = new ArrayList<>();

	public double getObjVal() {
		return objVal;
	}

	public LogiProblem(Instance para) {
		this.para = para;
	}
	
//	public static void main(String[] args) throws IloException {
//		Instance para = new Instance();
//		para.Initialization();
//		LogiProblem logiProblem = new LogiProblem(para);
//		logiProblem.buildModel();
//		logiProblem.getSolution();
//		System.out.println("run time: " + logiProblem.runTime);
//		System.out.println("optimal total completion time: " + logiProblem.objVal);
//	}

	public void buildModel() throws IloException {
		try {
			model = new IloCplex();
			model.setOut(null);
			y = new IloNumVar[para.logiTaskSet.size()+1][para.logiTaskSet.size()+1];//下标para.logiTaskSet.size()表示depot
			b = new IloNumVar[para.logiTaskSet.size()+1];
			r = new IloNumVar[para.logiTaskSet.size()+1];
			tau = new IloNumVar[para.logiTaskSet.size()+1];
			// 定义变量类型、范围、名称
			for (int i = 0; i < para.logiTaskSet.size()+1; i++) {
				b[i] = model.numVar(0, para.batteryCapacity, IloNumVarType.Float, "b[" + i + "]");
				r[i] = model.numVar(0, para.batteryCapacity, IloNumVarType.Float, "r[" + i + "]");
				tau[i] = model.numVar(0, 1e5, IloNumVarType.Float, "tau[" + i + "]");
				for (int j = 0; j < para.logiTaskSet.size()+1; j++) {
					if (i==j) {
						y[i][j] = model.numVar(0, 0, IloNumVarType.Int, "y[" + i + "," + j + "]");
					}else {
						y[i][j] = model.numVar(0, 1, IloNumVarType.Int, "y[" + i + "," + j + "]");
					}
				}
			}
			// 目标函数
			IloNumExpr objVal = model.numExpr();
			for (int i = 0; i < para.logiTaskSet.size(); i++) {
				objVal = model.sum(objVal, model.prod( para.ItoLogiTask.get(i).getLambda(),tau[i]) );
//				objVal = model.sum(objVal, model.prod( 1,tau[i]) );
			}
			model.addMinimize(objVal);

			// 约束
			for (int i = 0; i < para.logiTaskSet.size(); i++) {
				IloNumExpr expr21 = model.numExpr();
				IloNumExpr expr22 = model.numExpr();
				for (int j = 0; j < para.logiTaskSet.size()+1; j++) {
					expr21 = model.sum(expr21, y[j][i]);
					expr22 = model.sum(expr22, y[i][j]);
				}
				model.addEq(expr21, 1);
				model.addEq(expr22, 1);
			}
			
			IloNumExpr expr1 = model.numExpr();
			for (int i = 0; i < para.logiTaskSet.size()+1; i++) {
				expr1 = model.sum(expr1, y[para.logiTaskSet.size()][i]);
			}
			model.addLe(expr1, para.numVeh);

//			for (int i = 0; i < para.logiTaskSet.size(); i++) {
//				IloNumExpr expr3 = model.numExpr();
//				for (int j = 0; j <para.logiTaskSet.size()+1; j++) {
//					expr3 = model.sum(expr3, y[j][i]);
//				}
//				model.addEq(expr3, 1);
//			}
			
			for (int i = 0; i < para.logiTaskSet.size()+1; i++) {
				for (int j = 0; j < para.logiTaskSet.size(); j++) {
					IloNumExpr expr4 = model.numExpr();
					if (i==para.logiTaskSet.size() && j<para.logiTaskSet.size()) {//从depot到mach
						expr4 = model.diff(model.sum(tau[i], para.getDistance(para.depot,para.ItoLogiTask.get(j).getMach() )/para.vehSpeed ),
								model.prod(para.bigNum, model.diff(1, y[i][j])));
					}else if (i<para.logiTaskSet.size() && j==para.logiTaskSet.size()) {//从mach到depot
						expr4 = model.diff(model.sum(tau[i], para.getDistance(para.ItoLogiTask.get(i).getMach(),para.depot )/para.vehSpeed ),
								model.prod(para.bigNum, model.diff(1, y[i][j])));
					}else if(i==para.logiTaskSet.size() && j==para.logiTaskSet.size()) {//从depot到depot
						expr4 = model.diff(model.sum(tau[i], para.getDistance(para.depot,para.depot )/para.vehSpeed ),
								model.prod(para.bigNum, model.diff(1, y[i][j])));
					}else {//从mach到mach
						expr4 = model.diff(model.sum(tau[i], (para.getDistance(para.ItoLogiTask.get(i).getMach(),para.depot )+para.getDistance(para.depot,para.ItoLogiTask.get(j).getMach() ) ) /para.vehSpeed ),
								model.prod(para.bigNum, model.diff(1, y[i][j])));
					}
					model.addLe(expr4, tau[j]);
				}
			}


		} catch (IloException e) {
			System.err.println("Concert exception caught: " + e);
		}

	}

	public void getSolution() {
		try {
			double startTime = System.nanoTime();
			boolean isSolved = model.solve();
			double finishTime = System.nanoTime();
			runTime = (finishTime - startTime) / 1e9;// 求解时间，单位s
			if (isSolved == false) {
				// 模型不可解
				System.out.println("problem should not solve false!!!");
			} else {
				objVal=model.getObjValue();
				yVal=new double[para.logiTaskSet.size()+1][para.logiTaskSet.size()+1];
				for (int i = 0; i < para.logiTaskSet.size()+1; i++) {
					for (int j = 0; j < para.logiTaskSet.size()+1; j++) {
						yVal[i][j]=model.getValue(y[i][j]); 
					}
				}
				sol=new LogiSolution(para, yVal);
				sol.parse();
//				sol.solCheck();
				sol.print();
//				System.out.println(para.ItoLogiTask.get(0).getTransTime());
//				System.out.println(para.ItoLogiTask.get(1).getTransTime());
			}
		} catch (IloException e) {
			e.printStackTrace();
		}
	}

}
