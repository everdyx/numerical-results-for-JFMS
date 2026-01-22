package basics;

import java.util.ArrayList;

import ilog.concert.IloException;
import ilog.concert.IloNumExpr;
import ilog.concert.IloNumVar;
import ilog.concert.IloNumVarType;
import ilog.cplex.IloCplex;

public class LogiProblem_SumC {
	public Instance para;

	IloCplex model; // 定义cplex内部类的对象
	double[] tauVal;
	ArrayList<Integer> prodSeq = new ArrayList<>();

	IloNumVar[][] y;
	IloNumVar[] b;
	IloNumVar[] r;
	IloNumVar[] tau;
	double[][] yVal;
	LogiSolution_SumC sol;
	double runTime;
	ArrayList<ArrayList<Integer>> routes = new ArrayList<>();
	public double objVal; // 目标值object
//	ArrayList<ArrayList<Integer>> sequences = new ArrayList<>();

	public double getObjVal() {
		return objVal;
	}

	public LogiProblem_SumC(Instance para) {
		this.para = para;
	}
	
	public void buildModel() throws IloException {
		try {
			model = new IloCplex();
			model.setOut(null);

			int numNodes = para.logiTaskSet.size() + 1 + para.numCharVisit;
			y = new IloNumVar[numNodes][numNodes];// 下标para.logiTaskSet.size()表示depot
			b = new IloNumVar[numNodes];
			r = new IloNumVar[numNodes];
			tau = new IloNumVar[numNodes];
			// 定义变量类型、范围、名称
			for (int i = 0; i < numNodes; i++) {
				b[i] = model.numVar(0, para.batteryCapacity, IloNumVarType.Float, "b[" + i + "]");
				r[i] = model.numVar(0, para.batteryCapacity, IloNumVarType.Float, "r[" + i + "]");
				tau[i] = model.numVar(0, Float.MAX_VALUE, IloNumVarType.Float, "tau[" + i + "]");
				for (int j = 0; j < numNodes; j++) {
					if (i == j || (i > para.logiTaskSet.size() && j > para.logiTaskSet.size())) {
						y[i][j] = model.numVar(0, 0, IloNumVarType.Int, "y[" + i + "," + j + "]");
					} else {
						y[i][j] = model.numVar(0, 1, IloNumVarType.Int, "y[" + i + "," + j + "]");
					}
				}
			}
			// 目标函数
			IloNumExpr objVal = model.numExpr();
			for (int i = 0; i < para.logiTaskSet.size(); i++) {
				objVal = model.sum(objVal, model.prod( para.ItoLogiTask.get(i).getLambda(),tau[i]) );
			}
			model.addMinimize(objVal);

			// 约束
			for (int i = 0; i < numNodes; i++) {
				if (i == para.logiTaskSet.size()) {
					continue;
				} else {
					IloNumExpr expr21 = model.numExpr();
					IloNumExpr expr22 = model.numExpr();
					for (int j = 0; j < numNodes; j++) {
						expr21 = model.sum(expr21, y[j][i]);
						expr22 = model.sum(expr22, y[i][j]);
					}
					model.addEq(expr21, expr22);
				}
			}

			for (int i = 0; i < para.logiTaskSet.size(); i++) {
				IloNumExpr expr = model.numExpr();
				for (int j = 0; j < numNodes; j++) {
					expr = model.sum(expr, y[j][i]);
				}
				model.addEq(expr, 1);
			}

			for (int i = para.logiTaskSet.size() + 1; i < numNodes; i++) {
				IloNumExpr expr = model.numExpr();
				for (int j = 0; j < numNodes; j++) {
					expr = model.sum(expr, y[j][i]);
				}
				model.addLe(expr, 1);
			}

//			model.addEq(y[0][1], 1);
			IloNumExpr expr1 = model.numExpr();
			for (int i = 0; i < numNodes; i++) {
				expr1 = model.sum(expr1, y[para.logiTaskSet.size()][i]);
			}
			model.addEq(expr1, para.numVeh);

			for (int i = 0; i < numNodes; i++) {
				for (int j = 0; j < numNodes; j++) {
					if (j == para.logiTaskSet.size()) {// 目的地是depot的话跳过
						continue;
					} else {
						IloNumExpr expr4 = model.numExpr();
						double Tij;
						if (i == para.logiTaskSet.size() && j < para.logiTaskSet.size()) {// 从depot到mach
							Tij = para.getDistance(para.depot, para.ItoLogiTask.get(j).getMach()) / para.vehSpeed;
							expr4 = model.diff(model.sum(tau[i], Tij), model.prod(para.bigNum, model.diff(1, y[i][j])));
						} else if (i == para.logiTaskSet.size() && j > para.logiTaskSet.size()) {// 从depot到charging
							Tij = para.getDistance(para.depot, para.depot) / para.vehSpeed;
							expr4 = model.diff(model.sum(tau[i], model.sum(model.prod(r[j], para.chargingRate), Tij)),
									model.prod(para.bigNum, model.diff(1, y[i][j])));
						} else if (i < para.logiTaskSet.size() && j < para.logiTaskSet.size()) {// 从mach到mach
							Tij = (para.getDistance(para.ItoLogiTask.get(i).getMach(), para.depot)
									+ para.getDistance(para.depot, para.ItoLogiTask.get(j).getMach())) / para.vehSpeed;
							expr4 = model.diff(model.sum(tau[i], Tij), model.prod(para.bigNum, model.diff(1, y[i][j])));
						} else if (i < para.logiTaskSet.size() && j > para.logiTaskSet.size()) {// 从mach到charging
							Tij = para.getDistance(para.ItoLogiTask.get(i).getMach(), para.depot) / para.vehSpeed;
							expr4 = model.diff(model.sum(tau[i], model.sum(model.prod(r[j], para.chargingRate), Tij)),
									model.prod(para.bigNum, model.diff(1, y[i][j])));
						} else if (i > para.logiTaskSet.size() && j < para.logiTaskSet.size()) {// 从charing到mach
							Tij = para.getDistance(para.depot, para.ItoLogiTask.get(j).getMach()) / para.vehSpeed;
							expr4 = model.diff(model.sum(tau[i], Tij), model.prod(para.bigNum, model.diff(1, y[i][j])));
						} else {
							continue;
						}
						model.addLe(expr4, tau[j]);
					}
				}
			}

			for (int i = 0; i < numNodes; i++) {
				for (int j = 0; j < numNodes; j++) {
					if (j == para.logiTaskSet.size()) {// 目的地是depot的话跳过
						continue;
					} else {
						double Bij;
						IloNumExpr expr5 = model.numExpr();
						IloNumExpr expr6 = model.numExpr();
						if (i == para.logiTaskSet.size() && j < para.logiTaskSet.size()) {// 从depot到mach
							Bij = para.getDistance(para.depot, para.ItoLogiTask.get(j).getMach()) * para.consumingRate;
							expr5 = model.diff(model.diff(b[i], Bij), model.prod(para.bigNum, model.diff(1, y[i][j])));
							expr6 = model.sum(model.diff(b[i], Bij), model.prod(para.bigNum, model.diff(1, y[i][j])));
						} else if (i == para.logiTaskSet.size() && j > para.logiTaskSet.size()) {// 从depot到charging
							Bij = para.getDistance(para.depot, para.depot) * para.consumingRate;
							expr5 = model.diff(model.sum(model.diff(b[i], Bij), r[j]),
									model.prod(para.bigNum, model.diff(1, y[i][j])));
							expr6 = model.sum(model.sum(model.diff(b[i], Bij), r[j]),
									model.prod(para.bigNum, model.diff(1, y[i][j])));
						} else if (i < para.logiTaskSet.size() && j < para.logiTaskSet.size()) {// 从mach到mach
							Bij = (para.getDistance(para.ItoLogiTask.get(i).getMach(), para.depot)
									+ para.getDistance(para.depot, para.ItoLogiTask.get(j).getMach()))
									* para.consumingRate;
							expr5 = model.diff(model.diff(b[i], Bij), model.prod(para.bigNum, model.diff(1, y[i][j])));
							expr6 = model.sum(model.diff(b[i], Bij), model.prod(para.bigNum, model.diff(1, y[i][j])));
						} else if (i < para.logiTaskSet.size() && j > para.logiTaskSet.size()) {// 从mach到charging
							Bij = para.getDistance(para.ItoLogiTask.get(i).getMach(), para.depot) * para.consumingRate;
							expr5 = model.diff(model.sum(model.diff(b[i], Bij), r[j]),
									model.prod(para.bigNum, model.diff(1, y[i][j])));
							expr6 = model.sum(model.sum(model.diff(b[i], Bij), r[j]),
									model.prod(para.bigNum, model.diff(1, y[i][j])));
						} else if (i > para.logiTaskSet.size() && j < para.logiTaskSet.size()) {// 从charing到machine
							Bij = para.getDistance(para.depot, para.ItoLogiTask.get(j).getMach()) * para.consumingRate;
							expr5 = model.diff(model.diff(b[i], Bij), model.prod(para.bigNum, model.diff(1, y[i][j])));
							expr6 = model.sum(model.diff(b[i], Bij), model.prod(para.bigNum, model.diff(1, y[i][j])));
						} else {
							continue;
						}
						model.addLe(expr5, b[j]);
						model.addLe(b[j], expr6);
					}
				}
			}

			for (int i = 0; i < numNodes; i++) {
				if (i < para.logiTaskSet.size()) {
					model.addGe(b[i],
							para.getDistance(para.ItoLogiTask.get(i).getMach(), para.depot) * para.consumingRate);
				} else {
					model.addGe(b[i], para.getDistance(para.depot, para.depot) * para.consumingRate);
				}
				model.addLe(b[i], para.batteryCapacity);
			}

			model.addEq(b[para.logiTaskSet.size()], 0);

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
				objVal = model.getObjValue();

				int numNodes = para.logiTaskSet.size() + 1 + para.numCharVisit;
				yVal = new double[numNodes][numNodes];
				for (int i = 0; i < numNodes; i++) {
					for (int j = 0; j < numNodes; j++) {
						yVal[i][j] = model.getValue(y[i][j]);
					}
				}
				tauVal = new double[numNodes];
				for (int i = 0; i < numNodes; i++) {
					tauVal[i] = model.getValue(tau[i]);
				}
				double[] bVal = new double[numNodes];
				for (int i = 0; i < numNodes; i++) {
					bVal[i] = model.getValue(b[i]);
				}
				sol = new LogiSolution_SumC(para, yVal, tauVal);
				sol.parse();
				sol.print();
//				System.out.println("objVal="+objVal);
			}
		} catch (IloException e) {
			e.printStackTrace();
		}
	}

}
