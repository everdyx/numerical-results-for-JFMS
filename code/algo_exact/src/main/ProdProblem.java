package main;

import java.util.ArrayList;

import basics.Instance;
import ilog.concert.IloException;
import ilog.concert.IloNumExpr;
import ilog.concert.IloNumVar;
import ilog.concert.IloNumVarType;
import ilog.cplex.IloCplex;

public class ProdProblem {
	public Instance para;
	boolean show=false;
	public IloCplex model; // 定义cplex内部类的对象
	IloNumVar[][][] tss;
	IloNumVar[][][] q;
	IloNumVar[] taf;
	IloNumVar[][] x;
	double[][] xVal;
	public double[] tafVal;
	public double[][][] tssVal;
	double[] tauVal;
	public ArrayList<Integer> prodSeq = new ArrayList<>();

	IloNumVar[][] y;
	IloNumVar[] b;
	IloNumVar[] r;
	IloNumVar[] tau;
	double[][] yVal;
	Solution sol;
	Solution_1Prod sol_1Prod;
	double runTime;
	ArrayList<ArrayList<Integer>> routes = new ArrayList<>();
	public double objVal; // 目标值object
//	ArrayList<ArrayList<Integer>> sequences = new ArrayList<>();
	public int[] givenProdSeq;
	public boolean collaboration;
	public boolean given_q_pmk = true;

	private int pindex;

	private double benchTss0;

	private double benchTau;
	public double[][][] lambda;
	private double maxTaf;
	private double minTard;

//	private int preNumLogiTask;

	public double getObjVal() {
		return objVal;
	}
	
	public ProdProblem(Instance para2,double maxTaf,double minTard) {
		// TODO Auto-generated constructor stub
		this.para = para2;
		this.maxTaf=maxTaf;
		this.minTard=minTard;
//		this.lambda=lambda;
	}


	public void buildModel() throws IloException {
		try {
			model = new IloCplex();
			model.setOut(null);
			x = new IloNumVar[para.numProd + 1][para.numProd + 1];
			taf = new IloNumVar[para.numProd + 1];
			tss = new IloNumVar[para.numProd + 1][para.numMach][para.maxK];
			q = new IloNumVar[para.numProd + 1][para.numMach][para.maxK];

//			int numNodes = para.logiTaskSet.size() + 1 + para.numCharVisit;
//			y = new IloNumVar[numNodes][numNodes];// 下标para.logiTaskSet.size()表示depot
//			b = new IloNumVar[numNodes];
//			r = new IloNumVar[numNodes];
//			tau = new IloNumVar[numNodes];
			// 定义变量类型、范围、名称
			for (int p1 = 0; p1 < para.numProd + 1; p1++) {
				for (int p2 = 0; p2 < para.numProd + 1; p2++) {
					if (p1 == p2) {
						x[p1][p2] = model.numVar(0, 0, IloNumVarType.Int, "x[" + p1 + "," + p2 + "]");
					} else {
						x[p1][p2] = model.numVar(0, 1, IloNumVarType.Int, "x[" + p1 + "," + p2 + "]");
					}
				}
			}
			for (int p = 0; p < para.numProd + 1; p++) {
				taf[p] = model.numVar(0, maxTaf, IloNumVarType.Float, "taf[" + p + "]");
			}
			for (int p = 0; p < para.numProd + 1; p++) {
				for (int m = 0; m < para.numMach; m++) {
					for (int k = 0; k < para.K[p][m]; k++) {
						tss[p][m][k] = model.numVar(minTard, Float.MAX_VALUE, IloNumVarType.Float,
								"tss[" + p + "," + m + "," + k + "]");
						q[p][m][k] = model.numVar(0, para.quanProdMach[p][m], IloNumVarType.Int,
								"q[" + p + "," + m + "," + k + "]");
					}
				}
			}
			// 目标函数
			IloNumExpr objVal = model.numExpr();
			for (int p = 0; p < para.numProd; p++) {
				objVal = model.sum(objVal, taf[p]);
			}
			for (int p = 0; p < para.numProd; p++) {
				for (int m = 0; m < para.numMach; m++) {
					for (int k = 0; k < para.K[p][m]; k++) {
						String PMKstr = p + "_" + m + "_" + k;
						objVal=model.diff(objVal,model.prod(para.PMKtoLogiTask.get(PMKstr).getLambda(),tss[p][m][k] ) );
					}
				}
			}
			model.addMinimize(objVal);



			for (int p = 0; p < para.numProd + 1; p++) {
				IloNumExpr expr1 = model.numExpr();
				IloNumExpr expr2 = model.numExpr();
				for (int p1 = 0; p1 < para.numProd + 1; p1++) {
					expr1 = model.sum(expr1, x[p1][p]);
					expr2 = model.sum(expr2, x[p][p1]);
				}
				model.addEq(expr1, 1);
				model.addEq(expr2, 1);
			}


			if (given_q_pmk) {// 应用命题1，即已知最优的q_pmk，无需决策q_pmk的情况
				for (int p = 0; p < para.numProd; p++) {
					IloNumExpr expr1 = model.numExpr();
					expr1 = model.sum(tss[p][para.numMach - 1][para.K[p][para.numMach - 1] - 1],
							(para.demands[p]
									- (para.K[p][para.numMach - 1] - 1) * para.quanProdMach[p][para.numMach - 1])
									* para.cycTime[p]);// +(para.jobshopLength/para.numMach)/para.lineSpeed);
					model.addEq(taf[p], expr1);
				}

				for (int p2 = 0; p2 < para.numProd; p2++) {
					IloNumExpr expr1 = model.numExpr();
					for (int p1 = 0; p1 < para.numProd + 1; p1++) {
						expr1 = model.diff(1, x[p1][p2]);
						expr1 = model.prod(para.bigNum, expr1);
						expr1 = model.diff(taf[p1], expr1);
						expr1 = model.diff(expr1,
								(para.numMach - 1) * para.cycTime[p1]
										+ para.getDistance(para.machSet.get(0), para.machSet.get(para.numMach - 1))
												/ para.lineSpeed);// para.jobshopLength/para.lineSpeed);
						expr1 = model.sum(expr1, para.moldChangeTime);
						model.addGe(tss[p2][0][0], expr1);
					}
					model.addGe(tss[p2][0][0], 0);
				}

				for (int p = 0; p < para.numProd; p++) {
					for (int m = 1; m < para.numMach; m++) {
						IloNumExpr expr1 = model.numExpr();
						expr1 = model.sum(tss[p][m - 1][para.K[p][m - 1] - 1],
								(1 - (para.K[p][m - 1] - 1) * para.quanProdMach[p][m - 1]) * para.cycTime[p]
										+ para.getDistance(para.machSet.get(m - 1), para.machSet.get(m))
												/ para.lineSpeed);// (para.jobshopLength/para.numMach)/para.lineSpeed);
						model.addGe(tss[p][m][0], expr1);
					}
				}

				for (int p = 0; p < para.numProd; p++) {
					for (int m = 0; m < para.numMach; m++) {
						for (int k = 1; k < para.K[p][m]; k++) {
							IloNumExpr expr1 = model.numExpr();
							expr1 = model.sum(tss[p][m][k - 1], para.quanProdMach[p][m] * para.cycTime[p]);
							model.addGe(tss[p][m][k], expr1);
						}
					}
				}
			} else {// q_pmk未知，需要决策的情况
				for (int p = 0; p < para.numProd; p++) {
					IloNumExpr expr1 = model.numExpr();
					expr1 = model.sum(tss[p][para.numMach - 1][para.K[p][para.numMach - 1] - 1],
							model.prod(para.cycTime[p], q[p][para.numMach - 1][para.K[p][para.numMach - 1] - 1]));
					model.addEq(taf[p], expr1);
				}

				for (int p2 = 0; p2 < para.numProd; p2++) {
					IloNumExpr expr1 = model.numExpr();
					for (int p1 = 0; p1 < para.numProd + 1; p1++) {
						expr1 = model.diff(1, x[p1][p2]);
						expr1 = model.prod(para.bigNum, expr1);
						expr1 = model.diff(taf[p1], expr1);
						expr1 = model.diff(expr1, (para.demands[p1] + para.numMach - 1) * para.cycTime[p1]
								+ para.jobshopLength / para.lineSpeed);
						model.addGe(tss[p2][0][0], expr1);
					}
					model.addGe(tss[p2][0][0], 0);
				}

				for (int p = 0; p < para.numProd; p++) {
					for (int m = 1; m < para.numMach; m++) {
						IloNumExpr expr1 = model.numExpr();
						expr1 = model.diff(q[p][m - 1][para.K[p][m - 1] - 1], para.demands[p]);
						expr1 = model.sum(expr1, 1);
						expr1 = model.prod(para.cycTime[p], expr1);
						expr1 = model.sum(expr1, (para.jobshopLength / para.numMach) / para.lineSpeed);
						expr1 = model.sum(tss[p][m - 1][para.K[p][m - 1] - 1], expr1);
//						expr1 = model.sum(tss[p][m - 1][para.K[p][m - 1] - 1], 
//								model.prod(para.cycTime[p],
//										model.sum(
//												model.diff(q[p][m - 1][para.K[p][m - 1] - 1],para.demands[p])
//												,1)));
						model.addGe(tss[p][m][0], expr1);
					}
				}

				for (int p = 0; p < para.numProd; p++) {
					for (int m = 0; m < para.numMach; m++) {
						for (int k = 1; k < para.K[p][m]; k++) {
							IloNumExpr expr1 = model.numExpr();
							expr1 = model.sum(tss[p][m][k - 1], model.prod(q[p][m][k - 1], para.cycTime[p]));
							model.addGe(tss[p][m][k], expr1);
						}
					}
				}

				// 关于q_pmk的约束
				for (int p = 0; p < para.numProd; p++) {
					for (int m = 0; m < para.numMach; m++) {
						IloNumExpr expr = model.numExpr();
						for (int k = 0; k < para.K[p][m]; k++) {
							expr = model.sum(expr, q[p][m][k]);
						}
						model.addEq(expr, para.demands[p]);
					}
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
				System.out.println("problem should not solve, false!!!");
			} else {
				objVal = model.getObjValue();
				xVal = new double[para.numProd + 1][para.numProd + 1];
				for (int p = 0; p < para.numProd + 1; p++) {
					for (int p1 = 0; p1 < para.numProd + 1; p1++) {
						xVal[p][p1] = model.getValue(x[p][p1]);
					}
				}
				tafVal = new double[para.numProd + 1];
				tssVal = new double[para.numProd + 1][para.numMach][para.maxK];
				for (int p = 0; p < para.numProd; p++) {
					tafVal[p] = model.getValue(taf[p]);
					for (int m = 0; m < para.numMach; m++) {
						for (int k = 0; k < para.K[p][m]; k++) {
							tssVal[p][m][k] = model.getValue(tss[p][m][k]);
						}
					}
				}
				
				if (xVal!=null) {
					boolean flag1 = true;
					int cur = para.numProd;
					while (flag1) {
						for (int p = 0; p < para.numProd + 1; p++) {
							if (xVal[cur][p] >= 0.5) {
								prodSeq.add(p);
								cur = p;
								break;
							}
						}
						if (cur == para.numProd) {
							flag1 = false;
						}
					}
				}

//				sol = new Solution(para, xVal, yVal, tauVal, tssVal, tafVal);
//				if (show) {
//				sol.parse();
//				sol.print();
//				}
			}
		} catch (IloException e) {
			e.printStackTrace();
		}
	}

}
