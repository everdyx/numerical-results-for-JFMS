package basics;

import java.util.ArrayList;

import ilog.concert.IloException;
import ilog.concert.IloNumExpr;
import ilog.concert.IloNumVar;
import ilog.concert.IloNumVarType;
import ilog.cplex.IloCplex;
import main.Solution;

public class PrimalProblem_oneP_minTaf {
	public Instance para;

	IloCplex model; // 定义cplex内部类的对象
	IloNumVar[][][] tss;
	IloNumVar[][][] q;
	IloNumVar[] taf;
	IloNumVar[][] x;
	double[][] xVal;
	double[] tafVal;
	double[][][] tssVal;
	double[] tauVal;
	ArrayList<Integer> prodSeq = new ArrayList<>();

	IloNumVar[][] y;
	IloNumVar[] b;
	IloNumVar[] r;
	IloNumVar[] tau;
	double[][] yVal;
	double initialProdTime;
	double[] initialLogiTime;
	int pIndex;
	Solution sol;
	double runTime;
	ArrayList<ArrayList<Integer>> routes = new ArrayList<>();
	public double objVal; // 目标值object
//	ArrayList<ArrayList<Integer>> sequences = new ArrayList<>();
	public int[] givenProdSeq;
	public boolean collaboration;
//	ArrayList<Integer> forbidden=new ArrayList<Integer>();
	ArrayList<Integer> nodeSet=new ArrayList<Integer>();
	public boolean given_q_pmk = true;
	ArrayList<Integer> logiTaskIdSetForScheduling=new ArrayList<Integer>();
	public double getObjVal() {
		return objVal;
	}

	public PrimalProblem_oneP_minTaf(Instance para,int pIndex,double initialProdTime,double[] initialLogiTime,ArrayList<Integer> logiTaskIdSetForScheduling) {
		this.para = para;
		this.collaboration = false;
		this.givenProdSeq=new int[1];
		this.givenProdSeq[0]=pIndex;
		this.pIndex=pIndex;
		this.initialProdTime=initialProdTime;
		this.initialLogiTime=initialLogiTime;
		this.logiTaskIdSetForScheduling=logiTaskIdSetForScheduling;
	}

	public void buildModel() throws IloException {
		try {
			model = new IloCplex();
			model.setOut(null);
			x = new IloNumVar[para.numProd + 1][para.numProd + 1];
			taf = new IloNumVar[para.numProd + 1];
			tss = new IloNumVar[para.numProd + 1][para.numMach][para.maxK];
//			q = new IloNumVar[para.numProd + 1][para.numMach][para.maxK];

			int numNodes = para.logiTaskSet.size() + 1;// + logiTaskIdSetForScheduling.size();
			y = new IloNumVar[numNodes][numNodes];// 下标para.logiTaskSet.size()表示depot
//			b = new IloNumVar[numNodes];
//			r = new IloNumVar[numNodes];
			tau = new IloNumVar[numNodes];
			// 定义变量类型、范围、名称
			for (int p1 = 0; p1 < para.numProd + 1; p1++) {
				for (int p2 = 0; p2 < para.numProd + 1; p2++) {
					if ( (p1 ==para.logiTaskSet.size()  && p2==pIndex ) || (p2 ==para.logiTaskSet.size()  && p1==pIndex) ){
						x[p1][p2] = model.numVar(1, 1, IloNumVarType.Int, "x[" + p1 + "," + p2 + "]");
					} else {
						x[p1][p2] = model.numVar(0, 0, IloNumVarType.Int, "x[" + p1 + "," + p2 + "]");
					}
				}
			}
			for (int p = 0; p < para.numProd + 1; p++) {
				taf[p] = model.numVar(0, Float.MAX_VALUE, IloNumVarType.Float, "taf[" + p + "]");
			}
			for (int p = 0; p < para.numProd + 1; p++) {
				for (int m = 0; m < para.numMach; m++) {
					for (int k = 0; k < para.K[p][m]; k++) {
						tss[p][m][k] = model.numVar(0, Float.MAX_VALUE, IloNumVarType.Float,
								"tss[" + p + "," + m + "," + k + "]");
//						q[p][m][k] = model.numVar(0, para.quanProdMach[p][m], IloNumVarType.Int,
//								"q[" + p + "," + m + "," + k + "]");
					}
				}
			}
			
			for(int i = 0; i < numNodes; i++) {
				if (logiTaskIdSetForScheduling.contains(i)||i>=para.logiTaskSet.size()) {
					nodeSet.add(i);
				}
			}
			
			for (int i = 0; i < numNodes; i++) {
				if (nodeSet.contains(i)){
//					b[i] = model.numVar(0, para.batteryCapacity, IloNumVarType.Float, "b[" + i + "]");
//					r[i] = model.numVar(0, para.batteryCapacity, IloNumVarType.Float, "r[" + i + "]");
					tau[i] = model.numVar(0, Float.MAX_VALUE, IloNumVarType.Float, "tau[" + i + "]");
				}
				for (int j = 0; j < numNodes; j++) {
					if (nodeSet.contains(i) && nodeSet.contains(j)) {
						if (i == j || (i > para.logiTaskSet.size() && j > para.logiTaskSet.size())) {
							y[i][j] = model.numVar(0, 0, IloNumVarType.Int, "y[" + i + "," + j + "]");
						} else {
							y[i][j] = model.numVar(0, 1, IloNumVarType.Int, "y[" + i + "," + j + "]");
						}
					}
				}
			}
			// 目标函数
			IloNumExpr objVal = model.numExpr();
			objVal = model.sum(objVal, taf[pIndex]);
			model.addMinimize(objVal);

			// 生产与物流耦合约束
			for (int p = 0; p < para.numProd; p++) {
				for (int m = 0; m < para.numMach; m++) {
					for (int k = 0; k < para.K[p][m]; k++) {
						String PMKstr = p + "_" + m + "_" + k;
						if (p==pIndex) {
							model.addLe(tau[para.PMKtoLogiTask.get(PMKstr).getId()], tss[p][m][k]);
						}
					}
				}
			}

			// 生产子问题约束
			for (int p = 0; p < para.numProd; p++) {
				if (p==pIndex) {
					IloNumExpr expr1 = model.numExpr();
					expr1 = model.sum(tss[p][para.numMach - 1][para.K[p][para.numMach - 1] - 1],
						(para.demands[p] - (para.K[p][para.numMach - 1] - 1) * para.quanProdMach[p][para.numMach - 1])
								* para.cycTime[p]);
					model.addEq(taf[p], expr1);
				}
			}

//			for (int p2 = 0; p2 < para.numProd; p2++) {
//				IloNumExpr expr1 = model.numExpr();
//				for (int p1 = 0; p1 < para.numProd + 1; p1++) {
//					expr1 = model.diff(1, x[p1][p2]);
//					expr1 = model.prod(para.bigNum, expr1);
//					expr1 = model.diff(taf[p1], expr1);
//					model.addGe(tss[p2][0][0], expr1);
//				}
//				model.addGe(tss[p2][0][0], 0);
//			}
			
			model.addGe(tss[pIndex][0][0], initialProdTime);

			for (int p = 0; p < para.numProd; p++) {
				if (p==pIndex) {
					for (int m = 1; m < para.numMach; m++) {
						IloNumExpr expr1 = model.numExpr();
						expr1 = model.sum(tss[p][m - 1][para.K[p][m - 1] - 1],
								(1 - (para.K[p][m - 1] - 1) * para.quanProdMach[p][m - 1]) * para.cycTime[p]);
						model.addGe(tss[p][m][0], expr1);
					}
				}
			}

			for (int p = 0; p < para.numProd; p++) {
				if (p==pIndex) {
					for (int m = 0; m < para.numMach; m++) {
						for (int k = 1; k < para.K[p][m]; k++) {
							IloNumExpr expr1 = model.numExpr();
							expr1 = model.sum(tss[p][m][k - 1], para.quanProdMach[p][m] * para.cycTime[p]);
							model.addGe(tss[p][m][k], expr1);
						}
					}
				}
			}

			// 物流子问题约束
			for (int i:nodeSet) {
				if (i == para.logiTaskSet.size()) {
					continue;
				} else {
					IloNumExpr expr21 = model.numExpr();
					IloNumExpr expr22 = model.numExpr();
					for (int j:nodeSet) {
						expr21 = model.sum(expr21, y[j][i]);
						expr22 = model.sum(expr22, y[i][j]);
					}
					model.addEq(expr21, expr22);
				}
			}

				
			for (int i:nodeSet) {
				if (i<para.logiTaskSet.size()) {
					IloNumExpr expr = model.numExpr();
					for (int j:nodeSet) {
						expr = model.sum(expr, y[j][i]);
					}
					model.addEq(expr, 1);
				}
//				if (i>para.logiTaskSet.size()){
//					IloNumExpr expr = model.numExpr();
//					for (int j:nodeSet) {
//						expr = model.sum(expr, y[j][i]);
//					}
//					model.addLe(expr, 1);
//				}
			}
//			for (int i = 0; i < para.logiTaskSet.size(); i++) {
//				IloNumExpr expr = model.numExpr();
//				for (int j = 0; j < numNodes; j++) {
//					if (nodeSet.contains(i) && nodeSet.contains(j)) {
//						expr = model.sum(expr, y[j][i]);
//					}
//				}
//				model.addEq(expr, 1);
//			}
//
//			for (int i = para.logiTaskSet.size() + 1; i < numNodes; i++) {
//				IloNumExpr expr = model.numExpr();
//				for (int j = 0; j < numNodes; j++) {
//					if (nodeSet.contains(i) && nodeSet.contains(j)) {
//						expr = model.sum(expr, y[j][i]);
//					}
//				}
//				model.addLe(expr, 1);
//			}

			IloNumExpr expr1 = model.numExpr();
			for (int i = 0; i < numNodes; i++) {
				if (nodeSet.contains(i)) {
					expr1 = model.sum(expr1, y[para.logiTaskSet.size()][i]);
				}
			}
			model.addEq(expr1, para.numVeh);

			for (int i :nodeSet) {
				for (int j :nodeSet) {
					if (j == para.logiTaskSet.size()) {// 目的地是depot的话跳过
						continue;
					} else {
						IloNumExpr expr4 = model.numExpr();
						double Tij;
						if (i == para.logiTaskSet.size() && j < para.logiTaskSet.size()) {// 从depot到mach
							Tij = para.getDistance(para.depot, para.ItoLogiTask.get(j).getMach())*( 1/ para.vehSpeed+2*para.consumingRate*para.chargingRate);
							expr4 = model.diff(model.sum(tau[i], Tij), model.prod(para.bigNum, model.diff(1, y[i][j])));
//						} else if (i == para.logiTaskSet.size() && j > para.logiTaskSet.size()) {// 从depot到charging
//							Tij = para.getDistance(para.depot, para.depot) / para.vehSpeed;
//							expr4 = model.diff(model.sum(tau[i], model.sum(model.prod(r[j], para.chargingRate), Tij)),
//									model.prod(para.bigNum, model.diff(1, y[i][j])));
						} else if (i < para.logiTaskSet.size() && j < para.logiTaskSet.size()) {// 从mach到mach
							Tij = (para.getDistance(para.ItoLogiTask.get(i).getMach(), para.depot))/ para.vehSpeed
									+ para.getDistance(para.depot, para.ItoLogiTask.get(j).getMach())*( 1/ para.vehSpeed+2*para.consumingRate*para.chargingRate);
							expr4 = model.diff(model.sum(tau[i], Tij), model.prod(para.bigNum, model.diff(1, y[i][j])));
//						} else if (i < para.logiTaskSet.size() && j > para.logiTaskSet.size()) {// 从mach到charging
//							Tij = para.getDistance(para.ItoLogiTask.get(i).getMach(), para.depot) / para.vehSpeed;
//							expr4 = model.diff(model.sum(tau[i], model.sum(model.prod(r[j], para.chargingRate), Tij)),
//									model.prod(para.bigNum, model.diff(1, y[i][j])));
//						} else if (i > para.logiTaskSet.size() && j < para.logiTaskSet.size()) {// 从charing到mach
//							Tij = para.getDistance(para.depot, para.ItoLogiTask.get(j).getMach()) / para.vehSpeed;
//							expr4 = model.diff(model.sum(tau[i], Tij), model.prod(para.bigNum, model.diff(1, y[i][j])));
						} else {
							continue;
						}
						model.addLe(expr4, tau[j]);
					}
				}
			}

//			for (int i :nodeSet) {
//				for (int j :nodeSet) {
//					if (j == para.logiTaskSet.size()) {// 目的地是depot的话跳过
//						continue;
//					} else {
//						double Bij;
//						IloNumExpr expr5 = model.numExpr();
//						IloNumExpr expr6 = model.numExpr();
//						if (i == para.logiTaskSet.size() && j < para.logiTaskSet.size()) {// 从depot到mach
//							Bij = para.getDistance(para.depot, para.ItoLogiTask.get(j).getMach()) * para.consumingRate;
//							expr5 = model.diff(model.diff(b[i], Bij), model.prod(para.bigNum, model.diff(1, y[i][j])));
//							expr6 = model.sum(model.diff(b[i], Bij), model.prod(para.bigNum, model.diff(1, y[i][j])));
//						} else if (i == para.logiTaskSet.size() && j > para.logiTaskSet.size()) {// 从depot到charging
//							Bij = para.getDistance(para.depot, para.depot) * para.consumingRate;
//							expr5 = model.diff(model.sum(model.diff(b[i], Bij), r[j]),
//									model.prod(para.bigNum, model.diff(1, y[i][j])));
//							expr6 = model.sum(model.sum(model.diff(b[i], Bij), r[j]),
//									model.prod(para.bigNum, model.diff(1, y[i][j])));
//						} else if (i < para.logiTaskSet.size() && j < para.logiTaskSet.size()) {// 从mach到mach
//							Bij = (para.getDistance(para.ItoLogiTask.get(i).getMach(), para.depot)
//									+ para.getDistance(para.depot, para.ItoLogiTask.get(j).getMach()))
//									* para.consumingRate;
//							expr5 = model.diff(model.diff(b[i], Bij), model.prod(para.bigNum, model.diff(1, y[i][j])));
//							expr6 = model.sum(model.diff(b[i], Bij), model.prod(para.bigNum, model.diff(1, y[i][j])));
//						} else if (i < para.logiTaskSet.size() && j > para.logiTaskSet.size()) {// 从mach到charging
//							Bij = para.getDistance(para.ItoLogiTask.get(i).getMach(), para.depot) * para.consumingRate;
//							expr5 = model.diff(model.sum(model.diff(b[i], Bij), r[j]),
//									model.prod(para.bigNum, model.diff(1, y[i][j])));
//							expr6 = model.sum(model.sum(model.diff(b[i], Bij), r[j]),
//									model.prod(para.bigNum, model.diff(1, y[i][j])));
//						} else if (i > para.logiTaskSet.size() && j < para.logiTaskSet.size()) {// 从charing到machine
//							Bij = para.getDistance(para.depot, para.ItoLogiTask.get(j).getMach()) * para.consumingRate;
//							expr5 = model.diff(model.diff(b[i], Bij), model.prod(para.bigNum, model.diff(1, y[i][j])));
//							expr6 = model.sum(model.diff(b[i], Bij), model.prod(para.bigNum, model.diff(1, y[i][j])));
//						} else {
//							continue;
//						}
//						model.addLe(expr5, b[j]);
//						model.addLe(b[j], expr6);
//					}
//				}
//			}

//			for (int i :nodeSet) {
//				if (i < para.logiTaskSet.size()) {
//					model.addGe(b[i],
//							para.getDistance(para.ItoLogiTask.get(i).getMach(), para.depot) * para.consumingRate);
//				} else {
//					model.addGe(b[i], para.getDistance(para.depot, para.depot) * para.consumingRate);
//				}
//				model.addLe(b[i], para.batteryCapacity);
//			}
//
//			model.addEq(b[para.logiTaskSet.size()], 0);

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
//				xVal = new double[para.numProd + 1][para.numProd + 1];
//				for (int p = 0; p < para.numProd + 1; p++) {
//					for (int p1 = 0; p1 < para.numProd + 1; p1++) {
//						xVal[p][p1] = model.getValue(x[p][p1]);
//					}
//				}
				tafVal = new double[para.numProd + 1];
				tssVal = new double[para.numProd + 1][para.numMach][para.maxK];
				for (int p = 0; p < para.numProd; p++) {
					if (p==pIndex) {
					tafVal[p] = model.getValue(taf[p]);
					for (int m = 0; m < para.numMach; m++) {
						for (int k = 0; k < para.K[p][m]; k++) {
							tssVal[p][m][k] = model.getValue(tss[p][m][k]);
						}
					}
					}
				}

				int numNodes = para.logiTaskSet.size() + 1 + logiTaskIdSetForScheduling.size();
				yVal = new double[numNodes][numNodes];
				for (int i:nodeSet) {
					for (int j:nodeSet) {
						yVal[i][j] = model.getValue(y[i][j]);
					}
				}
				tauVal = new double[numNodes];
				for (int i :nodeSet) {
					tauVal[i] = model.getValue(tau[i]);
				}
//				double[] bVal = new double[numNodes];
//				for (int i :nodeSet) {
//					bVal[i] = model.getValue(b[i]);
//				}
				sol = new Solution(para, xVal, yVal, tauVal, tssVal, tafVal);
				sol.parse();
//				sol.solCheck();
				sol.print();
//				System.out.println("objVal="+objVal);
			}
		} catch (IloException e) {
			e.printStackTrace();
		}
	}

}
