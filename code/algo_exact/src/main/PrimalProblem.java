package main;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import basics.Instance;
import basics.Machine;
import ilog.concert.IloException;
import ilog.concert.IloNumExpr;
import ilog.concert.IloNumVar;
import ilog.concert.IloNumVarType;
import ilog.cplex.IloCplex;
import subgradient_method.ProdRatio;
import subgradient_method.ProdSortByWeightProcessTimeRatio;

public class PrimalProblem {
	public Instance para;
	public boolean show=false;
	public IloCplex model; // 定义cplex内部类的对象
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
	IloNumVar[] tauhat;
	double[][] yVal;
	Solution sol;
	Solution_1Prod sol_1Prod;
	double runTime;
	ArrayList<ArrayList<Integer>> routes = new ArrayList<>();
	public double objVal; // 目标值object
//	ArrayList<ArrayList<Integer>> sequences = new ArrayList<>();
	public int[] givenProdSeq;
	public boolean collaboration;
	public boolean given_q_pmk;
	private int pindex;

	private double benchTss0;

	private double benchTau;
	public boolean charge_given;
	public double z_f;

//	private int preNumLogiTask;

	public double getObjVal() {
		return objVal;
	}

	public PrimalProblem(Instance para,boolean given_q, boolean charge,double z_f) {
		this.para = para;
		this.collaboration = true;
		this.given_q_pmk=given_q;
		this.charge_given=charge;
		this.z_f=z_f;
	}

	public PrimalProblem(Instance para2, int[] prodSeq2,boolean given_q, boolean charge,double z_f) {
		// TODO Auto-generated constructor stub
		this.para = para2;
		this.givenProdSeq = prodSeq2;
		this.collaboration = false;
		this.given_q_pmk=given_q;
		this.charge_given=charge;
		this.z_f=z_f;
	}

	public PrimalProblem(Instance para2, int p) {
		// TODO Auto-generated constructor stub
		this.para = para2;
		this.collaboration = false;
		this.pindex = p;
		this.z_f=z_f;
	}
	
	public PrimalProblem(Instance para2, int p, double benchTss0,double benchTau,double z_f) {
		// TODO Auto-generated constructor stub
		this.para = para2;
		this.collaboration = false;
		this.pindex = p;
		this.benchTss0=benchTss0;
		this.benchTau=benchTau;
		this.z_f=z_f;
//		this.preNumLogiTask=preNumLogiTask;
	}

//	public static void main(String[] args) throws IloException {
//		Instance para = new Instance();
//		para.Initialization();
//		PrimalProblem primalProblem = new PrimalProblem(para);
//		primalProblem.buildModel();
//		primalProblem.getSolution();
//		System.out.println("run time: " + primalProblem.runTime);
//		System.out.println("optimal total completion time: " + primalProblem.objVal);
//	}
	public void buildModel() throws IloException {
		try {
			model = new IloCplex();
			model.setOut(null);
			// 定义变量类型、范围、名称
			x = new IloNumVar[para.numProd + 1][para.numProd + 1];
			taf = new IloNumVar[para.numProd + 1];
			tss = new IloNumVar[para.numProd + 1][para.numMach][para.maxK];
			q = new IloNumVar[para.numProd + 1][para.numMach][para.maxK];
			
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
				taf[p] = model.numVar(0, Float.MAX_VALUE, IloNumVarType.Float, "taf[" + p + "]");
			}
			for (int p = 0; p < para.numProd + 1; p++) {
				for (int m = 0; m < para.numMach; m++) {
					for (int k = 0; k < para.K[p][m]; k++) {
						tss[p][m][k] = model.numVar(0, Float.MAX_VALUE, IloNumVarType.Float,
								"tss[" + p + "," + m + "," + k + "]");
						q[p][m][k] = model.numVar(0, para.quanProdMach[p][m], IloNumVarType.Int,
								"q[" + p + "," + m + "," + k + "]");
					}
				}
			}
			
			int numNodes;
			if(!charge_given) {
				numNodes = para.logiTaskSet.size() + 1 + para.numCharVisit;
			}else {
				numNodes = para.logiTaskSet.size() + 1;
			}
			if(!charge_given) {
				y = new IloNumVar[numNodes][numNodes];// 下标para.logiTaskSet.size()表示depot
				b = new IloNumVar[numNodes];
				r = new IloNumVar[numNodes];
				tau = new IloNumVar[numNodes];
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
			}
			else {
				y = new IloNumVar[numNodes][numNodes];// 下标para.logiTaskSet.size()表示depot
				tau = new IloNumVar[numNodes];
				tauhat = new IloNumVar[numNodes];
				for (int i = 0; i < numNodes; i++) {
					tau[i] = model.numVar(0, Float.MAX_VALUE, IloNumVarType.Float, "tau[" + i + "]");
					tauhat[i] = model.numVar(0, Float.MAX_VALUE, IloNumVarType.Float, "tauhat[" + i + "]");
					for (int j = 0; j < numNodes; j++) {
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
			for (int p = 0; p < para.numProd; p++) {
				objVal = model.sum(objVal, taf[p]);
			}
//			objVal = model.prod(objVal, para.objValFactor);
			model.addMinimize(objVal);

			// 生产与物流耦合约束
			for (int p = 0; p < para.numProd; p++) {
				for (int m = 0; m < para.numMach; m++) {
					for (int k = 0; k < para.K[p][m]; k++) {
						String PMKstr = p + "_" + m + "_" + k;
						model.addLe(tau[para.PMKtoLogiTask.get(PMKstr).getId()], tss[p][m][k]);
					}
				}
			}

			// 生产子问题约束
			if (!collaboration) {// 如果不协同则要求将生产决策固定
				givenProdSeq = new int[para.numProd];
				ArrayList<ProdRatio> prodList1 = new ArrayList<ProdRatio>();
				for (int p = 0; p < para.numProd; p++) {
					ProdRatio prodRatio = new ProdRatio(p,
							1 / ((para.demands[p] + para.numMach - 1) * para.cycTime[p]));
					prodList1.add(prodRatio);
				}
				Collections.sort(prodList1, new ProdSortByWeightProcessTimeRatio());
				for (int p = 0; p < para.numProd; p++) {
					givenProdSeq[p] = prodList1.get(p).getId();
//					System.out.print("-> " + prodSeq[p] + " ");
				}

				for (int i = 0; i < givenProdSeq.length - 1; i++) {
					model.addEq(x[givenProdSeq[i]][givenProdSeq[i + 1]], 1);
				}
			}
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

//			model.addEq(x[1][0], 1);
//			model.addEq(x[0][2], 1);

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
						String key = p2 + "_" + 0 + "_" + 0;
						double X=(para.numMach - 1) * para.cycTime[p1]+ para.getDistance(para.machSet.get(0), para.machSet.get(para.numMach - 1))/ para.lineSpeed;
						double bigNum=z_f-X-para.getDistance(para.depot,para.PMKtoLogiTask.get(key).getMach())/ para.lineSpeed;
						expr1 = model.prod(bigNum, expr1);
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
//						expr1 = model.diff(expr1, (para.demands[p1] + para.numMach - 1) * para.cycTime[p1]
//								+ para.jobshopLength / para.lineSpeed);
						expr1=model.diff(expr1, ((para.numMach - 1) * para.cycTime[p1]
								+ para.getDistance(para.machSet.get(0), para.machSet.get(para.numMach - 1))
								/ para.lineSpeed));
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
						expr1 = model.sum(expr1, para.getDistance(para.machSet.get(m-1), para.machSet.get(m))
								/ para.lineSpeed);
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

			// 物流子问题约束
			
			if(charge_given){
				for (int i = 0; i < para.logiTaskSet.size(); i++) {
					IloNumExpr expr = model.numExpr();
					for (int j = 0; j < numNodes; j++) {
						expr = model.sum(expr, y[j][i]);
					}
					model.addEq(expr, 1);
				}

				for (int i = 0; i < para.logiTaskSet.size(); i++) {
					IloNumExpr expr = model.numExpr();
					for (int j = 0; j < numNodes; j++) {
						expr = model.sum(expr, y[i][j]);
					}
					model.addEq(expr, 1);
				}

				IloNumExpr expr1 = model.numExpr();
				for (int i = 0; i < numNodes; i++) {
					expr1 = model.sum(expr1, y[para.logiTaskSet.size()][i]);
				}
				model.addEq(expr1, para.numVeh);
				
				for (int i = 0; i < numNodes; i++) {
					if(i==para.logiTaskSet.size()) {
						model.addEq(model.sum(tau[i],para.getDistance(para.depot,para.depot) / para.vehSpeed), tauhat[i]);
					}else {
						model.addEq(model.sum(tau[i],para.getDistance(para.ItoLogiTask.get(i).getMach(),para.depot) / para.vehSpeed), tauhat[i]);
					}
				}
				
				for (int ip = 0; ip < numNodes; ip++) {
					for (int i = 0; i < para.logiTaskSet.size(); i++) {
						IloNumExpr expr = model.numExpr();
						expr=model.sum(tauhat[ip],(1+para.chargingRate*para.consumingRate)*2*para.getDistance(para.ItoLogiTask.get(i).getMach(),para.depot) / para.vehSpeed);
						double X=(1+para.chargingRate*para.consumingRate)*2*para.getDistance(para.ItoLogiTask.get(i).getMach(),para.depot) / para.vehSpeed;
						double bigNum=z_f+X-2*para.getDistance(para.ItoLogiTask.get(i).getMach(),para.depot) / para.vehSpeed;
						expr=model.diff(expr,model.prod(bigNum,model.diff(1,y[ip][i])));
						model.addGe(tauhat[i], expr);
					}
				}
				for (int i = 0; i < numNodes; i++) {
					if(i==para.logiTaskSet.size()) {
						model.addGe(tauhat[i], 2*para.getDistance(para.depot,para.depot) / para.vehSpeed);
					}else {
						model.addGe(tauhat[i], 2*para.getDistance(para.ItoLogiTask.get(i).getMach(),para.depot) / para.vehSpeed);
					}
				}
			}
			else {
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

//				model.addEq(y[0][1], 1);
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
								expr4 = model.diff(model.sum(tau[i], Tij),
										model.prod(para.bigNum, model.diff(1, y[i][j])));
							} else if (i == para.logiTaskSet.size() && j > para.logiTaskSet.size()) {// 从depot到charging
								Tij = para.getDistance(para.depot, para.depot) / para.vehSpeed;
								expr4 = model.diff(
										model.sum(tau[i], model.sum(model.prod(r[j], para.chargingRate), Tij)),
										model.prod(para.bigNum, model.diff(1, y[i][j])));
							} else if (i < para.logiTaskSet.size() && j < para.logiTaskSet.size()) {// 从mach到mach
								Tij = (para.getDistance(para.ItoLogiTask.get(i).getMach(), para.depot)
										+ para.getDistance(para.depot, para.ItoLogiTask.get(j).getMach()))
										/ para.vehSpeed;
								expr4 = model.diff(model.sum(tau[i], Tij),
										model.prod(para.bigNum, model.diff(1, y[i][j])));
							} else if (i < para.logiTaskSet.size() && j > para.logiTaskSet.size()) {// 从mach到charging
								Tij = para.getDistance(para.ItoLogiTask.get(i).getMach(), para.depot) / para.vehSpeed;
								expr4 = model.diff(
										model.sum(tau[i], model.sum(model.prod(r[j], para.chargingRate), Tij)),
										model.prod(para.bigNum, model.diff(1, y[i][j])));
							} else if (i > para.logiTaskSet.size() && j < para.logiTaskSet.size()) {// 从charing到mach
								Tij = para.getDistance(para.depot, para.ItoLogiTask.get(j).getMach()) / para.vehSpeed;
								expr4 = model.diff(model.sum(tau[i], Tij),
										model.prod(para.bigNum, model.diff(1, y[i][j])));
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
								Bij = para.getDistance(para.depot, para.ItoLogiTask.get(j).getMach())
										* para.consumingRate;
								expr5 = model.diff(model.diff(b[i], Bij),
										model.prod(para.bigNum, model.diff(1, y[i][j])));
								expr6 = model.sum(model.diff(b[i], Bij),
										model.prod(para.bigNum, model.diff(1, y[i][j])));
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
								expr5 = model.diff(model.diff(b[i], Bij),
										model.prod(para.bigNum, model.diff(1, y[i][j])));
								expr6 = model.sum(model.diff(b[i], Bij),
										model.prod(para.bigNum, model.diff(1, y[i][j])));
							} else if (i < para.logiTaskSet.size() && j > para.logiTaskSet.size()) {// 从mach到charging
								Bij = para.getDistance(para.ItoLogiTask.get(i).getMach(), para.depot)
										* para.consumingRate;
								expr5 = model.diff(model.sum(model.diff(b[i], Bij), r[j]),
										model.prod(para.bigNum, model.diff(1, y[i][j])));
								expr6 = model.sum(model.sum(model.diff(b[i], Bij), r[j]),
										model.prod(para.bigNum, model.diff(1, y[i][j])));
							} else if (i > para.logiTaskSet.size() && j < para.logiTaskSet.size()) {// 从charing到machine
								Bij = para.getDistance(para.depot, para.ItoLogiTask.get(j).getMach())
										* para.consumingRate;
								expr5 = model.diff(model.diff(b[i], Bij),
										model.prod(para.bigNum, model.diff(1, y[i][j])));
								expr6 = model.sum(model.diff(b[i], Bij),
										model.prod(para.bigNum, model.diff(1, y[i][j])));
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
			}

//			model.addEq(y[3][5],1);
//			model.addEq(y[5][0],1);
//			model.addEq(y[0][2],1);
//			model.addEq(y[8][3],1);
//			model.addEq(y[3][4],1);
//			model.addEq(y[9][8],1);

		} catch (IloException e) {
			System.err.println("Concert exception caught: " + e);
		}

	}

	public void buildModel_1Prod() throws IloException {
		try {
			int numLogiTask = 0;
			int maxK = 0;
			int numMach = para.numMach;
			int[] K = new int[numMach];
			for (int m = 0; m < numMach; m++) {
				K[m] = para.K[pindex][m];
			}
			for (int m = 0; m < numMach; m++) {
				numLogiTask = numLogiTask + K[m];
				if (K[m] > maxK) {
					maxK = K[m];
				}
			}
			int numCharVisit;
			if (para.chargingRate > 0) {
				numCharVisit = numLogiTask;
			} else {
				numCharVisit = 0;
			}

			HashMap<String, Integer> StoI = new HashMap<String, Integer>();
			HashMap<Integer, Machine> ItoM = new HashMap<Integer, Machine>();
			HashMap<Integer, Integer> ItoI = new HashMap<Integer, Integer>();
			int cnt = 0;
			for (int m = 0; m < numMach; m++) {
				for (int k = 0; k < K[m]; k++) {
					StoI.put(m + "_" + k, cnt);
					ItoM.put(cnt, para.PMKtoLogiTask.get(pindex + "_" + m + "_" + k).getMach());
					ItoI.put(cnt, para.PMKtoLogiTask.get(pindex + "_" + m + "_" + k).getId());
					cnt++;
				}
			}

			model = new IloCplex();
			model.setOut(null);
//			x = new IloNumVar[1 + 1][1 + 1];
			IloNumVar taf = model.numVar(0, Float.MAX_VALUE, IloNumVarType.Float, "taf");
			IloNumVar[][] tss = new IloNumVar[numMach][maxK];
			IloNumVar[][] q = new IloNumVar[numMach][maxK];

			int numNodes = numLogiTask + 1 + numCharVisit;
			IloNumVar[][] y = new IloNumVar[numNodes][numNodes];// 下标para.logiTaskSet.size()表示depot
			IloNumVar[] b = new IloNumVar[numNodes];
			IloNumVar[] r = new IloNumVar[numNodes];
			IloNumVar[] tau = new IloNumVar[numNodes];
			// 定义变量类型、范围、名称
//			for (int p1 = 0; p1 < 1 + 1; p1++) {
//				for (int p2 = 0; p2 < 1 + 1; p2++) {
//					if (p1 == p2) {
//						x[p1][p2] = model.numVar(0, 0, IloNumVarType.Int, "x[" + p1 + "," + p2 + "]");
//					} else {
//						x[p1][p2] = model.numVar(0, 1, IloNumVarType.Int, "x[" + p1 + "," + p2 + "]");
//					}
//				}
//			}
			for (int m = 0; m < numMach; m++) {
				for (int k = 0; k < K[m]; k++) {
					tss[m][k] = model.numVar(0, Float.MAX_VALUE, IloNumVarType.Float, "tss[" + m + "," + k + "]");
					q[m][k] = model.numVar(0, para.quanProdMach[pindex][m], IloNumVarType.Int,
							"q[" + m + "," + k + "]");
				}
			}
			for (int i = 0; i < numNodes; i++) {
				b[i] = model.numVar(0, para.batteryCapacity, IloNumVarType.Float, "b[" + i + "]");
				r[i] = model.numVar(0, para.batteryCapacity, IloNumVarType.Float, "r[" + i + "]");
				tau[i] = model.numVar(0, Float.MAX_VALUE, IloNumVarType.Float, "tau[" + i + "]");
				for (int j = 0; j < numNodes; j++) {
					if (i == j || (i > numLogiTask && j > numLogiTask)) {
						y[i][j] = model.numVar(0, 0, IloNumVarType.Int, "y[" + i + "," + j + "]");
					} else {
						y[i][j] = model.numVar(0, 1, IloNumVarType.Int, "y[" + i + "," + j + "]");
					}
				}
			}
			// 目标函数
			IloNumExpr objVal = model.numExpr();
			objVal = model.sum(objVal, taf);
			model.addMinimize(objVal);

			// 生产与物流耦合约束

			for (int m = 0; m < numMach; m++) {
				for (int k = 0; k < K[m]; k++) {
					String PMKstr = m + "_" + k;
					model.addLe(tau[StoI.get(PMKstr)], tss[m][k]);
				}
			}

			// 生产子问题约束
//			if (true) {// 如果不协同则要求将生产决策固定
//				model.addEq(x[0][1], 1);
//				model.addEq(x[1][0], 1);
//				model.addEq(x[0][0], 0);
//				model.addEq(x[1][1], 0);
//			}
//			for (int p = 0; p < 1 + 1; p++) {
//				IloNumExpr expr1 = model.numExpr();
//				IloNumExpr expr2 = model.numExpr();
//				for (int p1 = 0; p1 < 1 + 1; p1++) {
//					expr1 = model.sum(expr1, x[p1][p]);
//					expr2 = model.sum(expr2, x[p][p1]);
//				}
//				model.addEq(expr1, 1);
//				model.addEq(expr2, 1);
//			}

			// 应用命题1，即已知最优的q_pmk，无需决策q_pmk的情况
			IloNumExpr expr0 = model.numExpr();
			expr0 = model.sum(tss[numMach - 1][K[numMach - 1] - 1],
					(para.demands[pindex] - (K[numMach - 1] - 1) * para.quanProdMach[pindex][numMach - 1])
							* para.cycTime[pindex]);// +(para.jobshopLength/para.numMach)/para.lineSpeed);
			model.addEq(taf, expr0);

//			for (int p2 = 0; p2 < 1; p2++) {
//				IloNumExpr expr1 = model.numExpr();
//				for (int p1 = 0; p1 < 1 + 1; p1++) {
//					expr1 = model.diff(1, x[p1][p2]);
//					expr1 = model.prod(para.bigNum, expr1);
//					expr1 = model.diff(taf[p1], expr1);
//					expr1 = model.diff(expr1,
//							(para.numMach - 1) * para.cycTime[pindex]
//									+ para.getDistance(para.machSet.get(0), para.machSet.get(para.numMach - 1))
//											/ para.lineSpeed);// para.jobshopLength/para.lineSpeed);
//					expr1 = model.sum(expr1, para.moldChangeTime);
//					model.addGe(tss[p2][0][0], expr1);
//				}
//				model.addGe(tss[p2][0][0], 0);
//			}
			model.addGe(tss[0][0], para.moldChangeTime);
//			for (int m = 0; m < numMach; m++) {
//				for (int k = 0; k < K[m]; k++) {
//					double time=para.PMKtoLogiTask.get(pindex + "_" + m + "_" + k).getChargingTime()+para.PMKtoLogiTask.get(pindex + "_" + m + "_" + k).getTransTime()/2;
//					model.addGe(tss[m][k], time);
//				}
//			}

			for (int m = 1; m < numMach; m++) {
				IloNumExpr expr1 = model.numExpr();
				expr1 = model.sum(tss[m - 1][K[m - 1] - 1],
						(1 - (K[m - 1] - 1) * para.quanProdMach[pindex][m - 1]) * para.cycTime[pindex]
								+ para.getDistance(para.machSet.get(m - 1), para.machSet.get(m)) / para.lineSpeed);// (para.jobshopLength/para.numMach)/para.lineSpeed);
				model.addGe(tss[m][0], expr1);
			}

			for (int m = 0; m < numMach; m++) {
				for (int k = 1; k < K[m]; k++) {
					IloNumExpr expr1 = model.numExpr();
					expr1 = model.sum(tss[m][k - 1], para.quanProdMach[pindex][m] * para.cycTime[pindex]);
					model.addGe(tss[m][k], expr1);
				}
			}
			boolean considerLogi = true;
			if (considerLogi) {
				// 物流子问题约束
				for (int i = 0; i < numNodes; i++) {
					if (i == numLogiTask) {
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

				for (int i = 0; i < numLogiTask; i++) {
					IloNumExpr expr = model.numExpr();
					for (int j = 0; j < numNodes; j++) {
						expr = model.sum(expr, y[j][i]);
					}
					model.addEq(expr, 1);
				}

				for (int i = numLogiTask + 1; i < numNodes; i++) {
					IloNumExpr expr = model.numExpr();
					for (int j = 0; j < numNodes; j++) {
						expr = model.sum(expr, y[j][i]);
					}
					model.addLe(expr, 1);
				}

				IloNumExpr expr1 = model.numExpr();
				for (int i = 0; i < numNodes; i++) {
					expr1 = model.sum(expr1, y[numLogiTask][i]);
				}
				model.addEq(expr1, para.numVeh);

				for (int i = 0; i < numNodes; i++) {
					for (int j = 0; j < numNodes; j++) {
						if (j == numLogiTask) {// 目的地是depot的话跳过
							continue;
						} else {
							IloNumExpr expr4 = model.numExpr();
							double Tij;
							if (i == numLogiTask && j < numLogiTask) {// 从depot到mach
								Tij = para.getDistance(para.depot, ItoM.get(j)) / para.vehSpeed;
								expr4 = model.diff(model.sum(tau[i], Tij),
										model.prod(para.bigNum, model.diff(1, y[i][j])));
							} else if (i == numLogiTask && j > numLogiTask) {// 从depot到charging
								Tij = para.getDistance(para.depot, para.depot) / para.vehSpeed;
								expr4 = model.diff(
										model.sum(tau[i], model.sum(model.prod(r[j], para.chargingRate), Tij)),
										model.prod(para.bigNum, model.diff(1, y[i][j])));
							} else if (i < numLogiTask && j < numLogiTask) {// 从mach到mach
								Tij = (para.getDistance(ItoM.get(i), para.depot)
										+ para.getDistance(para.depot, ItoM.get(j))) / para.vehSpeed;
								expr4 = model.diff(model.sum(tau[i], Tij),
										model.prod(para.bigNum, model.diff(1, y[i][j])));
							} else if (i < numLogiTask && j > numLogiTask) {// 从mach到charging
								Tij = para.getDistance(ItoM.get(i), para.depot) / para.vehSpeed;
								expr4 = model.diff(
										model.sum(tau[i], model.sum(model.prod(r[j], para.chargingRate), Tij)),
										model.prod(para.bigNum, model.diff(1, y[i][j])));
							} else if (i > numLogiTask && j < numLogiTask) {// 从charing到mach
								Tij = para.getDistance(para.depot, ItoM.get(j)) / para.vehSpeed;
								expr4 = model.diff(model.sum(tau[i], Tij),
										model.prod(para.bigNum, model.diff(1, y[i][j])));
							} else {
								continue;
							}
							model.addLe(expr4, tau[j]);
						}
					}
				}

				for (int i = 0; i < numNodes; i++) {
					for (int j = 0; j < numNodes; j++) {
						if (j == numLogiTask) {// 目的地是depot的话跳过
							continue;
						} else {
							double Bij;
							IloNumExpr expr5 = model.numExpr();
							IloNumExpr expr6 = model.numExpr();
							if (i == numLogiTask && j < numLogiTask) {// 从depot到mach
								Bij = para.getDistance(para.depot, ItoM.get(j)) * para.consumingRate;
								expr5 = model.diff(model.diff(b[i], Bij),
										model.prod(para.bigNum, model.diff(1, y[i][j])));
								expr6 = model.sum(model.diff(b[i], Bij),
										model.prod(para.bigNum, model.diff(1, y[i][j])));
							} else if (i == numLogiTask && j > numLogiTask) {// 从depot到charging
								Bij = para.getDistance(para.depot, para.depot) * para.consumingRate;
								expr5 = model.diff(model.sum(model.diff(b[i], Bij), r[j]),
										model.prod(para.bigNum, model.diff(1, y[i][j])));
								expr6 = model.sum(model.sum(model.diff(b[i], Bij), r[j]),
										model.prod(para.bigNum, model.diff(1, y[i][j])));
							} else if (i < numLogiTask && j < numLogiTask) {// 从mach到mach
								Bij = (para.getDistance(ItoM.get(i), para.depot)
										+ para.getDistance(para.depot, ItoM.get(j))) * para.consumingRate;
								expr5 = model.diff(model.diff(b[i], Bij),
										model.prod(para.bigNum, model.diff(1, y[i][j])));
								expr6 = model.sum(model.diff(b[i], Bij),
										model.prod(para.bigNum, model.diff(1, y[i][j])));
							} else if (i < numLogiTask && j > numLogiTask) {// 从mach到charging
								Bij = para.getDistance(ItoM.get(i), para.depot) * para.consumingRate;
								expr5 = model.diff(model.sum(model.diff(b[i], Bij), r[j]),
										model.prod(para.bigNum, model.diff(1, y[i][j])));
								expr6 = model.sum(model.sum(model.diff(b[i], Bij), r[j]),
										model.prod(para.bigNum, model.diff(1, y[i][j])));
							} else if (i > numLogiTask && j < numLogiTask) {// 从charing到machine
								Bij = para.getDistance(para.depot, ItoM.get(j)) * para.consumingRate;
								expr5 = model.diff(model.diff(b[i], Bij),
										model.prod(para.bigNum, model.diff(1, y[i][j])));
								expr6 = model.sum(model.diff(b[i], Bij),
										model.prod(para.bigNum, model.diff(1, y[i][j])));
							} else {
								continue;
							}
							model.addLe(expr5, b[j]);
							model.addLe(b[j], expr6);
						}
					}
				}

				for (int i = 0; i < numNodes; i++) {
					if (i < numLogiTask) {
						model.addGe(b[i], para.getDistance(ItoM.get(i), para.depot) * para.consumingRate);
					} else {
						model.addGe(b[i], para.getDistance(para.depot, para.depot) * para.consumingRate);
					}
					model.addLe(b[i], para.batteryCapacity);
				}

				model.addEq(b[numLogiTask], 0);
			}
			
			boolean isSolved = model.solve();
			if (isSolved == false) {
				// 模型不可解
				System.out.println("problem should not solve false!!!");
			} else {
				double tafVal = model.getObjValue();
				double[][] tssVal = new double[numMach][maxK];
				for (int m = 0; m < numMach; m++) {
					for (int k = 0; k < K[m]; k++) {
						tssVal[m][k] = model.getValue(tss[m][k]);
					}
				}
				double[][] yVal = new double[numNodes][numNodes];
				for (int i = 0; i < numNodes; i++) {
					for (int j = 0; j < numNodes; j++) {
						yVal[i][j] = model.getValue(y[i][j]);
					}
				}
				double[] tauVal = new double[numNodes];
				for (int i = 0; i < numNodes; i++) {
					tauVal[i] = model.getValue(tau[i]);
				}
				double[] bVal = new double[numNodes];
				for (int i = 0; i < numNodes; i++) {
					bVal[i] = model.getValue(b[i]);
				}
//				HashMap<String, Integer> StoI = new HashMap<String, Integer>();
//				HashMap<Integer, Machine> ItoM = new HashMap<Integer, Machine>();
				sol_1Prod = new Solution_1Prod(para, yVal, tauVal, tssVal, tafVal,numLogiTask,numNodes,StoI,ItoM,ItoI,K);
//				sol_1Prod.parse_1Prod();
//				sol.solCheck();
//				sol_1Prod.print_1Prod();
			}
			

		} catch (IloException e) {
			System.err.println("Concert exception caught: " + e);
		}
		
		

	}

	public void buildModel_1Prod_v2() throws IloException {//考虑了因物流任务积累的提前时间
		try {
			int numLogiTask = 0;
			int maxK = 0;
			int numMach = para.numMach;
			int[] K = new int[numMach];
			for (int m = 0; m < numMach; m++) {
				K[m] = para.K[pindex][m];
			}
			for (int m = 0; m < numMach; m++) {
				numLogiTask = numLogiTask + K[m];
				if (K[m] > maxK) {
					maxK = K[m];
				}
			}
			int numCharVisit;
			if (para.chargingRate > 0) {
				numCharVisit = numLogiTask;
			} else {
				numCharVisit = 0;
			}

			HashMap<String, Integer> StoI = new HashMap<String, Integer>();
			HashMap<Integer, Machine> ItoM = new HashMap<Integer, Machine>();
			HashMap<Integer, Integer> ItoI = new HashMap<Integer, Integer>();
			int cnt = 0;
			for (int m = 0; m < numMach; m++) {
				for (int k = 0; k < K[m]; k++) {
					StoI.put(m + "_" + k, cnt);
					ItoM.put(cnt, para.PMKtoLogiTask.get(pindex + "_" + m + "_" + k).getMach());
					ItoI.put(cnt, para.PMKtoLogiTask.get(pindex + "_" + m + "_" + k).getId());
					cnt++;
				}
			}

			model = new IloCplex();
			model.setOut(null);
//			x = new IloNumVar[1 + 1][1 + 1];
			IloNumVar taf = model.numVar(0, Float.MAX_VALUE, IloNumVarType.Float, "taf");
			IloNumVar[][] tss = new IloNumVar[numMach][maxK];
			IloNumVar[][] q = new IloNumVar[numMach][maxK];

			int numNodes = numLogiTask + 1 + numCharVisit;
			IloNumVar[][] y = new IloNumVar[numNodes][numNodes];// 下标para.logiTaskSet.size()表示depot
			IloNumVar[] b = new IloNumVar[numNodes];
			IloNumVar[] r = new IloNumVar[numNodes];
			IloNumVar[] tau = new IloNumVar[numNodes];
			// 定义变量类型、范围、名称
//			for (int p1 = 0; p1 < 1 + 1; p1++) {
//				for (int p2 = 0; p2 < 1 + 1; p2++) {
//					if (p1 == p2) {
//						x[p1][p2] = model.numVar(0, 0, IloNumVarType.Int, "x[" + p1 + "," + p2 + "]");
//					} else {
//						x[p1][p2] = model.numVar(0, 1, IloNumVarType.Int, "x[" + p1 + "," + p2 + "]");
//					}
//				}
//			}
			for (int m = 0; m < numMach; m++) {
				for (int k = 0; k < K[m]; k++) {
					tss[m][k] = model.numVar(benchTss0, Float.MAX_VALUE, IloNumVarType.Float, "tss[" + m + "," + k + "]");
					q[m][k] = model.numVar(0, para.quanProdMach[pindex][m], IloNumVarType.Int,
							"q[" + m + "," + k + "]");
				}
			}
			for (int i = 0; i < numNodes; i++) {
				b[i] = model.numVar(0, para.batteryCapacity, IloNumVarType.Float, "b[" + i + "]");
				r[i] = model.numVar(0, para.batteryCapacity, IloNumVarType.Float, "r[" + i + "]");
				tau[i] = model.numVar(Float.MIN_VALUE, Float.MAX_VALUE, IloNumVarType.Float, "tau[" + i + "]");
				for (int j = 0; j < numNodes; j++) {
					if (i == j || (i > numLogiTask && j > numLogiTask)) {
						y[i][j] = model.numVar(0, 0, IloNumVarType.Int, "y[" + i + "," + j + "]");
					} else {
						y[i][j] = model.numVar(0, 1, IloNumVarType.Int, "y[" + i + "," + j + "]");
					}
				}
			}
			// 目标函数
			IloNumExpr objVal = model.numExpr();
			objVal = model.sum(objVal, taf);
			model.addMinimize(objVal);

			// 生产与物流耦合约束

			for (int m = 0; m < numMach; m++) {
				for (int k = 0; k < K[m]; k++) {
					String PMKstr = m + "_" + k;
					model.addLe(tau[StoI.get(PMKstr)], tss[m][k]);
				}
			}

			// 生产子问题约束
//			if (true) {// 如果不协同则要求将生产决策固定
//				model.addEq(x[0][1], 1);
//				model.addEq(x[1][0], 1);
//				model.addEq(x[0][0], 0);
//				model.addEq(x[1][1], 0);
//			}
//			for (int p = 0; p < 1 + 1; p++) {
//				IloNumExpr expr1 = model.numExpr();
//				IloNumExpr expr2 = model.numExpr();
//				for (int p1 = 0; p1 < 1 + 1; p1++) {
//					expr1 = model.sum(expr1, x[p1][p]);
//					expr2 = model.sum(expr2, x[p][p1]);
//				}
//				model.addEq(expr1, 1);
//				model.addEq(expr2, 1);
//			}

			// 应用命题1，即已知最优的q_pmk，无需决策q_pmk的情况
			IloNumExpr expr0 = model.numExpr();
			expr0 = model.sum(tss[numMach - 1][K[numMach - 1] - 1],
					(para.demands[pindex] - (K[numMach - 1] - 1) * para.quanProdMach[pindex][numMach - 1])
							* para.cycTime[pindex]);// +(para.jobshopLength/para.numMach)/para.lineSpeed);
			model.addEq(taf, expr0);

//			for (int p2 = 0; p2 < 1; p2++) {
//				IloNumExpr expr1 = model.numExpr();
//				for (int p1 = 0; p1 < 1 + 1; p1++) {
//					expr1 = model.diff(1, x[p1][p2]);
//					expr1 = model.prod(para.bigNum, expr1);
//					expr1 = model.diff(taf[p1], expr1);
//					expr1 = model.diff(expr1,
//							(para.numMach - 1) * para.cycTime[pindex]
//									+ para.getDistance(para.machSet.get(0), para.machSet.get(para.numMach - 1))
//											/ para.lineSpeed);// para.jobshopLength/para.lineSpeed);
//					expr1 = model.sum(expr1, para.moldChangeTime);
//					model.addGe(tss[p2][0][0], expr1);
//				}
//				model.addGe(tss[p2][0][0], 0);
//			}
			model.addGe(tss[0][0], para.moldChangeTime);
//			for (int m = 0; m < numMach; m++) {
//				for (int k = 0; k < K[m]; k++) {
//					double time=para.PMKtoLogiTask.get(pindex + "_" + m + "_" + k).getChargingTime()+para.PMKtoLogiTask.get(pindex + "_" + m + "_" + k).getTransTime()/2;
//					model.addGe(tss[m][k], time);
//				}
//			}

			for (int m = 1; m < numMach; m++) {
				IloNumExpr expr1 = model.numExpr();
				expr1 = model.sum(tss[m - 1][K[m - 1] - 1],
						(1 - (K[m - 1] - 1) * para.quanProdMach[pindex][m - 1]) * para.cycTime[pindex]
								+ para.getDistance(para.machSet.get(m - 1), para.machSet.get(m)) / para.lineSpeed);// (para.jobshopLength/para.numMach)/para.lineSpeed);
				model.addGe(tss[m][0], expr1);
			}

			for (int m = 0; m < numMach; m++) {
				for (int k = 1; k < K[m]; k++) {
					IloNumExpr expr1 = model.numExpr();
					expr1 = model.sum(tss[m][k - 1], para.quanProdMach[pindex][m] * para.cycTime[pindex]);
					model.addGe(tss[m][k], expr1);
				}
			}
			boolean considerLogi = true;
			if (considerLogi) {
				// 物流子问题约束
				for (int i = 0; i < numNodes; i++) {
					if (i == numLogiTask) {
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

				for (int i = 0; i < numLogiTask; i++) {
					IloNumExpr expr = model.numExpr();
					for (int j = 0; j < numNodes; j++) {
						expr = model.sum(expr, y[j][i]);
					}
					model.addEq(expr, 1);
				}

				for (int i = numLogiTask + 1; i < numNodes; i++) {
					IloNumExpr expr = model.numExpr();
					for (int j = 0; j < numNodes; j++) {
						expr = model.sum(expr, y[j][i]);
					}
					model.addLe(expr, 1);
				}

				IloNumExpr expr1 = model.numExpr();
				for (int i = 0; i < numNodes; i++) {
					expr1 = model.sum(expr1, y[numLogiTask][i]);
				}
				model.addEq(expr1, para.numVeh);
				
				// 新增的关于起始时间的约束
				model.addGe(tau[numLogiTask], benchTau);
				
				for (int i = 0; i < numNodes; i++) {
					for (int j = 0; j < numNodes; j++) {
						if (j == numLogiTask) {// 目的地是depot的话跳过
							continue;
						} else {
							IloNumExpr expr4 = model.numExpr();
							double Tij;
							if (i == numLogiTask && j < numLogiTask) {// 从depot到mach
								Tij = para.getDistance(para.depot, ItoM.get(j)) / para.vehSpeed;
								expr4 = model.diff(model.sum(tau[i], Tij),
										model.prod(para.bigNum, model.diff(1, y[i][j])));
							} else if (i == numLogiTask && j > numLogiTask) {// 从depot到charging
								Tij = para.getDistance(para.depot, para.depot) / para.vehSpeed;
								expr4 = model.diff(
										model.sum(tau[i], model.sum(model.prod(r[j], para.chargingRate), Tij)),
										model.prod(para.bigNum, model.diff(1, y[i][j])));
							} else if (i < numLogiTask && j < numLogiTask) {// 从mach到mach
								Tij = (para.getDistance(ItoM.get(i), para.depot)
										+ para.getDistance(para.depot, ItoM.get(j))) / para.vehSpeed;
								expr4 = model.diff(model.sum(tau[i], Tij),
										model.prod(para.bigNum, model.diff(1, y[i][j])));
							} else if (i < numLogiTask && j > numLogiTask) {// 从mach到charging
								Tij = para.getDistance(ItoM.get(i), para.depot) / para.vehSpeed;
								expr4 = model.diff(
										model.sum(tau[i], model.sum(model.prod(r[j], para.chargingRate), Tij)),
										model.prod(para.bigNum, model.diff(1, y[i][j])));
							} else if (i > numLogiTask && j < numLogiTask) {// 从charing到mach
								Tij = para.getDistance(para.depot, ItoM.get(j)) / para.vehSpeed;
								expr4 = model.diff(model.sum(tau[i], Tij),
										model.prod(para.bigNum, model.diff(1, y[i][j])));
							} else {
								continue;
							}
							model.addLe(expr4, tau[j]);
						}
					}
				}

				for (int i = 0; i < numNodes; i++) {
					for (int j = 0; j < numNodes; j++) {
						if (j == numLogiTask) {// 目的地是depot的话跳过
							continue;
						} else {
							double Bij;
							IloNumExpr expr5 = model.numExpr();
							IloNumExpr expr6 = model.numExpr();
							if (i == numLogiTask && j < numLogiTask) {// 从depot到mach
								Bij = para.getDistance(para.depot, ItoM.get(j)) * para.consumingRate;
								expr5 = model.diff(model.diff(b[i], Bij),
										model.prod(para.bigNum, model.diff(1, y[i][j])));
								expr6 = model.sum(model.diff(b[i], Bij),
										model.prod(para.bigNum, model.diff(1, y[i][j])));
							} else if (i == numLogiTask && j > numLogiTask) {// 从depot到charging
								Bij = para.getDistance(para.depot, para.depot) * para.consumingRate;
								expr5 = model.diff(model.sum(model.diff(b[i], Bij), r[j]),
										model.prod(para.bigNum, model.diff(1, y[i][j])));
								expr6 = model.sum(model.sum(model.diff(b[i], Bij), r[j]),
										model.prod(para.bigNum, model.diff(1, y[i][j])));
							} else if (i < numLogiTask && j < numLogiTask) {// 从mach到mach
								Bij = (para.getDistance(ItoM.get(i), para.depot)
										+ para.getDistance(para.depot, ItoM.get(j))) * para.consumingRate;
								expr5 = model.diff(model.diff(b[i], Bij),
										model.prod(para.bigNum, model.diff(1, y[i][j])));
								expr6 = model.sum(model.diff(b[i], Bij),
										model.prod(para.bigNum, model.diff(1, y[i][j])));
							} else if (i < numLogiTask && j > numLogiTask) {// 从mach到charging
								Bij = para.getDistance(ItoM.get(i), para.depot) * para.consumingRate;
								expr5 = model.diff(model.sum(model.diff(b[i], Bij), r[j]),
										model.prod(para.bigNum, model.diff(1, y[i][j])));
								expr6 = model.sum(model.sum(model.diff(b[i], Bij), r[j]),
										model.prod(para.bigNum, model.diff(1, y[i][j])));
							} else if (i > numLogiTask && j < numLogiTask) {// 从charing到machine
								Bij = para.getDistance(para.depot, ItoM.get(j)) * para.consumingRate;
								expr5 = model.diff(model.diff(b[i], Bij),
										model.prod(para.bigNum, model.diff(1, y[i][j])));
								expr6 = model.sum(model.diff(b[i], Bij),
										model.prod(para.bigNum, model.diff(1, y[i][j])));
							} else {
								continue;
							}
							model.addLe(expr5, b[j]);
							model.addLe(b[j], expr6);
						}
					}
				}

				for (int i = 0; i < numNodes; i++) {
					if (i < numLogiTask) {
						model.addGe(b[i], para.getDistance(ItoM.get(i), para.depot) * para.consumingRate);
					} else {
						model.addGe(b[i], para.getDistance(para.depot, para.depot) * para.consumingRate);
					}
					model.addLe(b[i], para.batteryCapacity);
				}

				model.addEq(b[numLogiTask], 0);
			}
			
			boolean isSolved = model.solve();
			if (isSolved == false) {
				// 模型不可解
				System.out.println("problem should not solve false!!!");
			} else {
				double tafVal = model.getObjValue();
				double[][] tssVal = new double[numMach][maxK];
				for (int m = 0; m < numMach; m++) {
					for (int k = 0; k < K[m]; k++) {
						tssVal[m][k] = model.getValue(tss[m][k]);
					}
				}
				double[][] yVal = new double[numNodes][numNodes];
				for (int i = 0; i < numNodes; i++) {
					for (int j = 0; j < numNodes; j++) {
						yVal[i][j] = model.getValue(y[i][j]);
					}
				}
				double[] tauVal = new double[numNodes];
				for (int i = 0; i < numNodes; i++) {
					tauVal[i] = model.getValue(tau[i]);
				}
				double[] bVal = new double[numNodes];
				for (int i = 0; i < numNodes; i++) {
					bVal[i] = model.getValue(b[i]);
				}
//				HashMap<String, Integer> StoI = new HashMap<String, Integer>();
//				HashMap<Integer, Machine> ItoM = new HashMap<Integer, Machine>();
				sol_1Prod = new Solution_1Prod(para, yVal, tauVal, tssVal, tafVal,numLogiTask,numNodes,StoI,ItoM,ItoI,K);
				if (show) {
				sol_1Prod.parse_1Prod();
//				sol.solCheck();
				sol_1Prod.print_1Prod();
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

				int numNodes;
				if(!charge_given) {
					numNodes= para.logiTaskSet.size() + 1 + para.numCharVisit;
				}
				else {
					numNodes= para.logiTaskSet.size() + 1;
				}
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
				if(!charge_given) {
					double[] bVal = new double[numNodes];
					for (int i = 0; i < numNodes; i++) {
						bVal[i] = model.getValue(b[i]);
					}
				}
				sol = new Solution(para, xVal, yVal, tauVal, tssVal, tafVal);
				sol.parse();
				if (show) {
				sol.parse();
//				sol.solCheck();
				sol.print();
				}
//				System.out.println("objVal="+objVal);
			}
		} catch (IloException e) {
			e.printStackTrace();
		}
	}

}
