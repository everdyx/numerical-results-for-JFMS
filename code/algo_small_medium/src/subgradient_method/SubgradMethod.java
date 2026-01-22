package subgradient_method;

import java.io.IOException;
import java.util.*;

import basics.AgvPath;
import basics.Instance;
import basics.LogiTask;
import basics.LogiTaskSortById;
import column_generation.ColGen;
import ilog.concert.IloException;
import main.PrimalProblem;

public class SubgradMethod {

    public double UBNoImprove=Float.MAX_VALUE;
    public double LBNoLift=Float.MIN_VALUE;
    boolean show = false;
    public Instance para;
    public double numerator_for_stepsize = 0.0002;
    public int iter = 0;
    public boolean infoOutput = false;
    public int iter_for_stepsize_factor = 0;
    public int iter_no_improve = 0;
    public int max_iter_no_improve = 10;
    public int max_iter_no_improve_for_stepsize_factor = 5;
    public double stepsize = Float.MAX_VALUE;
    public double stepsizeFactor;
    double[] stepsizeUpperBound_Prod;
    public double[] subgrad;
    //	public double[] lambda;
//	public double[] bestLambdaForUpperBound;
    public double[] bestLambdaForLowerBound;
    public double bestRelaxObj = Float.MIN_VALUE;
    public double prodObj;
    public double logiObj;
    public double earliestReadyTime = Float.MAX_VALUE;// 所有产品共同的最早开始时间
    public double makespan = 0;// 所有产品共同的最早开始时间
    double precision = 1e-4; // 0.0001 设置算法终止的精度
    public double objUpBound = Float.MAX_VALUE;
    double[] tad;
    double[] taf;
    public double bestObjUpBound = Float.MAX_VALUE;
    public int bestIter;
    private boolean improved;
    public double relaxObj;
    public ArrayList<AgvPath> optAgvPaths = new ArrayList<AgvPath>();
    public ArrayList<LogiTask> logiTaskSet = new ArrayList<LogiTask>();
    ArrayList<ArrayList<Integer>> logiSeq = new ArrayList<ArrayList<Integer>>();
    ArrayList<ArrayList<Integer>> logiSeqHeur3 = new ArrayList<ArrayList<Integer>>();
    //	public boolean compare=false;
    public int[] prodSeq;
    public int[] bestProdSeq;
    double[] tard;
    double minTard;
    double[] tardLB;
    double[] tard_j;
    double[] tard2;
    public HashMap<String, LogiTask> PMKtoLogiTask = new HashMap<String, LogiTask>();
    //	HashMap<String, Double> mapHeur2=new HashMap<String, Double>();
//	HashMap<String, ArrayList<ArrayList<Integer>> > mapHeurProdSeq_LogiSeq=new HashMap<String, ArrayList<ArrayList<Integer>> >();
    public ArrayList<ArrayList<Integer>> bestLogiSeq = new ArrayList<ArrayList<Integer>>();
    private boolean newLB;
    private double maxTaf;
    private double maxTard;
    private double[] tafUB;
    public int[] iniProdSeq;
    public double prob;
    private double[] tafLB;

    public SubgradMethod(Instance para, double stepsizeFactor, boolean newLB) {
        this.para = para;
        this.stepsizeFactor = stepsizeFactor;
        this.newLB = newLB;
    }

    public void initialization() throws IloException, ClassNotFoundException, IOException {
        subgrad = new double[para.numProd];
        prodSeq = new int[para.numProd];
        iniProdSeq = new int[para.numProd];
        taf = new double[para.numProd];

        tad = new double[para.logiTaskSet.size()];
        bestProdSeq = new int[para.numProd];
        this.PMKtoLogiTask = para.PMKtoLogiTask;
        this.logiTaskSet = para.logiTaskSet;
        Collections.sort(logiTaskSet, new LogiTaskSortById());

        stepsizeUpperBound_Prod = new double[para.numProd];// stepSize的上界，确保新的lambda求和不会超过1
        for (int p = 0; p < para.numProd; p++) {
            stepsizeUpperBound_Prod[p] = Float.MAX_VALUE;
        }

        // ###计算earliestReadyTime
        earliestReadyTime = 0;
        double[] denominator_wrt_prod = new double[para.numProd];

        for (int p = 0; p < para.numProd; p++) {
            denominator_wrt_prod[p] = (para.numMach - 1 + para.demands[p]) * para.cycTime[p] + para.moldChangeTime
                    + para.getDistance(para.machSet.get(0), para.machSet.get(para.numMach - 1)) / para.lineSpeed;
        }
        ArrayList<ProdRatio> prodList = new ArrayList<ProdRatio>();
        for (int p = 0; p < para.numProd; p++) {
            ProdRatio prodRatio = new ProdRatio(p, (para.objValFactor) / denominator_wrt_prod[p]);
            prodList.add(prodRatio);
        }
        Collections.sort(prodList, new ProdSortByWeightProcessTimeRatio());
        int[] prodSeq = new int[para.numProd];
        for (int p = 0; p < para.numProd; p++) {
            prodSeq[p] = prodList.get(p).getId();
            iniProdSeq[p] = prodList.get(p).getId();
//			System.out.print("-> " + prodSeq[p] + " ");
        }
        /// 计算一个基础时间/////////////////////////////////////////
        tard = new double[para.numProd];
        minTard = Float.MAX_VALUE;
        for (int p = 0; p < para.numProd; p++) {
            PrimalProblem primalProblem1 = new PrimalProblem(para, prodSeq[p]);
            primalProblem1.buildModel_1Prod();
            double productionTimeWithoutLogiConstriants = para.prodSet.get(prodSeq[p]).getCycleTime()
                    * (para.numMach - 1 + para.demands[prodSeq[p]])
                    + para.getDistance(para.machSet.get(0), para.machSet.get(para.numMach - 1)) / para.lineSpeed;
            tard[prodSeq[p]] = primalProblem1.model.getObjValue() - productionTimeWithoutLogiConstriants;
            if (show) {
                System.out.println(prodSeq[p] + "_" + tard[prodSeq[p]] + "=" + primalProblem1.model.getObjValue() + "-"
                        + productionTimeWithoutLogiConstriants);
            }
            if (tard[prodSeq[p]] < minTard) {
                minTard = tard[prodSeq[p]];
            }
        }

//		double minTimePerLogiTask = para.getDistance(para.depot, para.machSet.get(0)) * 2 * para.consumingRate
//				* para.chargingRate + para.getDistance(para.depot, para.machSet.get(0)) * 2 / para.vehSpeed;
//		double T_LO = 0;
//		tard_j = new double[para.numProd];
//		tard_j[0] = minTard;
//		double T_PR = tard_j[0]
//				+ para.prodSet.get(prodSeq[0]).getCycleTime() * (para.numMach - 1 + para.demands[prodSeq[0]])
//				+ para.getDistance(para.machSet.get(0), para.machSet.get(para.numMach - 1)) / para.lineSpeed;
//		for (int p = 1; p < para.numProd; p++)
//		{
//			PrimalProblem primalProblem1 = new PrimalProblem(para, prodSeq[p], T_PR, T_LO);
//			primalProblem1.buildModel_1Prod_v2();
//			double productionTimeWithoutLogiConstriants = para.prodSet.get(prodSeq[p]).getCycleTime()
//					* (para.numMach - 1 + para.demands[prodSeq[p]])
//					+ para.getDistance(para.machSet.get(0), para.machSet.get(para.numMach - 1)) / para.lineSpeed;
//			double tardxx = primalProblem1.model.getObjValue() - T_PR - productionTimeWithoutLogiConstriants;
////			if (tardxx>=0.0) {
//			T_LO = T_LO + Math.min(para.K[p][0] / para.numVeh, 1) * minTimePerLogiTask;
////			}
//			T_PR = primalProblem1.model.getObjValue();
//			tard_j[p] = tardxx;
//		}
///////////////////////////////////////////////////////////
        // 计算maxTaf
//		for (int p = 0; p < para.numProd; p++)
//		{
//			prodSeq[p] = prodList.get(para.numProd - 1 - p).getId();
//		}
        double timeCur = minTard;

        // 拿一个初始解，通过初始解拿到最大的taf
        for (int p = 0; p < para.numProd; p++) {
            prodSeq[p] = prodList.get(p).getId();
        }
        ArrayList<ArrayList<Integer>> logiSeq = new ArrayList<ArrayList<Integer>>();
        double[] timeCurVeh = new double[para.numVeh];
        double[] tad = new double[para.logiTaskSet.size()];
        for (int r = 0; r < para.numVeh; r++) {
            logiSeq.add(new ArrayList<Integer>());
        }
        for (int p = 0; p < para.numProd; p++) {// product从0开始编号
            int pindex = prodSeq[p];
            ArrayList<int[]> arrayofPMKI = para.PtoArrayofPMKI.get(pindex);
            ArrayList<Integer> arrayOfI = new ArrayList<Integer>();
            for (int index = 0; index < arrayofPMKI.size(); index++) {
                arrayOfI.add(arrayofPMKI.get(index)[3]);
            }

            for (int index = 0; index < arrayOfI.size(); index++) {
                int logiTaskId = arrayOfI.get(index);// 对应logiTask的ID
                int min_index = 0;
                double requiredDeliTime = para.ItoLogiTask.get(logiTaskId).getChargingTime()
                        + para.ItoLogiTask.get(logiTaskId).getTransTime() / 2;
                double min_val = Float.MAX_VALUE;
                for (int r = 0; r < para.numVeh; r++) {
                    if (timeCurVeh[r] + requiredDeliTime < min_val) {
                        min_val = timeCurVeh[r] + requiredDeliTime;
                        min_index = r;
                    }
                }
                logiSeq.get(min_index).add(logiTaskId);
                tad[logiTaskId] = timeCurVeh[min_index] + para.ItoLogiTask.get(logiTaskId).getChargingTime()
                        + para.ItoLogiTask.get(logiTaskId).getTransTime() / 2;
                timeCurVeh[min_index] = timeCurVeh[min_index] + para.ItoLogiTask.get(logiTaskId).getProcessTime();
            }
        }
        for (int i = 0; i < logiSeq.size(); i++) {// 第i个路径
            double cur = 0;
            for (int j = 0; j < logiSeq.get(i).size(); j++) {// 第i个路径上的第j个位置
                para.ItoLogiTask.get(logiSeq.get(i).get(j)).tadHeur = cur
                        + para.getDistance(para.depot, para.ItoLogiTask.get(logiSeq.get(i).get(j)).getMach())
                        * (1 / para.vehSpeed + 2 * para.consumingRate * para.chargingRate);
                cur = para.ItoLogiTask.get(logiSeq.get(i).get(j)).tadHeur;
                cur = cur + para.getDistance(para.depot, para.ItoLogiTask.get(logiSeq.get(i).get(j)).getMach())
                        / para.vehSpeed;
            }
        }

        String key;
        timeCur = 0;// 时间游标
        double sum_taf = 0;
        double[] taf = new double[para.numProd];
        key = prodSeq[0] + "_" + 0 + "_" + 0;
        timeCur = Math.max(para.moldChangeTime, PMKtoLogiTask.get(key).tadHeur);// 第1个产品第1个机器第1个批次的物料的开始时间
        // 必然是其对应的tadHeur(已经算好了)
        for (int p = 0; p < prodSeq.length; p++) {// product从0开始编号
            for (int m = 0; m < para.numMach; m++) {// machine从0开始编号
                if (m > 0) {// 自第一个机器后，最早开工时间需要做一次减法
                    timeCur = timeCur - para.prodSet.get(prodSeq[p]).getCycleTime() * (para.demands[prodSeq[p]] - 1)
                            + para.getDistance(para.machSet.get(m - 1), para.machSet.get(m)) / para.lineSpeed;// +(para.jobshopLength/para.numMach)/para.lineSpeed;
                }
                for (int k = 0; k < para.K[prodSeq[p]][m]; k++) {// batch从0开始编号
                    key = prodSeq[p] + "_" + m + "_" + k;
                    // 通过key获得PMKtoLogiTask.get(key)是p,m,k对应的运输任务
                    timeCur = Math.max(timeCur, PMKtoLogiTask.get(key).tadHeur); // 真正的开始时间是物料到达和最早开工时间取大；
                    timeCur = timeCur + PMKtoLogiTask.get(key).getAssemTime(); // 加上装配时间等于该批次的完工时间，也是下一批次的最早开工时间；
                }
            }
            taf[prodSeq[p]] = timeCur;
            sum_taf = sum_taf + timeCur * para.objValFactor;
            timeCur = timeCur - (para.numMach - 1) * para.cycTime[prodSeq[p]] + para.moldChangeTime
                    - para.getDistance(para.machSet.get(0), para.machSet.get(para.numMach - 1)) / para.lineSpeed;// -para.jobshopLength/para.lineSpeed;
        }
        maxTaf =Arrays.stream( taf).max().getAsDouble();
//		minTard=0;
        // 计算maxTard
//		double sum_processTime=0;
//		double max_processTime=-1;
//		for (LogiTask logiTask:logiTaskSet) {
//			sum_processTime=sum_processTime+logiTask.getProcessTime();
//			if (logiTask.getProcessTime()>max_processTime) {
//				max_processTime=logiTask.getProcessTime();
//			}
//		}
//		maxTard=sum_processTime/para.numVeh+(para.numVeh-1)*max_processTime/para.numVeh;

        // 计算maxTaf
//		prodObj=solveProdProblem();//子问题求解
//		logiObj=solveLogiProblem();//子问题求解
//		relaxObj = prodObj + logiObj;
//		String key;
//		double[] maxTaf=new double[para.numProd];
//		key = prodSeq[0] + "_" + 0 + "_" + 0;
//		timeCur = Math.max( PMKtoLogiTask.get(key).getTss(), para.moldChangeTime);// 第1个产品第1个机器第1个批次的物料的开始时间（生产问题求得），是整个系统的最早开始时间
//		for (int p = 0; p < para.numProd; p++) {// product从0开始编号
//			for (int m = 0; m < para.numMach; m++) {// machine从0开始编号
//				if (m > 0) {// 自第一个机器后，最早开工时间需要做一次减法
//					timeCur = timeCur - para.prodSet.get(prodSeq[p]).getCycleTime() * (para.demands[prodSeq[p]] - 1)+para.getDistance(para.machSet.get(m-1), para.machSet.get(m))/para.lineSpeed;//+(para.jobshopLength/para.numMach)/para.lineSpeed;
//				}
//				for (int k = 0; k < para.K[prodSeq[p]][m]; k++) {// batch从0开始编号
//					key = prodSeq[p] + "_" + m + "_" + k;
//					// 通过key获得PMKtoLogiTask.get(key)是p,m,k对应的运输任务
//					timeCur = Math.max(timeCur, PMKtoLogiTask.get(key).getTad()); // 真正的开始时间是物料到达和最早开工时间取大；
//					timeCur = timeCur + PMKtoLogiTask.get(key).getAssemTime(); // 加上装配时间等于该批次的完工时间，也是下一批次的最早开工时间；
//				}
//			}
//			maxTaf[p] = timeCur;
//			timeCur=timeCur-(para.numMach-1)*para.cycTime[prodSeq[p]]+para.moldChangeTime-para.getDistance(para.machSet.get(0), para.machSet.get(para.numMach-1))/para.lineSpeed;//-para.jobshopLength/para.lineSpeed;
//		}
//		tafUB[para.numProd-1]=maxTaf[para.numProd-1];
//		for (int p = 0; p < para.numProd; p++) {
//			tafUB[p]=maxTaf[p];
//		}
        //////////////////
        // 计算makespan
//		for (int p = 0; p < para.numProd; p++) {// product从0开始编号
//			makespan=makespan+para.prodSet.get(p).getCycleTime()*(para.numMach-1+para.demands[p]);
//		}
        tardLB = new double[para.numProd];
        double[] productionTimeWithoutLogiConstriants = new double[para.numProd];
        for (int p = 0; p < para.numProd; p++) {
            productionTimeWithoutLogiConstriants[p] = para.prodSet.get(prodSeq[p]).getCycleTime()
                    * (para.numMach - 1 + para.demands[prodSeq[p]])
                    + para.getDistance(para.machSet.get(0), para.machSet.get(para.numMach - 1)) / para.lineSpeed;
        }
        tafLB = new double[para.numProd];
        for (int p = 0; p < para.numProd; p++) {
            tafLB[p] = Float.MAX_VALUE;
        }
        getTafLB();
//        tardLB[prodSeq[0]] = Math.max(tafLB[prodSeq[0]] - productionTimeWithoutLogiConstriants[prodSeq[0]], 0);
//        double benchTime = 0;
//        for (int p = 1; p < para.numProd; p++) {
//            benchTime = tardLB[prodSeq[p - 1]] + productionTimeWithoutLogiConstriants[prodSeq[p - 1]];
//            benchTime += productionTimeWithoutLogiConstriants[prodSeq[p]];
//            tardLB[prodSeq[p]] = Math.max(tafLB[prodSeq[p]] - benchTime, 0);
//        }


    }

    public static void main(String[] args) throws IOException, ClassNotFoundException, IloException {

        Instance para = new Instance(2, 2, 2, false, 1, 0.5, 3,
                "data\\demand2.txt");
        para.Initialization(1);
        SubgradMethod sgd = new SubgradMethod(para, 1, false);
        sgd.initialization();
//        int[] b={1,0};
//        double a = sgd.LocalSearch(b);
//        ArrayList<Integer> a=new ArrayList<Integer>();
//        a.add(1);
//        a.add(2);
//        ArrayList<Integer> b=sgd.listCopy(a);
//        Collections.swap(a,0,1);
//        System.out.println(a.get(0));
//        System.out.println(b.get(0));
    }

    public void getTafLB() {

        double[] minBatchPerMach = new double[para.numMach];
        ArrayList<Integer[]> batchOrderAse = new ArrayList<Integer[]>();
//		int[][] K={{2,3},{1,1}};
        for (int m = 0; m < para.numMach; m++) {
            Integer[] batchOrder_at_m = new Integer[para.numProd];
            for (int p = 0; p < para.numProd; p++) {
                batchOrder_at_m[p] = para.K[p][m];
            }
            Arrays.sort(batchOrder_at_m);
            batchOrderAse.add(batchOrder_at_m);
        }

        ArrayList<Integer[]> accumBatchOrder = new ArrayList<Integer[]>();
        for (int m = 0; m < para.numMach; m++) {
            Integer[] accu = new Integer[para.numProd];
            accu[0] = batchOrderAse.get(m)[0];
            for (int p = 1; p < para.numProd; p++) {
                accu[p] = accu[p - 1] + batchOrderAse.get(m)[p];
            }
            accumBatchOrder.add(accu);
        }

        double[] logiLB = new double[para.numProd];
        for (int num_p = 0; num_p < para.numProd; num_p++) {
            logiLB[num_p] = getLB_ParallelMachine(para, accumBatchOrder, num_p);
        }


//		for (int m = 0; m < para.numMach; m++) {
//			minBatchPerMach[m]=Integer.MAX_VALUE;
//		}
//		for (int p = 0; p < para.numProd; p++) {
//			for (int m = 0; m < para.numMach; m++) {
//				if (para.K[p][m]<minBatchPerMach[m]){
//					minBatchPerMach[m]=para.K[p][m];
//				}
//			}
//		}
//		System.out.println(minBatchPerMach);


//		double[] tss_final_of_p=new double[para.numProd];
        // 先生成生产顺序
        ArrayList<ProdRatio> prodList = new ArrayList<ProdRatio>();
        for (int p = 0; p < para.numProd; p++) {
            ProdRatio prodRatio = new ProdRatio(p,
                    1 / ((para.numMach - 1 + para.demands[p]) * para.cycTime[p] + para.moldChangeTime
                            + para.getDistance(para.machSet.get(0), para.machSet.get(para.numMach - 1))
                            / para.lineSpeed));
            prodList.add(prodRatio);
        }
        Collections.sort(prodList, new ProdSortByWeightProcessTimeRatio());
        int[] prodSeq = new int[para.numProd];
        for (int p = 0; p < para.numProd; p++) {
            prodSeq[p] = prodList.get(p).getId();
        }
        // 后物流顺序
        // 先得到tss,并记录对应的pmk
//		ArrayList<TssWithoutLogiCons> tssWithoutLogiConsList = new ArrayList<TssWithoutLogiCons>();
        double[] tafWithoutLogiCons = new double[para.numProd];
        double[][][] tssWithoutLogiCons = new double[para.numProd][para.numMach][para.maxK];
        String key;
        double timeCur = 0;// earliestReadyTime;// 时间游标
        for (int p = 0; p < para.numProd; p++) {// product从0开始编号
            for (int m = 0; m < para.numMach; m++) {// machine从0开始编号
                if (m > 0) {// 自第一个机器后，最早开工时间需要做一次减法
                    timeCur = timeCur
                            - para.prodSet.get(prodSeq[p]).getCycleTime() * (para.demands[prodSeq[p]] - 1)
                            + para.getDistance(para.machSet.get(m - 1), para.machSet.get(m))
                            / para.lineSpeed;// (para.jobshopLength/para.numMach)/para.lineSpeed;
                }
                for (int k = 0; k < para.K[prodSeq[p]][m]; k++) {// batch从0开始编号
                    key = prodSeq[p] + "_" + m + "_" + k;
                    // 通过key获得PMKtoLogiTask.get(key)是p,m,k对应的运输任务
//					TssWithoutLogiCons tssWithoutLogiCons = new TssWithoutLogiCons(key, timeCur);
//					tssWithoutLogiConsList.add(tssWithoutLogiCons);
//								PMKtoLogiTask.get(key).setTss(timeCur);
                    tssWithoutLogiCons[prodSeq[p]][m][k] = timeCur;
                    timeCur = timeCur + para.PMKtoLogiTask.get(key).getAssemTime(); // 加上装配时间等于该批次的完工时间，也是下一批次的最早开工时间；
                    if (m == para.numMach - 1 && k == para.K[prodSeq[p]][m] - 1) {
                        tafWithoutLogiCons[prodSeq[p]] = timeCur;
                    }
                }
            }
            timeCur = timeCur - (para.numMach - 1) * para.cycTime[prodSeq[p]] + para.moldChangeTime
                    - para.getDistance(para.machSet.get(0), para.machSet.get(para.numMach - 1))
                    / para.lineSpeed;// -para.jobshopLength/para.lineSpeed;
        }
//		Collections.sort(tssWithoutLogiConsList, new LogiTaskSortByTss());
//		String[] pmk = new String[para.logiTaskSet.size()];
//		for (int s = 0; s < para.logiTaskSet.size(); s++) {
//			pmk[s] = tssWithoutLogiConsList.get(s).pmk;
//		}
//        double[] tafLB = new double[para.numProd];
//        for (int p = 0; p < para.numProd; p++) {
//            tafLB[p] = Float.MAX_VALUE;
//        }
//        double[] accTaf=new  double[para.numProd];
//        accTaf[0]=tafWithoutLogiCons[prodSeq[p]]
//        for (int p = 0; p < para.numProd; p++) {
//            accTaf[p]
//        }
        double[] minDelta= new double[para.numProd];
        int[] lastMach= new int[para.numProd];
        double minD=Float.MAX_VALUE;
        int lastM=0;
        for (int p = 0; p < para.numProd; p++) {
            minDelta[p]=Float.MAX_VALUE;
        }
        for (int p = 0; p < para.numProd; p++) {// product从0开始编号
            for (int m = 0; m < para.numMach; m++) {// machine从0开始编号
                for (int k = 0; k < para.K[prodSeq[p]][m]; k++) {// batch从0开始编号
                    double delta = tafWithoutLogiCons[prodSeq[p]] - tssWithoutLogiCons[prodSeq[p]][m][k];
                    if (delta<minDelta[prodSeq[p]]){
                        minDelta[prodSeq[p]]=delta;
                        lastMach[prodSeq[p]]=m;
                    }
                    if (delta<minD){
                        minD=delta;
                    }

//                    double tad = logiLB[p] - para.getDistance(para.depot, para.machSet.get(m)) / para.vehSpeed;
//                    double temp = tad + delta;
//                    if (temp < tafLB[p]) {
//                        tafLB[p] = temp;
//                    }
                }
            }
        }
        for (int p = 0; p < para.numProd; p++) {
            if (lastMach[p]>lastM){
                lastM=lastMach[p];
            }
        }

        double a=0;
//        double minDelta=0;//new double[para.numProd];
//        for (int p = 0; p < para.numProd; p++) {// product从0开始编号
////            minDelta[p]=Float.MAX_VALUE;
//            for (int m = 0; m < para.numMach; m++) {// machine从0开始编号
//                for (int k = 0; k < para.K[p][m]; k++) {// batch从0开始编号
//                    double delta[p][][] = tafWithoutLogiCons[p] - tssWithoutLogiCons[p][m][k];
//                    if (delta<minDelta) {
//                        minDelta = delta;
//                    }
//                }
//            }
//        }
        for (int p = 0; p < para.numProd; p++) {// product从0开始编号
//            double delta = tafWithoutLogiCons[prodSeq[p]] - tssWithoutLogiCons[prodSeq[p]][m][k];
            double tad = logiLB[p] - para.getDistance(para.depot, para.machSet.get(lastM)) / para.vehSpeed;
            tafLB[p] = tad + minD;
        }
    }

    public double getLB_ParallelMachine(Instance para, ArrayList<Integer[]> accumBatchOrder, int num_p) {
        double LB = Float.MIN_VALUE;
        int[] num_task = new int[para.numMach];

        double[] processTime = new double[para.numMach];
        double totalProcTime = 0;
        for (int m = 0; m < para.numMach; m++) {
            num_task[m] = accumBatchOrder.get(m)[num_p];
            processTime[m] = (para.getDistance(para.depot, para.machSet.get(m)) + para.getDistance(para.machSet.get(m), para.depot)) / para.vehSpeed +
                    (para.getDistance(para.depot, para.machSet.get(m)) + para.getDistance(para.machSet.get(m), para.depot)) * para.consumingRate * para.chargingRate;
            totalProcTime += num_task[m] * processTime[m];
        }
        ArrayList<Double> allProcessTimes = new ArrayList<Double>();
        for (int m = 0; m < para.numMach; m++) {
            for (int i = 0; i < num_task[m]; i++) {
                allProcessTimes.add(processTime[m]);
            }
        }
        allProcessTimes.sort(Comparator.naturalOrder());
        double[] candidate_LB = new double[3];
        candidate_LB[0] = allProcessTimes.get(allProcessTimes.size() - 1);
//		if (allProcessTimes.size()-1-para.numMach>=0) {
        candidate_LB[1] = 0;//allProcessTimes.get(allProcessTimes.size() - 1 - para.numMach) + allProcessTimes.get(allProcessTimes.size() - 1 - para.numMach + 1);
//		}
        candidate_LB[2] = totalProcTime / para.numMach;
        for (int i = 0; i < 3; i++) {
            if (candidate_LB[i] > LB) {
                LB = candidate_LB[i];
            }
        }
        for (int k = 1; k <= para.numMach; k++) {
//			int alphaUB=(allProcessTimes.size()-k)/para.numMach;
            for (int l = 1; l <= allProcessTimes.size(); l++) {
//				int l=alpha*para.numMach+k;
                int lambda = k * (l / para.numMach) + Math.min(k, l - (l / para.numMach) * para.numMach);
                ArrayList<Double> jobs = new ArrayList<Double>();
                for (int i = 0; i < lambda; i++) {
                    jobs.add(allProcessTimes.get(i));
                }
                double[] temp = new double[3];
                temp[0] = jobs.get(jobs.size() - 1);
                temp[1] = 0;
                double total = 0;
                for (int i = 0; i < jobs.size(); i++) {
                    total += jobs.get(i);
                }
                temp[2] = total / k;
                for (int i = 0; i < 3; i++) {
                    if (temp[i] > LB) {
                        LB = temp[i];
                    }
                }
            }
        }
        return LB;
    }

    public void run2() throws ClassNotFoundException, IOException, IloException {
        boolean flag = true;
        bestLambdaForLowerBound = new double[logiTaskSet.size()];
        double[] bestTssForLowerBound = new double[logiTaskSet.size()];
        double[] bestTadForLowerBound = new double[logiTaskSet.size()];
        // 为每一个logiTask设置其lambda（初始化拉格朗日乘子）
        for (int p = 0; p < para.numProd; p++) {// product从0开始编号
            for (int m = 0; m < para.numMach; m++) {// machine从0开始编号
                for (int k = 0; k < para.K[p][m]; k++) {// batch从0开始编号
                    PMKtoLogiTask.get(p + "_" + m + "_" + k).setLambda((double)Math.random()*0.1);// (double)Math.random()*0.1
                }
            }
        }
        for (LogiTask logiTask : logiTaskSet) {
            bestLambdaForLowerBound[logiTask.getId()] = 0.0;
        }
        iter = 0;
        while (flag) {// (flag && iter<1000) {
            iter++;
//			System.out.println("========iteration " + iter + ", incumbent bestRelaxObj="+bestRelaxObj+", bestObjUpBound="+bestObjUpBound +", bestIter="+bestIter+ " ============");
//			给定当前的lambda，求解lagrangian松弛问题 （最关键）
            prodObj = solveProdProblem();// 子问题求解
            logiObj = solveLogiProblem();// 子问题求解
            relaxObj = prodObj + logiObj;
            // 给定prodSeq得到两个算法下的可行目标值即上界
            double objUpBound1 = LagrangHeuristic();
            if (objUpBound1 < bestObjUpBound) {
                bestObjUpBound = objUpBound1;
            }
//			double objUpBound2 = LagrangHeuristic3(prodSeq);
//			if (objUpBound2<bestObjUpBound) {
//				bestObjUpBound=objUpBound2;
//			}

            // 另外两个启发式算法

            for (int i = 0; i < 3000; i++) {
                //对lambda扰动，让prodSeq有一定的随机性
//                double[] sum_lambda_wrt_prod = new double[para.numProd];
//                double[] denominator_wrt_prod = new double[para.numProd];
//                for (LogiTask logiTask : logiTaskSet) {// 遍历logiTask，按其对应的产品prod，将对应的乘子求和，注意prod的ID从0开始编号
//                    int corresProdId = logiTask.getProd().getId();
//                    sum_lambda_wrt_prod[corresProdId] = sum_lambda_wrt_prod[corresProdId] + logiTask.getLambda();//+(Math.random()*(1)-0.5);
//                }
//                for (int p = 0; p < para.numProd; p++) {
//                    denominator_wrt_prod[p] = (para.numMach - 1 + para.demands[p]) * para.cycTime[p] + para.moldChangeTime
//                            + para.getDistance(para.machSet.get(0), para.machSet.get(para.numMach - 1)) / para.lineSpeed;//+(Math.random()*(200)-100);//随机扰动
//                }
//                ArrayList<ProdRatio> prodList = new ArrayList<ProdRatio>();
//                for (int p = 0; p < para.numProd; p++) {
//                    ProdRatio prodRatio = new ProdRatio(p,
//                            (para.objValFactor - sum_lambda_wrt_prod[p]) / denominator_wrt_prod[p]);
//                    prodList.add(prodRatio);
//                }
//                Collections.sort(prodList, new ProdSortByWeightProcessTimeRatio());// 降序排列
//
//                int[] prodSeq_turbed=new int[para.numProd];
//                for (int p = 0; p < para.numProd; p++) {
//                    prodSeq_turbed[p] = prodList.get(p).getId();
//                }
                double temp = LagrangHeuristic3(prodSeq);
                if (temp < bestObjUpBound) {
                    bestObjUpBound = temp;
                    for (int ii = 0; ii < prodSeq.length; ii++) {
                        bestProdSeq[ii] = prodSeq[ii];
                    }

                    bestLogiSeq.clear();
                    for (int ii = 0; ii < logiSeq.size(); ii++) {
                        ArrayList<Integer> alist = new ArrayList<Integer>();
                        for (int jj = 0; jj < logiSeq.get(ii).size(); jj++) {
                            alist.add(logiSeq.get(ii).get(jj));
                        }
                        bestLogiSeq.add(alist);
                    }
                }
            }

//            double temp = LocalSearch(prodSeq);
//            if (temp < bestObjUpBound) {
//                bestObjUpBound = temp;
//            }

            if (relaxObj > bestRelaxObj) {// 本次迭代松弛问题的目标值有改进(提升)
                improved = true;
                bestRelaxObj = relaxObj;
                for (LogiTask logiTask : logiTaskSet) {
                    bestLambdaForLowerBound[logiTask.getId()] = logiTask.getLambda();
                    bestTssForLowerBound[logiTask.getId()] = logiTask.getTss();
                    bestTadForLowerBound[logiTask.getId()] = logiTask.getTad();
                }
                iter_no_improve = 0;
            } else {
                improved = false;
                iter_no_improve++;
//				if (iter>10) {
//				for (LogiTask logiTask:logiTaskSet) {
//					logiTask.setLambda(bestLambdaForLowerBound[logiTask.getId()]);
//					logiTask.setTss(bestTssForLowerBound[logiTask.getId()]);
//					logiTask.setTad(bestTadForLowerBound[logiTask.getId()]);
//				}
//				}
            }
//          获得梯度
            updateSubgradient();// update lambda with dualSolAndObj.solution
//          更新步长
            updateStepsize();
            // 输出信息
            if (show) {
                InfoOutput();
            }
//			更新乘子
            updateLambda();
            if (iter_no_improve > max_iter_no_improve) {
                flag = false;
            }
            if (iter > 1 && bestObjUpBound <= bestRelaxObj) {
//				bestRelaxObj = bestObjUpBound;
                flag = false;
            }
//			if (iter == 1 && bestRelaxObj - bestObjUpBound > 0.005)
            if (bestRelaxObj - bestObjUpBound > 0.005) {
                System.out.println("wrong");
            }
        }

    }

    public void run3(double obj_fcff, int[] prodSeq_fcff, ArrayList<ArrayList<Integer>> logiSeq_fcff) throws ClassNotFoundException, IOException, IloException {
    	UBNoImprove=obj_fcff;
    	bestObjUpBound=obj_fcff;
    	bestProdSeq=prodSeq_fcff.clone();
    	bestLogiSeq=listCopy2(logiSeq_fcff);
    	boolean flag = true;
        bestLambdaForLowerBound = new double[logiTaskSet.size()];
        double[] bestTssForLowerBound = new double[logiTaskSet.size()];
        double[] bestTadForLowerBound = new double[logiTaskSet.size()];
        // 为每一个logiTask设置其lambda（初始化拉格朗日乘子）
        for (int p = 0; p < para.numProd; p++) {// product从0开始编号
            for (int m = 0; m < para.numMach; m++) {// machine从0开始编号
                for (int k = 0; k < para.K[p][m]; k++) {// batch从0开始编号
                    PMKtoLogiTask.get(p + "_" + m + "_" + k).setLambda(Math.random()*0.001);// (double)Math.random()*0.1
                }
            }
        }
        for (LogiTask logiTask : logiTaskSet) {
            bestLambdaForLowerBound[logiTask.getId()] = 0.0;
        }
        iter = 0;
        while (flag) {// (flag && iter<1000) {
            iter++;
//			System.out.println("========iteration " + iter + ", incumbent bestRelaxObj="+bestRelaxObj+", bestObjUpBound="+bestObjUpBound +", bestIter="+bestIter+ " ============");
//			给定当前的lambda，求解lagrangian松弛问题 （最关键）
            prodObj = solveProdProblem();// 子问题求解
            logiObj = solveLogiProblem();// 子问题求解
            relaxObj = prodObj + logiObj;
            double prodObjNoImprove=solveProdProblemNoLift();
            LBNoLift=Math.max(LBNoLift,prodObjNoImprove+logiObj);
            // 给定prodSeq得到两个算法下的可行目标值即上界
            double objUpBound1 = LagrangHeuristic();
            if (objUpBound1 < bestObjUpBound) {
                bestObjUpBound = objUpBound1;
            }
            UBNoImprove=Math.min(objUpBound1,UBNoImprove);
//			double objUpBound2 = LagrangHeuristic3(prodSeq);
//			if (objUpBound2<bestObjUpBound) {
//				bestObjUpBound=objUpBound2;
//			}
            // 另外两个启发式算法
//            for (int i = 0; i < 3000; i++) {
//                double temp = LagrangHeuristic3(prodSeq);
//                if (temp < bestObjUpBound) {
//                    bestObjUpBound = temp;
//                    for (int ii = 0; ii < prodSeq.length; ii++) {
//                        bestProdSeq[ii] = prodSeq[ii];
//                    }
//
//                    bestLogiSeq.clear();
//                    for (int ii = 0; ii < logiSeq.size(); ii++) {
//                        ArrayList<Integer> alist = new ArrayList<Integer>();
//                        for (int jj = 0; jj < logiSeq.get(ii).size(); jj++) {
//                            alist.add(logiSeq.get(ii).get(jj));
//                        }
//                        bestLogiSeq.add(alist);
//                    }
//                }
//            }
            //对lambda扰动，让prodSeq有一定的随机性
//            double[] sum_lambda_wrt_prod = new double[para.numProd];
//            double[] denominator_wrt_prod = new double[para.numProd];
//            for (LogiTask logiTask : logiTaskSet) {// 遍历logiTask，按其对应的产品prod，将对应的乘子求和，注意prod的ID从0开始编号
//                int corresProdId = logiTask.getProd().getId();
//                sum_lambda_wrt_prod[corresProdId] = sum_lambda_wrt_prod[corresProdId] + logiTask.getLambda();//+(Math.random()*(1)-0.5);
//            }
//            for (int p = 0; p < para.numProd; p++) {
//                denominator_wrt_prod[p] = (para.numMach - 1 + para.demands[p]) * para.cycTime[p] + para.moldChangeTime
//                        + para.getDistance(para.machSet.get(0), para.machSet.get(para.numMach - 1)) / para.lineSpeed+(Math.random()*(100)-50);//随机扰动
//            }
//            ArrayList<ProdRatio> prodList = new ArrayList<ProdRatio>();
//            for (int p = 0; p < para.numProd; p++) {
//                ProdRatio prodRatio = new ProdRatio(p,
//                        (para.objValFactor - sum_lambda_wrt_prod[p]) / denominator_wrt_prod[p]);
//                prodList.add(prodRatio);
//            }
//            Collections.sort(prodList, new ProdSortByWeightProcessTimeRatio());// 降序排列
//
//            int[] prodSeq_turbed=new int[para.numProd];
//            for (int p = 0; p < para.numProd; p++) {
//                prodSeq_turbed[p] = prodList.get(p).getId();
//            }

            double temp = LocalSearch(prodSeq,logiSeq);
            if (temp < bestObjUpBound) {
                bestObjUpBound = temp;
//                for (int ii = 0; ii < prodSeq.length; ii++) {
//                    bestProdSeq[ii] = prodSeq[ii];
//                }
//
//                bestLogiSeq.clear();
//                for (int ii = 0; ii < logiSeq.size(); ii++) {
//                    ArrayList<Integer> alist = new ArrayList<Integer>();
//                    for (int jj = 0; jj < logiSeq.get(ii).size(); jj++) {
//                        alist.add(logiSeq.get(ii).get(jj));
//                    }
//                    bestLogiSeq.add(alist);
//                }
            }
            if (relaxObj > bestRelaxObj) {// 本次迭代松弛问题的目标值有改进(提升)
                improved = true;
                bestRelaxObj = relaxObj;
                for (LogiTask logiTask : logiTaskSet) {
                    bestLambdaForLowerBound[logiTask.getId()] = logiTask.getLambda();
                    bestTssForLowerBound[logiTask.getId()] = logiTask.getTss();
                    bestTadForLowerBound[logiTask.getId()] = logiTask.getTad();
                }
                iter_no_improve = 0;
            } else {
                improved = false;
                iter_no_improve++;
//				if (iter>10) {
//				for (LogiTask logiTask:logiTaskSet) {
//					logiTask.setLambda(bestLambdaForLowerBound[logiTask.getId()]);
//					logiTask.setTss(bestTssForLowerBound[logiTask.getId()]);
//					logiTask.setTad(bestTadForLowerBound[logiTask.getId()]);
//				}
//				}
            }
//          获得梯度
            updateSubgradient();// update lambda with dualSolAndObj.solution
//          更新步长
            updateStepsize();
            // 输出信息
            if (show) {
                InfoOutput();
            }
//			更新乘子
            updateLambda();
            if (iter_no_improve > max_iter_no_improve) {
                flag = false;
            }
            if (iter > 1 && bestObjUpBound <= bestRelaxObj) {
                bestRelaxObj = bestObjUpBound;
                flag = false;
            }
            if (iter == 1 && bestRelaxObj - bestObjUpBound > 0.005) {
                System.out.println("wrong");
            }
//            System.out.println("iter="+iter+", LB="+bestRelaxObj+", UB="+bestObjUpBound);
        }
//        System.out.println(tafLB[para.numProd-1]<para.prodSet.get(prodSeq[para.numProd-1]).getTaf());
//        System.out.println(para.prodSet.get(prodSeq[para.numProd-1]).getTaf());
//        System.out.println(" ");
    }

    public void InfoOutput() throws IloException {
        System.out.println(">>>>>>>>>>>>>生产子问题信息<<<<<<<<<<<<");
        // 子问题信息输出
        System.out.println("production sequence:");
        for (int p = 0; p < para.numProd; p++) {
            System.out.print("-> " + prodSeq[p] + " ");
        }
        System.out.println();
        System.out.printf("%-16s\t", "ProdID");
        for (int p = 0; p < para.prodSet.size(); p++) {
            System.out.printf("%10d", p);
        }
        System.out.println();

        System.out.printf("%-16s\t", "taf");
        for (int p = 0; p < para.prodSet.size(); p++) {
            System.out.printf("%10.2f", para.prodSet.get(p).getTaf());
        }
        System.out.println();

//		System.out.printf("%-16s\t", "UB_Prod");
//		for (int p = 0; p < para.prodSet.size(); p++) {
//			System.out.printf("%10.2e", stepsizeUpperBound_Prod[p]);
//		}
//		System.out.println();

        System.out.println("prodObj=" + prodObj);

//		ProdProblem prodProblem = new ProdProblem(para);
//		prodProblem.buildModel();
//		prodProblem.getSolution();
//		System.out.println("optimal total completion time derived from MIP model: " + prodProblem.objVal);

        System.out.println(">>>>>>>>>>>>>物流子问题信息<<<<<<<<<<<<");
        // 子问题信息输出
        System.out.println("logistics sequence:");
//		for (AgvPath agvPath : optAgvPaths) {
//			System.out.print(agvPath);
//        }

        System.out.printf("%-16s\t", "LogiID");
        for (int i = 0; i < para.logiTaskSet.size(); i++) {
            System.out.printf("%10d", i);
        }
        System.out.println();

        System.out.printf("%-16s\t", "lambda");
        for (int i = 0; i < para.logiTaskSet.size(); i++) {
            System.out.printf("%10.4f", para.ItoLogiTask.get(i).getLambda());
        }
        System.out.println();

        System.out.printf("%-16s\t", "tad");
        for (int i = 0; i < para.logiTaskSet.size(); i++) {
            System.out.printf("%10.2f", para.ItoLogiTask.get(i).getTad());
        }
        System.out.println();

        System.out.printf("%-16s\t", "tss");
        for (int i = 0; i < para.logiTaskSet.size(); i++) {
            System.out.printf("%10.2f", para.ItoLogiTask.get(i).getTss());
        }
        System.out.println();

        System.out.printf("%-16s\t", "subgrad");
        for (int i = 0; i < para.logiTaskSet.size(); i++) {
            System.out.printf("%10.2f", para.ItoLogiTask.get(i).getSubgrad());
        }
        System.out.println();

//		System.out.printf("%-16s\t", "UB_Logi");
//		for (int i = 0; i < para.logiTaskSet.size(); i++) {
//			if (para.ItoLogiTask.get(i).getTad() < para.ItoLogiTask.get(i).getTss()) {// 该下标集合内各元素对应的lambda会减小
//				// 较小的步长以确保下一次的lambda不能为负数
//				System.out.printf("%10.4f", para.ItoLogiTask.get(i).getLambda()
//						/ (para.ItoLogiTask.get(i).getTss() - para.ItoLogiTask.get(i).getTad()));
//			} else {
//				System.out.printf("%10s", "Inf");
//			}
//		}
        System.out.println();
        System.out.println("logiObj=" + logiObj);
        // 松弛问题信息输出
        System.out.println(">>>>>>>>>>>>>松弛问题信息<<<<<<<<<<<<");
        System.out.println("stepsizeFactor=" + stepsizeFactor + ", iter=" + iter + ", bestRelaxObj=" + bestRelaxObj
                + ", bestObjUpBound=" + bestObjUpBound + ", stepsize=" + stepsize);
    }

//	private void solveLagrangRelaxation() throws ClassNotFoundException, IOException, IloException{
////		次梯度算法的核心方法，给定每个logiTask对应的lagrangian乘子，原问题松弛后转变为一个生产子问题和一个物流子问题。
//		//原问题的求解转为了两个子问题的各自求解
//		prodObj=solveProdProblem();//子问题求解
//		
//		logiObj=solveLogiProblem();//子问题求解
//		relaxObj = prodObj + logiObj;
//		
//	}

    private double solveProdProblem() throws IloException {
//		解松弛后的生产子问题，更新logiTask的tss
//		System.out.println("<<<<<<<<<<< SWPT rule For Production Problem >>>>>>>>>>>\n");

        double[] sum_lambda_wrt_prod = new double[para.numProd];
        double[] denominator_wrt_prod = new double[para.numProd];
        double prodObj = 0;
        for (LogiTask logiTask : logiTaskSet) {// 遍历logiTask，按其对应的产品prod，将对应的乘子求和，注意prod的ID从0开始编号
            int corresProdId = logiTask.getProd().getId();
            sum_lambda_wrt_prod[corresProdId] = sum_lambda_wrt_prod[corresProdId] + logiTask.getLambda();
        }
        for (int p = 0; p < para.numProd; p++) {
            denominator_wrt_prod[p] = (para.numMach - 1 + para.demands[p]) * para.cycTime[p] + para.moldChangeTime
                    + para.getDistance(para.machSet.get(0), para.machSet.get(para.numMach - 1)) / para.lineSpeed;
        }
        ArrayList<ProdRatio> prodList = new ArrayList<ProdRatio>();
        for (int p = 0; p < para.numProd; p++) {
            ProdRatio prodRatio = new ProdRatio(p,
                    (para.objValFactor - sum_lambda_wrt_prod[p]) / denominator_wrt_prod[p]);
            prodList.add(prodRatio);
        }
        Collections.sort(prodList, new ProdSortByWeightProcessTimeRatio());// 降序排列


        for (int p = 0; p < para.numProd; p++) {
            prodSeq[p] = prodList.get(p).getId();
        }

        String key;
        double timeCur = minTard;//Math.max(minTard,tardLB[prodSeq[0]]);// minTard;//para.moldChangeTime;//earliestReadyTime;// 时间游标
//        prodObj = 0;
        double[] tempTaf = new double[para.numProd];
        for (int p = 0; p < para.numProd; p++) {// product从0开始编号
//			if (p>0) {
//				timeCur += tardLB[prodSeq[p]];
//			}
            for (int m = 0; m < para.numMach; m++) {// machine从0开始编号
                if (m > 0) {// 自第一个机器后，最早开工时间需要做一次减法
                    timeCur = timeCur - para.prodSet.get(prodSeq[p]).getCycleTime() * (para.demands[prodSeq[p]] - 1)
                            + para.getDistance(para.machSet.get(m - 1), para.machSet.get(m)) / para.lineSpeed;// (para.jobshopLength/para.numMach)/para.lineSpeed;
                }
                for (int k = 0; k < para.K[prodSeq[p]][m]; k++) {// batch从0开始编号
                    key = prodSeq[p] + "_" + m + "_" + k;
                    // 通过key获得PMKtoLogiTask.get(key)是p,m,k对应的运输任务
                    timeCur = timeCur + PMKtoLogiTask.get(key).getAssemTime(); // 加上装配时间等于该批次的完工时间，也是下一批次的最早开工时间；
                    if (m == para.numMach - 1 && k == para.K[prodSeq[p]][m] - 1) {
                        timeCur=Math.max( timeCur, tafLB[p] );
                        para.prodSet.get(prodSeq[p]).setTaf(timeCur);
                        tempTaf[prodSeq[p]] = timeCur;
                    }
                }
            }
//            prodObj = prodObj + timeCur * para.objValFactor;
            timeCur = timeCur - (para.numMach - 1) * para.cycTime[prodSeq[p]] + para.moldChangeTime
                    - para.getDistance(para.machSet.get(0), para.machSet.get(para.numMach - 1)) / para.lineSpeed;// -para.jobshopLength/para.lineSpeed;
        }

//        for (int p = 0; p < para.numProd; p++) {
//            double temp = para.prodSet.get(prodList.get(prodSeq[p]).getId()).getTaf();
//            if (temp < tafLB[p]) {
//                para.prodSet.get(prodList.get(prodSeq[p]).getId()).setTaf(tafLB[p]);
//            }
//            tempTaf[prodSeq[p]]=Math.max(temp,tafLB[p]);
//        }

        double[] latestTaf = new double[para.numProd];
        for (int p = 0; p < para.numProd; p++) {
            latestTaf[p] = tempTaf[p] + maxTaf - tempTaf[prodSeq[para.numProd - 1]];
        }

        for (int p = 0; p < para.numProd; p++) {
            if (prodList.get(p).getRatio() < 0)//只有总权重小于0才会令其Taf尽可能大
            {
                para.prodSet.get(prodList.get(p).getId()).setTaf(latestTaf[prodList.get(p).getId()]);
            }
        }


        // 基于新的taf更新tss
        for (int p = 0; p < para.numProd; p++) {// product从0开始编号
            for (int m = 0; m < para.numMach; m++) {// machine从0开始编号
                for (int k = 0; k < para.K[prodSeq[p]][m]; k++) {// batch从0开始编号
                    double tss = para.prodSet.get(prodSeq[p]).getTaf()
                            + ((k + 1 - 1) * para.quanProdMach[prodSeq[p]][m] - para.numMach + (m + 1)
                            - para.demands[prodSeq[p]]) * para.cycTime[prodSeq[p]]
                            - para.getDistance(para.machSet.get(m), para.machSet.get(para.numMach - 1))
                            / para.lineSpeed;
                    PMKtoLogiTask.get(prodSeq[p] + "_" + m + "_" + k).setTss(tss);
                }
            }
        }
        // 算最终的目标函数
        prodObj = 0;
        for (int p = 0; p < para.numProd; p++) {// product从0开始编号
            prodObj = prodObj + para.prodSet.get(prodSeq[p]).getTaf();
            for (int m = 0; m < para.numMach; m++) {// machine从0开始编号
                for (int k = 0; k < para.K[prodSeq[p]][m]; k++) {// batch从0开始编号
                    String key2 = prodSeq[p] + "_" + m + "_" + k;
                    // 通过key获得PMKtoLogiTask.get(key)是p,m,k对应的运输任务
                    prodObj = prodObj - (PMKtoLogiTask.get(key2).getTss()) * PMKtoLogiTask.get(key2).getLambda();
                }
            }
        }

        return prodObj;
    }

    private double solveProdProblemNoLift() throws IloException {
//		解松弛后的生产子问题，更新logiTask的tss
//		System.out.println("<<<<<<<<<<< SWPT rule For Production Problem >>>>>>>>>>>\n");

        double[] sum_lambda_wrt_prod = new double[para.numProd];
        double[] denominator_wrt_prod = new double[para.numProd];
        double prodObj = 0;
        for (LogiTask logiTask : logiTaskSet) {// 遍历logiTask，按其对应的产品prod，将对应的乘子求和，注意prod的ID从0开始编号
            int corresProdId = logiTask.getProd().getId();
            sum_lambda_wrt_prod[corresProdId] = sum_lambda_wrt_prod[corresProdId] + logiTask.getLambda();
        }
        for (int p = 0; p < para.numProd; p++) {
            denominator_wrt_prod[p] = (para.numMach - 1 + para.demands[p]) * para.cycTime[p] + para.moldChangeTime
                    + para.getDistance(para.machSet.get(0), para.machSet.get(para.numMach - 1)) / para.lineSpeed;
        }
        ArrayList<ProdRatio> prodList = new ArrayList<ProdRatio>();
        for (int p = 0; p < para.numProd; p++) {
            ProdRatio prodRatio = new ProdRatio(p,
                    (para.objValFactor - sum_lambda_wrt_prod[p]) / denominator_wrt_prod[p]);
            prodList.add(prodRatio);
        }
        Collections.sort(prodList, new ProdSortByWeightProcessTimeRatio());// 降序排列


        for (int p = 0; p < para.numProd; p++) {
            prodSeq[p] = prodList.get(p).getId();
        }

        String key;
        double timeCur = 0;//Math.max(minTard,tardLB[prodSeq[0]]);// minTard;//para.moldChangeTime;//earliestReadyTime;// 时间游标
//        prodObj = 0;
        double[] tempTaf = new double[para.numProd];
        for (int p = 0; p < para.numProd; p++) {// product从0开始编号
//			if (p>0) {
//				timeCur += tardLB[prodSeq[p]];
//			}
            for (int m = 0; m < para.numMach; m++) {// machine从0开始编号
                if (m > 0) {// 自第一个机器后，最早开工时间需要做一次减法
                    timeCur = timeCur - para.prodSet.get(prodSeq[p]).getCycleTime() * (para.demands[prodSeq[p]] - 1)
                            + para.getDistance(para.machSet.get(m - 1), para.machSet.get(m)) / para.lineSpeed;// (para.jobshopLength/para.numMach)/para.lineSpeed;
                }
                for (int k = 0; k < para.K[prodSeq[p]][m]; k++) {// batch从0开始编号
                    key = prodSeq[p] + "_" + m + "_" + k;
                    // 通过key获得PMKtoLogiTask.get(key)是p,m,k对应的运输任务
                    timeCur = timeCur + PMKtoLogiTask.get(key).getAssemTime(); // 加上装配时间等于该批次的完工时间，也是下一批次的最早开工时间；
                    if (m == para.numMach - 1 && k == para.K[prodSeq[p]][m] - 1) {
                        para.prodSet.get(prodSeq[p]).setTaf(timeCur);
                        tempTaf[prodSeq[p]] = timeCur;
                    }
                }
            }
//            prodObj = prodObj + timeCur * para.objValFactor;
            timeCur = timeCur - (para.numMach - 1) * para.cycTime[prodSeq[p]] + para.moldChangeTime
                    - para.getDistance(para.machSet.get(0), para.machSet.get(para.numMach - 1)) / para.lineSpeed;// -para.jobshopLength/para.lineSpeed;
        }

//        for (int p = 0; p < para.numProd; p++) {
//            double temp = para.prodSet.get(prodList.get(prodSeq[p]).getId()).getTaf();
//            if (temp < tafLB[p]) {
//                para.prodSet.get(prodList.get(prodSeq[p]).getId()).setTaf(tafLB[p]);
//            }
//            tempTaf[prodSeq[p]]=Math.max(temp,tafLB[p]);
//        }

        double[] latestTaf = new double[para.numProd];
        for (int p = 0; p < para.numProd; p++) {
            latestTaf[p] = tempTaf[p] + maxTaf - tempTaf[prodSeq[para.numProd - 1]];
        }

        for (int p = 0; p < para.numProd; p++) {
            if (prodList.get(p).getRatio() < 0)//只有总权重小于0才会令其Taf尽可能大
            {
                para.prodSet.get(prodList.get(p).getId()).setTaf(latestTaf[prodList.get(p).getId()]);
            }
        }


        // 基于新的taf更新tss
        for (int p = 0; p < para.numProd; p++) {// product从0开始编号
            for (int m = 0; m < para.numMach; m++) {// machine从0开始编号
                for (int k = 0; k < para.K[prodSeq[p]][m]; k++) {// batch从0开始编号
                    double tss = para.prodSet.get(prodSeq[p]).getTaf()
                            + ((k + 1 - 1) * para.quanProdMach[prodSeq[p]][m] - para.numMach + (m + 1)
                            - para.demands[prodSeq[p]]) * para.cycTime[prodSeq[p]]
                            - para.getDistance(para.machSet.get(m), para.machSet.get(para.numMach - 1))
                            / para.lineSpeed;
                    PMKtoLogiTask.get(prodSeq[p] + "_" + m + "_" + k).setTss(tss);
                }
            }
        }
        // 算最终的目标函数
        prodObj = 0;
        for (int p = 0; p < para.numProd; p++) {// product从0开始编号
            prodObj = prodObj + para.prodSet.get(prodSeq[p]).getTaf();
            for (int m = 0; m < para.numMach; m++) {// machine从0开始编号
                for (int k = 0; k < para.K[prodSeq[p]][m]; k++) {// batch从0开始编号
                    String key2 = prodSeq[p] + "_" + m + "_" + k;
                    // 通过key获得PMKtoLogiTask.get(key)是p,m,k对应的运输任务
                    prodObj = prodObj - (PMKtoLogiTask.get(key2).getTss()) * PMKtoLogiTask.get(key2).getLambda();
                }
            }
        }

        return prodObj;
    }


    private double solveLogiProblem() throws ClassNotFoundException, IOException, IloException {
//		解物流子问题，更新logiTask的tad
//		System.out.println("\n<<<<<<<<<<< Column Generation For Logistics Problem >>>>>>>>>>>");
        ColGen colGen = new ColGen(para);
        colGen.run();
        double logiObj = colGen.getObj();// 列生成获得的obj是不带 sum_i (lambda_i*TransTime_(0,i))的
//
        optAgvPaths = colGen.rstMasterProblem.getOptAgvPaths();// ArrayList<AgvPath>类型

        // 根据物流顺序返回tad的值
        double[] tad = new double[logiTaskSet.size()];
        for (int i = 0; i < logiTaskSet.size(); i++) {
            tad[i] = Float.MAX_VALUE;
        }
        for (AgvPath agvPath : optAgvPaths) {
            double baseTime = 0;
            for (LogiTask logiTask : agvPath.getExecuteSequence()) {
                double temp = baseTime + logiTask.getChargingTime() + logiTask.getTransTime() / 2;
                if (temp < tad[logiTask.getId()]) {
                    tad[logiTask.getId()] = temp;
                }
//				logiTask.setTad(baseTime +logiTask.getChargingTime()+ logiTask.getTransTime() / 2);
                baseTime = baseTime + logiTask.getProcessTime();
            }
        }
        for (LogiTask logiTask : logiTaskSet) {
            logiTask.setTad(tad[logiTask.getId()]);
        }
        for (LogiTask logiTask : logiTaskSet) {// 列生成获得的obj基础上减去sum_i (lambda_i*TransTime_(0,i))才是物流子问题的最终目标值
            logiObj = logiObj - logiTask.getLambda() * (logiTask.getTransTime() / 2);
        }

        logiSeq.clear();
        for (int r = 0; r < para.numVeh; r++) {
            logiSeq.add(new ArrayList<Integer>());
        }
        for (int r = 0; r < optAgvPaths.size(); r++) {
            AgvPath agvPath=optAgvPaths.get(r);
            for (LogiTask logiTask : agvPath.getExecuteSequence()) {
                logiSeq.get(r).add(logiTask.getId());
            }
        }

        return logiObj;

    }

    public double LagrangHeuristic() {
        String key;
        double timeCur;// 时间游标
        double sum_taf = 0;
        double iniMaxTaf = 0;
        key = prodSeq[0] + "_" + 0 + "_" + 0;
        timeCur = Math.max(PMKtoLogiTask.get(key).getTss(), para.moldChangeTime);// 第1个产品第1个机器第1个批次的物料的开始时间（生产问题求得），是整个系统的最早开始时间
        PMKtoLogiTask.get(key).tssHeur_Earl = Math.max(timeCur, PMKtoLogiTask.get(key).getTad());
        for (int p = 0; p < para.numProd; p++) {// product从0开始编号
            for (int m = 0; m < para.numMach; m++) {// machine从0开始编号
                if (m > 0) {// 自第一个机器后，最早开工时间需要做一次减法
                    timeCur = timeCur - para.prodSet.get(prodSeq[p]).getCycleTime() * (para.demands[prodSeq[p]] - 1)
                            + para.getDistance(para.machSet.get(m - 1), para.machSet.get(m)) / para.lineSpeed;// +(para.jobshopLength/para.numMach)/para.lineSpeed;
                }
                for (int k = 0; k < para.K[prodSeq[p]][m]; k++) {// batch从0开始编号
                    key = prodSeq[p] + "_" + m + "_" + k;
                    // 通过key获得PMKtoLogiTask.get(key)是p,m,k对应的运输任务
                    timeCur = Math.max(timeCur, PMKtoLogiTask.get(key).getTad()); // 真正的开始时间是物料到达和最早开工时间取大；
                    PMKtoLogiTask.get(key).tssHeur_Earl = timeCur;
                    timeCur = timeCur + PMKtoLogiTask.get(key).getAssemTime(); // 加上装配时间等于该批次的完工时间，也是下一批次的最早开工时间；
                    if (m == para.numMach - 1 && k == para.K[prodSeq[p]][m] - 1) {
                        para.prodSet.get(prodSeq[p]).tafHeur_Earl = (timeCur);
                    }
                }
            }
//			if (p==para.numProd-1) {
//				iniMaxTaf=timeCur;
//			}
//			timeCur=timeCur+(para.jobshopLength/para.numMach)/para.lineSpeed;
            sum_taf = sum_taf + timeCur * para.objValFactor;
            timeCur = timeCur - (para.numMach - 1) * para.cycTime[prodSeq[p]] + para.moldChangeTime
                    - para.getDistance(para.machSet.get(0), para.machSet.get(para.numMach - 1)) / para.lineSpeed;// -para.jobshopLength/para.lineSpeed;
        }
//		for (int p = 0; p < para.numProd; p++) {
//			tafUB[p]=tafUB[p]+iniMaxTaf-tafUB[para.numProd-1];
//		}
        return sum_taf;
    }

    public double LagrangHeuristic2(int[] prodSeq) {
//		if(Math.random()<0.1) {
//			Random rdm= new Random();
//			int idx1=rdm.nextInt(prodSeq.length);
//			int idx2=rdm.nextInt(prodSeq.length);
//			int temp=prodSeq[idx1];
//			prodSeq[idx1]=prodSeq[idx2];
//			prodSeq[idx2]=temp;
//		}

        logiSeq.clear();
        double[] timeCur = new double[para.numVeh];
        double v1 = 0;
        double v2 = 0;
        double v3 = 0;
//		ArrayList<ArrayList<Integer>> logiSeq= new ArrayList<ArrayList<Integer>>();
        ArrayList<ArrayList<Integer>> logiSeq1 = new ArrayList<ArrayList<Integer>>();
        ArrayList<ArrayList<Integer>> logiSeq2 = new ArrayList<ArrayList<Integer>>();
        ArrayList<ArrayList<Integer>> logiSeq3 = new ArrayList<ArrayList<Integer>>();
//		int[] preLogiTaskID=new int[para.numVeh];
        for (int r = 0; r < para.numVeh; r++) {
//			preLogiTaskID[r]=-1;
            logiSeq.add(new ArrayList<Integer>());
            logiSeq1.add(new ArrayList<Integer>());
            logiSeq2.add(new ArrayList<Integer>());
            logiSeq3.add(new ArrayList<Integer>());
        }
        for (int p = 0; p < prodSeq.length; p++) {// product从0开始编号
            int pindex = prodSeq[p];
            ArrayList<int[]> arrayofPMKI = para.PtoArrayofPMKI.get(pindex);
            ArrayList<Integer> arrayOfI = new ArrayList<Integer>();
            ArrayList<Integer> idSet = new ArrayList<Integer>();
            for (int index = 0; index < arrayofPMKI.size(); index++) {
                arrayOfI.add(arrayofPMKI.get(index)[3]);
            }

            while (!arrayOfI.isEmpty()) {
                Random random = new Random();
                int index = random.nextInt(arrayOfI.size());
                idSet.add(arrayOfI.get(index));
                arrayOfI.remove(index);
            }

            for (int index = 0; index < idSet.size(); index++) {
//			int machId = arrayofPMKI.get(index)[1];// 对应mach的ID
//			int batchId = arrayofPMKI.get(index)[2];// 对应batch的ID
                int logiTaskId = idSet.get(index);// 对应logiTask的ID
                int min_index = 0;
                if (Math.random() < prob) {
                    Random random = new Random();
                    min_index = random.nextInt(para.numVeh);

                    logiSeq1.get(min_index).add(logiTaskId);
                    tad[logiTaskId] = timeCur[min_index] + para.ItoLogiTask.get(logiTaskId).getChargingTime()
                            + para.ItoLogiTask.get(logiTaskId).getTransTime() / 2;
                    timeCur[min_index] = timeCur[min_index] + para.ItoLogiTask.get(logiTaskId).getProcessTime();
                } else {
                    double requiredDeliTime = para.ItoLogiTask.get(logiTaskId).getChargingTime()
                            + para.ItoLogiTask.get(logiTaskId).getTransTime() / 2;
                    double min_val = Float.MAX_VALUE;
                    for (int r = 0; r < para.numVeh; r++) {
                        if (timeCur[r] + requiredDeliTime < min_val) {
                            min_val = timeCur[r] + requiredDeliTime;
                            min_index = r;
                        }
                    }
                    logiSeq1.get(min_index).add(logiTaskId);
                    tad[logiTaskId] = timeCur[min_index] + para.ItoLogiTask.get(logiTaskId).getChargingTime()
                            + para.ItoLogiTask.get(logiTaskId).getTransTime() / 2;
                    timeCur[min_index] = timeCur[min_index] + para.ItoLogiTask.get(logiTaskId).getProcessTime();
                }
            }


            for (int index = 0; index < idSet.size(); index++) {
                int logiTaskId = idSet.get(index);// 对应logiTask的ID
                int min_index = 0;
                Random random = new Random();
                min_index = random.nextInt(para.numVeh);

                logiSeq2.get(min_index).add(logiTaskId);
                tad[logiTaskId] = timeCur[min_index] + para.ItoLogiTask.get(logiTaskId).getChargingTime()
                        + para.ItoLogiTask.get(logiTaskId).getTransTime() / 2;
                timeCur[min_index] = timeCur[min_index] + para.ItoLogiTask.get(logiTaskId).getProcessTime();
            }


            for (int index = 0; index < idSet.size(); index++) {
//					int machId = arrayofPMKI.get(index)[1];// 对应mach的ID
//					int batchId = arrayofPMKI.get(index)[2];// 对应batch的ID
                int logiTaskId = idSet.get(index);// 对应logiTask的ID
//						Random random2 = new Random();
//						int index2=random2.nextInt(para.numVeh);
//						logiSeq.get(index2).add(logiTaskId);
                int min_index = 0;
                double requiredDeliTime = para.ItoLogiTask.get(logiTaskId).getChargingTime()
                        + para.ItoLogiTask.get(logiTaskId).getTransTime() / 2;
                double min_val = Float.MAX_VALUE;
                for (int r = 0; r < para.numVeh; r++) {
                    if (timeCur[r] + requiredDeliTime < min_val) {
                        min_val = timeCur[r] + requiredDeliTime;
                        min_index = r;
                    }
                }
                logiSeq3.get(min_index).add(logiTaskId);
                tad[logiTaskId] = timeCur[min_index] + para.ItoLogiTask.get(logiTaskId).getChargingTime()
                        + para.ItoLogiTask.get(logiTaskId).getTransTime() / 2;
                timeCur[min_index] = timeCur[min_index] + para.ItoLogiTask.get(logiTaskId).getProcessTime();
            }

        }
        v1 = getObj(prodSeq, logiSeq1);
        v2 = getObj(prodSeq, logiSeq2);
        v3 = getObj(prodSeq, logiSeq3);

        double testObjVal = 0;
        if (v1 <= v2 && v1 <= v3) {
            logiSeq = logiSeq1;
        }
        if (v2 <= v1 && v2 <= v3) {
            logiSeq = logiSeq2;
        }
        if (v3 <= v1 && v3 <= v2) {
            logiSeq = logiSeq3;
        }

        testObjVal = getObj(prodSeq, logiSeq);
        return testObjVal;
    }

    public double LagrangHeuristic3(int[] prodSeq) {
        logiSeq.clear();
        double[] timeCur = new double[para.numVeh];
        double v1 = 0;
        double v2 = 0;
        double v3 = 0;
        int[][] posi1=new int[para.numProd][4];
        int[][] posi2=new int[para.numProd][4];
        int[][] posi3=new int[para.numProd][4];
//		ArrayList<ArrayList<Integer>> logiSeq= new ArrayList<ArrayList<Integer>>();
        ArrayList<ArrayList<Integer>> logiSeq1 = new ArrayList<ArrayList<Integer>>();
        ArrayList<ArrayList<Integer>> logiSeq2 = new ArrayList<ArrayList<Integer>>();
        ArrayList<ArrayList<Integer>> logiSeq3 = new ArrayList<ArrayList<Integer>>();
        ArrayList<Integer> first = new ArrayList<>();
        ArrayList<Integer> last = new ArrayList<>();
//		int[] preLogiTaskID=new int[para.numVeh];
        for (int r = 0; r < para.numVeh; r++) {
//			preLogiTaskID[r]=-1;
            logiSeq.add(new ArrayList<Integer>());
            logiSeq1.add(new ArrayList<Integer>());
            logiSeq2.add(new ArrayList<Integer>());
            logiSeq3.add(new ArrayList<Integer>());
        }
        for (int p = 0; p < prodSeq.length; p++) {// product从0开始编号
            int pindex = prodSeq[p];
            ArrayList<int[]> arrayofPMKI = para.PtoArrayofPMKI.get(pindex);
            ArrayList<Integer> arrayOfI = new ArrayList<Integer>();
            ArrayList<Integer> idSet = new ArrayList<Integer>();
            for (int index = 0; index < arrayofPMKI.size(); index++) {
                arrayOfI.add(arrayofPMKI.get(index)[3]);
            }
            first.add(arrayOfI.get(0));
            last.add(arrayOfI.get(arrayOfI.size()-1));

            while (!arrayOfI.isEmpty()) {
                Random random = new Random();
                int index = random.nextInt(arrayOfI.size());
                idSet.add(arrayOfI.get(index));
                arrayOfI.remove(index);
            }

            for (int index = 0; index < idSet.size(); index++) {
                int logiTaskId = idSet.get(index);// 对应logiTask的ID
                int min_index = 0;
                if (Math.random() < prob) {
                    Random random = new Random();
                    min_index = random.nextInt(para.numVeh);

                    logiSeq1.get(min_index).add(logiTaskId);
                    tad[logiTaskId] = timeCur[min_index] + para.ItoLogiTask.get(logiTaskId).getChargingTime()
                            + para.ItoLogiTask.get(logiTaskId).getTransTime() / 2;
                    timeCur[min_index] = timeCur[min_index] + para.ItoLogiTask.get(logiTaskId).getProcessTime();
                } else {
                    double requiredDeliTime = para.ItoLogiTask.get(logiTaskId).getChargingTime()
                            + para.ItoLogiTask.get(logiTaskId).getTransTime() / 2;
                    double min_val = Float.MAX_VALUE;
                    for (int r = 0; r < para.numVeh; r++) {
                        if (timeCur[r] + requiredDeliTime < min_val) {
                            min_val = timeCur[r] + requiredDeliTime;
                            min_index = r;
                        }
                    }
                    logiSeq1.get(min_index).add(logiTaskId);
                    if(first.contains(logiTaskId)){
                        posi1[pindex][0]=min_index;
                        posi1[pindex][1]=logiSeq1.get(min_index).size()-1;
                    }
                    if(last.contains(logiTaskId)){
                        posi1[pindex][2]=min_index;
                        posi1[pindex][3]=logiSeq1.get(min_index).size()-1;
                    }
                    tad[logiTaskId] = timeCur[min_index] + para.ItoLogiTask.get(logiTaskId).getChargingTime()
                            + para.ItoLogiTask.get(logiTaskId).getTransTime() / 2;
                    timeCur[min_index] = timeCur[min_index] + para.ItoLogiTask.get(logiTaskId).getProcessTime();
                }
            }

            for (int index = 0; index < idSet.size(); index++) {
                int logiTaskId = idSet.get(index);// 对应logiTask的ID
                int min_index = 0;
                Random random = new Random();
                min_index = random.nextInt(para.numVeh);

                logiSeq2.get(min_index).add(logiTaskId);
                if(first.contains(logiTaskId)){
                    posi2[pindex][0]=min_index;
                    posi2[pindex][1]=logiSeq2.get(min_index).size()-1;
                }
                if(last.contains(logiTaskId)){
                    posi2[pindex][2]=min_index;
                    posi2[pindex][3]=logiSeq2.get(min_index).size()-1;
                }
                tad[logiTaskId] = timeCur[min_index] + para.ItoLogiTask.get(logiTaskId).getChargingTime()
                        + para.ItoLogiTask.get(logiTaskId).getTransTime() / 2;
                timeCur[min_index] = timeCur[min_index] + para.ItoLogiTask.get(logiTaskId).getProcessTime();
            }

            for (int index = 0; index < idSet.size(); index++) {
                int logiTaskId = idSet.get(index);// 对应logiTask的ID
                int min_index = 0;
                double requiredDeliTime = para.ItoLogiTask.get(logiTaskId).getChargingTime()
                        + para.ItoLogiTask.get(logiTaskId).getTransTime() / 2;
                double min_val = Float.MAX_VALUE;
                for (int r = 0; r < para.numVeh; r++) {
                    if (timeCur[r] + requiredDeliTime < min_val) {
                        min_val = timeCur[r] + requiredDeliTime;
                        min_index = r;
                    }
                }
                logiSeq3.get(min_index).add(logiTaskId);
                if(first.contains(logiTaskId)){
                    posi3[pindex][0]=min_index;
                    posi3[pindex][1]=logiSeq3.get(min_index).size()-1;
                }
                if(last.contains(logiTaskId)){
                    posi3[pindex][2]=min_index;
                    posi3[pindex][3]=logiSeq3.get(min_index).size()-1;
                }
                tad[logiTaskId] = timeCur[min_index] + para.ItoLogiTask.get(logiTaskId).getChargingTime()
                        + para.ItoLogiTask.get(logiTaskId).getTransTime() / 2;
                timeCur[min_index] = timeCur[min_index] + para.ItoLogiTask.get(logiTaskId).getProcessTime();
            }
        }
        v1 = getObj(prodSeq, logiSeq1);
        v2 = getObj(prodSeq, logiSeq2);
        v3 = getObj(prodSeq, logiSeq3);

        int[][] posi=new int[para.numProd][4];
        double testObjVal = 0;
        if (v1 <= v2 && v1 <= v3) {
            logiSeq = logiSeq1;
            posi=posi1;
        }
        if (v2 <= v1 && v2 <= v3) {
            logiSeq = logiSeq2;
            posi=posi2;
        }
        if (v3 <= v1 && v3 <= v2) {
            logiSeq = logiSeq3;
            posi=posi3;
        }

        testObjVal = getObj(prodSeq, logiSeq);
        testObjVal=Math.min(testObjVal,LocalImprove(prodSeq, logiSeq,testObjVal,posi));
        return testObjVal;
    }

    public double LocalSearch_backup(int[] prodSeq) {
        double bestObjVal=Float.MAX_VALUE;

        ArrayList<ArrayList<Integer>> logiSeq=new ArrayList<ArrayList<Integer>>();
        for (int r = 0; r < para.numVeh; r++) {
            logiSeq.add(new ArrayList<Integer>());
        }
        double[] tad=new double[logiTaskSet.size()];
        double[] tss=new double[logiTaskSet.size()];
        double[] delta=new double[logiTaskSet.size()];
        double[] lambda=new double[logiTaskSet.size()];
        for (LogiTask logiTask:logiTaskSet) {
            lambda[logiTask.getId()]=logiTask.getLambda();
            tad[logiTask.getId()]=logiTask.getTad();
            tss[logiTask.getId()]=logiTask.getTss();
            delta[logiTask.getId()]=tad[logiTask.getId()]-tss[logiTask.getId()];
        }

        ArrayList<ArrayList<Integer>> iList=new ArrayList<ArrayList<Integer>>();
        for (int p = 0; p < prodSeq.length; p++) {
            iList.add(new ArrayList<Integer>());
        }

        for (int p = 0; p < prodSeq.length; p++) {
            ArrayList<TssWithoutLogiCons> tssWithoutLogiConsList = new ArrayList<TssWithoutLogiCons>();
            String key;
            double timeCur = 0;
            for (int m = 0; m < para.numMach; m++) {
                if (m > 0) {
                    timeCur = timeCur
                            - para.prodSet.get(prodSeq[p]).getCycleTime() * (para.demands[prodSeq[p]] - 1)
                            + para.getDistance(para.machSet.get(m - 1), para.machSet.get(m))
                            / para.lineSpeed;
                }
                for (int k = 0; k < para.K[prodSeq[p]][m]; k++) {
                    key = prodSeq[p] + "_" + m + "_" + k;
                    TssWithoutLogiCons tssWithoutLogiCons = new TssWithoutLogiCons(key, timeCur);
                    tssWithoutLogiConsList.add(tssWithoutLogiCons);
                    timeCur = timeCur + para.PMKtoLogiTask.get(key).getAssemTime(); // 加上装配时间等于该批次的完工时间，也是下一批次的最早开工时间；
                }
            }
            Collections.sort(tssWithoutLogiConsList, new LogiTaskSortByTss());
            for (int s = 0; s < tssWithoutLogiConsList.size(); s++) {
                String pmk = tssWithoutLogiConsList.get(s).pmk;
                iList.get(p).add(para.PMKtoLogiTask.get(pmk).getId());
            }
        }

        ArrayList<ArrayList<Integer>> cumb_iList=new ArrayList<ArrayList<Integer>>();
        ArrayList<ArrayList<Integer>> ini_iList=new ArrayList<ArrayList<Integer>>();
        for (int p = 0; p < prodSeq.length; p++) {
            cumb_iList.add(new ArrayList<Integer>());
            ini_iList.add(new ArrayList<Integer>());
        }
        for (int i = 0; i < iList.size(); i++) {
            for (int j = 0; j < iList.get(i).size(); j++) {
                int idx=iList.get(i).get(j);
                cumb_iList.get(i).add(idx);
                ini_iList.get(i).add(idx);
            }
        }
        int iter=0;
//        while (iter<1000){
//
//            ArrayList<ArrayList<Integer>> neigh_iList= neighbourOperator(cumb_iList,ini_iList);
//            logiSeq=routeGeneration(neigh_iList);
//            double testObjVal = getObj(prodSeq, logiSeq);
//            if(testObjVal<bestObjVal){
//                bestObjVal=testObjVal;
//                for (int i = 0; i < iList.size(); i++) {
//                    for (int j = 0; j < iList.get(i).size(); j++) {
//                        int idx=neigh_iList.get(i).get(j);
//                        cumb_iList.get(i).add(idx);
//                    }
//                }
//            }
//        }

        return bestObjVal;
    }

//    public double LocalSearch(int[] prodSeq) {
//
//        ArrayList<ArrayList<Integer>> logiPerProdList=new ArrayList<ArrayList<Integer>>();
//        for (int p = 0; p < prodSeq.length; p++) {
//            logiPerProdList.add(new ArrayList<Integer>());
//        }
//
//        for (int p = 0; p < prodSeq.length; p++) {
//            ArrayList<TssWithoutLogiCons> tssWithoutLogiConsList = new ArrayList<TssWithoutLogiCons>();
//            double timeCur = 0;
//            for (int m = 0; m < para.numMach; m++) {
//                if (m > 0) {
//                    timeCur = timeCur
//                            - para.prodSet.get(prodSeq[p]).getCycleTime() * (para.demands[prodSeq[p]] - 1)
//                            + para.getDistance(para.machSet.get(m - 1), para.machSet.get(m))
//                            / para.lineSpeed;
//                }
//                for (int k = 0; k < para.K[prodSeq[p]][m]; k++) {
//                    String key = prodSeq[p] + "_" + m + "_" + k;
//                    TssWithoutLogiCons tssWithoutLogiCons = new TssWithoutLogiCons(key, timeCur);
//                    tssWithoutLogiConsList.add(tssWithoutLogiCons);
//                    timeCur = timeCur + para.PMKtoLogiTask.get(key).getAssemTime(); // 加上装配时间等于该批次的完工时间，也是下一批次的最早开工时间；
//                }
//            }
//            Collections.sort(tssWithoutLogiConsList, new LogiTaskSortByTss());
//            for (int s = 0; s < tssWithoutLogiConsList.size(); s++) {
//                String pmk = tssWithoutLogiConsList.get(s).pmk;
//                logiPerProdList.get(p).add(para.PMKtoLogiTask.get(pmk).getId());
//            }
//        }
//
//        int[] firstP=new int[para.numProd];
//        int[] lastP=new int[para.numProd];
//        int[] preUP=new int[para.numProd];
//        int[] postUP=new int[para.numProd];
//        for (int p = 0; p < para.numProd; p++) {
//            firstP[p]=logiPerProdList.get(p).get(0);
//            lastP[p]=logiPerProdList.get(p).get( logiPerProdList.get(p).size()-1 );
//            if(p==0){
//                preUP[p]=0;
//                postUP[p]=logiPerProdList.get(p+1).size();
//            }else if(p==para.numProd-1){
//                preUP[p]=logiPerProdList.get(p-1).size();
//                postUP[p]=0;
//            }else {
//                postUP[p]=logiPerProdList.get(p+1).size();
//                preUP[p]=logiPerProdList.get(p-1).size();
//            }
//        }
//
//
//        ArrayList<ArrayList<Integer>> logiSeq=new ArrayList<ArrayList<Integer>>();
//        for (int r = 0; r < para.numVeh; r++) {
//            logiSeq.add(new ArrayList<Integer>());
//        }
//        double[] tad=new double[logiTaskSet.size()];
//        double[] tss=new double[logiTaskSet.size()];
//        double[] delta=new double[logiTaskSet.size()];
//        double[] lambda=new double[logiTaskSet.size()];
//        for (LogiTask logiTask:logiTaskSet) {
//            lambda[logiTask.getId()]=logiTask.getLambda();
//            tad[logiTask.getId()]=logiTask.getTad();
//            tss[logiTask.getId()]=logiTask.getTss();
//            delta[logiTask.getId()]=tad[logiTask.getId()]-tss[logiTask.getId()];
//        }
//
////        ArrayList<Integer> iList=new ArrayList<Integer>();
////        ArrayList<TssWithoutLogiCons> tssWithoutLogiConsList = new ArrayList<TssWithoutLogiCons>();
////        double timeCur = 0;
////        String key;
////        for (int p = 0; p < prodSeq.length; p++) {
////            for (int m = 0; m < para.numMach; m++) {
////                if (m > 0) {
////                    timeCur = timeCur
////                            - para.prodSet.get(prodSeq[p]).getCycleTime() * (para.demands[prodSeq[p]] - 1)
////                            + para.getDistance(para.machSet.get(m - 1), para.machSet.get(m))
////                            / para.lineSpeed;
////                }
////                for (int k = 0; k < para.K[prodSeq[p]][m]; k++) {
////                    key = prodSeq[p] + "_" + m + "_" + k;
////                    TssWithoutLogiCons tssWithoutLogiCons = new TssWithoutLogiCons(key, timeCur);
////                    tssWithoutLogiConsList.add(tssWithoutLogiCons);
////                    timeCur = timeCur + para.PMKtoLogiTask.get(key).getAssemTime(); // 加上装配时间等于该批次的完工时间，也是下一批次的最早开工时间；
////                }
////            }
////            timeCur = timeCur - (para.numMach - 1) * para.cycTime[prodSeq[p]] + para.moldChangeTime
////                    - para.getDistance(para.machSet.get(0), para.machSet.get(para.numMach - 1))
////                    / para.lineSpeed;// -para.jobshopLength/para.lineSpeed;
////        }
////        Collections.sort(tssWithoutLogiConsList, new LogiTaskSortByTss());
////        for (int s = 0; s < tssWithoutLogiConsList.size(); s++) {
////            String pmk = tssWithoutLogiConsList.get(s).pmk;
////            iList.add(para.PMKtoLogiTask.get(pmk).getId());
////        }
//
////        int[] iList2=new int[logiTaskSet.size()];
////        tssWithoutLogiConsList.clear();
////        for (LogiTask logiTask:logiTaskSet) {
////            key=logiTask.getProd().getId() + "_" + logiTask.getMach().getId() + "_" + logiTask.getBatchId();
////            tssWithoutLogiConsList.add(new TssWithoutLogiCons(key, logiTask.getLambda()));
////        }
////        Collections.sort(tssWithoutLogiConsList, new LogiTaskSortByTss());
////        for (int s = 0; s < tssWithoutLogiConsList.size(); s++) {
////            String pmk = tssWithoutLogiConsList.get(s).pmk;
////            iList2[s]=(para.PMKtoLogiTask.get(pmk).getId());
////        }
//
//        ArrayList<ArrayList<Integer>> idSetPerProd=new ArrayList<ArrayList<Integer>>();
//        for (int p = 0; p < para.numProd; p++) {
//            idSetPerProd.add(new ArrayList<Integer>());
//        }
//        for (LogiTask logiTask: para.logiTaskSet) {
//            idSetPerProd.get(logiTask.getProd().getId()).add(logiTask.getId());
//        }
//        ArrayList<Integer> iList=new ArrayList<Integer>();
//        for (int p = 0; p < para.numProd; p++) {
//            for (int i = 0; i < idSetPerProd.get(prodSeq[p]).size(); i++) {
//                iList.add(idSetPerProd.get(prodSeq[p]).get(i));
//            }
//        }
//
//        ArrayList<Integer> cumb_iList=listCopy(iList);
////        ArrayList<Integer> ini_iList=listCopy(iList);
//        double bestObjVal=Float.MAX_VALUE;
//        int iter=0;
//        while (iter<3000){
//            ArrayList<Integer> neigh_iList=listCopy(cumb_iList);
//            ArrayList<ArrayList<Integer>> idxPerProd=new ArrayList<ArrayList<Integer>>();
//            for (int i = 0; i < para.numProd; i++) {
//                idxPerProd.add(new ArrayList<Integer>());
//            }
//            for (int i = 0; i < neigh_iList.size(); i++) {
//                int p=para.ItoLogiTask.get(neigh_iList.get(i)).getProd().getId();
//                idxPerProd.get(p).add(i);
//            }
//            Random rnd = new Random();
//            int origPosition=rnd.nextInt(neigh_iList.size());
//            int rndLogiTaskId=neigh_iList.get(origPosition);
//            int rndP=para.ItoLogiTask.get(rndLogiTaskId).getProd().getId();
//            int posiLB,posiUB,newPosition;
////            if (Arrays.binarySearch(firstP,rndLogiTaskId)>=0){
////                posiLB=Math.max(origPosition-preUP[rndP],0);
////                ArrayList<Integer> candi=new ArrayList<Integer>();
////                for (int i = posiLB; i <origPosition ; i++) {
////                    candi.add(neigh_iList.get(i));
////                }
////                for (int i = 0; i < idxPerProd.get(rndP).size(); i++) {
////                    if (!candi.contains(idxPerProd.get(rndP).get(i))){
////                        candi.add(idxPerProd.get(rndP).get(i));
////                    }
////                }
////                newPosition=candi.get(rnd.nextInt(candi.size()));//rnd.nextInt(posiUB+1-posiLB)+posiLB;
////            } else if (Arrays.binarySearch(lastP,rndLogiTaskId)>=0) {
////                posiUB=Math.min(origPosition+postUP[rndP],neigh_iList.size()-1);
////                ArrayList<Integer> candi=new ArrayList<Integer>();
////                for (int i = origPosition; i <posiUB+1 ; i++) {
////                    candi.add(neigh_iList.get(i));
////                }
////                for (int i = 0; i < idxPerProd.get(rndP).size(); i++) {
////                    if (!candi.contains(idxPerProd.get(rndP).get(i))){
////                        candi.add(idxPerProd.get(rndP).get(i));
////                    }
////                }
////                newPosition=candi.get(rnd.nextInt(candi.size()));
////            }else {
////                newPosition=idxPerProd.get(rndP).get(rnd.nextInt(idxPerProd.get(rndP).size()));//rnd.nextInt(posiUB+1-posiLB)+posiLB;
////            }
////            Collections.swap(neigh_iList,origPosition,newPosition);
//            logiSeq=routeGeneration0(neigh_iList);
//            double testObjVal = getObj(prodSeq, logiSeq);
//            if(testObjVal<=bestObjVal){
//                bestObjVal=testObjVal;
//                cumb_iList=listCopy1(neigh_iList);
//            }
//            iter++;
//        }
//
//        return bestObjVal;
//    }

    public double LocalImprove(int[] prodSeq, ArrayList<ArrayList<Integer>> logiSeq, double obj, int[][] posi) {
        //posi[p][a]记录了产品p的第一个和最后一个logiTaskId所在的车辆a1和位置a2
        double new_obj=obj;

        ArrayList<Integer> prodPriority=new ArrayList<Integer>();
        for (int p = 0; p < prodSeq.length; p++) {
            prodPriority.add(prodSeq[p]);
        }
        for (int p = 0; p < prodSeq.length; p++) {
            int ovfirst = posi[p][0];
            int opfirst = posi[p][1];
            int ovlast = posi[p][2];
            int oplast = posi[p][3];
            if(prodSeq[0] == p){//产品p是第一个
                HashSet<int[]> insertSet = new HashSet<>();//相同优先级的插入（插入在指定位置之前）
                int o_LogiTaskId = logiSeq.get(ovlast).get(oplast);
                int o_ProdId = para.ItoLogiTask.get(o_LogiTaskId).getProd().getId();
                int o_Priority = prodPriority.indexOf(o_ProdId);
                for (int i = 0; i < logiSeq.size(); i++) {
                    for (int j = 0; j < logiSeq.get(i).size(); j++) {
                        int[] fea = {i, j};
                        int n_LogiTaskId = logiSeq.get(i).get(j);
                        int n_ProdId = para.ItoLogiTask.get(n_LogiTaskId).getProd().getId();
                        int n_Priority = prodPriority.indexOf(n_ProdId);
                        if (n_Priority == o_Priority+1) {
                            insertSet.add(fea);
                        }
                    }
                }
                Iterator<int[]> it = insertSet.iterator();
                while(it.hasNext()){
                    int[] fea=it.next();
                    ArrayList<ArrayList<Integer>> new_logiSeq = listCopy2(logiSeq);
                    insertOperator(new_logiSeq,ovlast,oplast,fea[0],fea[1]);
//                    it.remove();
                    double temp = getObj(prodSeq,new_logiSeq);
                    if (temp < new_obj){
                        new_obj = temp;
                        if(new_obj<bestObjUpBound) {
                			bestObjUpBound=new_obj;
                			bestProdSeq=prodSeq.clone();
                			bestLogiSeq=listCopy2(new_logiSeq);
                        }
                    }
                }
            } else if (prodSeq[prodSeq.length-1] == p) {//产品p是最后一个
                HashSet<int[]> insertSet = new HashSet<>();//相同优先级的插入（插入在指定位置之前）
                int o_LogiTaskId = logiSeq.get(ovfirst).get(opfirst);
                int o_ProdId = para.ItoLogiTask.get(o_LogiTaskId).getProd().getId();
                int o_Priority = prodPriority.indexOf(o_ProdId);
                for (int i = 0; i < logiSeq.size(); i++) {
                    for (int j = 0; j < logiSeq.get(i).size(); j++) {
                        int[] fea = {i, j};
                        int n_LogiTaskId = logiSeq.get(i).get(j);
                        int n_ProdId = para.ItoLogiTask.get(n_LogiTaskId).getProd().getId();
                        int n_Priority = prodPriority.indexOf(n_ProdId);
                        if (n_Priority == o_Priority-1) {
                            insertSet.add(fea);
                        }
                    }
                }
                Iterator<int[]> it = insertSet.iterator();
                while(it.hasNext()){
                    int[] fea=it.next();
                    ArrayList<ArrayList<Integer>> new_logiSeq = listCopy2(logiSeq);
                    insertOperator(new_logiSeq,ovfirst,opfirst,fea[0],fea[1]);
//                    it.remove();
                    double temp = getObj(prodSeq,new_logiSeq);
                    if (temp < new_obj){
                        new_obj = temp;
                        if(new_obj<bestObjUpBound) {
                			bestObjUpBound=new_obj;
                			bestProdSeq=prodSeq.clone();
                			bestLogiSeq=listCopy2(new_logiSeq);
                        }
                    }
                }
            } else {
                HashSet<int[]> insertSet1 = new HashSet<>();//相同优先级的插入（插入在指定位置之前）
                HashSet<int[]> insertSet2 = new HashSet<>();//相同优先级的插入（插入在指定位置之前）
                int o1_LogiTaskId = logiSeq.get(ovfirst).get(opfirst);
                int o1_ProdId = para.ItoLogiTask.get(o1_LogiTaskId).getProd().getId();
                int o1_Priority = prodPriority.indexOf(o1_ProdId);
                int o2_LogiTaskId = logiSeq.get(ovlast).get(oplast);
                int o2_ProdId = para.ItoLogiTask.get(o2_LogiTaskId).getProd().getId();
                int o2_Priority = prodPriority.indexOf(o2_ProdId);
                for (int i = 0; i < logiSeq.size(); i++) {
                    for (int j = 0; j < logiSeq.get(i).size(); j++) {
                        int[] fea = {i, j};
                        int n_LogiTaskId = logiSeq.get(i).get(j);
                        int n_ProdId = para.ItoLogiTask.get(n_LogiTaskId).getProd().getId();
                        int n_Priority = prodPriority.indexOf(n_ProdId);
                        if (n_Priority == o1_Priority-1) {
                            insertSet1.add(fea);
                        }
                        if (n_Priority == o2_Priority+1) {
                            insertSet2.add(fea);
                        }
                    }
                }
                Iterator<int[]> it1 = insertSet1.iterator();
                while(it1.hasNext()){
                    int[] fea=it1.next();
                    ArrayList<ArrayList<Integer>> new_logiSeq = listCopy2(logiSeq);
                    insertOperator(new_logiSeq,ovfirst,opfirst,fea[0],fea[1]);
//                    it1.remove();
                    double temp = getObj(prodSeq,new_logiSeq);
                    if (temp < new_obj){
                        new_obj = temp;
                        if(new_obj<bestObjUpBound) {
                			bestObjUpBound=new_obj;
                			bestProdSeq=prodSeq.clone();
                			bestLogiSeq=listCopy2(new_logiSeq);
                        }
                    }
                }
                Iterator<int[]> it2 = insertSet2.iterator();
                while(it2.hasNext()){
                    int[] fea=it2.next();
                    ArrayList<ArrayList<Integer>> new_logiSeq = listCopy2(logiSeq);
                    insertOperator(new_logiSeq,ovlast,oplast,fea[0],fea[1]);
//                    it2.remove();
                    double temp = getObj(prodSeq,new_logiSeq);
                    if (temp < new_obj){
                        new_obj = temp;
                        if(new_obj<bestObjUpBound) {
                			bestObjUpBound=new_obj;
                			bestProdSeq=prodSeq.clone();
                			bestLogiSeq=listCopy2(new_logiSeq);
                        }
                    }
                }
            }
        }

        return new_obj;
    }
//    public double LocalSearch(int[] prodSeq, ArrayList<ArrayList<Integer>> logiSeq) {
//
////        ArrayList<ArrayList<Integer>> logiSeq=new ArrayList<ArrayList<Integer>>();
////        for (int p = 0; p < para.numVeh; p++) {
////            logiSeq.add(new ArrayList<Integer>());
////        }
////
////        double[] time=new double[para.numVeh];
////        for (int p = 0; p < prodSeq.length; p++) {
////            ArrayList<TssWithoutLogiCons> tssWithoutLogiConsList = new ArrayList<TssWithoutLogiCons>();
////            String key;
////            double timeCur = 0;
////            for (int m = 0; m < para.numMach; m++) {
////                if (m > 0) {
////                    timeCur = timeCur
////                            - para.prodSet.get(prodSeq[p]).getCycleTime() * (para.demands[prodSeq[p]] - 1)
////                            + para.getDistance(para.machSet.get(m - 1), para.machSet.get(m))
////                            / para.lineSpeed;
////                }
////                for (int k = 0; k < para.K[prodSeq[p]][m]; k++) {
////                    key = prodSeq[p] + "_" + m + "_" + k;
////                    TssWithoutLogiCons tssWithoutLogiCons = new TssWithoutLogiCons(key, timeCur);
////                    tssWithoutLogiConsList.add(tssWithoutLogiCons);
////                    timeCur = timeCur + para.PMKtoLogiTask.get(key).getAssemTime(); // 加上装配时间等于该批次的完工时间，也是下一批次的最早开工时间；
////                }
////            }
////            Collections.sort(tssWithoutLogiConsList, new LogiTaskSortByTss());
////
////
////            for (int s = 0; s < tssWithoutLogiConsList.size(); s++) {
////                String pmk = tssWithoutLogiConsList.get(s).pmk;
//////                Random rnd=new Random();
//////                logiSeq.get(rnd.nextInt(para.numVeh)).add(para.PMKtoLogiTask.get(pmk).getId());
////                int logiTaskId=para.PMKtoLogiTask.get(pmk).getId();
////                int min_index = 0;
////                double requiredDeliTime = para.ItoLogiTask.get(logiTaskId).getChargingTime()
////                        + para.ItoLogiTask.get(logiTaskId).getTransTime() / 2;
////                double min_val = Float.MAX_VALUE;
////                for (int r = 0; r < para.numVeh; r++) {
////                    if (time[r] + requiredDeliTime < min_val) {
////                        min_val = time[r] + requiredDeliTime;
////                        min_index = r;
////                    }
////                }
////                logiSeq.get(min_index).add(logiTaskId);
////                time[min_index] = time[min_index] + para.ItoLogiTask.get(logiTaskId).getProcessTime();
////            }
////        }
//
//
//        double bestObjVal=getObj(prodSeq, logiSeq);
//        ArrayList<Integer> prodPriority=new ArrayList<Integer>();
//        for (int p = 0; p < prodSeq.length; p++) {
//            prodPriority.add(prodSeq[p]);
//        }
//
//        ArrayList<ArrayList<Integer>> best_logiSeq = listCopy2(logiSeq);
//
//        int iter=0;
//        Random rnd=new Random();
//        while (iter<1000) {
//            //随机找一个位置
//            //根据该位置返回可行的替换位置
//            //将两个位置做替换
//
//            int ov = rnd.nextInt(para.numVeh);
//            int op = rnd.nextInt(best_logiSeq.get(ov).size());
//            int o_LogiTaskId = best_logiSeq.get(ov).get(op);
//            int o_ProdId = para.ItoLogiTask.get(o_LogiTaskId).getProd().getId();
//            int o_Priority = prodPriority.indexOf(o_ProdId);
//
//            HashSet<int[]> feasInnerInsertPos = new HashSet<>();//相同优先级的插入（插入在指定位置之前）
//            HashSet<int[]> feasOuterInsertPos = new HashSet<>();//允许跨优先级的插入（插入在指定位置之前）
//            HashSet<int[]> feasAfterInsertPos = new HashSet<>();//插入在指定位置之后
//            HashSet<int[]> feasEndInsertPos = new HashSet<>();//插入在队列末尾
//            HashSet<int[]> feasSwapPos = new HashSet<>();//相同优先级的交换
//
//            for (int i = 0; i < best_logiSeq.size(); i++) {
//                for (int j = 0; j < best_logiSeq.get(i).size(); j++) {
//                    int[] fea = {i, j};
//                    int n_LogiTaskId = best_logiSeq.get(i).get(j);
//                    if (o_LogiTaskId != n_LogiTaskId) {
//                        int n_ProdId = para.ItoLogiTask.get(n_LogiTaskId).getProd().getId();
//                        int n_Priority = prodPriority.indexOf(n_ProdId);
//                        if (n_Priority == o_Priority) {
//                            feasSwapPos.add(fea);
//                            feasInnerInsertPos.add(fea);
//                        }
//                        if (n_Priority < o_Priority) {
//                            feasOuterInsertPos.add(fea);
//                        }
//                        if (j < best_logiSeq.get(i).size() - 1) {
//                            int n2_LogiTaskId = best_logiSeq.get(i).get(j + 1);
//                            int n2_ProdId = para.ItoLogiTask.get(n2_LogiTaskId).getProd().getId();
//                            int n2_Priority = prodPriority.indexOf(n2_ProdId);
//                            if (n_Priority <= o_Priority && n2_Priority > o_Priority) {
//                                feasAfterInsertPos.add(fea);
//                            }
//                        } else {
//                            if (n_Priority <= o_Priority) {
//                                feasEndInsertPos.add(fea);
//                            }
//                        }
//                    }
//                }
//            }
//
//            //找一个最好的领域
//            ArrayList<ArrayList<Integer>> bestNeigh_logiSeq=new ArrayList<ArrayList<Integer>>();
//            Random rndOperator = new Random();
//            boolean improve=false;
//            while(!feasSwapPos.isEmpty() || !feasInnerInsertPos.isEmpty() || !feasAfterInsertPos.isEmpty() || !feasEndInsertPos.isEmpty()){
//                switch (rndOperator.nextInt(4)) {
//                    case 0:
//                        if (!feasInnerInsertPos.isEmpty()) {
//                            Iterator<int[]> iterator = feasInnerInsertPos.iterator();
//                            int[] fea = iterator.next();
//                            feasInnerInsertPos.remove(fea);
//                            ArrayList<ArrayList<Integer>> candiNeigh_logiSeq = listCopy2(best_logiSeq);
//                            innerInsertOperator(candiNeigh_logiSeq, ov, op, fea[0], fea[1]);
//                            double testObjVal = getObj(prodSeq, candiNeigh_logiSeq);
//                            if (testObjVal <= bestObjVal) {
//                                bestObjVal = testObjVal;
//                                bestNeigh_logiSeq = listCopy2(candiNeigh_logiSeq);
//                                improve=true;
//                            }
//                        }
//                        break;
//                    case 1:
//                        if (!feasAfterInsertPos.isEmpty()) {
//                            Iterator<int[]> iterator = feasAfterInsertPos.iterator();
//                            int[] fea = iterator.next();
//                            feasAfterInsertPos.remove(fea);
//                            ArrayList<ArrayList<Integer>> candiNeigh_logiSeq = listCopy2(best_logiSeq);
//                            afterInsertOperator(candiNeigh_logiSeq, ov, op, fea[0], fea[1]);
//                            double testObjVal = getObj(prodSeq, candiNeigh_logiSeq);
//                            if (testObjVal <= bestObjVal) {
//                                bestObjVal = testObjVal;
//                                bestNeigh_logiSeq = listCopy2(candiNeigh_logiSeq);
//                                improve=true;
//                            }
//                        }
//                        break;
//                    case 2:
//                        if (!feasEndInsertPos.isEmpty()) {
//                            Iterator<int[]> iterator = feasEndInsertPos.iterator();
//                            int[] fea = iterator.next();
//                            feasEndInsertPos.remove(fea);
//                            ArrayList<ArrayList<Integer>> candiNeigh_logiSeq = listCopy2(best_logiSeq);
//                            endInsertOperator(candiNeigh_logiSeq, ov, op, fea[0], fea[1]);
//                            double testObjVal = getObj(prodSeq, candiNeigh_logiSeq);
//                            if (testObjVal <= bestObjVal) {
//                                bestObjVal = testObjVal;
//                                bestNeigh_logiSeq = listCopy2(candiNeigh_logiSeq);
//                                improve=true;
//                            }
//                        }
//                        break;
//                    case 3:
//                        if (!feasSwapPos.isEmpty()) {
//                            Iterator<int[]> iterator = feasSwapPos.iterator();
//                            int[] fea = iterator.next();
//                            feasSwapPos.remove(fea);
//                            ArrayList<ArrayList<Integer>> candiNeigh_logiSeq = listCopy2(best_logiSeq);
//                            swapOperator(candiNeigh_logiSeq, ov, op, fea[0], fea[1]);
//                            double testObjVal = getObj(prodSeq, candiNeigh_logiSeq);
//                            if (testObjVal <= bestObjVal) {
//                                bestObjVal = testObjVal;
//                                bestNeigh_logiSeq = listCopy2(candiNeigh_logiSeq);
//                                improve=true;
//                            }
//                        }
//                        break;
//                }
//
//            }
//            for (int[] fea : feasOuterInsertPos) {
//                ArrayList<ArrayList<Integer>> candiNeigh_logiSeq = listCopy2(best_logiSeq);
//                outerInsertOperator(candiNeigh_logiSeq, ov, op, fea[0], fea[1]);
//                double testObjVal = getObj(prodSeq, candiNeigh_logiSeq);
//                if (testObjVal <= bestObjVal) {
//                    bestObjVal = testObjVal;
//                    bestNeigh_logiSeq = listCopy2(candiNeigh_logiSeq);
//                    improve=true;
//                }
//            }
//
//
////            swapOperator(neigh_logiSeq,ov,op,nv,np);
////            afterInsertOperator(neigh_logiSeq,ov,op,nv,np);
////            endInsertOperator(neigh_logiSeq,ov,op,nv,np);
////            innerInsertOperator(neigh_logiSeq,ov,op,nv,np);
////            outerInsertOperator(neigh_logiSeq,ov,op,nv,np);
//            if(improve) {
//                best_logiSeq = bestNeigh_logiSeq;
//            }
//
//            iter++;
//        }
//
//        return bestObjVal;
//    }

    public double LocalSearch(int[] prodSeq,  ArrayList<ArrayList<Integer>> ini_logiSeq) {

        ArrayList<ArrayList<Integer>> cumb_logiSeq=new ArrayList<ArrayList<Integer>>();
        for (int r = 0; r < para.numVeh; r++) {
            cumb_logiSeq.add(new ArrayList<Integer>());
        }

        ArrayList<Integer> first = new ArrayList<>();
        ArrayList<Integer> last = new ArrayList<>();
        for (int p = 0; p < prodSeq.length; p++) {// product从0开始编号
            int pindex = prodSeq[p];
            ArrayList<int[]> arrayofPMKI = para.PtoArrayofPMKI.get(pindex);
            ArrayList<Integer> arrayOfI = new ArrayList<Integer>();
            for (int index = 0; index < arrayofPMKI.size(); index++) {
                arrayOfI.add(arrayofPMKI.get(index)[3]);
            }
            first.add(arrayOfI.get(0));
            last.add(arrayOfI.get(arrayOfI.size() - 1));
        }

        ArrayList<Integer> prodPriority=new ArrayList<Integer>();
        for (int p = 0; p < prodSeq.length; p++) {
            prodPriority.add(prodSeq[p]);
        }
        int[] priorityLogiTask=new int[logiTaskSet.size()];
        for (int i = 0; i < ini_logiSeq.size(); i++) {
            for (int j = 0; j < ini_logiSeq.get(i).size(); j++) {
                int logiTaskId=ini_logiSeq.get(i).get(j);
                int prodId=para.ItoLogiTask.get(logiTaskId).getProd().getId();
                priorityLogiTask[logiTaskId]=prodPriority.indexOf(prodId);
            }
        }
        for (int i = 0; i < ini_logiSeq.size(); i++) {
            ArrayList<LogiTaskPrior> priorSet=new ArrayList<LogiTaskPrior>();
            for (int j = 0; j < ini_logiSeq.get(i).size(); j++) {
                int logiTaskId = ini_logiSeq.get(i).get( j );
                priorSet.add(new LogiTaskPrior(logiTaskId,priorityLogiTask[logiTaskId]));
//                System.out.print(priorityLogiTask[logiTaskId]+" ");
            }
//            System.out.println(priorSet.size());
            Collections.sort(priorSet, new LogiTaskSortByPrior());

            for (int j = 0; j < priorSet.size(); j++) {
                cumb_logiSeq.get(i).add(priorSet.get(j).id);
            }

        }
        
        double temp1=getObj(prodSeq, cumb_logiSeq);
        double temp2=getObj(prodSeq, ini_logiSeq);

        double bestObjVal=Math.min(temp1,temp2);
        
        if (temp1<bestObjUpBound) {
			bestObjUpBound=temp1;
			bestProdSeq=prodSeq.clone();
			bestLogiSeq=listCopy2(cumb_logiSeq);
		}
        if (temp2<bestObjUpBound) {
			bestObjUpBound=temp2;
			bestProdSeq=prodSeq.clone();
			bestLogiSeq=listCopy2(ini_logiSeq);
		}

        
        

        ArrayList<ArrayList<Integer>> logiSeq = listCopy2(cumb_logiSeq);

        int iter=0;
        Random rnd=new Random();
        while (iter<2000) {
            int ov;
            int op ;
            while (true){
                int a = rnd.nextInt(para.numVeh);
                if (logiSeq.get(a).size()>1){
                    ov=a;
                    op=rnd.nextInt(logiSeq.get(a).size());
                    break;
                }
            }
            int o_LogiTaskId = logiSeq.get(ov).get(op);
            int o_ProdId = para.ItoLogiTask.get(o_LogiTaskId).getProd().getId();
            int o_Priority = prodPriority.indexOf(o_ProdId);

            HashSet<int[]> feasInnerInsertPos = new HashSet<>();//相同优先级的插入（插入在指定位置之前）
            HashSet<int[]> feasAfterInsertPos = new HashSet<>();//插入在指定位置之后
            HashSet<int[]> feasEndInsertPos = new HashSet<>();//插入在队列末尾
            HashSet<int[]> feasSwapPos = new HashSet<>();//相同优先级的交换

            for (int i = 0; i < logiSeq.size(); i++) {
                for (int j = 0; j < logiSeq.get(i).size(); j++) {
                    int[] fea = {i, j};
                    int n_LogiTaskId = logiSeq.get(i).get(j);
                    if (o_LogiTaskId != n_LogiTaskId) {
                        int n_ProdId = para.ItoLogiTask.get(n_LogiTaskId).getProd().getId();
                        int n_Priority = prodPriority.indexOf(n_ProdId);
                        if (n_Priority == o_Priority) {
                            feasSwapPos.add(fea);
                            feasInnerInsertPos.add(fea);
                        }
                        if (j < logiSeq.get(i).size() - 1) {
                            int n2_LogiTaskId = logiSeq.get(i).get(j + 1);
                            int n2_ProdId = para.ItoLogiTask.get(n2_LogiTaskId).getProd().getId();
                            int n2_Priority = prodPriority.indexOf(n2_ProdId);
                            if (n_Priority <= o_Priority && n2_Priority > o_Priority) {
                                feasAfterInsertPos.add(fea);
                            }
                        } else {
                            if (n_Priority <= o_Priority) {
                                feasEndInsertPos.add(fea);
                            }
                        }
                    }
                }
            }
            Random rndOperator = new Random();
            boolean improved=false;
            while( !improved && (!feasSwapPos.isEmpty() || !feasInnerInsertPos.isEmpty() || !feasAfterInsertPos.isEmpty() || !feasEndInsertPos.isEmpty() ) ){
                switch (rndOperator.nextInt(4)) {
                    case 0:
                        if (!feasInnerInsertPos.isEmpty()) {
                            Iterator<int[]> iterator = feasInnerInsertPos.iterator();
                            int[] fea = iterator.next();
                            feasInnerInsertPos.remove(fea);
                            ArrayList<ArrayList<Integer>> neigh_logiSeq = listCopy2(logiSeq);
                            innerInsertOperator(neigh_logiSeq, ov, op, fea[0], fea[1]);
                            double testObjVal = getObj(prodSeq, neigh_logiSeq);
                            if (testObjVal <= bestObjVal) {
                                bestObjVal = testObjVal;
                                if(bestObjVal<bestObjUpBound) {
                        			bestObjUpBound=bestObjVal;
                        			bestProdSeq=prodSeq.clone();
                        			bestLogiSeq=listCopy2(neigh_logiSeq);
                                }
                                logiSeq= neigh_logiSeq;
//                                if(!checkPriority(prodSeq,neigh_logiSeq)){
//                                    int a=0;
//                                }
                                improved = true;
                            }
                        }
                        break;
                    case 1:
                        if (!feasEndInsertPos.isEmpty()) {
                            Iterator<int[]> iterator = feasEndInsertPos.iterator();
                            int[] fea = iterator.next();
                            feasEndInsertPos.remove(fea);
                            ArrayList<ArrayList<Integer>> neigh_logiSeq = listCopy2(logiSeq);
                            endInsertOperator(neigh_logiSeq, ov, op, fea[0], fea[1]);
//                            if(!checkPriority(prodSeq,neigh_logiSeq)){
//                                int a=0;
//                            }
                            double testObjVal = getObj(prodSeq, neigh_logiSeq);
                            if (testObjVal <= bestObjVal) {
                                bestObjVal = testObjVal;
                                if(bestObjVal<bestObjUpBound) {
                        			bestObjUpBound=bestObjVal;
                        			bestProdSeq=prodSeq.clone();
                        			bestLogiSeq=listCopy2(neigh_logiSeq);
                                }
                                logiSeq= neigh_logiSeq;
                                improved = true;
                            }
                        }
                        break;
                    case 2:
                        if (!feasAfterInsertPos.isEmpty()) {
                            Iterator<int[]> iterator = feasAfterInsertPos.iterator();
                            int[] fea = iterator.next();
                            feasAfterInsertPos.remove(fea);
                            ArrayList<ArrayList<Integer>> neigh_logiSeq = listCopy2(logiSeq);
                            afterInsertOperator(neigh_logiSeq, ov, op, fea[0], fea[1]);
//                            if(!checkPriority(prodSeq,neigh_logiSeq)){
//                                int a=0;
//                            }
                            double testObjVal = getObj(prodSeq, neigh_logiSeq);
                            if (testObjVal <= bestObjVal) {
                                bestObjVal = testObjVal;
                                if(bestObjVal<bestObjUpBound) {
                        			bestObjUpBound=bestObjVal;
                        			bestProdSeq=prodSeq.clone();
                        			bestLogiSeq=listCopy2(neigh_logiSeq);
                                }
                                logiSeq= neigh_logiSeq;
                                improved = true;
                            }
                        }
                        break;
                    case 3:
                        if (!feasSwapPos.isEmpty()) {
                            Iterator<int[]> iterator = feasSwapPos.iterator();
                            int[] fea = iterator.next();
                            feasSwapPos.remove(fea);
                            ArrayList<ArrayList<Integer>> neigh_logiSeq = listCopy2(logiSeq);
                            swapOperator(neigh_logiSeq, ov, op, fea[0], fea[1]);
//                            if(!checkPriority(prodSeq,neigh_logiSeq)){
//                                int a=0;
//                            }
                            double testObjVal = getObj(prodSeq, neigh_logiSeq);
                            if (testObjVal <= bestObjVal) {
                                bestObjVal = testObjVal;
                                if(bestObjVal<bestObjUpBound) {
                        			bestObjUpBound=bestObjVal;
                        			bestProdSeq=prodSeq.clone();
                        			bestLogiSeq=listCopy2(neigh_logiSeq);
                                }
                                logiSeq= neigh_logiSeq;
                                improved = true;
                            }
                        }
                        break;
                }
                //做完简单的领域搜索后再试探一下跨优先级的更优解，但是这个解不会被记录下来
//                int[][] posi=new int[para.numProd][4];
//                for (int i = 0; i < logiSeq.size(); i++) {
//                    for (int j = 0; j < logiSeq.get(i).size(); j++) {
//                        int logiTaskId=logiSeq.get(i).get(j);
//                        int prodId=para.ItoLogiTask.get(logiTaskId).getProd().getId();
//                        if(first.contains(logiTaskId)){
//                            posi[prodId][0]=i;
//                            posi[prodId][1]=j;
//                        }
//                        if(last.contains(logiTaskId)){
//                            posi[prodId][2]=i;
//                            posi[prodId][3]=j;
//                        }
//                    }
//                }
//                ArrayList<ArrayList<Integer>> neigh_logiSeq=listCopy2(logiSeq);
//                double testObjVal=LocalImprove(prodSeq, neigh_logiSeq,bestObjVal,posi);
//                if (testObjVal <= bestObjVal) {
//                    bestObjVal = testObjVal;
//                    improved = true;
//                }
            }
            iter++;
        }

        int[][] posi=new int[para.numProd][4];
        for (int i = 0; i < logiSeq.size(); i++) {
            for (int j = 0; j < logiSeq.get(i).size(); j++) {
                int logiTaskId=logiSeq.get(i).get(j);
                int prodId=para.ItoLogiTask.get(logiTaskId).getProd().getId();
                if(first.contains(logiTaskId)){
                    posi[prodId][0]=i;
                    posi[prodId][1]=j;
                }
                if(last.contains(logiTaskId)){
                    posi[prodId][2]=i;
                    posi[prodId][3]=j;
                }
            }
        }
        ArrayList<ArrayList<Integer>> neigh_logiSeq=listCopy2(logiSeq);
        bestObjVal=Math.min(bestObjVal,LocalImprove(prodSeq, neigh_logiSeq,bestObjVal,posi));

        return bestObjVal;
    }

//    public double LocalSearch(int[] prodSeq) {
//
//
//        ArrayList<Integer> first = new ArrayList<>();
//        ArrayList<Integer> last = new ArrayList<>();
//        for (int p = 0; p < prodSeq.length; p++) {// product从0开始编号
//            int pindex = prodSeq[p];
//            ArrayList<int[]> arrayofPMKI = para.PtoArrayofPMKI.get(pindex);
//            ArrayList<Integer> arrayOfI = new ArrayList<Integer>();
//            for (int index = 0; index < arrayofPMKI.size(); index++) {
//                arrayOfI.add(arrayofPMKI.get(index)[3]);
//            }
//            first.add(arrayOfI.get(0));
//            last.add(arrayOfI.get(arrayOfI.size() - 1));
//        }
//
//        ArrayList<ArrayList<Integer>> ini_logiSeq= new ArrayList<ArrayList<Integer>>();
//        for (int i = 0; i < para.numVeh; i++) {
//            ini_logiSeq.add(new ArrayList<Integer>());
//        }
//        for (int p = 0; p < prodSeq.length; p++) {
//            int pindex = prodSeq[p];
//            ArrayList<int[]> arrayofPMKI = para.PtoArrayofPMKI.get(pindex);
//            ArrayList<Integer> arrayOfI = new ArrayList<Integer>();
//            ArrayList<Integer> idSet = new ArrayList<Integer>();
//            for (int index = 0; index < arrayofPMKI.size(); index++) {
//                arrayOfI.add(arrayofPMKI.get(index)[3]);
//            }
//            for (int s = 0; s < arrayOfI.size(); s++) {
//                Random rnd=new Random();
//                ini_logiSeq.get(rnd.nextInt(para.numVeh)).add(arrayOfI.get(s));
//            }
//        }
//
//        ArrayList<Integer> prodPriority=new ArrayList<Integer>();
//        for (int p = 0; p < prodSeq.length; p++) {
//            prodPriority.add(prodSeq[p]);
//        }
//
//        double bestObjVal=getObj(prodSeq, ini_logiSeq);
//
//        ArrayList<ArrayList<Integer>> logiSeq = listCopy2(ini_logiSeq);
//
//        int iter=0;
//        Random rnd=new Random();
//        while (iter<1000) {
//            int ov;
//            int op ;
//            while (true){
//                int a = rnd.nextInt(para.numVeh);
//                if (logiSeq.get(a).size()>1){
//                    ov=a;
//                    op=rnd.nextInt(logiSeq.get(a).size());
//                    break;
//                }
//            }
//            int o_LogiTaskId = logiSeq.get(ov).get(op);
//            int o_ProdId = para.ItoLogiTask.get(o_LogiTaskId).getProd().getId();
//            int o_Priority = prodPriority.indexOf(o_ProdId);
//
//            HashSet<int[]> feasInnerInsertPos = new HashSet<>();//相同优先级的插入（插入在指定位置之前）
//            HashSet<int[]> feasAfterInsertPos = new HashSet<>();//插入在指定位置之后
//            HashSet<int[]> feasEndInsertPos = new HashSet<>();//插入在队列末尾
//            HashSet<int[]> feasSwapPos = new HashSet<>();//相同优先级的交换
//
//            for (int i = 0; i < logiSeq.size(); i++) {
//                for (int j = 0; j < logiSeq.get(i).size(); j++) {
//                    int[] fea = {i, j};
//                    int n_LogiTaskId = logiSeq.get(i).get(j);
//                    if (o_LogiTaskId != n_LogiTaskId) {
//                        int n_ProdId = para.ItoLogiTask.get(n_LogiTaskId).getProd().getId();
//                        int n_Priority = prodPriority.indexOf(n_ProdId);
//                        if (n_Priority == o_Priority) {
//                            feasSwapPos.add(fea);
//                            feasInnerInsertPos.add(fea);
//                        }
//                        if (j < logiSeq.get(i).size() - 1) {
//                            int n2_LogiTaskId = logiSeq.get(i).get(j + 1);
//                            int n2_ProdId = para.ItoLogiTask.get(n2_LogiTaskId).getProd().getId();
//                            int n2_Priority = prodPriority.indexOf(n2_ProdId);
//                            if (n_Priority <= o_Priority && n2_Priority > o_Priority) {
//                                feasAfterInsertPos.add(fea);
//                            }
//                        } else {
//                            if (n_Priority <= o_Priority) {
//                                feasEndInsertPos.add(fea);
//                            }
//                        }
//                    }
//                }
//            }
//            Random rndOperator = new Random();
//            boolean improved=false;
//            while( !improved && (!feasSwapPos.isEmpty() || !feasInnerInsertPos.isEmpty() || !feasAfterInsertPos.isEmpty() || !feasEndInsertPos.isEmpty() ) ){
//                switch (rndOperator.nextInt(4)) {
//                    case 0:
//                        if (!feasInnerInsertPos.isEmpty()) {
//                            Iterator<int[]> iterator = feasInnerInsertPos.iterator();
//                            int[] fea = iterator.next();
//                            feasInnerInsertPos.remove(fea);
//                            ArrayList<ArrayList<Integer>> neigh_logiSeq = listCopy2(logiSeq);
//                            innerInsertOperator(neigh_logiSeq, ov, op, fea[0], fea[1]);
//                            double testObjVal = getObj(prodSeq, neigh_logiSeq);
//                            if (testObjVal <= bestObjVal) {
//                                bestObjVal = testObjVal;
//                                logiSeq= neigh_logiSeq;
//                                if(!checkPriority(prodSeq,neigh_logiSeq)){
//                                    int a=0;
//                                }
//                                improved = true;
//                            }
//                        }
//                        break;
//                    case 1:
//                        if (!feasEndInsertPos.isEmpty()) {
//                            Iterator<int[]> iterator = feasEndInsertPos.iterator();
//                            int[] fea = iterator.next();
//                            feasEndInsertPos.remove(fea);
//                            ArrayList<ArrayList<Integer>> neigh_logiSeq = listCopy2(logiSeq);
//                            endInsertOperator(neigh_logiSeq, ov, op, fea[0], fea[1]);
//                            if(!checkPriority(prodSeq,neigh_logiSeq)){
//                                int a=0;
//                            }
//                            double testObjVal = getObj(prodSeq, neigh_logiSeq);
//                            if (testObjVal <= bestObjVal) {
//                                bestObjVal = testObjVal;
//                                logiSeq= neigh_logiSeq;
//                                improved = true;
//                            }
//                        }
//                        break;
//                    case 2:
//                        if (!feasAfterInsertPos.isEmpty()) {
//                            Iterator<int[]> iterator = feasAfterInsertPos.iterator();
//                            int[] fea = iterator.next();
//                            feasAfterInsertPos.remove(fea);
//                            ArrayList<ArrayList<Integer>> neigh_logiSeq = listCopy2(logiSeq);
//                            afterInsertOperator(neigh_logiSeq, ov, op, fea[0], fea[1]);
//                            if(!checkPriority(prodSeq,neigh_logiSeq)){
//                                int a=0;
//                            }
//                            double testObjVal = getObj(prodSeq, neigh_logiSeq);
//                            if (testObjVal <= bestObjVal) {
//                                bestObjVal = testObjVal;
//                                logiSeq= neigh_logiSeq;
//                                improved = true;
//                            }
//                        }
//                        break;
//                    case 3:
//                        if (!feasSwapPos.isEmpty()) {
//                            Iterator<int[]> iterator = feasSwapPos.iterator();
//                            int[] fea = iterator.next();
//                            feasSwapPos.remove(fea);
//                            ArrayList<ArrayList<Integer>> neigh_logiSeq = listCopy2(logiSeq);
//                            swapOperator(neigh_logiSeq, ov, op, fea[0], fea[1]);
//                            if(!checkPriority(prodSeq,neigh_logiSeq)){
//                                int a=0;
//                            }
//                            double testObjVal = getObj(prodSeq, neigh_logiSeq);
//                            if (testObjVal <= bestObjVal) {
//                                bestObjVal = testObjVal;
//                                logiSeq= neigh_logiSeq;
//                                improved = true;
//                            }
//                        }
//                        break;
//                }
//                //做完简单的领域搜索后再试探一下跨优先级的更优解，但是这个解不会被记录下来
//                int[][] posi=new int[para.numProd][4];
//                for (int i = 0; i < logiSeq.size(); i++) {
//                    for (int j = 0; j < logiSeq.get(i).size(); j++) {
//                        int logiTaskId=logiSeq.get(i).get(j);
//                        int prodId=para.ItoLogiTask.get(logiTaskId).getProd().getId();
//                        if(first.contains(logiTaskId)){
//                            posi[prodId][0]=i;
//                            posi[prodId][1]=j;
//                        }
//                        if(last.contains(logiTaskId)){
//                            posi[prodId][2]=i;
//                            posi[prodId][3]=j;
//                        }
//                    }
//                }
//                ArrayList<ArrayList<Integer>> neigh_logiSeq=listCopy2(logiSeq);
//                double testObjVal=LocalImprove(prodSeq, neigh_logiSeq,bestObjVal,posi);
//                if (testObjVal <= bestObjVal) {
//                    bestObjVal = testObjVal;
//                    improved = true;
//                }
//            }
//            iter++;
//        }
//
//        int[][] posi=new int[para.numProd][4];
//        for (int i = 0; i < logiSeq.size(); i++) {
//            for (int j = 0; j < logiSeq.get(i).size(); j++) {
//                int logiTaskId=logiSeq.get(i).get(j);
//                int prodId=para.ItoLogiTask.get(logiTaskId).getProd().getId();
//                if(first.contains(logiTaskId)){
//                    posi[prodId][0]=i;
//                    posi[prodId][1]=j;
//                }
//                if(last.contains(logiTaskId)){
//                    posi[prodId][2]=i;
//                    posi[prodId][3]=j;
//                }
//            }
//        }
//        ArrayList<ArrayList<Integer>> neigh_logiSeq=listCopy2(logiSeq);
//        bestObjVal=Math.min(bestObjVal,LocalImprove(prodSeq, neigh_logiSeq,bestObjVal,posi));
//
//        return bestObjVal;
//    }

//    public boolean checkPriority(int[] prodSeq, ArrayList<ArrayList<Integer>> logiSeq){
//        ArrayList<Integer> prodPriority=new ArrayList<Integer>();
//        for (int p = 0; p < prodSeq.length; p++) {
//            prodPriority.add(prodSeq[p]);
//        }
//        boolean flag=true;
//        for (int i = 0; i < logiSeq.size(); i++) {
//            if (logiSeq.get(i).size()>1) {
//                for (int j = 0; j < logiSeq.get(i).size() - 1; j++) {
//                    int prior1=prodPriority.indexOf(  para.ItoLogiTask.get(  logiSeq.get(i).get(j)  ).getProd().getId() );
//                    int prior2=prodPriority.indexOf(  para.ItoLogiTask.get(  logiSeq.get(i).get(j+1)  ).getProd().getId() );
//                    if (prior1>prior2){
//                        flag=false;
//                        System.out.println(i+":  "+j+ " and " +(j+1));
//                    }
//                }
//            }
//        }
//        return flag;
//    }

    public void insertOperator(ArrayList<ArrayList<Integer>> logiSeq,int original_veh,int original_position,int new_veh,int new_position){
        int temp1=logiSeq.get(original_veh).get(original_position);
        if(original_position<=new_position) {
            logiSeq.get(new_veh).add(new_position, temp1);
            logiSeq.get(original_veh).remove(original_position);
        }else{
            logiSeq.get(original_veh).remove(original_position);
            logiSeq.get(new_veh).add(new_position, temp1);
        }
    }
    public void swapOperator(ArrayList<ArrayList<Integer>> logiSeq,int original_veh,int original_position,int new_veh,int new_position){
        int temp1=logiSeq.get(original_veh).get(original_position);
        int temp2=logiSeq.get(new_veh).get(new_position);
        logiSeq.get(original_veh).set(original_position,temp2);
        logiSeq.get(new_veh).set(new_position,temp1);
    }
    public void innerInsertOperator(ArrayList<ArrayList<Integer>> logiSeq,int original_veh,int original_position,int new_veh,int new_position){
        int temp1=logiSeq.get(original_veh).get(original_position);
        if(original_position<=new_position) {
            logiSeq.get(new_veh).add(new_position, temp1);
            logiSeq.get(original_veh).remove(original_position);
        }else{
            logiSeq.get(original_veh).remove(original_position);
            logiSeq.get(new_veh).add(new_position, temp1);
        }
    }
//    public void outerInsertOperator(ArrayList<ArrayList<Integer>> logiSeq,int original_veh,int original_position,int new_veh,int new_position){
//        int temp1=logiSeq.get(original_veh).get(original_position);
//        logiSeq.get(original_veh).remove(original_position);
//        logiSeq.get(new_veh).add(new_position,temp1);
//    }
    public void afterInsertOperator(ArrayList<ArrayList<Integer>> logiSeq,int original_veh,int original_position,int new_veh,int new_position){
        int temp1=logiSeq.get(original_veh).get(original_position);
        logiSeq.get(new_veh).add(new_position+1,temp1);
        logiSeq.get(original_veh).remove(original_position);
    }
    public void endInsertOperator(ArrayList<ArrayList<Integer>> logiSeq,int original_veh,int original_position,int new_veh,int new_position){
        int temp1=logiSeq.get(original_veh).get(original_position);
        logiSeq.get(new_veh).add(temp1);
        logiSeq.get(original_veh).remove(original_position);
    }


//    public ArrayList<Integer> listCopy1(ArrayList<Integer> a){
//        ArrayList<Integer> b=new ArrayList<Integer>();
//        for (int i = 0; i < a.size(); i++) {
//            b.add(a.get(i));
//        }
//        return b;
//    }
    public ArrayList<ArrayList<Integer>> listCopy2(ArrayList<ArrayList<Integer>> a){
        ArrayList<ArrayList<Integer>> b=new ArrayList<ArrayList<Integer>>();
        for (int i = 0; i < a.size(); i++) {
            b.add(new ArrayList<Integer>());
        }
        for (int i = 0; i < a.size(); i++) {
            for (int j = 0; j < a.get(i).size(); j++) {
                b.get(i).add(a.get(i).get(j));
            }
        }
        return b;
    }

//    public ArrayList<ArrayList<Integer>> swapOperator(ArrayList<ArrayList<Integer>> cumb_iList){
//        //initialize neigh_iList
//        ArrayList<ArrayList<Integer>> neigh_iList=new ArrayList<ArrayList<Integer>>();
//        for (int p = 0; p < prodSeq.length; p++) {
//            neigh_iList.add(new ArrayList<Integer>());
//        }
//        for (int i = 0; i < cumb_iList.size(); i++) {
//            for (int j = 0; j < cumb_iList.get(i).size(); j++) {
//                int idx=cumb_iList.get(i).get(j);
//                neigh_iList.get(i).add(idx);
//            }
//        }
//
//
//        Random rnd=new Random();
//        int rndP=rnd.nextInt(para.numProd);
//        int rnd1=rnd.nextInt(neigh_iList.get(rndP).size());
//        int rnd2=rnd.nextInt(neigh_iList.get(rndP).size());
//        Collections.swap(neigh_iList.get(rndP),rnd1,rnd2);
//
//        return neigh_iList;
//    }

//    public ArrayList<ArrayList<Integer>> listCopy(ArrayList<ArrayList<Integer>> list){
//
//    }

//    public ArrayList<ArrayList<Integer>> routeGeneration0(ArrayList<Integer> iList){
//        ArrayList<ArrayList<Integer>> logiSeq=new ArrayList<ArrayList<Integer>>();
//        for (int r = 0; r < para.numVeh; r++) {
//            logiSeq.add(new ArrayList<Integer>());
//        }
//        double[] timeCur = new double[para.numVeh];
//        for (int index = 0; index < iList.size(); index++) {
//            int logiTaskId = iList.get(index);// 对应logiTask的ID
//            int min_index = 0;
//            double requiredDeliTime = para.ItoLogiTask.get(logiTaskId).getChargingTime()
//                    + para.ItoLogiTask.get(logiTaskId).getTransTime() / 2;
//            double min_val = Float.MAX_VALUE;
//            for (int r = 0; r < para.numVeh; r++) {
//                if (timeCur[r] + requiredDeliTime < min_val) {
//                    min_val = timeCur[r] + requiredDeliTime;
//                    min_index = r;
//                }
//            }
//            logiSeq.get(min_index).add(logiTaskId);
//            timeCur[min_index] = timeCur[min_index] + para.ItoLogiTask.get(logiTaskId).getProcessTime();
//        }
//        return logiSeq;
//    }
//    public ArrayList<ArrayList<Integer>> routeGeneration(ArrayList<ArrayList<Integer>> iList){
//        ArrayList<ArrayList<Integer>> logiSeq=new ArrayList<ArrayList<Integer>>();
//        for (int r = 0; r < para.numVeh; r++) {
//            logiSeq.add(new ArrayList<Integer>());
//        }
//        double[] timeCur = new double[para.numVeh];
//        for (int p = 0; p < prodSeq.length; p++) {// product从0开始编号
//            ArrayList<Integer> idSet=iList.get(p);
//
//            for (int index = 0; index < idSet.size(); index++) {
//                int logiTaskId = idSet.get(index);// 对应logiTask的ID
//                int min_index = 0;
//                double requiredDeliTime = para.ItoLogiTask.get(logiTaskId).getChargingTime()
//                        + para.ItoLogiTask.get(logiTaskId).getTransTime() / 2;
//                double min_val = Float.MAX_VALUE;
//                for (int r = 0; r < para.numVeh; r++) {
//                    if (timeCur[r] + requiredDeliTime < min_val) {
//                        min_val = timeCur[r] + requiredDeliTime;
//                        min_index = r;
//                    }
//                }
//                logiSeq.get(min_index).add(logiTaskId);
//                timeCur[min_index] = timeCur[min_index] + para.ItoLogiTask.get(logiTaskId).getProcessTime();
//            }
//        }
//        return logiSeq;
//    }
//    public void getTad(int[] prodSeq, ArrayList<ArrayList<Integer>> logiSeq) {
//        for (int i = 0; i < logiSeq.size(); i++) {// 第i个路径
//            double cur = 0;
//            for (int j = 0; j < logiSeq.get(i).size(); j++) {// 第i个路径上的第j个位置
//                para.ItoLogiTask.get(logiSeq.get(i).get(j)).tadHeur = cur
//                        + para.getDistance(para.depot, para.ItoLogiTask.get(logiSeq.get(i).get(j)).getMach())
//                        * (1 / para.vehSpeed + 2 * para.consumingRate * para.chargingRate);
//                cur = para.ItoLogiTask.get(logiSeq.get(i).get(j)).tadHeur;
//                cur = cur + para.getDistance(para.depot, para.ItoLogiTask.get(logiSeq.get(i).get(j)).getMach())
//                        / para.vehSpeed;
//            }
//        }
//
//        String key;
//        double timeCur = 0;// 时间游标
//        double sum_taf = 0;
//        double[] taf = new double[para.numProd];
//        key = prodSeq[0] + "_" + 0 + "_" + 0;
//        timeCur = PMKtoLogiTask.get(key).tadHeur;// 第1个产品第1个机器第1个批次的物料的开始时间 必然是其对应的tadHeur(已经算好了)
//        for (int p = 0; p < prodSeq.length; p++) {// product从0开始编号
//            for (int m = 0; m < para.numMach; m++) {// machine从0开始编号
//                if (m > 0) {// 自第一个机器后，最早开工时间需要做一次减法
//                    timeCur = timeCur - para.prodSet.get(prodSeq[p]).getCycleTime() * (para.demands[prodSeq[p]] - 1);
//                }
//                for (int k = 0; k < para.K[prodSeq[p]][m]; k++) {// batch从0开始编号
//                    key = prodSeq[p] + "_" + m + "_" + k;
//                    // 通过key获得PMKtoLogiTask.get(key)是p,m,k对应的运输任务
//                    timeCur = Math.max(timeCur, PMKtoLogiTask.get(key).tadHeur); // 真正的开始时间是物料到达和最早开工时间取大；
//                    timeCur = timeCur + PMKtoLogiTask.get(key).getAssemTime(); // 加上装配时间等于该批次的完工时间，也是下一批次的最早开工时间；
//                }
//            }
//            taf[p] = timeCur;
//            sum_taf = sum_taf + timeCur;
//        }
//
//        for (int p = 0; p < para.numProd; p++) {
//            System.out.print("-> " + prodSeq[p] + "(" + taf[p] + ") ");
//        }
//        System.out.println();
//
//        for (int i = 0; i < logiSeq.size(); i++) {
//            for (int j = 0; j < logiSeq.get(i).size(); j++) {
//                System.out.print("-> " + logiSeq.get(i).get(j) + "("
//                        + para.ItoLogiTask.get(logiSeq.get(i).get(j)).tadHeur + ") ");
//            }
//            System.out.println();
//        }
//        System.out.println();
//
//    }

    public void updateSubgradient() {
//		update subgradient with the solution of lagrangian relaxation problem given the lagrangian multipliers
        for (LogiTask logiTask : logiTaskSet) {
//			logiTask.setSubgrad(logiTask.getTad() - (logiTask.getTss()+earliestReadyTime));
            logiTask.setSubgrad(logiTask.getTad() - logiTask.getTss());
//			logiTask.setSubgrad(logiTask.getTad() - (logiTask.getTss()-makespan));
//			logiTask.setSubgrad(logiTask.getTad() - (logiTask.getTss()+para.getDistance(para.depot, para.machSet.get(0))*2*para.consumingRate*para.chargingRate+para.getDistance(para.depot, para.machSet.get(0))/para.vehSpeed ));
        }

//		int cnt=0;
//		for(LogiTask logiTask : logiTaskSet) {
//			if (logiTask.getSubgrad()<=0) {
//				cnt++;
//			}
//		}
//		if (cnt==logiTaskSet.size()) {
//			for (LogiTask logiTask : logiTaskSet) {
//				logiTask.setSubgrad(logiTask.getTad() - (logiTask.getTss()));
//			}
//		}
    }

    public void updateStepsize() {
        boolean bb = true;
        if (bb) {
//			stepsizeFactor = 0.02;
            if (improved) {
                iter_for_stepsize_factor = 0;
            } else {
                iter_for_stepsize_factor++;
            }
            if (iter_for_stepsize_factor >= max_iter_no_improve_for_stepsize_factor) {// 若目标值在固定次数的迭代后仍无改进，stepsizeFactor值将减半
                stepsizeFactor = stepsizeFactor / 2;
                iter_for_stepsize_factor = 0;
            }
            double numerator = stepsizeFactor * (bestObjUpBound - relaxObj);
            double denominator = 0;
            for (LogiTask logiTask : logiTaskSet) {
                denominator = denominator
                        + ((logiTask.getTad() - logiTask.getTss()) * (logiTask.getTad() - logiTask.getTss()));
            }
            stepsize = numerator / denominator;
        } else {

            stepsize = numerator_for_stepsize / iter;
        }
////		计算stepSize的上界，确保对应各产品的所有新的lambda求和不会超过1
//		double[] sum1 = new double[para.numProd];
//		double[] sum2 = new double[para.numProd];
//		for (int p = 0; p < para.numProd; p++) {
//			sum1[p] = 0;
//			sum2[p] = para.objValFactor;
//			stepsizeUpperBound_Prod[p] = Float.MAX_VALUE;
//		}
//		for (LogiTask logiTask : logiTaskSet) {
//			sum1[logiTask.getProd().getId()] = sum1[logiTask.getProd().getId()] + logiTask.getTad() - logiTask.getTss();
//			sum2[logiTask.getProd().getId()] = sum2[logiTask.getProd().getId()] - logiTask.getLambda();
//		}
//		for (int p = 0; p < para.numProd; p++) {
//			if (sum1[p] > 0) {
//				stepsizeUpperBound_Prod[p] = sum2[p] / sum1[p];
//				stepsize = Math.min(stepsizeUpperBound_Prod[p], stepsize);
////				System.out.println("stepsizeUpperBound_Prod=(1-sum[lambda])/(sum[tad-tss])"+sum2[p]+"/"+sum1[p]+"="+stepsizeUpperBound_Prod[p]);
//			}
//		}
    }

    public void updateLambda() {
//		update lagrangian multiplers with current subgradient and stepsize  
        for (LogiTask logiTask : logiTaskSet) {
            logiTask.setLambda(Math.max(logiTask.getLambda() + stepsize * logiTask.getSubgrad(), 0));
        }
    }

    public double getObj(int[] prodSeq, ArrayList<ArrayList<Integer>> logiSeq) {

        for (int i = 0; i < logiSeq.size(); i++) {// 第i个路径
            double cur = 0;
            for (int j = 0; j < logiSeq.get(i).size(); j++) {// 第i个路径上的第j个位置
                para.ItoLogiTask.get(logiSeq.get(i).get(j)).tadHeur = cur
                        + para.getDistance(para.depot, para.ItoLogiTask.get(logiSeq.get(i).get(j)).getMach())
                        * (1 / para.vehSpeed + 2 * para.consumingRate * para.chargingRate);
                cur = para.ItoLogiTask.get(logiSeq.get(i).get(j)).tadHeur;
                cur = cur + para.getDistance(para.depot, para.ItoLogiTask.get(logiSeq.get(i).get(j)).getMach())
                        / para.vehSpeed;
            }
        }

        String key;
        double timeCur = 0;// 时间游标
        double sum_taf = 0;
        key = prodSeq[0] + "_" + 0 + "_" + 0;
        timeCur = Math.max(para.moldChangeTime, PMKtoLogiTask.get(key).tadHeur);// 第1个产品第1个机器第1个批次的物料的开始时间
        // 必然是其对应的tadHeur(已经算好了)
        for (int p = 0; p < prodSeq.length; p++) {// product从0开始编号
            for (int m = 0; m < para.numMach; m++) {// machine从0开始编号
                if (m > 0) {// 自第一个机器后，最早开工时间需要做一次减法
                    timeCur = timeCur - para.prodSet.get(prodSeq[p]).getCycleTime() * (para.demands[prodSeq[p]] - 1)
                            + para.getDistance(para.machSet.get(m - 1), para.machSet.get(m)) / para.lineSpeed;// +(para.jobshopLength/para.numMach)/para.lineSpeed;
                }
                for (int k = 0; k < para.K[prodSeq[p]][m]; k++) {// batch从0开始编号
                    key = prodSeq[p] + "_" + m + "_" + k;
                    // 通过key获得PMKtoLogiTask.get(key)是p,m,k对应的运输任务
                    timeCur = Math.max(timeCur, PMKtoLogiTask.get(key).tadHeur); // 真正的开始时间是物料到达和最早开工时间取大；
                    timeCur = timeCur + PMKtoLogiTask.get(key).getAssemTime(); // 加上装配时间等于该批次的完工时间，也是下一批次的最早开工时间；
                }
            }
//			timeCur=timeCur+(para.jobshopLength/para.numMach)/para.lineSpeed;
            taf[prodSeq[p]] = timeCur;
            sum_taf = sum_taf + timeCur * para.objValFactor;
            timeCur = timeCur - (para.numMach - 1) * para.cycTime[prodSeq[p]] + para.moldChangeTime
                    - para.getDistance(para.machSet.get(0), para.machSet.get(para.numMach - 1)) / para.lineSpeed;// -para.jobshopLength/para.lineSpeed;
        }

        maxTaf=Arrays.stream(taf).max().getAsDouble();
//        for (int p = 0; p < para.numProd ; p++) {
//            if(tafLB[p]>taf[prodSeq[p]]){
//                System.out.println("p"+p+", tafLB="+tafLB[p]+", taf="+taf[prodSeq[p]]);
//            }
//        }
        return sum_taf;
    }

    public void print3() {
        for (int i = 0; i < bestLogiSeq.size(); i++) {// 第i个路径
            double cur = 0;
            for (int j = 0; j < bestLogiSeq.get(i).size(); j++) {// 第i个路径上的第j个位置
                para.ItoLogiTask.get(bestLogiSeq.get(i).get(j)).tadHeur = cur
                        + para.getDistance(para.depot, para.ItoLogiTask.get(bestLogiSeq.get(i).get(j)).getMach())
                        * (1 / para.vehSpeed + 2 * para.consumingRate * para.chargingRate);
                cur = para.ItoLogiTask.get(bestLogiSeq.get(i).get(j)).tadHeur;
                cur = cur + para.getDistance(para.depot, para.ItoLogiTask.get(bestLogiSeq.get(i).get(j)).getMach())
                        / para.vehSpeed;
            }
        }

        String key;
        double timeCur = 0;// 时间游标
        double sum_taf = 0;
        key = bestProdSeq[0] + "_" + 0 + "_" + 0;
        timeCur = Math.max(para.moldChangeTime, PMKtoLogiTask.get(key).tadHeur);// 第1个产品第1个机器第1个批次的物料的开始时间
        // 必然是其对应的tadHeur(已经算好了)
        for (int p = 0; p < bestProdSeq.length; p++) {// product从0开始编号
            for (int m = 0; m < para.numMach; m++) {// machine从0开始编号
                if (m > 0) {// 自第一个机器后，最早开工时间需要做一次减法
                    timeCur = timeCur - para.prodSet.get(bestProdSeq[p]).getCycleTime() * (para.demands[bestProdSeq[p]] - 1)
                            + para.getDistance(para.machSet.get(m - 1), para.machSet.get(m)) / para.lineSpeed;// +(para.jobshopLength/para.numMach)/para.lineSpeed;
                }
                for (int k = 0; k < para.K[bestProdSeq[p]][m]; k++) {// batch从0开始编号
                    key = bestProdSeq[p] + "_" + m + "_" + k;
                    // 通过key获得PMKtoLogiTask.get(key)是p,m,k对应的运输任务
                    timeCur = Math.max(timeCur, PMKtoLogiTask.get(key).tadHeur); // 真正的开始时间是物料到达和最早开工时间取大；
                    timeCur = timeCur + PMKtoLogiTask.get(key).getAssemTime(); // 加上装配时间等于该批次的完工时间，也是下一批次的最早开工时间；
                }
            }
//			timeCur=timeCur+(para.jobshopLength/para.numMach)/para.lineSpeed;
            taf[bestProdSeq[p]] = timeCur;
            sum_taf = sum_taf + timeCur * para.objValFactor;
            timeCur = timeCur - (para.numMach - 1) * para.cycTime[bestProdSeq[p]] + para.moldChangeTime
                    - para.getDistance(para.machSet.get(0), para.machSet.get(para.numMach - 1)) / para.lineSpeed;// -para.jobshopLength/para.lineSpeed;
        }
        System.out.println("production sequence (Heur):");
        for (int i = 0; i < bestProdSeq.length; i++) {
            System.out.print(" -> " + bestProdSeq[i] + "(" + taf[bestProdSeq[i]] + ", P"
                    + para.prodSet.get(bestProdSeq[i]).getId() + ")");
        }
        System.out.println();

        System.out.println("logistics sequence (Heur): ");
        for (int i = 0; i < bestLogiSeq.size(); i++) {
            for (int j = 0; j < bestLogiSeq.get(i).size(); j++) {
                System.out.print(" -> " + bestLogiSeq.get(i).get(j) + "(" + tad[bestLogiSeq.get(i).get(j)] + ",P"
                        + para.ItoLogiTask.get(bestLogiSeq.get(i).get(j)).getProd().getId() + ",M"
                        + para.ItoLogiTask.get(bestLogiSeq.get(i).get(j)).getMach().getId() + ")");
            }
            System.out.println();
        }

    }

    public void print(int[] bestProdSeq2, ArrayList<ArrayList<Integer>> bestLogiSeq2) {

        double[] tad = new double[para.logiTaskSet.size()];
//		double[] tss=new double[para.logiTaskSet.size()];
        double[] taf = new double[para.numProd];
        for (int i = 0; i < bestLogiSeq2.size(); i++) {// 第i个路径
            double cur = 0;
            for (int j = 0; j < bestLogiSeq2.get(i).size(); j++) {// 第i个路径上的第j个位置
                tad[bestLogiSeq2.get(i).get(j)] = cur
                        + para.getDistance(para.depot, para.ItoLogiTask.get(bestLogiSeq2.get(i).get(j)).getMach())
                        * (1 / para.vehSpeed + 2 * para.consumingRate * para.chargingRate);
                cur = tad[bestLogiSeq2.get(i).get(j)];
                cur = cur + para.getDistance(para.depot, para.ItoLogiTask.get(bestLogiSeq2.get(i).get(j)).getMach())
                        / para.vehSpeed;
            }
        }

        String key;
        double timeCur = 0;// 时间游标
//        key = bestProdSeq2[0] + "_" + 0 + "_" + 0;
        timeCur = para.moldChangeTime;// 第1个产品第1个机器第1个批次的物料的开始时间 必然是其对应的tadHeur(已经算好了)
        for (int p = 0; p < bestProdSeq2.length; p++) {// product从0开始编号
            for (int m = 0; m < para.numMach; m++) {// machine从0开始编号
                if (m > 0) {// 自第一个机器后，最早开工时间需要做一次减法
//                    timeCur = timeCur - para.prodSet.get(prodSeq[p]).getCycleTime() * (para.demands[prodSeq[p]] - 1)
//                            + para.getDistance(para.machSet.get(m - 1), para.machSet.get(m)) / para.lineSpeed;// +(para.jobshopLength/para.numMach)/para.lineSpeed;
                    timeCur = timeCur
                            - para.prodSet.get(bestProdSeq2[p]).getCycleTime() * (para.demands[bestProdSeq2[p]] - 1)
                            + para.getDistance(para.machSet.get(m - 1), para.machSet.get(m)) / para.lineSpeed;
                }
                for (int k = 0; k < para.K[bestProdSeq2[p]][m]; k++) {// batch从0开始编号
                    key = bestProdSeq2[p] + "_" + m + "_" + k;
                    // 通过key获得PMKtoLogiTask.get(key)是p,m,k对应的运输任务
                    timeCur = Math.max(timeCur, tad[PMKtoLogiTask.get(key).getId()]); // 真正的开始时间是物料到达和最早开工时间取大；
                    timeCur = timeCur + PMKtoLogiTask.get(key).getAssemTime(); // 加上装配时间等于该批次的完工时间，也是下一批次的最早开工时间；
                }
            }
            taf[bestProdSeq2[p]] = timeCur;
        }

        System.out.println("production sequence (Heur):");
        for (int i = 0; i < bestProdSeq2.length; i++) {
            System.out.print(" -> " + bestProdSeq2[i] + "(" + taf[bestProdSeq2[i]] + ", P"
                    + para.prodSet.get(bestProdSeq2[i]).getId() + ")");
        }
        System.out.println();

        System.out.println("logistics sequence (Heur): ");
        for (int i = 0; i < bestLogiSeq2.size(); i++) {
            for (int j = 0; j < bestLogiSeq2.get(i).size(); j++) {
                System.out.print(" -> " + bestLogiSeq2.get(i).get(j) + "(" + tad[bestLogiSeq2.get(i).get(j)] + ",P"
                        + para.ItoLogiTask.get(bestLogiSeq2.get(i).get(j)).getProd().getId() + ",M"
                        + para.ItoLogiTask.get(bestLogiSeq2.get(i).get(j)).getMach().getId() + ")");
            }
            System.out.println();
        }

    }

    public void print2(int[] bestProdSeq2, ArrayList<ArrayList<Integer>> bestLogiSeq2) {// 输出tss和tad，看哪里的充电可以提前
        double[] tad = new double[para.logiTaskSet.size()];
        double[] tss = new double[para.logiTaskSet.size()];
        double[] tssNoLogi = new double[para.logiTaskSet.size()];
        double[] tssLate = new double[para.logiTaskSet.size()];
        double[] taf = new double[para.numProd];
        for (int i = 0; i < bestLogiSeq2.size(); i++) {// 第i个路径
            double cur = 0;
            for (int j = 0; j < bestLogiSeq2.get(i).size(); j++) {// 第i个路径上的第j个位置
                tad[bestLogiSeq2.get(i).get(j)] = cur
                        + para.getDistance(para.depot, para.ItoLogiTask.get(bestLogiSeq2.get(i).get(j)).getMach())
                        * (1 / para.vehSpeed + 2 * para.consumingRate * para.chargingRate);
                cur = tad[bestLogiSeq2.get(i).get(j)];
                cur = cur + para.getDistance(para.depot, para.ItoLogiTask.get(bestLogiSeq2.get(i).get(j)).getMach())
                        / para.vehSpeed;
            }
        }

        String key;
        double timeCur = 0;// 时间游标
        key = bestProdSeq2[0] + "_" + 0 + "_" + 0;
        timeCur = tad[PMKtoLogiTask.get(key).getId()];// 第1个产品第1个机器第1个批次的物料的开始时间 必然是其对应的tadHeur(已经算好了)
        for (int p = 0; p < bestProdSeq2.length; p++) {// product从0开始编号
            for (int m = 0; m < para.numMach; m++) {// machine从0开始编号
                if (m > 0) {// 自第一个机器后，最早开工时间需要做一次减法
                    timeCur = timeCur
                            - para.prodSet.get(bestProdSeq2[p]).getCycleTime() * (para.demands[bestProdSeq2[p]] - 1);
                }
                for (int k = 0; k < para.K[bestProdSeq2[p]][m]; k++) {// batch从0开始编号
                    key = bestProdSeq2[p] + "_" + m + "_" + k;
                    // 通过key获得PMKtoLogiTask.get(key)是p,m,k对应的运输任务
                    timeCur = Math.max(timeCur, tad[PMKtoLogiTask.get(key).getId()]); // 真正的开始时间是物料到达和最早开工时间取大；
                    tss[PMKtoLogiTask.get(key).getId()] = timeCur;
                    timeCur = timeCur + PMKtoLogiTask.get(key).getAssemTime(); // 加上装配时间等于该批次的完工时间，也是下一批次的最早开工时间；
                }
            }
            taf[bestProdSeq2[p]] = timeCur;
        }

        double[][] tss_PM = new double[para.numProd][para.numMach];// 产品p在机器m上的基准时间
        double[] tss_P = new double[para.numProd];// 产品p的基准时间
        tss_P[bestProdSeq2[0]] = para.getDistance(para.depot, para.machSet.get(0)) * 2 * para.consumingRate
                * para.chargingRate + para.getDistance(para.depot, para.machSet.get(0)) / para.vehSpeed;
        para.prodSet.get(bestProdSeq2[0]).setTaf(tss_P[bestProdSeq2[0]]
                + (para.demands[bestProdSeq2[0]] + para.numMach - 1) * para.cycTime[bestProdSeq2[0]]);
        for (int i = 1; i < para.numProd; i++) { // 下一个产品的开始时间由上一个产品的开始时间决定
            tss_P[bestProdSeq2[i]] = tss_P[bestProdSeq2[i - 1]]
                    + (para.demands[bestProdSeq2[i - 1]] + para.numMach - 1) * para.cycTime[bestProdSeq2[i - 1]];
            para.prodSet.get(bestProdSeq2[i]).setTaf(tss_P[bestProdSeq2[i]]
                    + (para.demands[bestProdSeq2[i]] + para.numMach - 1) * para.cycTime[bestProdSeq2[i]]);
        }
        for (int i = 0; i < para.numProd; i++) {
            for (int j = 0; j < para.numMach; j++) {// machine从0开始编号
                tss_PM[bestProdSeq2[i]][j] = tss_P[bestProdSeq2[i]] + j * para.cycTime[bestProdSeq2[i]];
            }
        }

        // 更新tss
        double tss_PMK = 0;
        for (int i = 0; i < para.numProd; i++) {// product从0开始编号
            for (int j = 0; j < para.numMach; j++) {// machine从0开始编号
                for (int k = 0; k < para.K[bestProdSeq2[i]][j]; k++) {// batch从0开始编号
                    tss_PMK = tss_PM[bestProdSeq2[i]][j]
                            + k * para.quanProdMach[bestProdSeq2[i]][j] * para.cycTime[bestProdSeq2[i]];
                    tssNoLogi[PMKtoLogiTask.get(bestProdSeq2[i] + "_" + j + "_" + k).getId()] = tss_PMK;
                    tssLate[PMKtoLogiTask.get(bestProdSeq2[i] + "_" + j + "_" + k).getId()] = taf[bestProdSeq2[i]]
                            + ((k + 1 - 1) * para.quanProdMach[bestProdSeq2[i]][j] - para.numMach + (j + 1)
                            - para.demands[bestProdSeq2[i]]) * para.cycTime[bestProdSeq2[i]];
                }
            }
        }

        System.out.println("production sequence:");
        for (int i = 0; i < bestProdSeq2.length; i++) {
            System.out.print(" -> " + bestProdSeq2[i] + "(" + taf[bestProdSeq2[i]] + ", P"
                    + para.prodSet.get(bestProdSeq2[i]).getId() + ")");
        }
        System.out.println();

        System.out.println("logistics sequence: ");
        for (int i = 0; i < bestLogiSeq2.size(); i++) {
            for (int j = 0; j < bestLogiSeq2.get(i).size(); j++) {
                System.out.print(" -> " + bestLogiSeq2.get(i).get(j) + "(" + tad[bestLogiSeq2.get(i).get(j)] + ","
                        + tss[bestLogiSeq2.get(i).get(j)] + ",P"
                        + para.ItoLogiTask.get(bestLogiSeq2.get(i).get(j)).getProd().getId() + "M"
                        + para.ItoLogiTask.get(bestLogiSeq2.get(i).get(j)).getMach().getId() + ")");
            }
            System.out.println();
        }

        // 生产信息输出
        System.out.printf("%-16s", "ProdID");
        for (int p = 0; p < para.prodSet.size(); p++) {
            System.out.printf("%10d", bestProdSeq2[p]);
        }
        System.out.println();
        System.out.printf("%-16s", "taf");
        for (int p = 0; p < para.prodSet.size(); p++) {
            System.out.printf("%10.2f", taf[bestProdSeq2[p]]);
        }
        System.out.println();
        System.out.println();

        // 物流信息输出

        for (int i = 0; i < bestLogiSeq2.size(); i++) {
            System.out.printf("%-16s", "LogiID");
            for (int j = 0; j < bestLogiSeq2.get(i).size(); j++) {
                System.out.printf("%10s",
                        bestLogiSeq2.get(i).get(j) + "(P"
                                + para.ItoLogiTask.get(bestLogiSeq2.get(i).get(j)).getProd().getId() + "M"
                                + para.ItoLogiTask.get(bestLogiSeq2.get(i).get(j)).getMach().getId() + ")");
            }
            System.out.println();
            System.out.printf("%-16s", "tad");
            for (int j = 0; j < bestLogiSeq2.get(i).size(); j++) {
                System.out.printf("%10.2f", tad[bestLogiSeq2.get(i).get(j)]);
            }
            System.out.println();
            System.out.printf("%-16s", "tss");
            for (int j = 0; j < bestLogiSeq2.get(i).size(); j++) {
                System.out.printf("%10.2f", tss[bestLogiSeq2.get(i).get(j)]);
            }
            System.out.println();
            System.out.printf("%-16s", "tssNoLogi");
            for (int j = 0; j < bestLogiSeq2.get(i).size(); j++) {
                System.out.printf("%10.2f", tssNoLogi[bestLogiSeq2.get(i).get(j)]);
            }
            System.out.println();
            System.out.printf("%-16s", "tssLate");
            for (int j = 0; j < bestLogiSeq2.get(i).size(); j++) {
                System.out.printf("%10.2f", tssLate[bestLogiSeq2.get(i).get(j)]);
            }
            System.out.println();
            System.out.println();
        }

    }
}
